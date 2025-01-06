// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Transporter drives the time integration of transport equations
  \details   Transporter drives the time integration of transport equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/transporter.ci.
*/
// *****************************************************************************

#include <string>
#include <iomanip>
#include <cstddef>
#include <unordered_set>

#include "Transporter.hpp"
#include "Fields.hpp"
#include "UniPDF.hpp"
#include "PDFWriter.hpp"
#include "ContainerUtil.hpp"
#include "LoadDistributor.hpp"
#include "ExodusIIMeshReader.hpp"
#include "InciterConfig.hpp"
#include "DiagWriter.hpp"
#include "Diagnostics.hpp"
#include "Integrals.hpp"
#include "Callback.hpp"
#include "Problems.hpp"

#include "NoWarning/inciter.decl.h"
#include "NoWarning/partitioner.decl.h"

extern CProxy_Main mainProxy;

namespace inciter {

extern ctr::Config g_cfg;
extern int g_nrestart;

}

using inciter::Transporter;

Transporter::Transporter() :
  m_input{ g_cfg.get< tag::input >() },
  m_nchare( m_input.size() ),
  m_ncit( m_nchare.size(), 0 ),
  m_nload( 0 ),
  m_npart( 0 ),
  m_nstat( 0 ),
  m_ndisc( 0 ),
  m_nchk( 0 ),
  m_ncom( 0 ),
  m_nt0refit( m_nchare.size(), 0 ),
  m_ndtrefit( m_nchare.size(), 0 ),
  m_noutrefit( m_nchare.size(), 0 ),
  m_noutderefit( m_nchare.size(), 0 ),
  m_nelem( m_nchare.size() ),
  m_finished( m_nchare.size(), 0 ),
  m_meshvol( m_nchare.size() ),
  m_minstat( m_nchare.size() ),
  m_maxstat( m_nchare.size() ),
  m_avgstat( m_nchare.size() ),
  m_progMesh( g_cfg.get< tag::feedback >(), ProgMeshPrefix, ProgMeshLegend ),
  m_progWork( g_cfg.get< tag::feedback >(), ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  const auto nstep = g_cfg.get< tag::nstep >();
  const auto t0 = g_cfg.get< tag::t0 >();
  const auto term = g_cfg.get< tag::term >();
  const auto constdt = g_cfg.get< tag::dt >();

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the constant time
  // step size (if that is used) is smaller than the duration of the time to be
  // simulated, we have work to do, otherwise, finish right away. If a constant
  // dt is not used, that part of the logic is always true as the default
  // constdt is zero.
  if ( nstep != 0 && term > t0 && constdt < term-t0 ) {

    // Enable SDAG waits for collecting mesh statistics
    thisProxy.wait4stat();

    // Configure and write diagnostics file header
    diagHeader();

    // Configure and write integrals file header
    integralsHeader();

    // Create mesh partitioner AND boundary condition object group
    createPartitioner();

  } else finish();      // stop if no time stepping requested
}

Transporter::Transporter( CkMigrateMessage* m ) :
  CBase_Transporter( m ),
  m_progMesh( g_cfg.get< tag::feedback >(), ProgMeshPrefix, ProgMeshLegend ),
  m_progWork( g_cfg.get< tag::feedback >(), ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Migrate constructor: returning from a checkpoint
//! \param[in] m Charm++ migrate message
// *****************************************************************************
{
  auto print = tk::Print();
  print << "\nXyst> Restarted from checkpoint\n";
  inthead( print );
}

bool
Transporter::matchBCs( std::map< int, std::vector< std::size_t > >& bnd )
// *****************************************************************************
// Verify that side sets specified in the control file exist in mesh file
//! \details This function does two things: (1) it verifies that the side
//!   sets used in the input file (either to which boundary conditions (BC)
//!   are assigned or listed as field output by the user in the
//!   input file) all exist among the side sets read from the input mesh
//!   file and errors out if at least one does not, and (2) it matches the
//!   side set ids at which the user has configured BCs (or listed as an output
//!   surface) to side set ids read from the mesh file and removes those face
//!   and node lists associated to side sets that the user did not set BCs or
//!   listed as field output on (as they will not need processing further since
//!   they will not be used).
//! \param[in,out] bnd Node or face lists mapped to side set ids
//! \return True if sidesets have been used and found in mesh
// *****************************************************************************
{
  std::unordered_set< int > usersets;

  // Collect side sets at which BCs are set

  for (const auto& s : g_cfg.get< tag::bc_dir >()) {
    if (!s.empty()) usersets.insert( s[0] );
  }

  for (auto s : g_cfg.get< tag::bc_sym >()) usersets.insert( s );

  for (auto s : g_cfg.get< tag::bc_far >()) usersets.insert( s );
 
  for (const auto& s : g_cfg.get< tag::bc_pre >()) {
    if (!s.empty()) usersets.insert( s[0] );
  }
 
  for (const auto& s : g_cfg.get< tag::pre_bc_dir >()) {
    if (!s.empty()) usersets.insert( s[0] );
  }

  for (auto s : g_cfg.get< tag::pre_bc_sym >()) usersets.insert( s );

  // Add sidesets requested for field output
  for (auto s : g_cfg.get< tag::fieldout >()) usersets.insert( s );
  // Add sidesets requested for integral output
  for (auto s : g_cfg.get< tag::integout >()) usersets.insert( s );

  // Find user-configured side set ids among side sets read from mesh file
  std::unordered_set< int > sidesets_used;
  for (auto i : usersets) {       // for all side sets used in control file
    if (bnd.find(i) != end(bnd))  // used set found among side sets in file
      sidesets_used.insert( i );  // store side set id configured as BC
    else {
      Throw( "Boundary conditions specified on side set " + std::to_string(i) +
             " which does not exist in mesh file" );
    }
  }
 
  // Remove sidesets not used (will not process those further)
  tk::erase_if( bnd, [&]( auto& item ) {
    return sidesets_used.find( item.first ) == end(sidesets_used);
  });

  return not bnd.empty();
}

void
Transporter::createPartitioner()
// *****************************************************************************
// Create mesh partitioner AND boundary conditions group
// *****************************************************************************
{
  // cppcheck-suppress unreadVariable
  auto print = tk::Print();

  // Create partitioner callbacks (order important)
  tk::PartitionerCallback cbp {{
      CkCallback( CkReductionTarget(Transporter,queriedPart), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,respondedPart), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,partitioned), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refinserted), thisProxy )
  }};

  // Create refiner callbacks (order important)
  tk::RefinerCallback cbr {{
      CkCallback( CkReductionTarget(Transporter,queriedRef), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,respondedRef), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,compatibility), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,bndint), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,matched), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create sorter callbacks (order important)
  tk::SorterCallback cbs {{
      CkCallback( CkReductionTarget(Transporter,queried), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,responded), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,discinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,workinserted), thisProxy )
  }};

  // Start timer measuring preparation of mesh(es) for partitioning
  m_timer[ TimerTag::MESH_READ ];

  ErrChk( !m_input.empty(), "No input mesh" );

  // Start preparing mesh(es)
  print.section( "Reading mesh" + std::string(m_input.size()>1?"es":"") );

  // Read boundary (side set) data from a list of input mesh files
  std::size_t meshid = 0;
  for (const auto& filename : m_input) {
    // Create mesh reader for reading side sets from file
    tk::ExodusIIMeshReader mr( filename );

    // Read out total number of mesh points from mesh file
    m_npoin.push_back( mr.npoin() );

    std::map< int, std::vector< std::size_t > > bface;
    std::map< int, std::vector< std::size_t > > faces;
    std::map< int, std::vector< std::size_t > > bnode;

    // Read boundary-face connectivity on side sets
    mr.readSidesetFaces( bface, faces );

    bool bcs_set = false;
    // Read node lists on side sets
    bnode = mr.readSidesetNodes();
    // Verify boundarty condition (BC) side sets used exist in mesh file
    bcs_set = matchBCs( bnode );
    bcs_set = bcs_set || matchBCs( bface );

    // Warn on no BCs
    if (!bcs_set) print << "\n>>> WARNING: No boundary conditions set\n\n";

    // Create empty discretization chare array
    m_discretization.push_back( CProxy_Discretization::ckNew() );
    CkArrayOptions opt;
    opt.bindTo( m_discretization.back() );

    // Create empty discretization scheme chare array (bound to discretization)
    CProxy_RieCG riecg;
    CProxy_LaxCG laxcg;
    CProxy_ZalCG zalcg;
    CProxy_KozCG kozcg;
    CProxy_ChoCG chocg;
    CProxy_LohCG lohcg;
    tk::CProxy_ConjugateGradients cgpre, cgmom;
    const auto& solver = g_cfg.get< tag::solver >();
    if (solver == "riecg") {
      m_riecg.push_back( CProxy_RieCG::ckNew(opt) );
      riecg = m_riecg.back();
    }
    else if (solver == "laxcg") {
      m_laxcg.push_back( CProxy_LaxCG::ckNew(opt) );
      laxcg = m_laxcg.back();
    }
    else if (solver == "zalcg") {
      m_zalcg.push_back( CProxy_ZalCG::ckNew(opt) );
      zalcg = m_zalcg.back();
    }
    else if (solver == "kozcg") {
      m_kozcg.push_back( CProxy_KozCG::ckNew(opt) );
      kozcg = m_kozcg.back();
    }
    else if (solver == "chocg") {
      m_chocg.push_back( CProxy_ChoCG::ckNew(opt) );
      chocg = m_chocg.back();
      m_cgpre.push_back( tk::CProxy_ConjugateGradients::ckNew(opt) );
      cgpre = m_cgpre.back();
      m_cgmom.push_back( tk::CProxy_ConjugateGradients::ckNew(opt) );
      cgmom = m_cgmom.back();
    }
    else if (solver == "lohcg") {
      m_lohcg.push_back( CProxy_LohCG::ckNew(opt) );
      lohcg = m_lohcg.back();
      m_cgpre.push_back( tk::CProxy_ConjugateGradients::ckNew(opt) );
      cgpre = m_cgpre.back();
    }
    else {
      Throw( "Unknown solver: " + solver );
    }

    // Create empty mesh refiner chare array (bound to discretization)
    m_refiner.push_back( CProxy_Refiner::ckNew(opt) );
    // Create empty mesh sorter Charm++ chare array (bound to discretization)
    m_sorter.push_back( CProxy_Sorter::ckNew(opt) );

    // Create MeshWriter chare group for mesh
    m_meshwriter.push_back(
      tk::CProxy_MeshWriter::ckNew(
        g_cfg.get< tag::benchmark >(), m_input.size() ) );

    // Create mesh partitioner Charm++ chare nodegroup for all meshes
    m_partitioner.push_back(
      CProxy_Partitioner::ckNew( meshid, filename, cbp, cbr, cbs,
        thisProxy, m_refiner.back(), m_sorter.back(), m_meshwriter.back(),
        m_discretization.back(), riecg, laxcg, zalcg, kozcg, chocg, lohcg,
        cgpre, cgmom, bface, faces, bnode ) );

    ++meshid;
  }
}

void
Transporter::load( std::size_t meshid, std::size_t nelem )
// *****************************************************************************
// Reduction target: the mesh has been read from file on all PEs
//! \param[in] meshid Mesh id (summed accross all compute nodes)
//! \param[in] nelem Number of mesh elements per mesh (summed across all
//!    compute nodes)
// *****************************************************************************
{
  meshid /= static_cast< std::size_t >( CkNumNodes() );
  Assert( meshid < m_nelem.size(), "MeshId indexing out" );
  m_nelem[meshid] = nelem;

  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare[meshid] = static_cast<int>(
    tk::linearLoadDistributor(
       g_cfg.get< tag::virt >(),
       m_nelem[meshid], CkNumPes(), chunksize, remainder ) );

  // Store sum of meshids (across all chares, key) for each meshid (value).
  // This is used to look up the mesh id after collectives that sum their data.
  m_meshid[ static_cast<std::size_t>(m_nchare[meshid])*meshid ] = meshid;
  Assert( meshid < m_nelem.size(), "MeshId indexing out" );

  // Partition first mesh
  if (meshid == 0) {
    m_timer[ TimerTag::MESH_PART ];  // start timer measuring mesh partitioning
    m_partitioner[0].partition( m_nchare[0] );
  }

  if (++m_nload == m_nelem.size()) {     // all meshes have been loaded
    m_nload = 0;
    auto print = tk::Print();

    auto& timer = tk::ref_find( m_timer, TimerTag::MESH_READ );
    timer.second = timer.first.dsec();
    print << "Mesh read time: " + std::to_string( timer.second ) + " sec\n";

    // Print out mesh partitioning configuration
    print.section( "Partitioning mesh" );
    print.item( "Partitioner", g_cfg.get< tag::part >() );
    print.item( "Virtualization", g_cfg.get< tag::virt >() );
    // Print out initial mesh statistics
    meshstat( "Mesh read from file" );

    // Tell meshwriter the total number of chares
    m_meshwriter[meshid].nchare( m_nchare[meshid] );

    // Query number of initial mesh refinement steps
    int nref = 0;
    if (g_cfg.get< tag::href_t0 >()) {
      nref = static_cast<int>(g_cfg.get< tag::href_init >().size());
    }

    // Query if PE-local reorder is configured
    int nreord = 0;
    if (g_cfg.get< tag::reorder >()) nreord = m_nchare[0];

    print << '\n';
    m_progMesh.start( print, "Preparing mesh", {{ CkNumPes(), CkNumPes(), nref,
      m_nchare[0], m_nchare[0], nreord, nreord }} );
  }
}

void
Transporter::partitioned( std::size_t meshid )
// *****************************************************************************
// Reduction target: a mesh has been partitioned
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  if (++m_npart == m_nelem.size()) {     // all meshes have been partitioned
    m_npart = 0;
    auto& timer = tk::ref_find( m_timer, TimerTag::MESH_PART );
    timer.second = timer.first.dsec();
    m_timer[ TimerTag::MESH_DIST ];  // start timer measuring mesh distribution
  } else { // partition next mesh
    m_partitioner[meshid+1].partition( m_nchare[meshid+1] );
  }
}

void
Transporter::distributed( std::size_t meshid )
// *****************************************************************************
// Reduction target: all compute nodes have distributed their mesh after
// partitioning
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_partitioner[meshid].refine();
  auto& timer = tk::ref_find( m_timer, TimerTag::MESH_DIST );
  timer.second = timer.first.dsec();
}

void
Transporter::refinserted( std::size_t meshid, std::size_t error )
// *****************************************************************************
// Reduction target: all compute nodes have created the mesh refiners
//! \param[in] meshid Mesh id (aggregated across all compute nodes with operator
//!   max)
//! \param[in] error Error code (aggregated across all compute nodes with
//!   operator max)
// *****************************************************************************
{
  if (error) {

    tk::Print() <<
        "\n>>> ERROR: A worker chare was not assigned any mesh "
        "elements after distributing mesh " + std::to_string(meshid) +
        ". This can happen in SMP-mode with a large +ppn "
        "parameter (number of worker threads per logical node) and is "
        "most likely the fault of the mesh partitioning algorithm not "
        "tolerating the case when it is asked to divide the "
        "computational domain into a number of partitions different "
        "than the number of ranks it is called on, i.e., in case of "
        "overdecomposition and/or calling the partitioner in SMP mode "
        "with +ppn larger than 1. Solution 1: Try a different "
        "partitioning algorithm. Solution 2: Decrease +ppn.";
    finish( meshid );

  } else {

    m_refiner[meshid].doneInserting();

  }
}

void
Transporter::queriedRef( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have queried their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_refiner[meshid].response();
}

void
Transporter::respondedRef( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have setup their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_refiner[meshid].refine();
}

void
Transporter::compatibility( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have received a round of edges,
// and have run their compatibility algorithm
//! \param[in] meshid Mesh id (aggregated across all chares using operator max)
//! \details This is called iteratively, until convergence by Refiner. At this
//!   point all Refiner chares have received a round of edge data (tags whether
//!   an edge needs to be refined, etc.), and applied the compatibility
//!   algorithm independent of other Refiner chares. We keep going until the
//!   mesh is no longer modified by the compatibility algorithm, based on a new
//!   round of edge data communication started in Refiner::comExtra().
// *****************************************************************************
{
  m_refiner[meshid].correctref();
}

void
Transporter::matched( std::size_t summeshid,
                      std::size_t nextra,
                      std::size_t nref,
                      std::size_t nderef,
                      std::size_t sumrefmode )
// *****************************************************************************
// Reduction target: all Refiner chares have matched/corrected the tagging
// of chare-boundary edges, all chares are ready to perform refinement.
//! \param[in] summeshid Mesh id (summed across all chares)
//! \param[in] nextra Sum (across all chares) of the number of edges on each
//!   chare that need correction along chare boundaries
//! \param[in] nref Sum of number of refined tetrahedra across all chares.
//! \param[in] nderef Sum of number of derefined tetrahedra across all chares.
//! \param[in] sumrefmode Sum of contributions from all chares, encoding
//!   refinement mode of operation.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );

  // If at least a single edge on a chare still needs correction, do correction,
  // otherwise, this mesh refinement step is complete
  if (nextra > 0) {

    ++m_ncit[meshid];
    m_refiner[meshid].comExtra();

  } else {

    tk::Print print;

    // decode refmode
    auto refmode = static_cast< Refiner::RefMode >(
                     sumrefmode / static_cast<std::size_t>(m_nchare[meshid]) );

    if (refmode == Refiner::RefMode::T0REF) {

      if (!g_cfg.get< tag::feedback >()) {
        const auto& initref = g_cfg.get< tag::href_init >();
        print << '\n';
        print.diag( { "meshid", "t0ref", "type", "nref", "nderef", "ncorr" },
                    { std::to_string(meshid),
                      std::to_string(m_nt0refit[meshid]),
                      initref[ m_nt0refit[ meshid ] ],
                      std::to_string(nref),
                      std::to_string(nderef),
                      std::to_string(m_ncit[meshid]) } );
        ++m_nt0refit[meshid];
        if (m_nt0refit[meshid] == initref.size()) print << '\n';
      }
      m_progMesh.inc< REFINE >( print );

    } else if (refmode == Refiner::RefMode::DTREF) {

      print.diag( { "meshid", "dtref", "type", "nref", "nderef", "ncorr" },
                  { std::to_string(meshid),
                    std::to_string(++m_ndtrefit[meshid]),
                    "error",
                    std::to_string(nref),
                    std::to_string(nderef),
                    std::to_string(m_ncit[meshid]) } );

    } else Throw( "RefMode not implemented" );

    m_ncit[meshid] = 0;
    m_refiner[meshid].perform();

  }
}

void
Transporter::bndint( tk::real sx, tk::real sy, tk::real sz, tk::real cb,
                     tk::real summeshid )
// *****************************************************************************
// Compute surface integral across the whole problem and perform leak-test
//! \param[in] sx X component of vector summed
//! \param[in] sy Y component of vector summed
//! \param[in] sz Z component of vector summed
//! \param[in] cb Invoke callback if positive
//! \param[in] summeshid Mesh id (summed accross all chares)
//! \details This function aggregates partial surface integrals across the
//!   boundary faces of the whole problem. After this global sum a
//!   non-zero vector result indicates a leak, e.g., a hole in the boundary,
//!   which indicates an error in the boundary face data structures used to
//!   compute the partial surface integrals.
// *****************************************************************************
{
  /*auto meshid =*/tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  std::stringstream err;
  if (cb < 0.0) {
    err << "Mesh boundary leaky after mesh refinement step; this is due to a "
     "problem with updating the side sets used to specify boundary conditions "
     "on faces: ";
  } else if (cb > 0.0) {
    err << "Mesh boundary leaky during initialization; this is due to "
    "incorrect or incompletely specified boundary conditions for a given input "
    "mesh: ";
  }

  auto eps = 1.0e-10;
  if (std::abs(sx) > eps || std::abs(sy) > eps || std::abs(sz) > eps) {
    err << "Integral result must be a zero vector: " << std::setprecision(12) <<
           std::abs(sx) << ", " << std::abs(sy) << ", " << std::abs(sz) <<
           ", eps = " << eps;
    Throw( err.str() );
  }
}

void
Transporter::refined( std::size_t summeshid,
                      std::size_t nelem,
                      std::size_t npoin )
// *****************************************************************************
// Reduction target: all chares have refined their mesh
//! \param[in] summeshid Mesh id (summed accross all Refiner chares)
//! \param[in] nelem Total number of elements in mesh summed across the
//!   distributed mesh
//! \param[in] npoin Total number of mesh points summed across the distributed
//!   mesh. Note that in parallel this is larger than the number of points in
//!   the mesh, because the boundary nodes are multi-counted. But we only need
//!   an equal or larger than npoin for Sorter::setup, so this is okay.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );

  // Store new number of elements for initially refined mesh
  m_nelem[meshid] = nelem;

  m_sorter[meshid].doneInserting();
  m_sorter[meshid].setup( npoin );
}

void
Transporter::queriedPart( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Partitioner nodes have queried their mesh graphs
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_partitioner[meshid].response();
}

void
Transporter::respondedPart( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Partitioner nodes have responded with their mesh graphs
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_partitioner[meshid].load();
}

void
Transporter::queried( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Sorter chares have queried their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_sorter[meshid].response();
}

void
Transporter::responded( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Sorter chares have responded with their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_sorter[meshid].start();
}

void
Transporter::resized( std::size_t meshid )
// *****************************************************************************
// Reduction target: all worker chares have resized their own mesh data after
//! \param[in] meshid Mesh id
//! \note Only used for nodal schemes
// *****************************************************************************
{
  m_discretization[ meshid ].vol();

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    m_riecg[ meshid ].feop();
  }
  else if (solver == "laxcg") {
    m_laxcg[ meshid ].feop();
  }
  else if (solver == "zalcg") {
    m_zalcg[ meshid ].feop();
  }
  else if (solver == "kozcg") {
    m_kozcg[ meshid ].feop();
  }
  else if (solver == "chocg") {
    m_chocg[ meshid ].feop();
  }
  else if (solver == "lohcg") {
    m_lohcg[ meshid ].feop();
  }
  else {
    Throw( "Unknown solver: " + solver  );
  }
}

void
Transporter::discinserted( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Discretization chares have been inserted
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_discretization[ meshid ].doneInserting();
}

void
Transporter::meshstat( const std::string& header ) const
// *****************************************************************************
// Print out mesh statistics
//! \param[in] header Section header
// *****************************************************************************
{
  tk::Print print;

  print.section( header );

  if (m_nelem.size() > 1) {
    print.item( "Number of tetrahedra (per mesh)",tk::parameters(m_nelem) );
    print.item( "Number of points (per mesh)", tk::parameters(m_npoin) );
    print.item( "Number of work units (per mesh)", tk::parameters(m_nchare) );
  }

  print.item( "Total number of tetrahedra",
              std::accumulate( begin(m_nelem), end(m_nelem), 0UL ) );
  print.item( "Total number of points",
              std::accumulate( begin(m_npoin), end(m_npoin), 0UL ) );
  print.item( "Total number of work units",
              std::accumulate( begin(m_nchare), end(m_nchare), 0 ) );
}

void
Transporter::disccreated( std::size_t summeshid, std::size_t npoin )
// *****************************************************************************
// Reduction target: all Discretization constructors have been called
//! \param[in] summeshid Mesh id (summed accross all chares)
//! \param[in] npoin Total number of mesh points (summed across all chares)
//!  Note that as opposed to npoin in refined(), this npoin is not
//!  multi-counted, and thus should be correct in parallel.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );

  // Update number of mesh points for mesh, since it may have been refined
  if (g_cfg.get< tag::href_t0 >()) m_npoin[meshid] = npoin;

  if (++m_ndisc == m_nelem.size()) { // all Disc arrays have been created
    m_ndisc = 0;
    tk::Print print;
    m_progMesh.end( print );
    if (g_cfg.get< tag::href_t0 >()) {
      meshstat( "Mesh initially refined" );
    }
  }

  m_refiner[ meshid ].sendProxy();
  m_discretization[ meshid ].vol();

  m_discretization[0][0].npoin(
    std::accumulate( begin(m_npoin), end(m_npoin), 0UL ) );
}

void
Transporter::workinserted( std::size_t meshid )
// *****************************************************************************
// Reduction target: all worker (derived discretization) chares have been
// inserted
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    m_riecg[ meshid ].doneInserting();
  }
  else if (solver == "laxcg") {
    m_laxcg[ meshid ].doneInserting();
  }
  else if (solver == "zalcg") {
    m_zalcg[ meshid ].doneInserting();
  }
  else if (solver == "kozcg") {
    m_kozcg[ meshid ].doneInserting();
  }
  else if (solver == "chocg") {
    m_chocg[ meshid ].doneInserting();
    m_cgpre[ meshid ].doneInserting();
    m_cgmom[ meshid ].doneInserting();
  }
  else if (solver == "lohcg") {
    m_lohcg[ meshid ].doneInserting();
    m_cgpre[ meshid ].doneInserting();
  }
  else {
    Throw( "Unknown solver: " + solver );
  }
}

void
Transporter::diagHeader()
// *****************************************************************************
// Configure and write diagnostics file header
// *****************************************************************************
{
  // Output header for diagnostics output file
  tk::DiagWriter dw( g_cfg.get< tag::diag >(),
                     g_cfg.get< tag::diag_format >(),
                     g_cfg.get< tag::diag_precision >() );

  std::vector< std::string > d;

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg" ||
      solver == "laxcg" ||
      solver == "zalcg" ||
      solver == "kozcg")
  {

    // Collect variables names for integral/diagnostics output
    std::vector< std::string > var{ "r", "ru", "rv", "rw", "rE" };
    auto ncomp = g_cfg.get< tag::problem_ncomp >();
    for (std::size_t c=5; c<ncomp; ++c)
      var.push_back( "c" + std::to_string(c-5) );

    auto nv = var.size();

    // Add 'L2(var)' for all variables
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(" + var[i] + ')' );

    // Add L2-norm of the residuals
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(d" + var[i] + ')' );

    // Add total energy
    d.push_back( "mE" );

    // Augment diagnostics variables by error norms (if computed)
    if (problems::SOL()) {
      d.push_back( "L2(err:r)" );
      d.push_back( "L2(err:u)" );
      d.push_back( "L2(err:v)" );
      d.push_back( "L2(err:w)" );
      d.push_back( "L2(err:e)" );
      for (std::size_t i=5; i<nv; ++i) d.push_back( "L2(err:" + var[i] + ')' );
      d.push_back( "L1(err:r)" );
      d.push_back( "L1(err:u)" );
      d.push_back( "L1(err:v)" );
      d.push_back( "L1(err:w)" );
      d.push_back( "L1(err:e)" );
      for (std::size_t i=5; i<nv; ++i) d.push_back( "L1(err:" + var[i] + ')' );
    }

  }
  else if (solver == "chocg") {

    // query function to evaluate analytic solution (if defined)
    auto pressure_sol = problems::PRESSURE_SOL();

    // Collect variables names for integral/diagnostics output
    std::vector< std::string > var{ "p" };
    if (!pressure_sol) {
      var.push_back( "u" );
      var.push_back( "v" );
      var.push_back( "w" );
    }

    auto nv = var.size();

    // Add 'L2(var)' for all variables
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(" + var[i] + ')' );

    // Add L2-norm of the residuals
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(d" + var[i] + ')' );

    // Augment diagnostics variables by error norms (if computed)
    if (pressure_sol) {
      d.push_back( "L2(err:p)" );
      d.push_back( "L1(err:p)" );
    }

  }
  else if (solver == "lohcg") {

    // Collect variables names for integral/diagnostics output
    std::vector< std::string > var{ "p" };
    var.push_back( "u" );
    var.push_back( "v" );
    var.push_back( "w" );

    auto nv = var.size();

    // Add 'L2(var)' for all variables
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(" + var[i] + ')' );

    // Add L2-norm of the residuals
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(d" + var[i] + ')' );

  }
  else {
    Throw( "Unknown solver: " + solver );
  }

  // Write diagnostics header
  dw.header( d );
}

void
Transporter::integralsHeader()
// *****************************************************************************
// Configure and write integrals file header
// *****************************************************************************
{
  const auto& sidesets_integral = g_cfg.get< tag::integout >();

  if (sidesets_integral.empty()) return;

  auto filename = g_cfg.get< tag::output >() + ".int";
  tk::DiagWriter dw( filename,
                     g_cfg.get< tag::integout_format >(),
                     g_cfg.get< tag::integout_precision >() );

  // Collect variables names for integral output
  std::vector< std::string > var;
  const auto& reqv = g_cfg.get< tag::integout_integrals >();
  std::unordered_set< std::string > req( begin(reqv), end(reqv) );
  for (auto s : sidesets_integral) {
    if (req.count( "mass_flow_rate" )) {
      var.push_back( "mass_flow_rate:" + std::to_string(s) );
    }
    if (req.count( "force" )) {
      auto si = std::to_string( s );
      var.push_back( "force_x:" + si );
      var.push_back( "force_y:" + si );
      var.push_back( "force_z:" + si );
    }
  }

  // Write integrals header
  dw.header( var );
}

void
Transporter::totalvol( tk::real v, tk::real initial, tk::real summeshid )
// *****************************************************************************
// Reduction target summing total mesh volume across all workers
//! \param[in] v Mesh volume summed across the distributed mesh
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during setup, if zero, during time stepping.
//! \param[in] summeshid Mesh id (summed accross the distributed mesh)
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  m_meshvol[meshid] = v;

  if (initial > 0.0) {   // during initialization

    m_discretization[ meshid ].stat( v );

  } else {               // during AMR

    const auto& solver = g_cfg.get< tag::solver >();
    if (solver == "riecg") {
      m_riecg[ meshid ].resize_complete();
    }
    else if (solver == "laxcg") {
      m_laxcg[ meshid ].resize_complete();
    }
    else if (solver == "zalcg") {
      m_zalcg[ meshid ].resize_complete();
    }
    else if (solver == "kozcg") {
      m_kozcg[ meshid ].resize_complete();
    }
    else if (solver == "chocg") {
      m_chocg[ meshid ].resize_complete();
    }
    else if (solver == "lohcg") {
      m_lohcg[ meshid ].resize_complete();
    }
    else {
      Throw( "Unknown solver: " + solver );
    }

  }
}

void
Transporter::minstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3,
                      tk::real d4, tk::real d5, tk::real rmeshid )
// *****************************************************************************
// Reduction target yielding minimum mesh statistcs across all workers
//! \param[in] d0 Minimum mesh statistics collected over all chares
//! \param[in] d1 Minimum mesh statistics collected over all chares
//! \param[in] d2 Minimum mesh statistics collected over all chares
//! \param[in] d3 Minimum mesh statistics collected over all chares
//! \param[in] d4 Minimum mesh statistics collected over all chares
//! \param[in] d5 Minimum mesh statistics collected over all chares
//! \param[in] rmeshid Mesh id as a real
// *****************************************************************************
{
  auto meshid = static_cast<std::size_t>(rmeshid);

  m_minstat[meshid][0] = d0;  // minimum edge length
  m_minstat[meshid][1] = d1;  // minimum cell volume cubic root
  m_minstat[meshid][2] = d2;  // minimum number of elements on chare
  m_minstat[meshid][3] = d3;  // minimum number of points on chare
  m_minstat[meshid][4] = d4;  // minimum number of edges on chare
  m_minstat[meshid][5] = d5;  // minimum number of comm/total points on chare

  minstat_complete(meshid);
}

void
Transporter::maxstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3,
                      tk::real d4, tk::real d5, tk::real rmeshid )
// *****************************************************************************
// Reduction target yielding the maximum mesh statistics across all workers
//! \param[in] d0 Maximum mesh statistics collected over all chares
//! \param[in] d1 Maximum mesh statistics collected over all chares
//! \param[in] d2 Maximum mesh statistics collected over all chares
//! \param[in] d3 Maximum mesh statistics collected over all chares
//! \param[in] d4 Minimum mesh statistics collected over all chares
//! \param[in] d5 Minimum mesh statistics collected over all chares
//! \param[in] rmeshid Mesh id as a real
// *****************************************************************************
{
  auto meshid = static_cast<std::size_t>(rmeshid);

  m_maxstat[meshid][0] = d0;  // maximum edge length
  m_maxstat[meshid][1] = d1;  // maximum cell volume cubic root
  m_maxstat[meshid][2] = d2;  // maximum number of elements on chare
  m_maxstat[meshid][3] = d3;  // maximum number of points on chare
  m_maxstat[meshid][4] = d4;  // maximum number of edges on chare
  m_maxstat[meshid][5] = d5;  // maximum number of comm/total points on chare

  maxstat_complete(meshid);
}

void
Transporter::sumstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3,
                      tk::real d4, tk::real d5, tk::real d6, tk::real d7,
                      tk::real d8, tk::real summeshid )
// *****************************************************************************
// Reduction target yielding the sum mesh statistics across all workers
//! \param[in] d0 Sum mesh statistics collected over all chares
//! \param[in] d1 Sum mesh statistics collected over all chares
//! \param[in] d2 Sum mesh statistics collected over all chares
//! \param[in] d3 Sum mesh statistics collected over all chares
//! \param[in] d4 Sum mesh statistics collected over all chares
//! \param[in] d5 Sum mesh statistics collected over all chares
//! \param[in] d6 Sum mesh statistics collected over all chares
//! \param[in] d7 Sum mesh statistics collected over all chares
//! \param[in] d8 Sum mesh statistics collected over all chares
//! \param[in] summeshid Mesh id (summed accross the distributed mesh)
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  m_avgstat[meshid][0] = d1 / d0;  // avg edge length
  m_avgstat[meshid][1] = d3 / d2;  // avg cell volume cubic root
  m_avgstat[meshid][2] = d5 / d4;  // avg number of elements per chare
  m_avgstat[meshid][3] = d6 / d4;  // avg number of points per chare
  m_avgstat[meshid][4] = d7 / d4;  // avg number of edges per chare
  m_avgstat[meshid][5] = d8 / d4;  // avg number of comm/total points per chare

  sumstat_complete(meshid);
}

void
Transporter::pdfstat( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target yielding PDF of mesh statistics across all workers
//! \param[in] msg Serialized PDF
// *****************************************************************************
{
  std::size_t meshid;
  std::vector< tk::UniPDF > pdf;

  // Deserialize final PDF
  PUP::fromMem creator( msg->getData() );
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | pdf;
  delete msg;

  // cppcheck-suppress uninitvar
  auto id = std::to_string(meshid);

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfe( "mesh_edge_pdf." + id + ".txt" );
  // Output edgelength PDF
  // cppcheck-suppress containerOutOfBounds
  pdfe.writeTxt( pdf[0],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"edgelength"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfv( "mesh_vol_pdf." + id + ".txt" );
  // Output cell volume cubic root PDF
  // cppcheck-suppress containerOutOfBounds
  pdfv.writeTxt( pdf[1],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"V^{1/3}"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfn( "mesh_ntet_pdf." + id + ".txt" );
  // Output number of cells PDF
  // cppcheck-suppress containerOutOfBounds
  pdfn.writeTxt( pdf[2],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"ntets"}, 0, 0.0 } );

  pdfstat_complete(meshid);
}

void
Transporter::stat()
// *****************************************************************************
// Echo diagnostics on mesh statistics
// *****************************************************************************
{
  tk::Print print;

  if (++m_nstat == m_nelem.size()) {     // stats from all meshes have arrived
    m_nstat = 0;

    auto& t = tk::ref_find( m_timer, TimerTag::MESH_PART );
    print << '\n';
    print << "Mesh partitioning time: " + std::to_string(t.second) + " sec\n";
    t = tk::ref_find( m_timer, TimerTag::MESH_DIST );
    print << "Mesh distribution time: " + std::to_string(t.second) + " sec\n";


    for (std::size_t i=0; i<m_nelem.size(); ++i) {
      if (m_nelem.size() > 1) {
        print.section("Mesh " + std::to_string(i) + " distribution statistics");
      } else {
        print.section( "Mesh distribution statistics" );
      }
      print <<
        "min/max/avg edgelength = " +
        std::to_string( m_minstat[i][0] ) + " / " +
        std::to_string( m_maxstat[i][0] ) + " / " +
        std::to_string( m_avgstat[i][0] ) + "\n" +
        "min/max/avg V^(1/3) = " +
        std::to_string( m_minstat[i][1] ) + " / " +
        std::to_string( m_maxstat[i][1] ) + " / " +
        std::to_string( m_avgstat[i][1] ) + "\n" +
        "min/max/avg nelem = " +
        std::to_string( static_cast<std::size_t>(m_minstat[i][2]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_maxstat[i][2]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_avgstat[i][2]) ) + "\n" +
        "min/max/avg npoin = " +
        std::to_string( static_cast<std::size_t>(m_minstat[i][3]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_maxstat[i][3]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_avgstat[i][3]) ) + "\n" +
        "min/max/avg nedge = " +
        std::to_string( static_cast<std::size_t>(m_minstat[i][4]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_maxstat[i][4]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_avgstat[i][4]) ) + '\n' +
        "min/max/avg ncompoin/npoin = " +
        std::to_string( m_minstat[i][5] ) + " / " +
        std::to_string( m_maxstat[i][5] ) + " / " +
        std::to_string( m_avgstat[i][5] ) + '\n';
    }

    // Print out time integration header to screen
    inthead( print );

    m_progWork.start( print, "Preparing workers", {{ m_nchare[0] }} );

    // Create "derived-class" workers
    for (std::size_t i=0; i<m_nelem.size(); ++i) m_sorter[i].createWorkers();
  }
}

void
Transporter::boxvol( tk::real v, tk::real summeshid )
// *****************************************************************************
// Reduction target computing total volume of IC box(es)
//! \param[in] v Total volume within user-specified IC box(es)
//! \param[in] summeshid Mesh id as a real (summed accross the distributed mesh)
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );
  if (v > 0.0) tk::Print() << "IC-box-volume sum: " + std::to_string(v) << '\n';

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    m_riecg[ meshid ].setup( v );
  }
  else if (solver == "laxcg") {
    m_laxcg[ meshid ].setup( v );
  }
  else if (solver == "zalcg") {
    m_zalcg[ meshid ].setup( v );
  }
  else if (solver == "kozcg") {
    m_kozcg[ meshid ].setup( v );
  }
  else if (solver == "chocg") {
    m_chocg[ meshid ].setup( v );
  }
  else if (solver == "lohcg") {
    m_lohcg[ meshid ].setup( v );
  }
  else {
    Throw( "Unknown solver: " + solver );
  }

  // Turn on automatic load balancing
  if (++m_ncom == m_nelem.size()) { // all worker arrays have finished
    m_ncom = 0;
    tk::Print print;
    m_progWork.end( print );
    tk::CProxy_LBSwitch::ckNew();
  }
}

void
Transporter::inthead( const tk::Print& print )
// *****************************************************************************
// Print out time integration header to screen
//! \param[in] print Pretty printer object to use for printing
// *****************************************************************************
{
  const auto dea = g_cfg.get< tag::deactivate >();
  const auto solver = g_cfg.get< tag::solver >();
  const auto pre = solver == "chocg" ? 1 : 0;
  const auto theta = g_cfg.get< tag::theta >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto mom = solver == "chocg" and theta > eps ? 1 : 0;

  print.section( "Time integration" );
  print <<
  "Legend: it - iteration count\n"
  "         t - physics time\n"
  "        dt - physics time step size\n"
  "       ETE - estimated wall-clock time elapsed (h:m:s)\n"
  "       ETA - estimated wall-clock time for accomplishment (h:m:s)\n"
  "       EGT - estimated grind wall-clock time (1e-6sec/timestep)\n"
  "       EGP - estimated grind performance: wall-clock time "
                "(1e-6sec/DOF/timestep)\n"
  "       flg - status flags, legend:\n"
  "             f - field (volume and surface) output\n"
  "             i - integral output\n"
  "             d - diagnostics output\n"
  "             t - physics time history output\n"
  "             h - h-refinement\n"
  "             l - load balancing\n"
  "             c - checkpoint\n" << (dea ?
  "             e:x/y - x of y work units deactivated\n" : "") << (pre ?
  "             p:it - pressure linear solve iterations\n" : "") << (mom ?
  "             m:it - momentum/transport linear solve iterations\n" : "") <<
  "\n      it             t            dt        ETE        ETA        EGT"
  "           EGP  flg\n"
  "-----------------------------------------------------------------------"
  "-----------------\n";
}

void
Transporter::rhodiagnostics( CkReductionMsg* msg )
// *****************************************************************************
//  Reduction target collecting diagnostics from density-based solvers
//! \param[in] msg Serialized diagnostics vector aggregated across all PEs
// *****************************************************************************
{
  using namespace diagnostics;

  std::size_t meshid;
  std::size_t ncomp;
  std::vector< std::vector< tk::real > > d;

  // Deserialize diagnostics vector
  PUP::fromMem creator( msg->getData() );
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | ncomp;
  creator | d;
  delete msg;

  // cppcheck-suppress uninitvar
  // cppcheck-suppress unreadVariable
  auto id = std::to_string(meshid);

  Assert( ncomp > 0, "Number of scalar components must be positive");
  Assert( d.size() == NUMDIAG, "Diagnostics vector size mismatch" );

  // cppcheck-suppress unsignedLessThanZero
  for (std::size_t i=0; i<d.size(); ++i) {
     Assert( d[i].size() == ncomp, "Size mismatch at final stage of "
             "diagnostics aggregation for mesh " + id );
  }

  // Allocate storage for those diagnostics that are always computed
  std::vector< tk::real > diag( ncomp, 0.0 );

  // Finish computing the L2 norm of conserved variables
  for (std::size_t i=0; i<d[L2SOL].size(); ++i) {
    // cppcheck-suppress uninitvar
    diag[i] = sqrt( d[L2SOL][i] / m_meshvol[meshid] );
  }
 
  // Finish computing the L2 norm of the residuals
  std::vector< tk::real > l2res( d[L2RES].size(), 0.0 );
  for (std::size_t i=0; i<d[L2RES].size(); ++i) {
    // cppcheck-suppress uninitvar
    l2res[i] = std::sqrt( d[L2RES][i] / m_meshvol[meshid] );
    diag.push_back( l2res[i] );
  }

  // Append total energy
  diag.push_back( d[TOTALEN][0] );

  // Finish computing norms of the numerical - analytical solution
  if (problems::SOL()) {
    for (std::size_t i=0; i<d[L2ERR].size(); ++i) {
      // cppcheck-suppress uninitvar
      diag.push_back( std::sqrt( d[L2ERR][i] / m_meshvol[meshid] ) );
    }
    for (std::size_t i=0; i<d[L1ERR].size(); ++i) {
      // cppcheck-suppress uninitvar
      diag.push_back( d[L1ERR][i] / m_meshvol[meshid] );
    }
  }
 
  // Append diagnostics file at selected times
  auto filename = g_cfg.get< tag::diag >();
  if (m_nelem.size() > 1) filename += '.' + id;
  tk::DiagWriter dw( filename,
                     g_cfg.get< tag::diag_format >(),
                     g_cfg.get< tag::diag_precision >(),
                     std::ios_base::app );
  dw.write( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    // cppcheck-suppress uninitvar
    m_riecg[ meshid ].evalres( l2res );
  }
  else if (solver == "laxcg") {
    // cppcheck-suppress uninitvar
    m_laxcg[ meshid ].evalres( l2res );
  }
  else if (solver == "zalcg") {
    // cppcheck-suppress uninitvar
    m_zalcg[ meshid ].evalres( l2res );
  }
  else if (solver == "kozcg") {
    // cppcheck-suppress uninitvar
    m_kozcg[ meshid ].evalres( l2res );
  }
  else {
    Throw( "Unknown solver: " + solver );
  }
}

void
Transporter::prediagnostics( CkReductionMsg* msg )
// *****************************************************************************
//  Reduction target collecting diagnostics from pressure-based solvers
//! \param[in] msg Serialized diagnostics vector aggregated across all PEs
// *****************************************************************************
{
  using namespace diagnostics;

  std::size_t meshid;
  std::size_t ncomp;
  std::vector< std::vector< tk::real > > d;

  // Deserialize diagnostics vector
  PUP::fromMem creator( msg->getData() );
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | ncomp;
  creator | d;
  delete msg;

  // cppcheck-suppress uninitvar
  // cppcheck-suppress unreadVariable
  auto id = std::to_string(meshid);

  Assert( ncomp > 0, "Number of scalar components must be positive");
  Assert( d.size() == NUMDIAG, "Diagnostics vector size mismatch" );

  // cppcheck-suppress unsignedLessThanZero
  for (std::size_t i=0; i<d.size(); ++i) {
     Assert( d[i].size() == ncomp, "Size mismatch at final stage of "
             "diagnostics aggregation for mesh " + id );
  }

  // Allocate storage for those diagnostics that are always computed
  std::vector< tk::real > diag( ncomp, 0.0 );

  // Finish computing the L2 norm of conserved variables
  for (std::size_t i=0; i<d[L2SOL].size(); ++i) {
    // cppcheck-suppress uninitvar
    diag[i] = sqrt( d[L2SOL][i] / m_meshvol[meshid] );
  }

  // Finish computing the L2 norm of the residuals
  std::vector< tk::real > l2res( d[L2RES].size(), 0.0 );
  for (std::size_t i=0; i<d[L2RES].size(); ++i) {
    // cppcheck-suppress uninitvar
    l2res[i] = std::sqrt( d[L2RES][i] / m_meshvol[meshid] );
    diag.push_back( l2res[i] );
  }

  // Finish computing norms of the numerical - analytical solution
  if (problems::PRESSURE_SOL()) {
    for (std::size_t i=0; i<d[L2ERR].size(); ++i) {
      // cppcheck-suppress uninitvar
      diag.push_back( std::sqrt( d[L2ERR][i] / m_meshvol[meshid] ) );
    }
    for (std::size_t i=0; i<d[L1ERR].size(); ++i) {
      // cppcheck-suppress uninitvar
      diag.push_back( d[L1ERR][i] / m_meshvol[meshid] );
    }
  }

  // Append diagnostics file at selected times
  auto filename = g_cfg.get< tag::diag >();
  if (m_nelem.size() > 1) filename += '.' + id;
  tk::DiagWriter dw( filename,
                     g_cfg.get< tag::diag_format >(),
                     g_cfg.get< tag::diag_precision >(),
                     std::ios_base::app );
  dw.write( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "chocg") {
    // cppcheck-suppress uninitvar
    m_chocg[ meshid ].evalres( l2res );
  }
  else {
    Throw( "Unknown solver: " + solver );
  }
}

void
Transporter::acdiagnostics( CkReductionMsg* msg )
// *****************************************************************************
//  Reduction target collecting diagnostics from artificial compressibility
//  solvers
//! \param[in] msg Serialized diagnostics vector aggregated across all PEs
// *****************************************************************************
{
  using namespace diagnostics;

  std::size_t meshid;
  std::size_t ncomp;
  std::vector< std::vector< tk::real > > d;

  // Deserialize diagnostics vector
  PUP::fromMem creator( msg->getData() );
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | ncomp;
  creator | d;
  delete msg;

  // cppcheck-suppress uninitvar
  // cppcheck-suppress unreadVariable
  auto id = std::to_string(meshid);

  Assert( ncomp > 0, "Number of scalar components must be positive");
  Assert( d.size() == NUMDIAG, "Diagnostics vector size mismatch" );

  // cppcheck-suppress unsignedLessThanZero
  for (std::size_t i=0; i<d.size(); ++i) {
     Assert( d[i].size() == ncomp, "Size mismatch at final stage of "
             "diagnostics aggregation for mesh " + id );
  }

  // Allocate storage for those diagnostics that are always computed
  std::vector< tk::real > diag( ncomp, 0.0 );

  // Finish computing the L2 norm of conserved variables
  for (std::size_t i=0; i<d[L2SOL].size(); ++i) {
    // cppcheck-suppress uninitvar
    diag[i] = sqrt( d[L2SOL][i] / m_meshvol[meshid] );
  }

  // Finish computing the L2 norm of the residuals
  std::vector< tk::real > l2res( d[L2RES].size(), 0.0 );
  for (std::size_t i=0; i<d[L2RES].size(); ++i) {
    // cppcheck-suppress uninitvar
    l2res[i] = std::sqrt( d[L2RES][i] / m_meshvol[meshid] );
    diag.push_back( l2res[i] );
  }

  // Append diagnostics file at selected times
  auto filename = g_cfg.get< tag::diag >();
  if (m_nelem.size() > 1) filename += '.' + id;
  tk::DiagWriter dw( filename,
                     g_cfg.get< tag::diag_format >(),
                     g_cfg.get< tag::diag_precision >(),
                     std::ios_base::app );
  dw.write( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "lohcg") {
    // cppcheck-suppress uninitvar
    m_lohcg[ meshid ].evalres( l2res );
  }
  else {
    Throw( "Unknown solver: " + solver );
  }
}

void
Transporter::integrals( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target optionally collecting integrals
//! \param[in] msg Serialized integrals aggregated across all PEs
// *****************************************************************************
{
  using namespace integrals;

  // cppcheck-suppress unassignedVariable
  std::size_t meshid;
  std::vector< std::map< int, tk::real > > d;

  // Deserialize integrals vector
  PUP::fromMem creator( msg->getData() );
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | d;
  delete msg;

  const auto& sidesets_integral = g_cfg.get< tag::integout >();
  // cppcheck-suppress
  if (not sidesets_integral.empty()) {

    Assert( d.size() == NUMINT, "Integrals vector size mismatch" );

    // Collect integrals for output
    std::vector< tk::real > ints;
    for (const auto& [s,m] : d[MASS_FLOW_RATE]) ints.push_back( m );
    for (const auto& [s,m] : d[FORCE_X]) ints.push_back( m );
    for (const auto& [s,m] : d[FORCE_Y]) ints.push_back( m );
    for (const auto& [s,m] : d[FORCE_Z]) ints.push_back( m );

    // Append integrals file at selected times
    auto filename = g_cfg.get< tag::output >() + ".int";
    tk::DiagWriter dw( filename,
                       g_cfg.get< tag::integout_format >(),
                       g_cfg.get< tag::integout_precision >(),
                       std::ios_base::app );
    // cppcheck-suppress containerOutOfBounds
    dw.write( static_cast<uint64_t>(tk::cref_find( d[ITER], 0 )),
              // cppcheck-suppress containerOutOfBounds
              tk::cref_find( d[TIME], 0 ),
              // cppcheck-suppress containerOutOfBounds
              tk::cref_find( d[DT], 0 ),
              ints );
  }

  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    // cppcheck-suppress uninitvar
    m_riecg[ meshid ].step();
  }
  else if (solver == "laxcg") {
    // cppcheck-suppress uninitvar
    m_laxcg[ meshid ].step();
  }
  else if (solver == "zalcg") {
    // cppcheck-suppress uninitvar
    m_zalcg[ meshid ].step();
  }
  else if (solver == "kozcg") {
    // cppcheck-suppress uninitvar
    m_kozcg[ meshid ].step();
  }
  else if (solver == "chocg") {
    // cppcheck-suppress uninitvar
    m_chocg[ meshid ].step();
  }
  else if (solver == "lohcg") {
    // cppcheck-suppress uninitvar
    m_lohcg[ meshid ].step();
  }
  else
    Throw( "Unknown solver: " + solver );
}

void
Transporter::resume()
// *****************************************************************************
// Resume execution from checkpoint/restart files
//! \details This is invoked by Charm++ after the checkpoint is done, as well as
//!   when the restart (returning from a checkpoint) is complete
// *****************************************************************************
{
  if (std::any_of(begin(m_finished), end(m_finished), [](auto f){return !f;})) {

    // If just restarted from a checkpoint, Main( CkMigrateMessage* msg ) has
    // increased g_nrestart, but only on PE 0, so broadcast.

    const auto& solver = g_cfg.get< tag::solver >();
    if (solver == "riecg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_riecg[i].evalLB( g_nrestart );
      }
    }
    else if (solver == "laxcg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_laxcg[i].evalLB( g_nrestart );
      }
    }
    else if (solver == "zalcg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_zalcg[i].evalLB( g_nrestart );
      }
    }
    else if ( solver == "kozcg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_kozcg[i].evalLB( g_nrestart );
      }
    }
    else if ( solver == "chocg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_chocg[i].evalLB( g_nrestart );
      }
    }
    else if ( solver == "lohcg") {
      for (std::size_t i=0; i<m_nelem.size(); ++i) {
        m_lohcg[i].evalLB( g_nrestart );
      }
    }
    else {
      Throw( "Unknown solver: " + solver );
    }


  } else {

    mainProxy.finalize();

  }
}

void
Transporter::checkpoint( std::size_t finished, std::size_t meshid )
// *****************************************************************************
// Save checkpoint/restart files
//! \param[in] finished Nonzero if finished with time stepping
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_finished[meshid] = finished;

  if (++m_nchk == m_nelem.size()) { // all worker arrays have checkpointed
    m_nchk = 0;
    if (not g_cfg.get< tag::benchmark >()) {
      const auto& ckptdir = g_cfg.get< tag::checkpoint >();
      CkCallback res( CkIndex_Transporter::resume(), thisProxy );
      CkStartCheckpoint( ckptdir.c_str(), res );
      //CkStartMemCheckpoint( res );
    } else {
      resume();
    }
  }
}

void
Transporter::finish( std::size_t meshid )
// *****************************************************************************
// Normal finish of time stepping
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  checkpoint( /* finished = */ 1, meshid );
}

#include "NoWarning/transporter.def.h"
