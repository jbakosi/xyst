// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \details   Data and functionality common to all discretization schemes
  \see       Discretization.h and Discretization.C for more info.
*/
// *****************************************************************************

#include <iomanip>

#include "Reorder.hpp"
#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Discretization.hpp"
#include "MeshWriter.hpp"
#include "DiagWriter.hpp"
#include "InciterConfig.hpp"
#include "Print.hpp"
#include "Around.hpp"
#include "PDFReducer.hpp"
#include "XystBuildConfig.hpp"
#include "Box.hpp"
#include "Transfer.hpp"

namespace inciter {

static CkReduction::reducerType PDFMerger;
extern ctr::Config g_cfg;

} // inciter::

using inciter::Discretization;

Discretization::Discretization(
  std::size_t meshid,
  const std::vector< CProxy_Discretization >& disc,
  const CProxy_Transporter& transporter,
  const tk::CProxy_MeshWriter& meshwriter,
  const tk::UnsMesh::CoordMap& coordmap,
  const tk::UnsMesh::Chunk& el,
  const std::map< int, std::unordered_set< std::size_t > >& nodeCommMap,
  int nc ) :
  m_meshid( meshid ),
  m_nchare( nc ),
  m_it( 0 ),
  m_itr( 0 ),
  m_itf( 0 ),
  m_initial( 1 ),
  m_t( g_cfg.get< tag::t0 >() ),
  m_lastDumpTime( -std::numeric_limits< tk::real >::max() ),
  m_physFieldFloor( 0.0 ),
  m_physHistFloor( 0.0 ),
  m_physIntegFloor( 0.0 ),
  m_rangeFieldFloor( g_cfg.get< tag::fieldout, tag::range >().size(), 0.0 ),
  m_rangeHistFloor( g_cfg.get< tag::histout, tag::range >().size(), 0.0 ),
  m_rangeIntegFloor( g_cfg.get< tag::integout, tag::range >().size(), 0.0 ),
  m_dt( g_cfg.get< tag::dt >() ),
  m_dtn( m_dt ),
  m_nvol( 0 ),
  m_disc( disc ),
  m_transporter( transporter ),
  m_meshwriter( meshwriter ),
  m_el( el ),     // fills m_inpoel, m_gid, m_lid
  m_coord( setCoord( coordmap ) ),
  m_meshvol( 0.0 ),
  m_v( m_gid.size(), 0.0 ),
  m_vol( m_gid.size(), 0.0 ),
  m_refined( 0 ),
  m_prevstatus( std::chrono::high_resolution_clock::now() ),
  m_nrestart( 0 ),
  m_res( 0.0 ),
  m_res0( 0.0 ),
  m_res1( 0.0 ),
  m_dea( 0 ),
  m_deastarted( 0 ),
  m_pit( 0 ),
  m_mit( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] meshid Mesh ID
//! \param[in] disc Discretization proxy for all meshes
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] meshwriter Mesh writer proxy
//! \param[in] coordmap Coordinates of mesh nodes and their global IDs
//! \param[in] el Elements of the mesh chunk we operate on
//! \param[in] nodeCommMap Node lists associated to chare IDs bordering the
//!   mesh chunk we operate on
//! \param[in] nc Total number of Discretization chares
// *****************************************************************************
{
  Assert( !m_inpoel.empty(), "No elements assigned to Discretization chare" );
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Jacobian in input mesh to Discretization non-positive" );
  #if not defined(__INTEL_COMPILER) || defined(NDEBUG)
  // The above ifdef skips running the conformity test with the intel compiler
  // in debug mode only. This is necessary because in tk::conforming(), filling
  // up the map can fail with some meshes (only in parallel), e.g., tube.exo,
  // used by some regression tests, due to the intel compiler generating some
  // garbage incorrect code - only in debug, only in parallel, only with that
  // mesh.
  Assert( tk::conforming( m_inpoel, m_coord ),
          "Input mesh to Discretization not conforming" );
  #endif

  // Store node communication map
  for (const auto& [c,map] : nodeCommMap) m_nodeCommMap[c] = map;

  // Get ready for computing/communicating nodal volumes
  startvol();

  // Find host elements of user-specified points where time histories are
  // saved, and save the shape functions evaluated at the point locations
  const auto& pt = g_cfg.get< tag::histout, tag::points >();
  for (std::size_t p=0; p<pt.size(); ++p) {
    std::array< tk::real, 4 > N;
    const auto& l = pt[p];
    for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
      if (tk::intet( m_coord, m_inpoel, l, e, N )) {
        m_histdata.push_back( HistData{{ "p"+std::to_string(p+1), e, N }} );
        break;
      }
    }
  }

  // Register with mesh-transfer (if coupled)
  if (m_disc.size() == 1) {
    transfer_initialized();
  } else {
    if (thisIndex == 0) {
      transfer::addMesh( thisProxy, m_nchare,
        CkCallback(CkIndex_Discretization::transfer_initialized(), thisProxy) );
    }
  }
}

void
Discretization::transfer_initialized()
// *****************************************************************************
// Our mesh has been registered with the mesh-to-mesh transfer (if coupled)
// *****************************************************************************
{
  // Compute number of mesh points owned
  std::size_t npoin = m_gid.size();
  for (auto g : m_gid) if (tk::slave( m_nodeCommMap, g, thisIndex ) ) --npoin;

  // Tell the RTS that the Discretization chares have been created and compute
  // the total number of mesh points across the distributed mesh
  std::vector< std::size_t > meshdata{ m_meshid, npoin };
  contribute( meshdata, CkReduction::sum_ulong,
    CkCallback( CkReductionTarget(Transporter,disccreated), m_transporter ) );
}

void
Discretization::transfer( tk::Fields& u, CkCallback c )
// *****************************************************************************
// Initiate solution transfer (if coupled) in 'to' direction
// *****************************************************************************
{
  if (m_disc.size() == 1) {     // not coupled

    c.send();

  }
  else {

    m_transfer_complete = c;
    m_transfer_sol = static_cast< tk::Fields* >( &u );

    // Initiate transfer in 'to' direction
    if (m_meshid == 0) {
      transfer::setSourceTets( thisProxy, thisIndex, m_inpoel, m_coord, u,
        CkCallback( CkIndex_Discretization::transfer_from(),
                    thisProxy[thisIndex] ) );
    }
    else {
      transfer::setDestPoints( thisProxy, thisIndex, m_coord, u,
        CkCallback( CkIndex_Discretization::transfer_from(),
                    thisProxy[thisIndex] ) );
    }

  }
}

void
Discretization::transfer_from()
// *****************************************************************************
// Initiate solution transfer from overset to background mesh
// *****************************************************************************
{
  if (g_cfg.get< tag::overset, tag::oneway >()) {

    m_transfer_complete.send();

  }
  else {

    if (m_meshid == 0) {
      transfer::setDestPoints( thisProxy, thisIndex, m_coord, *m_transfer_sol,
                               m_transfer_complete );
    }
    else {
      transfer::setSourceTets( thisProxy, thisIndex, m_inpoel, m_coord,
                               *m_transfer_sol, m_transfer_complete );
    }

  }
}

void
Discretization::intergrid(
  const std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
//  Prepare integrid-boundary data structures (if coupled)
//! \param[in] bnode Boundary-node lists mapped to side sets used in input file
// *****************************************************************************
{
  bool multi = g_cfg.get< tag::input >().size() > 1;
  if (!multi) return;
  m_transfer_flag.resize( m_coord[0].size(), -1 );
  if (m_meshid == 0) return;

  // Access intergrid-boundary side set ids for this mesh
  const auto& setids = g_cfg.get< tag::overset, tag::intergrid_ >()[ m_meshid ];

  // Compile unique set of intergrid-boundary side set ids for this mesh
  std::unordered_set< std::size_t > ibs( begin(setids), end(setids) );
  if (ibs.empty()) return;

  // Flag points on intergrid boundary
  std::unordered_set< std::size_t > bp; // points flagged
  for (const auto& [setid,n] : bnode) {
    if (ibs.count( static_cast<std::size_t>(setid) )) {
      for (auto g : n) {
        auto i = tk::cref_find(m_lid,g);
        m_transfer_flag[i] = 1;
        bp.insert(i);
      }
    }
  }

  // Add a some layers to intergrid boundary
  auto psup = tk::genPsup( m_inpoel, 4, tk::genEsup( m_inpoel, 4 ) );
  auto layers = g_cfg.get< tag::overset, tag::layers_ >()[ m_meshid ];
  for (int n=0; n<layers; ++n) {
    std::unordered_set< std::size_t > add;
    for (auto p : bp) {
      for (auto q : tk::Around(psup,p)) {
        m_transfer_flag[q] = 1;
        add.insert(q);
      }
    }
    bp.merge( add );
    add.clear();
  }

  // Mark next outer layers for transfer in opposite direction
  for (int n=0; n<layers; ++n) {
    std::unordered_set< std::size_t > add;
    for (auto p : bp) {
      for (auto q : tk::Around(psup,p)) {
        if (m_transfer_flag[q] == -1) m_transfer_flag[q] = 0;
        add.insert(q);
      }
    }
    bp.merge( add );
    add.clear();
  }
}

void
Discretization::resizePostAMR(
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< int, std::unordered_set< std::size_t > >&
    nodeCommMap,
  const std::set< std::size_t >& /*removedNodes*/ )
// *****************************************************************************
//  Resize mesh data structures after mesh refinement
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] nodeCommMap New node communication map
//! \param[in] removedNodes Newly removed mesh node local ids
// *****************************************************************************
{
  m_el = chunk;                 // updates m_inpoel, m_gid, m_lid
  m_nodeCommMap = nodeCommMap;

  // Update mesh volume container size
  m_vol.resize( m_gid.size(), 0.0 );

  // update mesh node coordinates
  m_coord = coord;
}

void
Discretization::startvol()
// *****************************************************************************
//  Get ready for (re-)computing/communicating nodal volumes
// *****************************************************************************
{
  m_nvol = 0;
  thisProxy[ thisIndex ].wait4vol();

  // Zero out mesh volume container
  std::fill( begin(m_vol), end(m_vol), 0.0 );

  // Clear receive buffer that will be used for collecting nodal volumes
  m_volc.clear();
}

void
Discretization::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//!  \details Since this is a [initnode] routine, see the .ci file, the
//!   Charm++ runtime system executes the routine exactly once on every
//!   logical node early on in the Charm++ init sequence. Must be static as
//!   it is called without an object. See also: Section "Initializations at
//!   Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  PDFMerger = CkReduction::addReducer( tk::mergeUniPDFs );
}

tk::UnsMesh::Coords
Discretization::setCoord( const tk::UnsMesh::CoordMap& coordmap )
// *****************************************************************************
// Set mesh coordinates based on coordinates map
// *****************************************************************************
{
  Assert( coordmap.size() == m_gid.size(), "Size mismatch" );
  Assert( coordmap.size() == m_lid.size(), "Size mismatch" );

  tk::UnsMesh::Coords coord;
  coord[0].resize( coordmap.size() );
  coord[1].resize( coordmap.size() );
  coord[2].resize( coordmap.size() );

  for (const auto& [ gid, coords ] : coordmap) {
    auto i = tk::cref_find( m_lid, gid );
    coord[0][i] = coords[0];
    coord[1][i] = coords[1];
    coord[2][i] = coords[2];
  }

  return coord;
}

void
Discretization::remap(
  const std::unordered_map< std::size_t, std::size_t >& map )
// *****************************************************************************
//  Remap mesh data based on new local ids
//! \param[in] map Mapping of old->new local ids
// *****************************************************************************
{
  // Remap connectivity containing local IDs
  for (auto& l : m_inpoel) l = tk::cref_find(map,l);

  // Remap global->local id map
  for (auto& [g,l] : m_lid) l = tk::cref_find(map,l);

  // Remap global->local id map
  auto maxid = std::numeric_limits< std::size_t >::max();
  std::vector< std::size_t > newgid( m_gid.size(), maxid );
  for (const auto& [o,n] : map) newgid[n] = m_gid[o];
  m_gid = std::move( newgid );

  Assert( std::all_of( m_gid.cbegin(), m_gid.cend(),
            [=](std::size_t i){ return i < maxid; } ),
          "Not all gid have been remapped" );

  // Remap nodal volumes (with contributions along chare-boundaries)
  std::vector< tk::real > newvol( m_vol.size(), 0.0 );
  for (const auto& [o,n] : map) newvol[n] = m_vol[o];
  m_vol = std::move( newvol );

  // Remap nodal volumes (without contributions along chare-boundaries)
  std::vector< tk::real > newv( m_v.size(), 0.0 );
  for (const auto& [o,n] : map) newv[n] = m_v[o];
  m_v = std::move( newv );

  // Remap locations of node coordinates
  tk::UnsMesh::Coords newcoord;
  auto npoin = m_coord[0].size();
  newcoord[0].resize( npoin );
  newcoord[1].resize( npoin );
  newcoord[2].resize( npoin );
  for (const auto& [o,n] : map) {
    newcoord[0][n] = m_coord[0][o];
    newcoord[1][n] = m_coord[1][o];
    newcoord[2][n] = m_coord[2][o];
  }
  m_coord = std::move( newcoord );
}

void
Discretization::setRefiner( const CProxy_Refiner& ref )
// *****************************************************************************
//  Set Refiner Charm++ proxy
//! \param[in] ref Incoming refiner proxy to store
// *****************************************************************************
{
  m_refiner = ref;
}

void
Discretization::vol()
// *****************************************************************************
// Sum mesh volumes to nodes, start communicating them on chare-boundaries
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Compute nodal volumes on our chunk of the mesh
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const auto N = m_inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) / 24.0;
    ErrChk( J > 0, "Element Jacobian non-positive: PE:" +
                   std::to_string(CkMyPe()) + ", node IDs: " +
                   std::to_string(m_gid[N[0]]) + ',' +
                   std::to_string(m_gid[N[1]]) + ',' +
                   std::to_string(m_gid[N[2]]) + ',' +
                   std::to_string(m_gid[N[3]]) + ", coords: (" +
                   std::to_string(x[N[0]]) + ", " +
                   std::to_string(y[N[0]]) + ", " +
                   std::to_string(z[N[0]]) + "), (" +
                   std::to_string(x[N[1]]) + ", " +
                   std::to_string(y[N[1]]) + ", " +
                   std::to_string(z[N[1]]) + "), (" +
                   std::to_string(x[N[2]]) + ", " +
                   std::to_string(y[N[2]]) + ", " +
                   std::to_string(z[N[2]]) + "), (" +
                   std::to_string(x[N[3]]) + ", " +
                   std::to_string(y[N[3]]) + ", " +
                   std::to_string(z[N[3]]) + ')' );
    // scatter add V/4 to nodes
    for (std::size_t j=0; j<4; ++j) m_vol[N[j]] += J;
  }

  // Store nodal volumes without contributions from other chares on
  // chare-boundaries
  m_v = m_vol;

  // Send our nodal volume contributions to neighbor chares
  if (m_nodeCommMap.empty()) {
    comvol_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< tk::real > v( n.size() );
      std::size_t j = 0;
      for (auto i : n) v[ j++ ] = m_vol[ tk::cref_find(m_lid,i) ];
      thisProxy[c].comvol( thisIndex,
                           std::vector<std::size_t>(begin(n), end(n)), v );
    }
  }

  ownvol_complete();
}

void
Discretization::comvol( int c,
                        const std::vector< std::size_t >& gid,
                        const std::vector< tk::real >& nodevol )
// *****************************************************************************
//  Receive nodal volumes on chare-boundaries
//! \param[in] c Sender chare id
//! \param[in] gid Global mesh node IDs at which we receive volume contributions
//! \param[in] nodevol Partial sums of nodal volume contributions to
//!    chare-boundary nodes
// *****************************************************************************
{
  Assert( nodevol.size() == gid.size(), "Size mismatch" );

  auto& cvolc = m_cvolc[c];
  for (std::size_t i=0; i<gid.size(); ++i) {
    m_volc[ gid[i] ] += nodevol[i];
    cvolc[ tk::cref_find(m_lid,gid[i]) ] = nodevol[i];
  }

  if (++m_nvol == m_nodeCommMap.size()) {
    m_nvol = 0;
    comvol_complete();
  }
}

void
Discretization::totalvol()
// *****************************************************************************
// Sum mesh volumes and contribute own mesh volume to total volume
// *****************************************************************************
{
  // Add received contributions to nodal volumes
  for (const auto& [gid, vol] : m_volc)
    m_vol[ tk::cref_find(m_lid,gid) ] += vol;

  // Clear receive buffer
  tk::destroy(m_volc);

  // Sum mesh volume to host
  std::vector< tk::real > tvol{ 0.0,
                                static_cast<tk::real>(m_initial),
                                static_cast<tk::real>(m_meshid) };
  for (auto v : m_v) tvol[0] += v;
  contribute( tvol, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,totalvol), m_transporter) );
}

void
Discretization::stat( tk::real mesh_volume )
// *****************************************************************************
// Compute mesh cell statistics
//! \param[in] mesh_volume Total mesh volume
// *****************************************************************************
{
  // Store total mesh volume
  m_meshvol = mesh_volume;

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  auto MIN = -std::numeric_limits< tk::real >::max();
  auto MAX = std::numeric_limits< tk::real >::max();
  std::vector< tk::real > min( 6, MAX );
  std::vector< tk::real > max( 6, MIN );
  std::vector< tk::real > sum( 9, 0.0 );
  tk::UniPDF edgePDF( 1e-4 );
  tk::UniPDF volPDF( 1e-4 );
  tk::UniPDF ntetPDF( 1e-4 );

  // Compute points surrounding points
  auto psup = tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) );
  Assert( psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Compute edge length statistics
  // Note that while the min and max edge lengths are independent of the number
  // of CPUs (by the time they are aggregated across all chares), the sum of
  // the edge lengths and the edge length PDF are not. This is because the
  // edges on the chare-boundary are counted multiple times and we
  // conscientiously do not make an effort to precisely compute this, because
  // that would require communication and more complex logic. Since these
  // statistics are intended as simple average diagnostics, we ignore these
  // small differences. For reproducible average edge lengths and edge length
  // PDFs, run the mesh in serial.
  tk::UnsMesh::EdgeSet edges;
  for (std::size_t p=0; p<m_gid.size(); ++p) {
    for (auto i : tk::Around(psup,p)) {
       const auto dx = x[ i ] - x[ p ];
       const auto dy = y[ i ] - y[ p ];
       const auto dz = z[ i ] - z[ p ];
       const auto length = std::sqrt( dx*dx + dy*dy + dz*dz );
       if (length < min[0]) min[0] = length;
       if (length > max[0]) max[0] = length;
       sum[0] += 1.0;
       sum[1] += length;
       edgePDF.add( length );
       edges.insert( { m_gid[i], m_gid[p] } );
    }
  }

  // Compute mesh cell volume statistics
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
    if (L < min[1]) min[1] = L;
    if (L > max[1]) max[1] = L;
    sum[2] += 1.0;
    sum[3] += L;
    volPDF.add( L );
  }

  // Contribute statistics
  sum[4] = 1.0;
  min[2] = max[2] = sum[5] = static_cast< tk::real >( m_inpoel.size() / 4 );
  min[3] = max[3] = sum[6] = static_cast< tk::real >( m_gid.size() );
  min[4] = max[4] = sum[7] = static_cast< tk::real >( edges.size() );
  min[5] = max[5] = sum[8] =
    static_cast< tk::real >( tk::sumvalsize(m_nodeCommMap) ) /
    static_cast< tk::real >( m_gid.size() );
  ntetPDF.add( min[2] );

  min.push_back( static_cast<tk::real>(m_meshid) );
  max.push_back( static_cast<tk::real>(m_meshid) );
  sum.push_back( static_cast<tk::real>(m_meshid) );

  // Contribute to mesh statistics across all Discretization chares
  contribute( min, CkReduction::min_double,
    CkCallback(CkReductionTarget(Transporter,minstat), m_transporter) );
  contribute( max, CkReduction::max_double,
    CkCallback(CkReductionTarget(Transporter,maxstat), m_transporter) );
  contribute( sum, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,sumstat), m_transporter) );

  // Serialize PDFs to raw stream
  auto stream = tk::serialize( m_meshid, { edgePDF, volPDF, ntetPDF } );
  // Create Charm++ callback function for reduction of PDFs with
  // Transporter::pdfstat() as the final target where the results will appear.
  CkCallback cb( CkIndex_Transporter::pdfstat(nullptr), m_transporter );
  // Contribute serialized PDF of partial sums to host via Charm++ reduction
  contribute( stream.first, stream.second.get(), PDFMerger, cb );
}

void
Discretization::boxvol()
// *****************************************************************************
// Compute total box IC volume
// *****************************************************************************
{
  // Determine which nodes reside in user-defined IC box(es) if any
  m_boxnodes = problems::boxnodes( m_coord );

  // Compute partial box IC volume (just add up all boxes)
  tk::real boxvol = 0.0;
  // cppcheck-suppress useStlAlgorithm
  for (const auto& b : m_boxnodes) for (auto i : b) boxvol += m_v[i];

  // Sum up box IC volume across all chares
  std::vector< tk::real > meshdata{ boxvol, static_cast<tk::real>(m_meshid) };
  contribute( meshdata, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,boxvol), m_transporter) );
}

void
Discretization::write(
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel,
  const std::vector< std::string>& elemfieldnames,
  const std::vector< std::string>& nodefieldnames,
  const std::vector< std::string>& elemsurfnames,
  const std::vector< std::string>& nodesurfnames,
  const std::vector< std::vector< tk::real > >& elemfields,
  const std::vector< std::vector< tk::real > >& nodefields,
  const std::vector< std::vector< tk::real > >& elemsurfs,
  const std::vector< std::vector< tk::real > >& nodesurfs,
  CkCallback c )
// *****************************************************************************
//  Output mesh and fields data (solution dump) to file(s)
//! \param[in] inpoel Mesh connectivity for the mesh chunk to be written
//! \param[in] coord Node coordinates of the mesh chunk to be written
//! \param[in] bface Map of boundary-face lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \param[in] elemfieldnames Names of element fields to be output to file
//! \param[in] nodefieldnames Names of node fields to be output to file
//! \param[in] elemsurfnames Names of elem surface fields to be output to file
//! \param[in] nodesurfnames Names of node surface fields to be output to file
//! \param[in] elemfields Field data in mesh elements to output to file
//! \param[in] nodefields Field data in mesh nodes to output to file
//! \param[in] elemsurfs Surface field data in mesh elements to output to file
//! \param[in] nodesurfs Surface field data in mesh nodes to output to file
//! \param[in] c Function to continue with after the write
//! \details Since m_meshwriter is a Charm++ chare group, it never migrates and
//!   an instance is guaranteed on every PE. We index the first PE on every
//!   logical compute node. In Charm++'s non-SMP mode, a node is the same as a
//!   PE, so the index is the same as CkMyPe(). In SMP mode the index is the
//!   first PE on every logical node. In non-SMP mode this yields one or more
//!   output files per PE with zero or non-zero virtualization, respectively. If
//!   there are multiple chares on a PE, the writes are serialized per PE, since
//!   only a single entry method call can be executed at any given time. In SMP
//!   mode, still the same number of files are output (one per chare), but the
//!   output is serialized through the first PE of each compute node. In SMP
//!   mode, channeling multiple files via a single PE on each node is required
//!   by NetCDF and HDF5, as well as ExodusII, since none of these libraries are
//!   thread-safe.
// *****************************************************************************
{
  // If the previous iteration refined (or moved) the mesh or this is called
  // before the first time step, we also output the mesh.
  bool meshoutput = m_itf == 0 ? true : false;

  auto eps = std::numeric_limits< tk::real >::epsilon();
  bool fieldoutput = false;

  // Output field data only if there is no dump at this physical time yet
  if (std::abs(m_lastDumpTime - m_t) > eps ) {
    m_lastDumpTime = m_t;
    ++m_itf;
    fieldoutput = true;
  }

  const auto& f = g_cfg.get< tag::fieldout, tag::sidesets >();
  std::set< int > outsets( begin(f), end(f) );

  m_meshwriter[ CkNodeFirst( CkMyNode() ) ].
    write( m_meshid, meshoutput, fieldoutput, m_itr, m_itf, m_t, thisIndex,
           g_cfg.get< tag::output >(),
           inpoel, coord, bface, bnode, triinpoel, elemfieldnames,
           nodefieldnames, elemsurfnames, nodesurfnames, elemfields, nodefields,
           elemsurfs, nodesurfs, outsets, c );
}

void
Discretization::setdt( tk::real newdt )
// *****************************************************************************
// Set time step size
//! \param[in] newdt Size of the new time step
// *****************************************************************************
{
  m_dtn = m_dt;
  m_dt = newdt;

  // Truncate the size of last time step
  const auto term = g_cfg.get< tag::term >();
  if (m_t+m_dt > term) m_dt = term - m_t;
}

void
Discretization::next()
// *****************************************************************************
// Prepare for next step
// *****************************************************************************
{
  // Update floor of physics time divided by output interval times
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto ft = g_cfg.get< tag::fieldout, tag::time >();
  if (ft > eps) m_physFieldFloor = std::floor( m_t / ft );
  const auto ht = g_cfg.get< tag::histout, tag::time >();
  if (ht > eps) m_physHistFloor = std::floor( m_t / ht );
  const auto it = g_cfg.get< tag::integout, tag::time >();
  if (it > eps) m_physIntegFloor = std::floor( m_t / it );

  // Update floors of physics time divided by output interval times for ranges
  const auto& rf = g_cfg.get< tag::fieldout, tag::range >();
  for (std::size_t i=0; i<rf.size(); ++i) {
    if (m_t > rf[i][0] and m_t < rf[i][1]) {
      m_rangeFieldFloor[i] = std::floor( m_t / rf[i][2] );
    }
  }

  const auto& rh = g_cfg.get< tag::histout, tag::range >();
  for (std::size_t i=0; i<rh.size(); ++i) {
    if (m_t > rh[i][0] and m_t < rh[i][1]) {
      m_rangeHistFloor[i] = std::floor( m_t / rh[i][2] );
    }
  }

  const auto& ri = g_cfg.get< tag::integout, tag::range >();
  for (std::size_t i=0; i<ri.size(); ++i) {
    if (m_t > ri[i][0] and m_t < ri[i][1]) {
      m_rangeIntegFloor[i] = std::floor( m_t / ri[i][2] );
    }
  }

  ++m_it;
  m_t += m_dt;
}

void
Discretization::grindZero()
// *****************************************************************************
//  Zero grind-time
// *****************************************************************************
{
  m_prevstatus = std::chrono::high_resolution_clock::now();

  if (thisIndex == 0 && m_meshid == 0) {
    tk::Print() << "Starting time stepping ...\n";
  }
}

bool
Discretization::restarted( int nrestart )
// *****************************************************************************
//  Detect if just returned from a checkpoint and if so, zero timers
//! \param[in] nrestart Number of times restarted
//! \return True if restart detected
// *****************************************************************************
{
  // Detect if just restarted from checkpoint:
  //   nrestart == -1 if there was no checkpoint this step
  //   m_nrestart == nrestart if there was a checkpoint this step
  //   if both false, just restarted from a checkpoint
  bool restarted = nrestart != -1 and m_nrestart != nrestart;

   // If just restarted from checkpoint
  if (restarted) {
    // Update number of restarts
    m_nrestart = nrestart;
    // Start timer measuring time stepping wall clock time
    m_timer.zero();
    // Zero grind-timer
    grindZero();
  }

  return restarted;
}

std::string
Discretization::histfilename( const std::string& id, std::streamsize precision )
// *****************************************************************************
//  Construct history output filename
//! \param[in] id History point id
//! \param[in] precision Floating point precision to use for output
//! \return History file name
// *****************************************************************************
{
  auto of = g_cfg.get< tag::output >();
  std::stringstream ss;

  ss << std::setprecision(static_cast<int>(precision)) << of << ".hist." << id;

  return ss.str();
}

void
Discretization::histheader( std::vector< std::string >&& names )
// *****************************************************************************
//  Output headers for time history files (one for each point)
//! \param[in] names History output variable names
// *****************************************************************************
{
  for (const auto& h : m_histdata) {
    auto prec = g_cfg.get< tag::histout, tag:: precision >();
    tk::DiagWriter hw( histfilename( h.get< tag::id >(), prec ),
                       g_cfg.get< tag::histout, tag::format >(),
                       prec );
    hw.header( names );
  }
}

void
Discretization::history( std::vector< std::vector< tk::real > >&& data )
// *****************************************************************************
//  Output time history for a time step
//! \param[in] data Time history data for all variables and equations integrated
// *****************************************************************************
{
  Assert( data.size() == m_histdata.size(), "Size mismatch" );

  std::size_t i = 0;
  for (const auto& h : m_histdata) {
    auto prec = g_cfg.get< tag::histout, tag::precision >();
    tk::DiagWriter hw( histfilename( h.get< tag::id >(), prec ),
                       g_cfg.get< tag::histout, tag::format >(),
                       prec,
                       std::ios_base::app );
    hw.write( m_it, m_t, m_dt, data[i] );
    ++i;
  }
}

bool
Discretization::fielditer() const
// *****************************************************************************
//  Decide if field output iteration count interval is hit
//! \return True if field output iteration count interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  return m_it % g_cfg.get< tag::fieldout, tag::iter >() == 0;
}

bool
Discretization::fieldtime() const
// *****************************************************************************
//  Decide if field output physics time interval is hit
//! \return True if field output physics time interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto ft = g_cfg.get< tag::fieldout, tag::time >();

  if (ft < eps) return false;

  return std::floor(m_t/ft) - m_physFieldFloor > eps;
}

bool
Discretization::fieldrange() const
// *****************************************************************************
//  Decide if physics time falls into a field output time range
//! \return True if physics time falls into a field output time range
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  const auto eps = std::numeric_limits< tk::real >::epsilon();

  bool output = false;

  const auto& rf = g_cfg.get< tag::fieldout, tag::range >();
  for (std::size_t i=0; i<rf.size(); ++i) {
    if (m_t > rf[i][0] and m_t < rf[i][1]) {
      output |= std::floor(m_t/rf[i][2]) - m_rangeFieldFloor[i] > eps;
    }
  }

  return output;
}

bool
Discretization::histiter() const
// *****************************************************************************
//  Decide if history output iteration count interval is hit
//! \return True if history output iteration count interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto hist = g_cfg.get< tag::histout, tag::iter >();
  const auto& hist_points = g_cfg.get< tag::histout, tag::points >();

  return m_it % hist == 0 and not hist_points.empty();
}

bool
Discretization::histtime() const
// *****************************************************************************
//  Decide if history output physics time interval is hit
//! \return True if history output physics time interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto ht = g_cfg.get< tag::histout, tag::time >();

  if (ht < eps) return false;

  return std::floor(m_t/ht) - m_physHistFloor > eps;
}

bool
Discretization::histrange() const
// *****************************************************************************
//  Decide if physics time falls into a history output time range
//! \return True if physics time falls into a history output time range
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto eps = std::numeric_limits< tk::real >::epsilon();

  bool output = false;

  const auto& rh = g_cfg.get< tag::histout, tag::range >();
  for (std::size_t i=0; i<rh.size(); ++i) {
    if (m_t > rh[i][0] and m_t < rh[i][1]) {
      output |= std::floor(m_t/rh[i][2]) - m_rangeHistFloor[i] > eps;
    }
  }

  return output;
}

bool
Discretization::integiter() const
// *****************************************************************************
//  Decide if integral output iteration count interval is hit
//! \return True if integral output iteration count interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto integ = g_cfg.get< tag::integout, tag::iter >();
  const auto& sidesets_integral = g_cfg.get< tag::integout, tag::sidesets >();

  return m_it % integ == 0 and not sidesets_integral.empty();
}

bool
Discretization::integtime() const
// *****************************************************************************
//  Decide if integral output physics time interval is hit
//! \return True if integral output physics time interval is hit
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto it = g_cfg.get< tag::integout, tag::time >();

  if (it < eps) return false;

  return std::floor(m_t/it) - m_physIntegFloor > eps;
}

bool
Discretization::integrange() const
// *****************************************************************************
//  Decide if physics time falls into a integral output time range
//! \return True if physics time falls into a integral output time range
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) return false;

  auto eps = std::numeric_limits< tk::real >::epsilon();

  bool output = false;

  const auto& ri = g_cfg.get< tag::integout, tag::range >();
  for (std::size_t i=0; i<ri.size(); ++i) {
    if (m_t > ri[i][0] and m_t < ri[i][1]) {
      output |= std::floor(m_t/ri[i][2]) - m_rangeIntegFloor[i] > eps;
    }
  }

  return output;
}

bool
Discretization::finished() const
// *****************************************************************************
//  Decide if this is the last time step
//! \return True if this is the last time step
// *****************************************************************************
{
  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto nstep = g_cfg.get< tag::nstep >();
  auto term = g_cfg.get< tag::term >();
  auto residual = g_cfg.get< tag::residual >();

  return std::abs(m_t-term) < eps or m_it >= nstep or
         (m_res > 0.0 and m_res < residual);
}

void
Discretization::residual( tk::real r )
// *****************************************************************************
//  Update residual (during convergence to steady state)
//! \param[in] r Current residual
// *****************************************************************************
{
  auto ttyi = g_cfg.get< tag::ttyi >();

  if (m_it % ttyi == 0) {
    m_res0 = m_res1;
    m_res1 = r;
  }

  m_res = r;
}

bool
Discretization::deastart()
// *****************************************************************************
//  Decide whether to start the deactivation procedure in this time step
//! \return True to start the deactivation prcedure in this time step
// *****************************************************************************
{
  if (not m_deastarted and m_t > g_cfg.get<tag::deatime>() + m_dt) {
    m_deastarted = 1;
    return true;
  }

  return false;
}

bool
Discretization::deactivate() const
// *****************************************************************************
//  Decide whether to run deactivation procedure in this time step
//! \return True to run deactivation prcedure in this time step
// *****************************************************************************
{
  auto dea = g_cfg.get< tag::deactivate >();
  auto deafreq = g_cfg.get< tag::deafreq >();

  if (dea and !m_nodeCommMap.empty() and m_it % deafreq == 0)
    return true;
  else
    return false;
}

void
Discretization::deactivated( int n )
// *****************************************************************************
//  Receive deactivation report
//! \param[in] n Sum of deactivated chares
// *****************************************************************************
{
  m_dea = n;
}

bool
Discretization::lb() const
// *****************************************************************************
//  Decide whether to do load balancing this time step
//! \return True to do load balancing in this time step
// *****************************************************************************
{
  auto lbfreq = g_cfg.get< tag::lbfreq >();
  auto lbtime = g_cfg.get< tag::lbtime >();

  if ((m_t > lbtime and m_it % lbfreq == 0) or m_it == 2)
    return true;
  else
    return false;
}

void
Discretization::pit( std::size_t it )
// *****************************************************************************
//  Update number of pressure linear solve iterations taken
//! \param[in] it Number of pressure linear solve iterations taken
// *****************************************************************************
{
  m_pit = it;
}

void
Discretization::mit( std::size_t it )
// *****************************************************************************
//  Update number of momentum/transport linear solve iterations taken
//! \param[in] it Number of momentum/transport linear solve iterations taken
// *****************************************************************************
{
  m_mit = it;
}

void
Discretization::npoin( std::size_t n )
// *****************************************************************************
//  Set number of mesh points (across all meshes)
//! \param[in] n Number of mesh points
// *****************************************************************************
{
  m_npoin = n;
}

void
Discretization::status()
// *****************************************************************************
// Output one-liner status report
// *****************************************************************************
{
  const auto ttyi = g_cfg.get< tag::ttyi >();

  if (thisIndex == 0 and m_meshid == 0 and m_it % ttyi == 0) {

    // estimate grind time (taken between this and the previous time step)
    using std::chrono::duration_cast;
    using us = std::chrono::microseconds;
    using clock = std::chrono::high_resolution_clock;
    auto grind_time = duration_cast< us >(clock::now() - m_prevstatus).count()
                      / static_cast< long >( ttyi );
    auto grind_perf = static_cast<tk::real>(grind_time)
                      / static_cast<tk::real>(m_npoin);
    m_prevstatus = clock::now();

    const auto term = g_cfg.get< tag::term >();
    const auto t0 = g_cfg.get< tag::t0 >();
    const auto nstep = g_cfg.get< tag::nstep >();
    const auto diag = g_cfg.get< tag::diag_iter >();
    const auto rsfreq = g_cfg.get< tag::rsfreq >();
    const auto benchmark = g_cfg.get< tag::benchmark >();
    const auto residual = g_cfg.get< tag::residual >();
    const auto solver = g_cfg.get< tag::solver >();
    const auto pre = solver == "chocg" ? 1 : 0;
    const auto theta = g_cfg.get< tag::theta >();
    const auto eps = std::numeric_limits< tk::real >::epsilon();
    const auto mom = solver == "chocg" and theta > eps ? 1 : 0;

    // estimate time elapsed and time for accomplishment
    tk::Timer::Watch ete, eta;
    m_timer.eta( term-t0, m_t-t0, nstep, m_it, m_res0, m_res1, residual,
                 ete, eta );

    tk::Print print;

    // Output one-liner
    print << std::setfill(' ') << std::setw(8) << m_it << "  "
          << std::scientific << std::setprecision(6)
          << std::setw(12) << m_t << "  "
          << m_dt << "  "
          << std::setfill('0')
          << std::setw(3) << ete.hrs.count() << ":"
          << std::setw(2) << ete.min.count() << ":"
          << std::setw(2) << ete.sec.count() << "  "
          << std::setw(3) << eta.hrs.count() << ":"
          << std::setw(2) << eta.min.count() << ":"
          << std::setw(2) << eta.sec.count() << "  "
          << std::scientific << std::setprecision(6) << std::setfill(' ')
          << std::setw(9) << grind_time << "  "
          << std::setw(9) << grind_perf << "  ";

    // Augment one-liner status with output indicators
    if (fielditer() or fieldtime() or fieldrange()) print << 'f';
    if (not (m_it % diag)) print << 'd';
    if (histiter() or histtime() or histrange()) print << 't';
    if (integiter() or integtime() or integrange()) print << 'i';
    if (m_refined) print << 'h';
    if (lb() and not finished()) print << 'l';
    if (not benchmark && (not (m_it % rsfreq) || finished())) print << 'c';
    if (m_deastarted and deactivate()) {
      print << "\te:" << m_dea << '/' << m_nchare;
    }
    if (pre) print << "\tp:" << m_pit;
    if (mom) print << "\tm:" << m_mit;

    print << '\n';
  }
}

#include "NoWarning/discretization.def.h"
