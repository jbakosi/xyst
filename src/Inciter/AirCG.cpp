// *****************************************************************************
/*!
  \file      src/Inciter/AirCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     AirCG for a PDE system with continuous Galerkin + ALE + RK
  \details   AirCG advances a system of partial differential equations (PDEs)
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme..
*/
// *****************************************************************************

#include "XystBuildConfig.hpp"
#include "AirCG.hpp"
#include "Vector.hpp"
#include "Reader.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "ExodusIIMeshWriter.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "DerivedData.hpp"
#include "Discretization.hpp"
#include "DiagReducer.hpp"
#include "Refiner.hpp"
#include "Reorder.hpp"
#include "Around.hpp"
#include "Operators.hpp"
#include "EOS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;

//! Runge-Kutta coefficients
static const std::array< tk::real, 3 > rkcoef{{ 1.0/3.0, 1.0/2.0, 1.0 }};

} // inciter::

using inciter::AirCG;

AirCG::AirCG( const CProxy_Discretization& disc,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ngrad( 0 ),
  m_nrhs( 0 ),
  m_nbnorm( 0 ),
  m_ndfnorm( 0 ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_bndel( Disc()->bndel() ),
  m_dfnorm(),
  m_dfnormc(),
  m_dfn(),
  m_esup( tk::genEsup( Disc()->Inpoel(), 4 ) ),
  m_psup( tk::genPsup( Disc()->Inpoel(), 4, m_esup ) ),
  m_u( Disc()->Gid().size(), 5 ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_rhsc(),
  m_chBndGrad( Disc()->Bid().size(), m_u.nprop()*3 ),
  m_chBndGradc(),
  m_diag(),
  m_bnorm(),
  m_bnormc(),
  m_dirbc(),
  m_symbcnodes(),
  m_farfieldbcnodes(),
  m_symbctri(),
  m_spongenodes(),
  m_timedepbcnodes(),
  m_timedepbcFn(),
  m_stage( 0 ),
  m_boxnodes(),
  m_edgenode(),
  m_edgeid(),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_finished( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side sets used in the input file
//! \param[in] bnode Boundary-node lists mapped to side sets used in input file
//! \param[in] triinpoel Boundary-face connectivity where BCs set (global ids)
// *****************************************************************************
//! [Constructor]
{
  usesAtSync = true;    // enable migration at AtSync

  auto d = Disc();

  // Perform optional operator-access-pattern mesh node reordering
  if (g_inputdeck.get< tag::discr, tag::operator_reorder >()) {

    // Create new local ids based on access pattern of PDE operators
    std::unordered_map< std::size_t, std::size_t > map;
    std::size_t n = 0;

    for (std::size_t p=0; p<m_u.nunk(); ++p) {  // for each point p
      if (map.find(p) == end(map)) map[p] = n++;
      for (auto q : tk::Around(m_psup,p)) {     // for each edge p-q
        if (map.find(q) == end(map)) map[q] = n++;
      }
    }

    Assert( map.size() == d->Gid().size(), "Map size mismatch" );

    // Remap data in bound Discretization object
    d->remap( map );
    // Recompute elements surrounding points
    m_esup = tk::genEsup( d->Inpoel(), 4 );
    // Recompute points surrounding points
    m_psup = tk::genPsup( d->Inpoel(), 4, m_esup );
    // Remap boundary triangle face connectivity
    tk::remap( m_triinpoel, map );
  }

  // Query/update boundary-conditions-related data structures
  queryBnd();

  // Activate SDAG wait for initially computing normals
  thisProxy[ thisIndex ].wait4norm();

  // Signal the runtime system that the workers have been created
  auto meshid = d->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::sum_ulong,
    CkCallback(CkReductionTarget(Transporter,comfinal), d->Tr()) );
}
//! [Constructor]

std::unordered_map< std::size_t, std::vector< std::pair< bool, tk::real > > >
AirCG::match( [[maybe_unused]] tk::ctr::ncomp_t nprop,
              tk::real t,
              tk::real dt,
              const std::vector< tk::real >& tp,
              const std::vector< tk::real >& dtp,
              const tk::UnsMesh::Coords& coord,
              const std::unordered_map< std::size_t, std::size_t >& lid,
              const std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
//  Match user-specified boundary conditions at nodes for side sets
//! \param[in] nprop Number of scalar components in all PDE systems solved
//! \param[in] t Physical time at which to query boundary conditions
//! \param[in] dt Time step size (for querying BC increments in time)
//! \param[in] tp Physical time for each mesh node
//! \param[in] dtp Time step size for each mesh node
//! \param[in] coord Mesh node coordinates
//! \param[in] lid Local node IDs associated to local node IDs
//! \param[in] bnode Map storing global mesh node IDs mapped to side set ids
//! \return Vector of pairs of bool and boundary condition value associated to
//!   local mesh node IDs at which the user has set Dirichlet boundary
//!   conditions for all systems of PDEs integrated. The bool indicates whether
//!   the BC is set at the node for that component: if true, the real value is
//!   the increment (from t to dt) in (or the value of) the BC specified for a
//!   component.
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets (i.e., a
//!   list of elements). This function queries Dirichlet boundary condition
//!   values from all PDEs in the multiple systems of PDEs integrated at the
//!   node lists associated to side set IDs, given by bnode. Each
//!   PDE system returns a BC data structure. Note that the BC mesh nodes that
//!   this function results in (stored in dirbc) only contains those nodes that
//!   are supplied via bnode. i.e., in parallel only a part of the mesh is
//!   worked on.
// *****************************************************************************
{
  // Vector of pairs of bool and boundary condition value associated to mesh
  // node IDs at which the user has set Dirichlet BCs for all PDEs integrated.
  std::unordered_map< std::size_t,
    std::vector< std::pair< bool, tk::real > > > dirbc;

  // Details for the algorithm below: PDE::dirbc() returns a new map that
  // associates a vector of pairs associated to local node IDs. (The pair is a
  // pair of bool and real value, the former is the fact that the BC is to be
  // set while the latter is the value if it is to be set). The length of this
  // NodeBC vector, returning from each system of PDEs equals to the number of
  // scalar components the given PDE integrates. Here we contatenate this map
  // for all PDEs being integrated. If there are multiple BCs set at a mesh node
  // (dirbc::key), either because (1) in the same PDE system the user prescribed
  // BCs on side sets that share nodes or (2) because more than a single PDE
  // system assigns BCs to a given node (on different variables), the NodeBC
  // vector must be correctly stored. "Correctly" here means that the size of
  // the NodeBC vectors must all be the same and equal to the sum of all scalar
  // components integrated by all PDE systems integrated. Example: single-phase
  // compressible flow (density, momentum, energy = 5) + transported scalars of
  // 10 variables -> NodeBC vector length = 15. Note that in case (1) above a
  // new node encountered must "overwrite" the already existing space for the
  // NodeBC vector. "Overwrite" here means that it should keep the existing BCs
  // and add the new ones yielding the union the two prescription for BCs but in
  // the same space that already exist in the NodeBC vector. In case (2),
  // however, the NodeBC pairs must go to the location in the vector assigned to
  // the given PDE system, i.e., using the above example BCs for the 10 (or
  // less) scalars should go in the positions starting at 5, leaving the first 5
  // false, indicating no BCs for the flow variables.
  //
  // When a particular node belongs to two or more side sets with different BCs,
  // there is an ambiguity as to which of the multiple BCs should be applied to
  // the node. This issue is described in case (1) above. In the current
  // implementation, every side set applies the BC to the common node in
  // question, successively overwriting the BC applied by the previous side set.
  // Effectively, the BC corresponding to the last side set ID is applied to the
  // common node. Since bnode is an ordered map, the side set with a larger
  // id wins if a node belongs to multiple side sets.

  // Lambda to convert global to local node ids of a list of nodes
  auto local = [ &lid ]( const std::vector< std::size_t >& gnodes ){
    std::vector< std::size_t > lnodes( gnodes.size() );
    for (std::size_t i=0; i<gnodes.size(); ++i)
      lnodes[i] = tk::cref_find( lid, gnodes[i] );
    return lnodes;
  };

  // Query Dirichlet BCs for all PDEs integrated and assign to nodes
  for (const auto& s : bnode) { // for all side sets passed in
    auto l = local(s.second);   // generate local node ids on side set
    // query Dirichlet BCs at nodes of this side set
    auto eqbc = inciter::dirbc( t, dt, tp, dtp, {s.first,l}, coord );
    std::size_t offset = 0;
    for (const auto& [node,bcval] : eqbc) {
      auto& nodebc = dirbc[node];
      nodebc.clear(); // multiple BCs at a node: last one wins
      for (std::size_t c=0; c<nprop; ++c) nodebc.push_back( {false,0.0} );
      for (std::size_t c=0; c<5; ++c) nodebc[offset+c] = bcval[c];
    }
  }

  // Verify the size of each NodeBC vectors. They must have the same lengths and
  // equal to the total number of scalar components for all systems of PDEs
  // integrated.
  Assert( std::all_of( begin(dirbc), end(dirbc),
            [ nprop ]( const auto& n ){ return n.second.size() == nprop; } ),
          "Size of NodeBC vector incorrect" );

  return dirbc;
}

void
AirCG::queryBnd()
// *****************************************************************************
// Query/update boundary-conditions-related data structures
// *****************************************************************************
{
  auto d = Disc();

  // Query and match user-specified Dirichlet boundary conditions to side sets
  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  if (steady) for (auto& deltat : m_dtp) deltat *= rkcoef[m_stage];
  m_dirbc = match( m_u.nprop(), d->T(), rkcoef[m_stage] * d->Dt(),
                   m_tp, m_dtp, d->Coord(), d->Lid(), m_bnode );
  if (steady) for (auto& deltat : m_dtp) deltat /= rkcoef[m_stage];

  // Prepare unique set of symmetry BC nodes
  auto sym = d->bcnodes< tag::bc, tag::bcsym >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : sym)
    m_symbcnodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of farfield BC nodes
  auto far = d->bcnodes< tag::bc, tag::bcfarfield >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : far)
    m_farfieldbcnodes.insert( begin(nodes), end(nodes) );

  // If farfield BC is set on a node, will not also set symmetry BC
  for (auto fn : m_farfieldbcnodes) m_symbcnodes.erase(fn);

  // Prepare boundary nodes contiguously accessible from a triangle-face loop
  m_symbctri.resize( m_triinpoel.size()/3, 0 );
  for (std::size_t e=0; e<m_triinpoel.size()/3; ++e)
    if (m_symbcnodes.find(m_triinpoel[e*3+0]) != end(m_symbcnodes))
      m_symbctri[e] = 1;

  // Prepare unique set of sponge nodes
  auto sponge = d->bcnodes< tag::sponge, tag::sideset >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : sponge)
    m_spongenodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of time dependent BC nodes
  const auto& timedep =
    g_inputdeck.template get< tag::param, tag::compflow, tag::bctimedep >();
  if (!timedep.empty()) {
    m_timedepbcnodes.resize(timedep[0].size());
    m_timedepbcFn.resize(timedep[0].size());
    std::size_t ib=0;
    for (const auto& bndry : timedep[0]) {
      std::unordered_set< std::size_t > nodes;
      for (const auto& s : bndry.template get< tag::sideset >()) {
        auto k = m_bnode.find( std::stoi(s) );
        if (k != end(m_bnode)) {
          for (auto g : k->second) {      // global node ids on side set
            nodes.insert( tk::cref_find(d->Lid(),g) );
          }
        }
      }
      m_timedepbcnodes[ib].insert( begin(nodes), end(nodes) );

      // Store user defined discrete function in time. This is done in the same
      // loop as the BC nodes, so that the indices for the two vectors
      // m_timedepbcnodes and m_timedepbcFn are consistent with each other
      auto fn = bndry.template get< tag::fn >();
      for (std::size_t ir=0; ir<fn.size()/6; ++ir) {
        m_timedepbcFn[ib].push_back({{ fn[ir*6+0], fn[ir*6+1], fn[ir*6+2],
          fn[ir*6+3], fn[ir*6+4], fn[ir*6+5] }});
      }
      ++ib;
    }
  }

  Assert(m_timedepbcFn.size() == m_timedepbcnodes.size(), "Incorrect number of "
    "time dependent functions.");
}

void
AirCG::norm()
// *****************************************************************************
// Start (re-)computing boundary point-, and dual-face normals
// *****************************************************************************
{
  auto d = Disc();

  // Query nodes at which symmetry BCs are specified
  auto bn = d->bcnodes< tag::bc, tag::bcsym >( m_bface, m_triinpoel );

  // Query nodes at which farfield BCs are specified
  auto far = d->bcnodes< tag::bc, tag::bcfarfield >( m_bface, m_triinpoel );
  // Merge BC data where boundary-point normals are required
  for (const auto& [s,n] : far) bn[s].insert( begin(n), end(n) );

  // Compute boundary point normals
  bnorm( bn );

  // Compute dual-face normals associated to edges
  dfnorm();
}

std::array< tk::real, 3 >
AirCG::edfnorm( const tk::UnsMesh::Edge& edge,
                const std::unordered_map< tk::UnsMesh::Edge,
                        std::vector< std::size_t >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& esued )
const
// *****************************************************************************
//  Compute normal of dual-mesh associated to edge
//! \param[in] edge Edge whose dual-face normal to compute given by local ids
//! \param[in] esued Elements surrounding edges
//! \return Dual-face normal for edge
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  std::array< tk::real, 3 > n{ 0.0, 0.0, 0.0 };

  for (auto e : tk::cref_find(esued,edge)) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );
    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    // sum normal contributions
    // The constant 1/48: Eq (12) from Waltz et al. Computers & fluids (92) 2014
    // The result of the integral of shape function N on a tet is V/4.
    // This can be written as J/(6*4). Eq (12) has a 1/2 multiplier.
    // This leads to J/48.
    auto J48 = J/48.0;
    for (const auto& [a,b] : tk::lpoed) {
      auto s = tk::orient( {N[a],N[b]}, edge );
      for (std::size_t j=0; j<3; ++j)
        n[j] += J48 * s * (grad[a][j] - grad[b][j]);
    }
  }

  return n;
}

void
AirCG::dfnorm()
// *****************************************************************************
// Compute dual-face normals associated to edges
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& gid = d->Gid();

  // compute derived data structures
  auto esued = tk::genEsued( inpoel, 4, m_esup );

  // Compute dual-face normals for domain edges
  for (std::size_t p=0; p<gid.size(); ++p)    // for each point p
    for (auto q : tk::Around(m_psup,p))       // for each edge p-q
      if (gid[p] < gid[q])
        m_dfnorm[{gid[p],gid[q]}] = edfnorm( {p,q}, esued );

  // Send our dual-face normal contributions to neighbor chares
  if (d->EdgeCommMap().empty())
    comdfnorm_complete();
  else {
    for (const auto& [c,edges] : d->EdgeCommMap()) {
      decltype(m_dfnorm) exp;
      for (const auto& e : edges) exp[e] = tk::cref_find(m_dfnorm,e);
      thisProxy[c].comdfnorm( exp );
    }
  }

  owndfnorm_complete();
}

void
AirCG::comdfnorm( const std::unordered_map< tk::UnsMesh::Edge,
                    std::array< tk::real, 3 >,
                    tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& dfnorm )
// *****************************************************************************
// Receive contributions to dual-face normals on chare-boundaries
//! \param[in] dfnorm Incoming partial sums of dual-face normals associated to
//!   chare-boundary edges
// *****************************************************************************
{
  // Buffer up inccoming contributions to dual-face normals
  for (const auto& [e,n] : dfnorm) {
    auto& dfn = m_dfnormc[e];
    dfn[0] += n[0];
    dfn[1] += n[1];
    dfn[2] += n[2];
  }

  if (++m_ndfnorm == Disc()->EdgeCommMap().size()) {
    m_ndfnorm = 0;
    comdfnorm_complete();
  }
}

void
AirCG::bnorm( const std::unordered_map< int,
                std::unordered_set< std::size_t > >& bcnodes )
// *****************************************************************************
//  Compute boundary point normals
//! \param[in] bcnodes Local node ids associated to side set ids at which BCs
//!    are set that require normals
//*****************************************************************************
{
  auto d = Disc();

  m_bnorm = tk::bnorm( m_bface, m_triinpoel, d->Coord(), d->Gid(), bcnodes );

  // Send our nodal normal contributions to neighbor chares
  if (d->NodeCommMap().empty())
    comnorm_complete();
  else
    for (const auto& [ neighborchare, sharednodes ] : d->NodeCommMap()) {
      std::unordered_map< int,
        std::unordered_map< std::size_t, std::array< tk::real, 4 > > > exp;
      for (auto i : sharednodes) {
        for (const auto& [s,norms] : m_bnorm) {
          auto j = norms.find(i);
          if (j != end(norms)) exp[s][i] = j->second;
        }
      }
      thisProxy[ neighborchare ].comnorm( exp );
    }

  ownnorm_complete();
}

void
AirCG::comnorm( const std::unordered_map< int,
  std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& innorm )
// *****************************************************************************
// Receive boundary point normals on chare-boundaries
//! \param[in] innorm Incoming partial sums of boundary point normal
//!   contributions to normals (first 3 components), inverse distance squared
//!   (4th component), associated to side set ids
// *****************************************************************************
{
  // Buffer up incoming boundary-point normal vector contributions
  for (const auto& [s,norms] : innorm) {
    auto& bnorms = m_bnormc[s];
    for (const auto& [p,n] : norms) {
      auto& bnorm = bnorms[p];
      bnorm[0] += n[0];
      bnorm[1] += n[1];
      bnorm[2] += n[2];
      bnorm[3] += n[3];
    }
  }

  if (++m_nbnorm == Disc()->NodeCommMap().size()) {
    m_nbnorm = 0;
    comnorm_complete();
  }
}

void
AirCG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types initiated from this chare array
//! \details Since this is a [initnode] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  NodeDiagnostics::registerReducers();
}

void
AirCG::ResumeFromSync()
// *****************************************************************************
//  Return from migration
//! \details This is called when load balancing (LB) completes. The presence of
//!   this function does not affect whether or not we block on LB.
// *****************************************************************************
{
  if (Disc()->It() == 0) Throw( "it = 0 in ResumeFromSync()" );

  if (!g_inputdeck.get< tag::cmd, tag::nonblocking >()) next();
}

//! [setup]
void
AirCG::setup()
// *****************************************************************************
// Start setup for solution
// *****************************************************************************
{
  auto d = Disc();

  // Determine nodes inside user-defined IC box
  ICBoxNodes( d->Coord(), m_boxnodes );

  // Compute volume of user-defined box IC
  d->boxvol( m_boxnodes );

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    //for (const auto& eq : g_cgpde) {
    //  auto n = eq.histNames();
    //  histnames.insert( end(histnames), begin(n), end(n) );
    //}
    d->histheader( std::move(histnames) );
  }
}
//! [setup]

void
AirCG::box( tk::real v )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions for all PDEs
  initialize( d->Coord(), m_u, d->T() );

  // Compute boundary point-, and dual-face normals
  norm();
}

//! [start]
void
AirCG::start()
// *****************************************************************************
// Start time stepping
// *****************************************************************************
{
  // Set flag that indicates that we are now during time stepping
  Disc()->Initial( 0 );
  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();
  // Zero grind-timer
  Disc()->grindZero();
  // Continue to first time step
  next();
}
//! [start]

//! [Merge normals and continue]
void
AirCG::mergenorm()
// *****************************************************************************
// The own and communication portion of the left-hand side is complete
// *****************************************************************************
{
  // Combine own and communicated contributions of normals
  normfinal();

  if (Disc()->Initial()) {
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_AirCG::start(), thisProxy[thisIndex]) );
  } else {
    norm_complete();
  }
}
//! [Merge normals and continue]

void
AirCG::normfinal()
// *****************************************************************************
//  Finish computing dual-face and boundary point normals
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

  // Combine own and communicated contributions to boundary point normals
  for (const auto& [s,norms] : m_bnormc) {
    auto& bnorms = m_bnorm[s];
    for (const auto& [p,n] : norms) {
      auto& norm = bnorms[p];
      norm[0] += n[0];
      norm[1] += n[1];
      norm[2] += n[2];
      norm[3] += n[3];
    }
  }
  tk::destroy( m_bnormc );

  // Divide summed point normals by the sum of inverse distance squared
  for (auto& [s,norms] : m_bnorm)
    for (auto& [p,n] : norms) {
      n[0] /= n[3];
      n[1] /= n[3];
      n[2] /= n[3];
      Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
              1.0e+3*std::numeric_limits< tk::real >::epsilon(),
              "Non-unit normal" );
    }

  // Replace global->local ids associated to boundary point normals
  decltype(m_bnorm) bnorm;
  for (auto& [s,norms] : m_bnorm) {
    auto& bnorms = bnorm[s];
    for (auto&& [g,n] : norms)
      bnorms[ tk::cref_find(lid,g) ] = std::move(n);
  }
  m_bnorm = std::move(bnorm);

  // Count contributions to chare-boundary edges
  std::unordered_map< tk::UnsMesh::Edge, std::size_t,
    tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > edge_node_count;
  for (const auto& [c,edges] : d->EdgeCommMap())
    for (const auto& e : edges)
      ++edge_node_count[e];

  // Combine and weigh communicated contributions to dual-face normals
  for (auto& [e,n] : m_dfnormc) {
    const auto& dfn = tk::cref_find( m_dfnorm, e );
    n[0] += dfn[0];
    n[1] += dfn[1];
    n[2] += dfn[2];
    auto count = static_cast< tk::real >( tk::cref_find( edge_node_count, e ) );
    auto factor = 1.0/(count + 1.0);
    for (auto & x : n) x *= factor;
  }

  // Generate list of unique edges
  tk::UnsMesh::EdgeSet uedge;
  for (std::size_t p=0; p<m_u.nunk(); ++p)
    for (auto q : tk::Around(m_psup,p))
      uedge.insert( {p,q} );

  // Flatten edge list
  m_edgenode.resize( uedge.size() * 2 );
  std::size_t f = 0;
  const auto& gid = d->Gid();
  for (auto&& [p,q] : uedge) {
    if (gid[p] > gid[q]) {
      m_edgenode[f+0] = std::move(q);
      m_edgenode[f+1] = std::move(p);
    } else {
      m_edgenode[f+0] = std::move(p);
      m_edgenode[f+1] = std::move(q);
    }
    f += 2;
  }
  tk::destroy(uedge);

  // Convert dual-face normals to streamable (and vectorizable) data structure
  m_dfn.resize( m_edgenode.size() * 3 );      // 2 vectors per access
  std::unordered_map< tk::UnsMesh::Edge, std::size_t,
                      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > eid;
  for (std::size_t e=0; e<m_edgenode.size()/2; ++e) {
    auto p = m_edgenode[e*2+0];
    auto q = m_edgenode[e*2+1];
    eid[{p,q}] = e;
    std::array< std::size_t, 2 > g{ gid[p], gid[q] };
    auto n = tk::cref_find( m_dfnorm, g );
    // figure out if this is an edge on the parallel boundary
    auto nit = m_dfnormc.find( g );
    auto m = ( nit != m_dfnormc.end() ) ? nit->second : n;
    m_dfn[e*6+0] = n[0];
    m_dfn[e*6+1] = n[1];
    m_dfn[e*6+2] = n[2];
    m_dfn[e*6+3] = m[0];
    m_dfn[e*6+4] = m[1];
    m_dfn[e*6+5] = m[2];
  }

  tk::destroy( m_dfnorm );
  tk::destroy( m_dfnormc );

  // Flatten edge id data structure
  m_edgeid.resize( m_psup.first.size() );
  for (std::size_t p=0,k=0; p<m_u.nunk(); ++p)
    for (auto q : tk::Around(m_psup,p))
      m_edgeid[k++] = tk::cref_find( eid, {p,q} );
}

void
AirCG::BC()
// *****************************************************************************
// Apply boundary conditions
// \details The following BC enforcement changes the initial condition or
//!   updated solution (dependending on when it is called) to ensure strong
//!   imposition of the BCs. This is a matter of choice. Another alternative is
//!   to only apply BCs when computing fluxes at boundary faces, thereby only
//!   weakly enforcing the BCs. The former is conventionally used in continunous
//!   Galerkin finite element methods (such as AirCG implements), whereas the
//!   latter, in finite volume methods.
// *****************************************************************************
{
  auto d = Disc();

  // Apply Dirichlet BCs
  for (const auto& [b,bc] : m_dirbc)
    for (ncomp_t c=0; c<m_u.nprop(); ++c)
      if (bc[c].first) m_u(b,c,0) = bc[c].second;

  const auto& coord = d->Coord();

  // Apply symmetry BCs
  symbc( m_u, coord, m_bnorm, m_symbcnodes );

  // Apply farfield BCs
  farfieldbc( m_u, coord, m_bnorm, m_farfieldbcnodes );

  // Apply sponge conditions
  sponge( m_u, coord, m_spongenodes );

  // Apply user defined time dependent BCs
  timedepbc( d->T(), m_u, m_timedepbcnodes, m_timedepbcFn );
}

void
AirCG::next()
// *****************************************************************************
// Continue to next time step
// *****************************************************************************
{
  dt();
}

void
AirCG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  auto d = Disc();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    //! [Find the minimum dt across all PDEs integrated]
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

      // compute new dt for each mesh point
      inciter::dt( d->It(), d->Vol(), m_u, m_dtp );

      // find the smallest dt of all nodes on this chare
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );

    } else {    // compute new dt for this chare

      // find the smallest dt of all equations on this chare
      mindt = inciter::dt( d->Coord(), d->Inpoel(), d->T(), m_u );

    }
    //! [Find the minimum dt across all PDEs integrated]

  }

  //! [Advance]
  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();
  thisProxy[ thisIndex ].wait4stage();

  // Contribute to minimum dt across all chares and advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(AirCG,advance), thisProxy) );
  //! [Advance]
}

void
AirCG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  // Compute gradients for next time step
  chBndGrad();
}

void
AirCG::chBndGrad()
// *****************************************************************************
// Compute nodal gradients at chare-boundary nodes. Gradients at internal nodes
// are calculated locally as needed and are not stored.
// *****************************************************************************
{
  auto d = Disc();

  // Compute own portion of gradients for all equations
  bndgrad( d->Coord(), d->Inpoel(), m_bndel, d->Gid(), d->Bid(), m_u,
           m_chBndGrad );

  // Communicate gradients to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comgrad_complete();
  else // send gradients contributions to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > g( n.size() );
      std::size_t j = 0;
      for (auto i : n) g[ j++ ] = m_chBndGrad[ tk::cref_find(d->Bid(),i) ];
      thisProxy[c].comChBndGrad( std::vector<std::size_t>(begin(n),end(n)), g );
    }

  owngrad_complete();
}

void
AirCG::comChBndGrad( const std::vector< std::size_t >& gid,
                     const std::vector< std::vector< tk::real > >& G )
// *****************************************************************************
//  Receive contributions to nodal gradients on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive grad contributions
//! \param[in] G Partial contributions of gradients to chare-boundary nodes
//! \details This function receives contributions to m_chBndGrad, which stores
//!   nodal gradients at mesh chare-boundary nodes. While m_chBndGrad stores
//!   own contributions, m_chBndGradc collects the neighbor chare
//!   contributions during communication. This way work on m_chBndGrad and
//!   m_chBndGradc is overlapped. The two are combined in rhs().
// *****************************************************************************
{
  Assert( G.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) m_chBndGradc[ gid[i] ] += G[i];

  if (++m_ngrad == Disc()->NodeCommMap().size()) {
    m_ngrad = 0;
    comgrad_complete();
  }
}

void
AirCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions to nodal gradients
  for (const auto& [gid,g] : m_chBndGradc) {
    auto bid = tk::cref_find( d->Bid(), gid );
    for (ncomp_t c=0; c<m_chBndGrad.nprop(); ++c)
      m_chBndGrad(bid,c,0) += g[c];
  }

  // clear gradients receive buffer
  tk::destroy(m_chBndGradc);

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  inciter::rhs( d->T() + prev_rkcoef * d->Dt(), d->Coord(), d->Inpoel(),
       m_triinpoel, d->Gid(), d->Bid(), d->Lid(), m_dfn, m_esup, m_psup,
       m_symbctri, m_spongenodes, d->Vol(), m_edgenode, m_edgeid, m_tp,
       m_chBndGrad, m_u, m_rhs );

  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];

  // Query/update boundary-conditions-related data structures from user input
  queryBnd();

  // Communicate rhs to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comrhs_complete();
  else // send contributions of rhs to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > r( n.size() );
      std::size_t j = 0;
      for (auto i : n) r[ j++ ] = m_rhs[ tk::cref_find(d->Lid(),i) ];
      thisProxy[c].comrhs( std::vector<std::size_t>(begin(n),end(n)), r );
    }

  ownrhs_complete();
}

void
AirCG::comrhs( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& R )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] R Partial contributions of RHS to chare-boundary nodes
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped. The two
//!   are combined in solve().
// *****************************************************************************
{
  Assert( R.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) m_rhsc[ gid[i] ] += R[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nrhs == Disc()->NodeCommMap().size()) {
    m_nrhs = 0;
    comrhs_complete();
  }
}

void
AirCG::solve()
// *****************************************************************************
//  Advance systems of equations
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions to rhs
  for (const auto& b : m_rhsc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<m_rhs.nprop(); ++c) m_rhs(lid,c,0) += b.second[c];
  }

  // clear receive buffer
  tk::destroy(m_rhsc);

  // Update state at time n
  if (m_stage == 0) m_un = m_u;

  // Solve the sytem
  if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

    // Advance solution, converging to steady state
    for (std::size_t i=0; i<m_u.nunk(); ++i)
      for (ncomp_t c=0; c<m_u.nprop(); ++c)
        m_u(i,c,0) = m_un(i,c,0) + rkcoef[m_stage] * m_dtp[i] * m_rhs(i,c,0);

  } else {

    // Advance unsteady solution
    const auto& V = d->Vol();
    auto adt = rkcoef[m_stage] * d->Dt();
    for (std::size_t i=0; i<m_u.nunk(); ++i)
      for (ncomp_t c=0; c<5; ++c)
        m_u(i,c,0) = m_un(i,c,0) + adt * m_rhs(i,c,0) / V[i];

  }

  // Enforce boundary conditions
  BC();

  if (m_stage < 2) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4grad();
    thisProxy[ thisIndex ].wait4rhs();

    // start next time step stage
    stage();

  } else {

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u, m_un, m_bnorm,
                                         m_symbcnodes, m_farfieldbcnodes );
    // Increase number of iterations and physical time
    d->next();
    // Advance physical time for local time stepping
    if (g_inputdeck.get< tag::discr, tag::steady_state >())
      for (std::size_t i=0; i<m_u.nunk(); ++i) m_tp[i] += m_dtp[i];
    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 1.0 ) );

  }
}

//! [Refine]
void
AirCG::refine( const std::vector< tk::real >& l2res )
// *****************************************************************************
// Optionally refine/derefine mesh
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  const auto residual = g_inputdeck.get< tag::discr, tag::residual >();
  const auto rc = g_inputdeck.get< tag::discr, tag::rescomp >() - 1;

  if (steady) {

    // this is the last time step if max time of max number of time steps
    // reached or the residual has reached its convergence criterion
    if (d->finished() or l2res[rc] < residual) m_finished = 1;

  } else {

    // this is the last time step if max time or max iterations reached
    if (d->finished()) m_finished = 1;

  }

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // Activate SDAG waits for re-computing the normals
  thisProxy[ thisIndex ].wait4norm();

  // if t>0 refinement enabled and we hit the frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    d->startvol();
    d->Ref()->dtref( m_bface, m_bnode, m_triinpoel );
    d->refined() = 1;

  } else {      // do not refine

    d->refined() = 0;
    norm_complete();
    resized();

  }
}
//! [Refine]

//! [Resize]
void
AirCG::resizePostAMR(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
  const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
  const std::set< std::size_t >& removedNodes,
  const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] removedNodes Newly removed mesh node local ids
//! \param[in] amrNodeMap Node id map after amr (local ids)
//! \param[in] nodeCommMap New node communication map
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] bnode Boundary-node lists mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  auto d = Disc();

  d->Itf() = 0;  // Zero field output iteration count if AMR
  ++d->Itr();    // Increase number of iterations with a change in the mesh

  // Resize mesh data structures after mesh refinement
  d->resizePostAMR( chunk, coord, amrNodeMap, nodeCommMap, removedNodes );

  // Remove newly removed nodes from solution vectors
  m_u.rm(removedNodes);
  m_un.rm(removedNodes);
  m_rhs.rm(removedNodes);

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  m_u.resize( npoin );
  m_un.resize( npoin );
  m_rhs.resize( npoin );
  m_chBndGrad.resize( d->Bid().size() );
  tk::destroy(m_psup);
  m_esup = tk::genEsup( d->Inpoel(), 4 );
  m_psup = tk::genPsup( d->Inpoel(), 4, m_esup );

  // Update solution on new mesh
  for (const auto& n : addedNodes)
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      Assert(n.first < m_u.nunk(), "Added node index out of bounds post-AMR");
      Assert(n.second[0] < m_u.nunk() && n.second[1] < m_u.nunk(),
        "Indices of parent-edge nodes out of bounds post-AMR");
      m_u(n.first,c,0) = (m_u(n.second[0],c,0) + m_u(n.second[1],c,0))/2.0;
    }

  // Update physical-boundary node-, face-, and element lists
  m_bnode = bnode;
  m_bface = bface;
  m_triinpoel = tk::remap( triinpoel, d->Lid() );

  auto meshid = d->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
              CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );
}
//! [Resize]

void
AirCG::resized()
// *****************************************************************************
// Resizing data sutrctures after mesh refinement has been completed
// *****************************************************************************
{
  resize_complete();
}

//! [stage]
void
AirCG::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  // Increment Runge-Kutta stage counter
  ++m_stage;

  // if not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise output field data to file(s)
  if (m_stage < 3) chBndGrad(); else out();
}
//! [stage]

void
AirCG::writeFields( CkCallback c )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    c.send();

  } else {

    auto d = Disc();
    const auto& coord = d->Coord();

    std::vector< std::string > nodefieldnames{ "r", "u", "v", "w", "E", "p" };

    using tk::operator/=;
    auto r = m_u.extract( 0, 0 );
    auto u = m_u.extract( 1, 0 );  u /= r;
    auto v = m_u.extract( 2, 0 );  v /= r;
    auto w = m_u.extract( 3, 0 );  w /= r;
    auto E = m_u.extract( 4, 0 );
    auto p = r;
    for (std::size_t i=0; i<m_u.nunk(); ++i)
      p[i] = eos_pressure( r[i], u[i], v[i], w[i], E[i] );
    E /= r;

    std::vector< std::vector< tk::real > > nodefields{
      std::move(r), std::move(u), std::move(v), std::move(w), std::move(E),
      std::move(p) };

    ////! Lambda to put in a field for output if not empty
    //auto add_node_field = [&]( const auto& name, const auto& field ){
    //  if (not field.empty()) {
    //    nodefieldnames.push_back( name );
    //    nodefields.push_back( field );
    //  }
    //};

    //// Collect field output names for analytical solutions
    //for (const auto& eq : g_cgpde) analyticFieldNames( eq, nodefieldnames );

    //// Collect field output from analytical solutions (if exist)
    //for (const auto& eq : g_cgpde)
    //  analyticFieldOutput( eq, coord[0], coord[1],
    //                       coord[2], d->T(), nodefields );

    //// Query and collect block and surface field names from PDEs integrated
    std::vector< std::string > nodesurfnames;
    //for (const auto& eq : g_cgpde) {
    //  auto s = eq.surfNames();
    //  nodesurfnames.insert( end(nodesurfnames), begin(s), end(s) );
    //}

    //// Collect node block and surface field solution
    std::vector< std::vector< tk::real > > nodesurfs;
    //for (const auto& eq : g_cgpde) {
    //  auto s = eq.surfOutput( tk::bfacenodes(m_bface,m_triinpoel), m_u );
    //  nodesurfs.insert( end(nodesurfs), begin(s), end(s) );
    //}

    Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

    // Send mesh and fields data (solution dump) for output to file
    d->write( d->Inpoel(), coord, m_bface, tk::remap(m_bnode,d->Lid()),
              m_triinpoel, {}, nodefieldnames, nodesurfnames, {}, nodefields,
              nodesurfs, c );

  }
}

void
AirCG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    //for (const auto& eq : g_cgpde) {
    //  auto h = eq.histOutput( d->Hist(), d->Inpoel(), m_u );
    //  hist.insert( end(hist), begin(h), end(h) );
    //}
    d->history( std::move(hist) );
  }

  // Output field data
  if (d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished)
    writeFields( CkCallback(CkIndex_AirCG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
AirCG::evalLB( int nrestart )
// *****************************************************************************
// Evaluate whether to do load balancing
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  auto d = Disc();

  // Detect if just returned from a checkpoint and if so, zero timers and
  // finished flag
  if (d->restarted( nrestart )) m_finished = 0;

  const auto lbfreq = g_inputdeck.get< tag::cmd, tag::lbfreq >();
  const auto nonblocking = g_inputdeck.get< tag::cmd, tag::nonblocking >();

  // Load balancing if user frequency is reached or after the second time-step
  if ( (d->It()) % lbfreq == 0 || d->It() == 2 ) {

    AtSync();
    if (nonblocking) next();

  } else {

    next();

  }
}

void
AirCG::evalRestart()
// *****************************************************************************
// Evaluate whether to save checkpoint/restart
// *****************************************************************************
{
  auto d = Disc();

  const auto rsfreq = g_inputdeck.get< tag::cmd, tag::rsfreq >();
  const auto benchmark = g_inputdeck.get< tag::cmd, tag::benchmark >();

  if ( !benchmark && (d->It()) % rsfreq == 0 ) {

    std::vector< std::size_t > meshdata{ /* finished = */ 0, d->MeshId() };
    contribute( meshdata, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,checkpoint), d->Tr()) );

  } else {

    evalLB( /* nrestart = */ -1 );

  }
}

void
AirCG::step()
// *****************************************************************************
// Evaluate whether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  d->status();
  // Reset Runge-Kutta stage counter
  m_stage = 0;

  if (not m_finished) {

    evalRestart();

  } else {

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/aircg.def.h"
