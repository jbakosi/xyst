// *****************************************************************************
/*!
  \file      src/Inciter/AirCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     AirCG: continuous Galerkin finite elements + Runge Kutta
  \details   AirCG solves the compressible Euler or Navier-Stokes equations
    coupled to a number of scalars usng a continuous Galerkin (CG) finite
    element (FE) spatial discretization (using linear shapefunctions on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping scheme.
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
#include "Problems.hpp"
#include "EOS.hpp"
#include "BC.hpp"

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
  m_nrhs( 0 ),
  m_nnorm( 0 ),
  m_nbpint( 0 ),
  m_nbeint( 0 ),
  m_ndeint( 0 ),
  m_ngrad( 0 ),
  m_ncomp( 5 +
     (g_inputdeck.get<tag::component>().get<tag::transport>().empty() ? 0 :
      g_inputdeck.get<tag::component>().get<tag::transport>()[0] ) ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_u( Disc()->Gid().size(), m_ncomp ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_grad( m_u.nunk(), m_u.nprop()*3 ),
  m_stage( 0 ),
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

  // Create new local ids based on mesh locality
  std::unordered_map< std::size_t, std::size_t > map;
  std::size_t n = 0;

  auto psup = tk::genPsup( d->Inpoel(), 4, tk::genEsup( d->Inpoel(), 4 ) );
  for (std::size_t p=0; p<m_u.nunk(); ++p) {  // for each point p
    if (!map.count(p)) map[p] = n++;
    for (auto q : tk::Around(psup,p)) {       // for each edge p-q
      if (!map.count(q)) map[q] = n++;
    }
  }

  Assert( map.size() == d->Gid().size(),
          "Mesh-locality reorder map size mismatch" );

  // Remap data in bound Discretization object
  d->remap( map );
  // Remap boundary triangle face connectivity
  tk::remap( m_triinpoel, map );

  // Activate SDAG wait for initially computing integrals
  thisProxy[ thisIndex ].wait4int();

  // Signal the runtime system that the workers have been created
  auto meshid = d->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::sum_ulong,
    CkCallback(CkReductionTarget(Transporter,comfinal), d->Tr()) );
}
//! [Constructor]

void
AirCG::integrals()
// *****************************************************************************
// Start (re-)computing domain and boundary integrals
// *****************************************************************************
{
  auto d = Disc();

  // Query BC nodes associated to side sets
  auto dir = d->bcnodes< tag::bc, tag::bcdir >( m_bface, m_triinpoel );
  auto sym = d->bcnodes< tag::bc, tag::bcsym >( m_bface, m_triinpoel );
  auto far = d->bcnodes< tag::bc, tag::bcfarfield >( m_bface, m_triinpoel );

  // Compile unique set of BC nodes
  for (const auto& [s,n] : dir)
    m_dirbcnodes.insert( end(m_dirbcnodes), begin(n), end(n) );
  tk::unique( m_dirbcnodes );

  for (const auto& [s,n] : sym) m_symbcnodeset.insert( begin(n), end(n) );
  for (const auto& [s,n] : far) m_farbcnodeset.insert( begin(n), end(n) );

  // If farfield BC is set on a node, will not also set symmetry BC
  for (auto i : m_farbcnodeset) m_symbcnodeset.erase(i);

  // Compute local contributions to boundary normals and integrals
  bndint();
  // Compute local contributions to domain edge integrals
  domint();

  // Send boundary point normal contributions to neighbor chares
  if (d->NodeCommMap().empty()) {
    comnorm_complete();
  } else {
    for (const auto& [c,nodes] : d->NodeCommMap()) {
      decltype(m_bnorm) exp;
      for (auto i : nodes) {
        for (const auto& [s,b] : m_bnorm) {
          auto k = b.find(i);
          if (k != end(b)) exp[s][i] = k->second;
        }
      }
      thisProxy[c].comnorm( exp );
    }
  }
  ownnorm_complete();
}

void
AirCG::bndint()
// *****************************************************************************
//! Compute local contributions to boundary normals and integrals
// *****************************************************************************
{
  auto d = Disc();
  const auto& coord = d->Coord();
  const auto& gid = d->Gid();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Lambda to compute the inverse distance squared between boundary face
  // centroid and boundary point. Here p is the global node id and c is the
  // the boundary face centroid.
  auto invdistsq = [&]( const std::array< tk::real, 3 >& c, std::size_t p ){
    return 1.0 / ( (c[0] - x[p]) * (c[0] - x[p]) +
                   (c[1] - y[p]) * (c[1] - y[p]) +
                   (c[2] - z[p]) * (c[2] - z[p]) );
  };

  tk::destroy( m_bnorm );
  tk::destroy( m_bndpoinint );
  tk::destroy( m_bndedgeint );

  // Compute boundary point normals and boundary point-, and edge-integrals.
  // The boundary point normals are computed by summing
  // inverse-distance-weighted boundary face normals to boundary points. Note
  // that these are only partial sums at shared boundary points and edges in
  // parallel.

   for (const auto& [ setid, faceids ] : m_bface) { // for all side sets
     for (auto f : faceids) { // for all side set triangles
 
       const std::array< std::size_t, 3 >
         N{ m_triinpoel[f*3+0], m_triinpoel[f*3+1], m_triinpoel[f*3+2] };
       const std::array< tk::real, 3 >
         ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
         ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
       auto n = tk::cross( ba, ca );
       auto A = tk::length( n );
       n[0] /= A;
       n[1] /= A;
       n[2] /= A;
       A /= 2.0;
       const std::array< tk::real, 3 > centroid{
         (x[N[0]] + x[N[1]] + x[N[2]]) / 3.0,
         (y[N[0]] + y[N[1]] + y[N[2]]) / 3.0,
         (z[N[0]] + z[N[1]] + z[N[2]]) / 3.0 };

       for (std::size_t j=0; j<3; ++j) {
         auto p = N[j];
         tk::real r = invdistsq( centroid, p );
         auto& v = m_bnorm[setid];      // associate side set id
         auto& bpn = v[gid[p]];         // associate global node id of bnd pnt
         bpn[0] += r * n[0];            // inv.dist.sq-weighted normal
         bpn[1] += r * n[1];
         bpn[2] += r * n[2];
         bpn[3] += r;                   // inv.dist.sq of node from centroid
         auto& b = m_bndpoinint[gid[p]];// assoc global id of bnd point
         b[0] += n[0] * A / 3.0;        // bnd-point integral
         b[1] += n[1] * A / 3.0;
         b[2] += n[2] * A / 3.0;
         auto q = N[ tk::lpoet[j][1] ]; // the other node of bnd edge
         tk::UnsMesh::Edge ed{ gid[p], gid[q] };
         tk::real sig = ed[0] < ed[1] ? 1.0 : -1.0;
         if (ed[0] > ed[1]) std::swap( ed[0], ed[1] );
         auto& e = m_bndedgeint[ ed ];
         e[0] += sig * n[0] * A / 12.0; // bnd-edge integral
         e[1] += sig * n[1] * A / 12.0;
         e[2] += sig * n[2] * A / 12.0;
       }

    }
  }
}

void
AirCG::domint()
// *****************************************************************************
//! Compute local contributions to domain edge integrals
// *****************************************************************************
{
  auto d = Disc();

  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();

  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  tk::destroy( m_domedgeint );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
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
    auto J48 = J/48.0;
    for (const auto& [p,q] : tk::lpoed) {
      tk::UnsMesh::Edge ed{ gid[N[p]], gid[N[q]] };
      tk::real sig = ed[0] < ed[1] ? 1.0 : -1.0;
      if (ed[0] > ed[1]) std::swap( ed[0], ed[1] );
      auto& n = m_domedgeint[ ed ];
      for (std::size_t j=0; j<3; ++j)
        n[j] += J48 * sig * (grad[p][j] - grad[q][j]);
    }
  }
}

void
AirCG::comnorm( const decltype(m_bnorm)& inbnd )
// *****************************************************************************
// Receive contributions to boundary point normals on chare-boundaries
//! \param[in] inbnd Incoming partial sums of boundary point normals
// *****************************************************************************
{
  // Buffer up incoming boundary point normal vector contributions
  for (const auto& [s,b] : inbnd) {
    auto& bndnorm = m_bnormc[s];
    for (const auto& [p,n] : b) {
      auto& norm = bndnorm[p];
      norm[0] += n[0];
      norm[1] += n[1];
      norm[2] += n[2];
      norm[3] += n[3];
    }
  }

  if (++m_nnorm == Disc()->NodeCommMap().size()) {
    m_nnorm = 0;
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

  // Set initial conditions
  problems::initialize( d->Coord(), m_u, d->T() );

  // Compute volume of user-defined box IC
  d->boxvol( m_boxnodes );

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    d->histheader( { "density", "xvelocity", "yvelocity", "zvelocity",
                     "energy", "pressure" } );
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
  // Store user-defined box IC volume
  Disc()->Boxvol() = v;

  // Compute edge integrals
  integrals();
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
AirCG::merge()
// *****************************************************************************
// Combine own and communicated portions of the integrals
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();

  // Combine own and communicated contributions to boundary point normals
  for (const auto& [s,b] : m_bnormc) {
    auto& bndnorm = m_bnorm[s];
    for (const auto& [g,n] : b) {
      auto& norm = bndnorm[g];
      norm[0] += n[0];
      norm[1] += n[1];
      norm[2] += n[2];
      norm[3] += n[3];
    }
  }
  tk::destroy( m_bnormc );

  // Divide summed point normals by the sum of the inverse distance squared
  for (auto& [s,b] : m_bnorm) {
    for (auto& [g,n] : b) {
      n[0] /= n[3];
      n[1] /= n[3];
      n[2] /= n[3];
      Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
              1.0e+3*std::numeric_limits< tk::real >::epsilon(),
              "Non-unit normal" );
    }
  }

  // Replace global->local ids associated to boundary point normals
  decltype(m_bnorm) loc;
  for (auto& [s,b] : m_bnorm) {
    auto& bnd = loc[s];
    for (auto&& [g,n] : b) {
      bnd[ tk::cref_find(lid,g) ] = std::move(n);
    }
  }
  m_bnorm = std::move(loc);

  // Convert boundary point integrals into streamable data structures
  m_bpoin.resize( m_bndpoinint.size() );
  m_bpsym.resize( m_bndpoinint.size() );
  m_bpint.resize( m_bndpoinint.size() * 3 );
  std::size_t i = 0;
  for (const auto& [g,b] : m_bndpoinint) {
    m_bpoin[i] = tk::cref_find( lid, g );
    m_bpsym[i] = static_cast< std::uint8_t >(m_symbcnodeset.count(m_bpoin[i]));
    m_bpint[i*3+0] = b[0];
    m_bpint[i*3+1] = b[1];
    m_bpint[i*3+2] = b[2];
    ++i;
  }
  tk::destroy( m_bndpoinint );

  // Convert boundary edge integrals into streamable data structures
  m_bedge.resize( m_bndedgeint.size() * 2 );
  m_besym.resize( m_bndedgeint.size() * 2 );
  m_beint.resize( m_bndedgeint.size() * 3 );
  std::size_t j = 0;
  for (const auto& [ed,b] : m_bndedgeint) {
    auto p = tk::cref_find( lid, ed[0] );
    auto q = tk::cref_find( lid, ed[1] );
    m_bedge[j*2+0] = p;
    m_bedge[j*2+1] = q;
    m_besym[j*2+0] = static_cast< std::uint8_t >( m_symbcnodeset.count(p) );
    m_besym[j*2+1] = static_cast< std::uint8_t >( m_symbcnodeset.count(q) );
    m_beint[j*3+0] = b[0];
    m_beint[j*3+1] = b[1];
    m_beint[j*3+2] = b[2];
    ++j;
  }
  tk::destroy( m_bndedgeint );

  // Convert domain edge integrals into streamable data structures
  m_dedge.resize( m_domedgeint.size() * 2 );
  m_deint.resize( m_domedgeint.size() * 3 );
  std::size_t k = 0;
  for (const auto& [ed,d] : m_domedgeint) {
    auto p = tk::cref_find( lid, ed[0] );
    auto q = tk::cref_find( lid, ed[1] );
    m_dedge[k*2+0] = p;
    m_dedge[k*2+1] = q;
    m_deint[k*3+0] = d[0];
    m_deint[k*3+1] = d[1];
    m_deint[k*3+2] = d[2];
    ++k;
  }
  tk::destroy( m_domedgeint );

  // Convert symmetry BC data to streamable data structures
  const auto& sbc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, tag::bcsym >();
  for (auto p : m_symbcnodeset) {
    for (const auto& s : sbc[0]) {
      auto m = m_bnorm.find(std::stoi(s));
      if (m != end(m_bnorm)) {
        auto r = m->second.find(p);
        if (r != end(m->second)) {
          m_symbcnodes.push_back( p );
          m_symbcnorms.push_back( r->second[0] );
          m_symbcnorms.push_back( r->second[1] );
          m_symbcnorms.push_back( r->second[2] );
        }
      }
    }
  }
  tk::destroy( m_symbcnodeset );

  // Convert farfield BC data to streamable data structures
  const auto& fbc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, tag::bcfarfield >();
  for (auto p : m_farbcnodeset) {
    for (const auto& s : fbc[0]) {
      auto n = m_bnorm.find(std::stoi(s));
      if (n != end(m_bnorm)) {
        auto a = n->second.find(p);
        if (a != end(n->second)) {
          m_farbcnodes.push_back( p );
          m_farbcnorms.push_back( a->second[0] );
          m_farbcnorms.push_back( a->second[1] );
          m_farbcnorms.push_back( a->second[2] );
        }
      }
    }
  }
  tk::destroy( m_farbcnodeset );

  if (Disc()->Initial()) {
    // Enforce boundary conditions on initial conditions
    BC();
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_AirCG::start(), thisProxy[thisIndex]) );
  } else {
    //integrals_complete();
  }
}
//! [Merge normals and continue]

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
  const auto& coord = d->Coord();

  // Apply Dirichlet BCs
  auto t = d->T() + rkcoef[m_stage] * d->Dt();
  physics::dirbc( m_u, t, coord, m_dirbcnodes );

  // Apply symmetry BCs
  physics::symbc( m_u, m_symbcnodes, m_symbcnorms );

  // Apply farfield BCs
  physics::farbc( m_u, m_farbcnodes, m_farbcnorms );
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
      physics::dt( d->Vol(), m_u, m_dtp );

      // find the smallest dt of all nodes on this chare
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );

    } else {    // compute new dt for this chare

      // find the smallest dt of all equations on this chare
      mindt = physics::dt( d->Vol(), m_u );

    }
    //! [Find the minimum dt across all PDEs integrated]

  }

  //! [Advance]
  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

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
  // Set new time step size
  if (m_stage == 0) Disc()->setdt( newdt );

  // Compute gradients for next time step stage
  grad();
}

void
AirCG::grad()
// *****************************************************************************
// Compute gradients for next time step
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

  physics::grad( m_dedge, m_deint, m_bpoin, m_bpint, m_bedge, m_beint, m_u,
                 m_grad );

  // Send gradient contributions to neighbor chares
  if (d->NodeCommMap().empty()) {
    comgrad_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_gradc) exp;
      for (auto g : n) exp[g] = m_grad[ tk::cref_find(lid,g) ];
      thisProxy[c].comgrad( exp );
    }
  }

  owngrad_complete();
}

void
AirCG::comgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to node gradients on chare-boundaries
//! \param[in] ingrad Partial contributions to chare-boundary nodes. Key: 
//!   global mesh node IDs, value: contributions for all scalar components.
//! \param[in] inrhs Partial contributions of gradients at chare-boundary nodes
//! \details This function receives contributions to m_grad, which stores the
//!   gradients at mesh nodes. While m_grad stores own contributions, m_gradc
//!   collects the neighbor chare contributions during communication. This way
//!   work on m_grad and m_gradc is overlapped. The two are combined in rhs().
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : ingrad) m_gradc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
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
  const auto& lid = d->Lid();

  // Combine own and communicated contributions to gradients
  for (const auto& [g,r] : m_gradc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_grad(i,c,0) += r[c];
  }
  tk::destroy(m_gradc);

  // divide weak result in gradients by nodal volume
  const auto& vol = d->Vol();
  for (std::size_t p=0; p<m_grad.nunk(); ++p)
    for (std::size_t c=0; c<m_grad.nprop(); ++c)
      m_grad(p,c,0) /= vol[p];

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];

  if (steady) {
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  }

  physics::rhs( m_dedge, m_deint, m_bpoin, m_bpint, m_bedge, m_beint, m_bpsym,
    m_besym, d->Coord(), m_grad, m_u, d->V(), d->T(), m_tp, m_rhs );

  if (steady) {
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];
  }

  // Communicate rhs to other chares on chare-boundary
  if (d->NodeCommMap().empty()) {
    comrhs_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_rhsc) exp;
      for (auto g : n) exp[g] = m_rhs[ tk::cref_find(lid,g) ];
      thisProxy[c].comrhs( exp );
    }
  }
  ownrhs_complete();
}

void
AirCG::comrhs(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& inrhs )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] inrhs Partial contributions of RHS to chare-boundary nodes. Key: 
//!   global mesh node IDs, value: contributions for all scalar components.
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped. The two
//!   are combined in solve().
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : inrhs) m_rhsc[g] += r;

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
  const auto lid = d->Lid();

  // Combine own and communicated contributions to rhs
  for (const auto& [g,r] : m_rhsc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_rhs(i,c,0) += r[c];
  }
  tk::destroy(m_rhsc);

  // divide weak result in rhs by nodal volume
  const auto& vol = d->Vol();
  for (std::size_t p=0; p<m_rhs.nunk(); ++p)
    for (std::size_t c=0; c<m_rhs.nprop(); ++c)
      m_rhs(p,c,0) /= vol[p];

  // Update state at time n
  if (m_stage == 0) m_un = m_u;

  // Solve the sytem
  if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

    // Advance solution, converging to steady state
    for (std::size_t i=0; i<m_u.nunk(); ++i) {
      for (ncomp_t c=0; c<m_u.nprop(); ++c) {
        m_u(i,c,0) = m_un(i,c,0) - rkcoef[m_stage] * m_dtp[i] * m_rhs(i,c,0);
      }
    }

  } else {

    // Advance unsteady solution
    m_u = m_un - rkcoef[m_stage] * d->Dt() * m_rhs;

  }

  // Enforce boundary conditions
  BC();

  // Activate SDAG wait for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  if (m_stage < 2) {

    // start next time step stage
    stage();

  } else {

    // Activate SDAG waits for finishing a this time step stage
    thisProxy[ thisIndex ].wait4stage();
    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u, m_un );
    // Increase number of iterations and physical time
    d->next();
    // Advance physical time for local time stepping
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) {
      using tk::operator+=;
      m_tp += m_dtp;
    }
    // Continue to mesh refinement
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

  if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

    const auto residual = g_inputdeck.get< tag::discr, tag::residual >();
    const auto rc = g_inputdeck.get< tag::discr, tag::rescomp >() - 1;

    // this is the last time step if max time of max number of time steps
    // reached or the residual has reached its convergence criterion
    if (d->finished() or l2res[rc] < residual) m_finished = 1;

  } else {

    // this is the last time step if max time or max iterations reached
    if (d->finished()) m_finished = 1;

  }

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // if t>0 refinement enabled and we hit the frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    d->startvol();
    d->Ref()->dtref( m_bface, m_bnode, m_triinpoel );
    d->refined() = 1;

    // Activate SDAG waits for re-computing the integrals
    thisProxy[ thisIndex ].wait4int();

  } else {      // do not refine

    d->refined() = 0;
    integrals_complete();
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
  m_grad.resize( npoin );

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
  m_triinpoel = triinpoel;

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

void
AirCG::writeFields( CkCallback cb )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] cb Function to continue with after the write
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    cb.send();

  } else {

    auto d = Disc();
    auto ncomp = m_u.nprop();

    // Field output

    std::vector< std::string > nodefieldnames
      {"density", "xvelocity", "yvelocity", "zvelocity", "energy", "pressure"};

    using tk::operator/=;
    auto r = m_u.extract( 0, 0 );
    auto u = m_u.extract( 1, 0 );  u /= r;
    auto v = m_u.extract( 2, 0 );  v /= r;
    auto w = m_u.extract( 3, 0 );  w /= r;
    auto e = m_u.extract( 4, 0 );  e /= r;
    std::vector< tk::real > p( m_u.nunk() );
    for (std::size_t i=0; i<p.size(); ++i) {
      auto ei = e[i] - 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
      p[i] = eos::pressure( r[i], ei );
    }

    std::vector< std::vector< tk::real > > nodefields{
      std::move(r), std::move(u), std::move(v), std::move(w), std::move(e),
      std::move(p) };

    for (std::size_t c=0; c<ncomp-5; ++c) {
      nodefieldnames.push_back( "c" + std::to_string(c) );
      nodefields.push_back( m_u.extract( 5+c, 0 ) );
    }

    Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

    // Surface output

    std::vector< std::string > nodesurfnames
      {"density", "xvelocity", "yvelocity", "zvelocity", "energy", "pressure"};

    for (std::size_t c=1; c<ncomp-5; ++c) {
      nodesurfnames.push_back( "c" + std::to_string(c) );
    }

    std::vector< std::vector< tk::real > > nodesurfs;

    const auto& lid = d->Lid();
    auto bnode = tk::bfacenodes( m_bface, m_triinpoel );
    for (auto sideset : g_inputdeck.outsets()) {
      auto b = bnode.find(sideset);
      if (b == end(bnode)) continue;
      const auto& nodes = b->second;
      auto i = nodesurfs.size();
      nodesurfs.insert( end(nodesurfs), ncomp + 1,
                        std::vector< tk::real >( nodes.size() ) );
      std::size_t j = 0;
      for (auto n : nodes) {
        const auto s = m_u[n];
        nodesurfs[i+0][j] = s[0];
        nodesurfs[i+1][j] = s[1]/s[0];
        nodesurfs[i+2][j] = s[2]/s[0];
        nodesurfs[i+3][j] = s[3]/s[0];
        nodesurfs[i+4][j] = s[4]/s[0];
        auto ei = s[4]/s[0] - 0.5*(s[1]*s[1] + s[2]*s[2] + s[3]*s[3])/s[0]/s[0];
        nodesurfs[i+5][j] = eos::pressure( s[0], ei );
        for (std::size_t c=0; c<ncomp-5; ++c) nodesurfs[i+1+c][j] = s[5+c];
        ++j;
      } 
    }   

    // Send mesh and fields data (solution dump) for output to file
    d->write( d->Inpoel(), d->Coord(), m_bface, tk::remap(m_bnode,lid),
              m_triinpoel, {}, nodefieldnames, nodesurfnames, {}, nodefields,
              nodesurfs, cb );

  }
}

void
AirCG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Time history

  if (d->histiter() or d->histtime() or d->histrange()) {

    const auto& inpoel = d->Inpoel();
    std::vector< std::vector< tk::real > > hist( d->Hist().size() );
    std::size_t j = 0;
    for (const auto& p : d->Hist()) {
      auto e = p.get< tag::elem >();        // host element id
      const auto& n = p.get< tag::fn >();   // shapefunctions evaluated at point
      hist[j].resize( 6, 0.0 );
      for (std::size_t i=0; i<4; ++i) {
        const auto u = m_u[ inpoel[e*4+i] ];
        hist[j][0] += n[i] * u[0];
        hist[j][1] += n[i] * u[1]/u[0];
        hist[j][2] += n[i] * u[2]/u[0];
        hist[j][3] += n[i] * u[3]/u[0];
        hist[j][4] += n[i] * u[4]/u[0];
        auto ei = u[4]/u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0]/u[0];
        hist[j][5] += n[i] * eos::pressure( u[0], ei );
      }
      ++j;
    }
    d->history( std::move(hist) );

  }

  // Field data
  if (d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished)
    writeFields( CkCallback(CkIndex_AirCG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
AirCG::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  // Increment Runge-Kutta stage counter
  ++m_stage;

  // If not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise output field data to file(s)
  if (m_stage < 3) grad(); else out();
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
