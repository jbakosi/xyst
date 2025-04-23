// *****************************************************************************
/*!
  \file      src/Inciter/KozCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     KozCG: Taylor-Galerkin, FCT, element-based continuous Galerkin
*/
// *****************************************************************************

#include "XystBuildConfig.hpp"
#include "KozCG.hpp"
#include "Vector.hpp"
#include "Reader.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "ExodusIIMeshWriter.hpp"
#include "InciterConfig.hpp"
#include "DerivedData.hpp"
#include "Discretization.hpp"
#include "DiagReducer.hpp"
#include "IntegralReducer.hpp"
#include "Integrals.hpp"
#include "Refiner.hpp"
#include "Reorder.hpp"
#include "Around.hpp"
#include "Kozak.hpp"
#include "Problems.hpp"
#include "EOS.hpp"
#include "BC.hpp"

namespace inciter {

extern ctr::Config g_cfg;

static CkReduction::reducerType IntegralsMerger;

} // inciter::

using inciter::g_cfg;
using inciter::KozCG;

KozCG::KozCG( const CProxy_Discretization& disc,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_nrhs( 0 ),
  m_nnorm( 0 ),
  m_naec( 0 ),
  m_nalw( 0 ),
  m_nlim( 0 ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_u( Disc()->Gid().size(), g_cfg.get< tag::problem_ncomp >() ),
  m_p( m_u.nunk(), m_u.nprop()*2 ),
  m_q( m_u.nunk(), m_u.nprop()*2 ),
  m_a( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_cfg.get< tag::t0 >() ),
  m_finished( 0 ),
  m_freezeflow( 1.0 )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side sets used in the input file
//! \param[in] bnode Boundary-node lists mapped to side sets used in input file
//! \param[in] triinpoel Boundary-face connectivity where BCs set (global ids)
// *****************************************************************************
{
  usesAtSync = true;    // enable migration at AtSync

  auto d = Disc();

  // Compute total box IC volume
  d->boxvol();

  // Activate SDAG wait for initially computing integrals
  thisProxy[ thisIndex ].wait4int();
}

void
KozCG::setupBC()
// *****************************************************************************
// Prepare boundary condition data structures
// *****************************************************************************
{
  // Query Dirichlet BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > dir;
  for (const auto& s : g_cfg.get< tag::bc_dir >()) {
    auto k = m_bface.find(s[0]);
    if (k != end(m_bface)) {
      auto& n = dir[ k->first ];
      for (auto f : k->second) {
        const auto t = m_triinpoel.data() + f*3;
        n.insert( t[0] );
        n.insert( t[1] );
        n.insert( t[2] );
      }
    }
  }

  // Augment Dirichlet BC nodes with nodes not necessarily part of faces
  const auto& lid = Disc()->Lid();
  for (const auto& s : g_cfg.get< tag::bc_dir >()) {
    auto k = m_bnode.find(s[0]);
    if (k != end(m_bnode)) {
      auto& n = dir[ k->first ];
      for (auto g : k->second) {
        n.insert( tk::cref_find(lid,g) );
      }
    }
  }

  // Collect unique set of nodes + Dirichlet BC components mask
  auto ncomp = m_u.nprop();
  auto nmask = ncomp + 1;
  const auto& dbc = g_cfg.get< tag::bc_dir >();
  std::unordered_map< std::size_t, std::vector< int > > dirbcset;
  for (const auto& mask : dbc) {
    ErrChk( mask.size() == nmask, "Incorrect Dirichlet BC mask ncomp" );
    auto n = dir.find( mask[0] );
    if (n != end(dir)) {
      for (auto p : n->second) {
        auto& m = dirbcset[p];
        if (m.empty()) m.resize( ncomp, 0 );
        for (std::size_t c=0; c<ncomp; ++c) {
          if (!m[c]) m[c] = mask[c+1];  // overwrite mask if 0 -> 1
        }
      }
    }
  }
  // Compile streamable list of nodes + Dirichlet BC components mask
  tk::destroy( m_dirbcmasks );
  for (const auto& [p,mask] : dirbcset) {
    m_dirbcmasks.push_back( p );
    m_dirbcmasks.insert( end(m_dirbcmasks), begin(mask), end(mask) );
  }
  ErrChk( m_dirbcmasks.size() % nmask == 0, "Dirichlet BC masks incomplete" );

  // Query pressure BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > pre;
  for (const auto& ss : g_cfg.get< tag::bc_pre >()) {
    for (const auto& s : ss) {
      auto k = m_bface.find(s);
      if (k != end(m_bface)) {
        auto& n = pre[ k->first ];
        for (auto f : k->second) {
          const auto t = m_triinpoel.data() + f*3;
          n.insert( t[0] );
          n.insert( t[1] );
          n.insert( t[2] );
        }
      }
    }
  }

  // Augment Pressure BC nodes with nodes not necessarily part of faces
  for (const auto& s : g_cfg.get< tag::bc_pre >()) {
    auto k = m_bnode.find(s[0]);
    if (k != end(m_bnode)) {
      auto& n = pre[ k->first ];
      for (auto g : k->second) {
        n.insert( tk::cref_find(lid,g) );
      }
    }
  }

  // Prepare density and pressure values for pressure BC nodes
  const auto& pbc_set = g_cfg.get< tag::bc_pre >();
  if (!pbc_set.empty()) {
    const auto& pbc_r = g_cfg.get< tag::bc_pre_density >();
    ErrChk( pbc_r.size() == pbc_set.size(), "Pressure BC density unspecified" );
    const auto& pbc_p = g_cfg.get< tag::bc_pre_pressure >();
    ErrChk( pbc_p.size() == pbc_set.size(), "Pressure BC pressure unspecified" );
    tk::destroy( m_prebcnodes );
    tk::destroy( m_prebcvals );
    for (const auto& [s,n] : pre) {
      m_prebcnodes.insert( end(m_prebcnodes), begin(n), end(n) );
      for (std::size_t p=0; p<pbc_set.size(); ++p) {
        for (auto u : pbc_set[p]) {
          if (s == u) {
            for (std::size_t i=0; i<n.size(); ++i) {
              m_prebcvals.push_back( pbc_r[p] );
              m_prebcvals.push_back( pbc_p[p] );
            }
          }
        }
      }
    }
    ErrChk( m_prebcnodes.size()*2 == m_prebcvals.size(),
            "Pressure BC data incomplete" );
  }

  // Query symmetry BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > sym;
  for (auto s : g_cfg.get< tag::bc_sym >()) {
    auto k = m_bface.find(s);
    if (k != end(m_bface)) {
      auto& n = sym[ k->first ];
      for (auto f : k->second) {
        const auto t = m_triinpoel.data() + f*3;
        n.insert( t[0] );
        n.insert( t[1] );
        n.insert( t[2] );
      }
    }
  }

  // Query farfield BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > far;
  for (auto s : g_cfg.get< tag::bc_far >()) {
    auto k = m_bface.find(s);
    if (k != end(m_bface)) {
      auto& n = far[ k->first ];
      for (auto f : k->second) {
        const auto t = m_triinpoel.data() + f*3;
        n.insert( t[0] );
        n.insert( t[1] );
        n.insert( t[2] );
      }
    }
  }

  // Generate unique set of symmetry BC nodes
  tk::destroy( m_symbcnodeset );
  for (const auto& [s,n] : sym) m_symbcnodeset.insert( begin(n), end(n) );
  // Generate unique set of farfield BC nodes
  tk::destroy( m_farbcnodeset );
  for (const auto& [s,n] : far) m_farbcnodeset.insert( begin(n), end(n) );

  // If farfield BC is set on a node, will not also set symmetry BC
  for (auto i : m_farbcnodeset) m_symbcnodeset.erase(i);
  // If pressure BC is set on a node, will not also set symmetry BC
  for (auto i : m_prebcnodes) m_symbcnodeset.erase(i);
}

void
KozCG::feop()
// *****************************************************************************
// Start (re-)computing finite element domain and boundary operators
// *****************************************************************************
{
  auto d = Disc();

  // Prepare boundary conditions data structures
  setupBC();

  // Compute local contributions to boundary normals and integrals
  bndint();

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
KozCG::bndint()
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
  auto invdistsq = [&]( const tk::real c[], std::size_t p ){
    return 1.0 / ( (c[0] - x[p]) * (c[0] - x[p]) +
                   (c[1] - y[p]) * (c[1] - y[p]) +
                   (c[2] - z[p]) * (c[2] - z[p]) );
  };

  tk::destroy( m_bnorm );
  tk::destroy( m_bndpoinint );

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
       const tk::real centroid[3] = {
         (x[N[0]] + x[N[1]] + x[N[2]]) / 3.0,
         (y[N[0]] + y[N[1]] + y[N[2]]) / 3.0,
         (z[N[0]] + z[N[1]] + z[N[2]]) / 3.0 };

       // contribute all edges of triangle
       for (const auto& [i,j] : tk::lpoet) {
         auto p = N[i];
         tk::real r = invdistsq( centroid, p );
         auto& v = m_bnorm[setid];      // associate side set id
         auto& bn = v[gid[p]];          // associate global node id of bnd pnt
         bn[0] += r * n[0];             // inv.dist.sq-weighted normal
         bn[1] += r * n[1];
         bn[2] += r * n[2];
         bn[3] += r;                    // inv.dist.sq of node from centroid
         auto& b = m_bndpoinint[gid[p]];// assoc global id of bnd point
         b[0] += n[0] * A / 3.0;        // bnd-point integral
         b[1] += n[1] * A / 3.0;
         b[2] += n[2] * A / 3.0;
       }
    }
  }
}

void
KozCG::comnorm( const decltype(m_bnorm)& inbnd )
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
KozCG::registerReducers()
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
  IntegralsMerger = CkReduction::addReducer( integrals::mergeIntegrals );
}

void
// cppcheck-suppress unusedFunction
KozCG::ResumeFromSync()
// *****************************************************************************
//  Return from migration
//! \details This is called when load balancing (LB) completes. The presence of
//!   this function does not affect whether or not we block on LB.
// *****************************************************************************
{
  if (Disc()->It() == 0) Throw( "it = 0 in ResumeFromSync()" );

  if (!g_cfg.get< tag::nonblocking >()) dt();
}

void
KozCG::setup( tk::real v )
// *****************************************************************************
// Start setup for solution
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions
  problems::initialize( d->Coord(), m_u, d->T(), d->BoxNodes() );

  // Query time history field output labels from all PDEs integrated
  if (!g_cfg.get< tag::histout >().empty()) {
    std::vector< std::string > var
      {"density", "xvelocity", "yvelocity", "zvelocity", "energy", "pressure"};
    auto ncomp = m_u.nprop();
    for (std::size_t c=5; c<ncomp; ++c)
      var.push_back( "c" + std::to_string(c-5) );
    d->histheader( std::move(var) );
  }

  // Compute finite element operators
  feop();
}

void
KozCG::start()
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
  dt();
}

void
KozCG::bnorm()
// *****************************************************************************
// Combine own and communicated portions of the boundary point normals
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
}

void
KozCG::streamable()
// *****************************************************************************
// Convert integrals into streamable data structures
// *****************************************************************************
{
  // Query surface integral output nodes
  std::unordered_map< int, std::vector< std::size_t > > surfintnodes;
  const auto& is = g_cfg.get< tag::integout >();
  std::set< int > outsets( begin(is), end(is) );
  for (auto s : outsets) {
    auto m = m_bface.find(s);
    if (m != end(m_bface)) {
      auto& n = surfintnodes[ m->first ];       // associate set id
      for (auto f : m->second) {                // face ids on side set
        n.push_back( m_triinpoel[f*3+0] );      // nodes on side set
        n.push_back( m_triinpoel[f*3+1] );
        n.push_back( m_triinpoel[f*3+2] );
      }
    }
  }
  for (auto& [s,n] : surfintnodes) tk::unique( n );
  // Prepare surface integral data
  tk::destroy( m_surfint );
  const auto& gid = Disc()->Gid();
  for (auto&& [s,n] : surfintnodes) {
    auto& sint = m_surfint[s];  // associate set id
    auto& nodes = sint.first;
    auto& ndA = sint.second;
    nodes = std::move(n);
    ndA.resize( nodes.size()*3 );
    std::size_t a = 0;
    for (auto p : nodes) {
      const auto& b = tk::cref_find( m_bndpoinint, gid[p] );
      ndA[a*3+0] = b[0];        // store ni * dA
      ndA[a*3+1] = b[1];
      ndA[a*3+2] = b[2];
      ++a;
    }
  }

  // Convert symmetry BC data to streamable data structures
  tk::destroy( m_symbcnodes );
  tk::destroy( m_symbcnorms );
  for (auto p : m_symbcnodeset) {
    for (const auto& s : g_cfg.get< tag::bc_sym >()) {
      auto m = m_bnorm.find(s);
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
  tk::destroy( m_farbcnodes );
  tk::destroy( m_farbcnorms );
  for (auto p : m_farbcnodeset) {
    for (const auto& s : g_cfg.get< tag::bc_far >()) {
      auto n = m_bnorm.find(s);
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
  tk::destroy( m_bnorm );
}


void
// cppcheck-suppress unusedFunction
KozCG::merge()
// *****************************************************************************
// Combine own and communicated portions of the integrals
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions to boundary point normals
  bnorm();

  // Convert integrals into streamable data structures
  streamable();

  // Enforce boundary conditions using (re-)computed boundary data
  BC( m_u, d->T() );

  if (d->Initial()) {
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_KozCG::start(), thisProxy[thisIndex]) );
  } else {
    feop_complete();
  }
}

void
KozCG::BC( tk::Fields& u, tk::real t )
// *****************************************************************************
// Apply boundary conditions
//! \param[in,out] u Solution to apply BCs to
//! \param[in] t Physical time
// *****************************************************************************
{
  auto d = Disc();

  // Apply Dirichlet BCs
  physics::dirbc( u, t, d->Coord(), d->BoxNodes(), m_dirbcmasks );

  // Apply symmetry BCs
  physics::symbc( u, m_symbcnodes, m_symbcnorms, /*pos=*/1 );

  // Apply farfield BCs
  physics::farbc( u, m_farbcnodes, m_farbcnorms );

  // Apply pressure BCs
  physics::prebc( u, m_prebcnodes, m_prebcvals );
}

void
KozCG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_cfg.get< tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto d = Disc();

  // use constant dt if configured
  if (std::abs(const_dt) > eps) {

    // cppcheck-suppress redundantInitialization
    mindt = const_dt;

  } else {

    const auto& vol = d->Vol();
    auto cfl = g_cfg.get< tag::cfl >();

    if (g_cfg.get< tag::steady >()) {

      for (std::size_t p=0; p<m_u.nunk(); ++p) {
        auto r = m_u(p,0);
        auto u = m_u(p,1)/r;
        auto v = m_u(p,2)/r;
        auto w = m_u(p,3)/r;
        auto pr = eos::pressure( m_u(p,4) - 0.5*r*(u*u + v*v + w*w) );
        auto c = eos::soundspeed( r, std::max(pr,0.0) );
        auto L = std::cbrt( vol[p] );
        auto vel = std::sqrt( u*u + v*v + w*w );
        m_dtp[p] = L / std::max( vel+c, 1.0e-8 ) * cfl;
      }
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );

    } else {

      for (std::size_t p=0; p<m_u.nunk(); ++p) {
        auto r = m_u(p,0);
        auto u = m_u(p,1)/r;
        auto v = m_u(p,2)/r;
        auto w = m_u(p,3)/r;
        auto pr = eos::pressure( m_u(p,4) - 0.5*r*(u*u + v*v + w*w) );
        auto c = eos::soundspeed( r, std::max(pr,0.0) );
        auto L = std::cbrt( vol[p] );
        auto vel = std::sqrt( u*u + v*v + w*w );
        auto euler_dt = L / std::max( vel+c, 1.0e-8 );
        mindt = std::min( mindt, euler_dt );
      }
      mindt *= cfl;

    }

    // Freeze flow if configured and apply multiplier on scalar(s)
    if (d->T() > g_cfg.get< tag::freezetime >()) {
      m_freezeflow = g_cfg.get< tag::freezeflow >();
    }

    mindt *= m_freezeflow;

  }

  // Actiavate SDAG waits for next time step
  thisProxy[ thisIndex ].wait4rhs();
  thisProxy[ thisIndex ].wait4aec();
  thisProxy[ thisIndex ].wait4alw();
  thisProxy[ thisIndex ].wait4sol();
  thisProxy[ thisIndex ].wait4step();

  // Contribute to minimum dt across all chares and advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(KozCG,advance), thisProxy) );
}

void
KozCG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
// *****************************************************************************
{
  // Detect blowup
  auto eps = std::numeric_limits< tk::real >::epsilon();
  if (newdt < eps) m_finished = 1;

  // Set new time step size
  Disc()->setdt( newdt );

  // Compute rhs
  rhs();
}

void
KozCG::rhs()
// *****************************************************************************
// Compute right-hand side
// *****************************************************************************
{
  auto d = Disc();

  // Compute own portion of right-hand side for all equations
  kozak::rhs( d->Inpoel(), d->Coord(), d->T(), d->Dt(), m_tp, m_dtp, m_u,
              m_rhs );

  // Communicate rhs to other chares on chare-boundary
  if (d->NodeCommMap().empty()) {
    comrhs_complete();
  } else {
    const auto& lid = d->Lid();
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_rhsc) exp;
      for (auto g : n) exp[g] = m_rhs[ tk::cref_find(lid,g) ];
      thisProxy[c].comrhs( exp );
    }
  }
  ownrhs_complete();
}

void
KozCG::comrhs(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& inrhs )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] inrhs Partial contributions of RHS to chare-boundary nodes. Key: 
//!   global mesh node IDs, value: contributions for all scalar components.
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
KozCG::fct()
// *****************************************************************************
// Continue with flux-corrected transport if enabled
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

  // Combine own and communicated contributions to rhs
  for (const auto& [g,r] : m_rhsc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_rhs(i,c) += r[c];
  }
  tk::destroy(m_rhsc);

  if (g_cfg.get< tag::fct >()) aec(); else solve();
}

void
// cppcheck-suppress unusedFunction
KozCG::aec()
// *****************************************************************************
// Compute antidiffusive contributions: P+/-
// *****************************************************************************
{
  auto d = Disc();
  const auto ncomp = m_u.nprop();
  const auto& lid = d->Lid();

  // Antidiffusive contributions: P+/-

  auto ctau = g_cfg.get< tag::fctdif >();
  m_p.fill( 0.0 );

  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );
    for (std::size_t c=0; c<ncomp; ++c) {
      auto p = c*2;
      auto n = p+1;
      tk::real aec[4] = { 0.0, 0.0, 0.0, 0.0 };
      for (std::size_t a=0; a<4; ++a) {
        for (std::size_t b=0; b<4; ++b) {
          auto m = J/120.0 * ((a == b) ? 3.0 : -1.0);
          aec[a] += m * ctau * m_u(N[b],c);
        }
        m_p(N[a],p) += std::max(0.0,aec[a]);
        m_p(N[a],n) += std::min(0.0,aec[a]);
      }
    }
  }

  // Apply symmetry BCs on AEC
  for (std::size_t i=0; i<m_symbcnodes.size(); ++i) {
    auto p = m_symbcnodes[i];
    auto nx = m_symbcnorms[i*3+0];
    auto ny = m_symbcnorms[i*3+1];
    auto nz = m_symbcnorms[i*3+2];
    auto rvnp = m_p(p,2)*nx + m_p(p,4)*ny + m_p(p,6)*nz;
    auto rvnn = m_p(p,3)*nx + m_p(p,5)*ny + m_p(p,7)*nz;
    m_p(p,2) -= rvnp * nx;
    m_p(p,3) -= rvnn * nx;
    m_p(p,4) -= rvnp * ny;
    m_p(p,5) -= rvnn * ny;
    m_p(p,6) -= rvnp * nz;
    m_p(p,7) -= rvnn * nz;
  }

  // Communicate antidiffusive edge and low-order solution contributions
  if (d->NodeCommMap().empty()) {
    comaec_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_pc) exp;
      for (auto g : n) exp[g] = m_p[ tk::cref_find(lid,g) ];
      thisProxy[c].comaec( exp );
    }
  }
  ownaec_complete();
}

void
KozCG::comaec( const std::unordered_map< std::size_t,
                       std::vector< tk::real > >& inaec )
// *****************************************************************************
//  Receive antidiffusive and low-order contributions on chare-boundaries
//! \param[in] inaec Partial contributions of antidiffusive edge and low-order
//!   solution contributions on chare-boundary nodes. Key: global mesh node IDs,
//!   value: 0: antidiffusive contributions, 1: low-order solution.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,a] : inaec) m_pc[g] += a;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_naec == Disc()->NodeCommMap().size()) {
    m_naec = 0;
    comaec_complete();
  }
}

void
KozCG::alw()
// *****************************************************************************
// Compute allowed limits, Q+/-
// *****************************************************************************
{
  auto d = Disc();
  const auto steady = g_cfg.get< tag::steady >();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& lid = d->Lid();
  const auto& vol = d->Vol();
  const auto& inpoel = d->Inpoel();

  // Combine own and communicated contributions to antidiffusive contributions
  // and low-order solution
  for (const auto& [g,p] : m_pc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<p.size(); ++c) m_p(i,c) += p[c];
  }
  tk::destroy(m_pc);

  // Finish computing antidiffusive contributions and low-order solution
  auto dt = d->Dt();
  for (std::size_t i=0; i<npoin; ++i) {
    if (steady) dt = m_dtp[i];
    for (std::size_t c=0; c<ncomp; ++c) {
      auto p = c*2;
      auto n = p+1;
      m_p(i,p) /= vol[i];
      m_p(i,n) /= vol[i];
      // low-order solution
      m_rhs(i,c) = m_u(i,c) + dt*m_rhs(i,c)/vol[i] - m_p(i,p) - m_p(i,n);
    }
  }

  // Allowed limits: Q+/-

  using std::max;
  using std::min;

  auto large = std::numeric_limits< tk::real >::max();
  for (std::size_t i=0; i<m_q.nunk(); ++i) {
    for (std::size_t c=0; c<m_q.nprop()/2; ++c) {
      m_q(i,c*2+0) = -large;
      m_q(i,c*2+1) = +large;
    }
  }

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto N = inpoel.data() + e*4;
    for (std::size_t c=0; c<ncomp; ++c) {
      auto alwp = -large;
      auto alwn = +large;
      for (std::size_t a=0; a<4; ++a) {
        if (g_cfg.get< tag::fctclip >()) {
          alwp = max( alwp, m_rhs(N[a],c) );
          alwn = min( alwn, m_rhs(N[a],c) );
        } else {
          alwp = max( alwp, max(m_rhs(N[a],c), m_u(N[a],c)) );
          alwn = min( alwn, min(m_rhs(N[a],c), m_u(N[a],c)) );
        }
      }
      auto p = c*2;
      auto n = p+1;
      for (std::size_t a=0; a<4; ++a) {
        m_q(N[a],p) = max(m_q(N[a],p), alwp);
        m_q(N[a],n) = min(m_q(N[a],n), alwn);
      }
    }
  }

  // Communicate allowed limits contributions
  if (d->NodeCommMap().empty()) {
    comalw_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_qc) exp;
      for (auto g : n) exp[g] = m_q[ tk::cref_find(lid,g) ];
      thisProxy[c].comalw( exp );
    }
  }
  ownalw_complete();
}

void
KozCG::comalw( const std::unordered_map< std::size_t,
                       std::vector< tk::real > >& inalw )
// *****************************************************************************
//  Receive allowed limits contributions on chare-boundaries
//! \param[in] inalw Partial contributions of allowed limits contributions on
//!   chare-boundary nodes. Key: global mesh node IDs, value: allowed limit
//!   contributions.
// *****************************************************************************
{
  for (const auto& [g,alw] : inalw) {
    auto& q = m_qc[g];
    q.resize( alw.size() );
    for (std::size_t c=0; c<alw.size()/2; ++c) {
      auto p = c*2;
      auto n = p+1;
      q[p] = std::max( q[p], alw[p] );
      q[n] = std::min( q[n], alw[n] );
    }
  }

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nalw == Disc()->NodeCommMap().size()) {
    m_nalw = 0;
    comalw_complete();
  }
}

void
KozCG::lim()
// *****************************************************************************
// Compute limit coefficients
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& lid = d->Lid();

  using std::max;
  using std::min;

  // Combine own and communicated contributions to allowed limits
  for (const auto& [g,alw] : m_qc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<alw.size()/2; ++c) {
      auto p = c*2;
      auto n = p+1;
      m_q(i,p) = max( m_q(i,p), alw[p] );
      m_q(i,n) = min( m_q(i,n), alw[n] );
    }
  }
  tk::destroy(m_qc);

  // Finish computing allowed limits
  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      auto p = c*2;
      auto n = p+1;
      m_q(i,p) -= m_rhs(i,c);
      m_q(i,n) -= m_rhs(i,c);
    }
  }

  // Limit coefficients, C

  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      auto p = c*2;
      auto n = p+1;
      auto eps = std::numeric_limits< tk::real >::epsilon();
      m_q(i,p) = m_p(i,p) <  eps ? 0.0 : min(1.0, m_q(i,p)/m_p(i,p));
      m_q(i,n) = m_p(i,n) > -eps ? 0.0 : min(1.0, m_q(i,n)/m_p(i,n));
    }
  }

  // Limited antidiffusive contributions

  auto ctau = g_cfg.get< tag::fctdif >();
  m_a.fill( 0.0 );

  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto fctsys = g_cfg.get< tag::fctsys >();
  for (auto& c : fctsys) --c;

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );
    tk::real coef[ncomp], aec[ncomp][4];
    for (std::size_t c=0; c<ncomp; ++c) {
      auto p = c*2;
      auto n = p+1;
      coef[c] = 1.0;
      for (std::size_t a=0; a<4; ++a) {
        aec[c][a] = 0.0;
        for (std::size_t b=0; b<4; ++b) {
          auto m = J/120.0 * ((a == b) ? 3.0 : -1.0);
          aec[c][a] += m * ctau * m_u(N[b],c);
        }
        coef[c] = min(coef[c], aec[c][a] > 0.0 ? m_q(N[a],p) : m_q(N[a],n));
      }
    }
    tk::real cs = 1.0;
    for (auto c : fctsys) cs = min( cs, coef[c] );
    for (auto c : fctsys) coef[c] = cs;
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t a=0; a<4; ++a) {
        m_a(N[a],c) += coef[c] * aec[c][a];
      }
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif

  // Communicate limited antidiffusive contributions
  if (d->NodeCommMap().empty()) {
    comlim_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_ac) exp;
      for (auto g : n) exp[g] = m_a[ tk::cref_find(lid,g) ];
      thisProxy[c].comlim( exp );
    }
  }
  ownlim_complete();
}

void
KozCG::comlim( const std::unordered_map< std::size_t,
                       std::vector< tk::real > >& inlim )
// *****************************************************************************
//  Receive limited antidiffusive contributions on chare-boundaries
//! \param[in] inlim Partial contributions of limited contributions on
//!   chare-boundary nodes. Key: global mesh node IDs, value: limited
//!   contributions.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,a] : inlim) m_ac[g] += a;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nlim == Disc()->NodeCommMap().size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
KozCG::solve()
// *****************************************************************************
// Compute limit coefficients
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& lid = d->Lid();
  const auto& vol = d->Vol();
  const auto steady = g_cfg.get< tag::steady >();

  // Combine own and communicated contributions to limited antidiffusive
  // contributions
  for (const auto& [g,a] : m_ac) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<a.size(); ++c) m_a(i,c) += a[c];
  }
  tk::destroy(m_ac);

  tk::Fields u;
  std::size_t cstart = m_freezeflow > 1.0 ? 5 : 0;
  if (cstart) u = m_u;

  if (g_cfg.get< tag::fct >()) {
    // Apply limited antidiffusive contributions to low-order solution
    for (std::size_t i=0; i<npoin; ++i) {
      for (std::size_t c=0; c<ncomp; ++c) {
        m_a(i,c) = m_rhs(i,c) + m_a(i,c)/vol[i];
      }
    }
  } else {
    // Apply rhs
    auto dt = d->Dt();
    for (std::size_t i=0; i<npoin; ++i) {
      if (steady) dt = m_dtp[i];
      for (std::size_t c=0; c<ncomp; ++c) {
        m_a(i,c) = m_u(i,c) + dt*m_rhs(i,c)/vol[i];
      }
    }
  }

  // Apply scalar source to solution (if defined)
  auto src = problems::PHYS_SRC();
  if (src) src( d->Coord(), d->T(), m_a );

  // Enforce boundary conditions
  BC( m_a, d->T() + d->Dt() );

  // Explicitly zero out flow for freezeflow
  if (cstart) {
    for (std::size_t i=0; i<npoin; ++i) {
      for (std::size_t c=0; c<cstart; ++c) {
        m_a(i,c) = u(i,c);
      }
    }
  }

  // Compute diagnostics, e.g., residuals
  auto diag_iter = g_cfg.get< tag::diag_iter >();
  auto diag = m_diag.rhocompute( *d, m_a, m_u, diag_iter );

  // Update solution
  m_u = m_a;
  m_a.fill( 0.0 );

  // Increase number of iterations and physical time
  d->next();

  // Advance physical time for local time stepping
  if (steady) {
    using tk::operator+=;
    m_tp += m_dtp;
  }

  // Evaluate residuals
  if (!diag) evalres( std::vector< tk::real >( ncomp, 1.0 ) );
}

void
KozCG::evalres( const std::vector< tk::real >& l2res )
// *****************************************************************************
//  Evaluate residuals
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
// *****************************************************************************
{
  if (g_cfg.get< tag::steady >()) {
    const auto rc = g_cfg.get< tag::rescomp >() - 1;
    Disc()->residual( l2res[rc] );
  }

  refine();
}

void
KozCG::refine()
// *****************************************************************************
// Optionally refine/derefine mesh
// *****************************************************************************
{
  auto d = Disc();

  // See if this is the last time step
  if (d->finished()) m_finished = 1;

  auto dtref = g_cfg.get< tag::href_dt >();
  auto dtfreq = g_cfg.get< tag::href_dtfreq >();

  // if t>0 refinement enabled and we hit the frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    d->refined() = 1;
    d->startvol();
    d->Ref()->dtref( m_bface, m_bnode, m_triinpoel );

    // Activate SDAG waits for re-computing the integrals
    thisProxy[ thisIndex ].wait4int();

  } else {      // do not refine

    d->refined() = 0;
    feop_complete();
    resize_complete();

  }
}

void
KozCG::resizePostAMR(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
  const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
  const std::set< std::size_t >& removedNodes,
  const std::unordered_map< int, std::unordered_set< std::size_t > >&
    nodeCommMap,
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
  d->resizePostAMR( chunk, coord, nodeCommMap, removedNodes );

  Assert(coord[0].size() == m_u.nunk()-removedNodes.size()+addedNodes.size(),
    "Incorrect vector length post-AMR: expected length after resizing = " +
    std::to_string(coord[0].size()) + ", actual unknown vector length = " +
    std::to_string(m_u.nunk()-removedNodes.size()+addedNodes.size()));

  // Remove newly removed nodes from solution vectors
  m_u.rm( removedNodes );
  m_rhs.rm( removedNodes );

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  m_u.resize( npoin );
  m_rhs.resize( npoin );

  // Update solution on new mesh
  for (const auto& n : addedNodes)
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      Assert(n.first < m_u.nunk(), "Added node index out of bounds post-AMR");
      Assert(n.second[0] < m_u.nunk() && n.second[1] < m_u.nunk(),
        "Indices of parent-edge nodes out of bounds post-AMR");
      m_u(n.first,c) = (m_u(n.second[0],c) + m_u(n.second[1],c))/2.0;
    }

  // Update physical-boundary node-, face-, and element lists
  m_bnode = bnode;
  m_bface = bface;
  m_triinpoel = tk::remap( triinpoel, d->Lid() );

  auto meshid = d->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
              CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );
}

void
KozCG::writeFields( CkCallback cb )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] cb Function to continue with after the write
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) { cb.send(); return; }

  auto d = Disc();
  auto ncomp = m_u.nprop();

  // Field output

  std::vector< std::string > nodefieldnames
    {"density", "velocityx", "velocityy", "velocityz", "energy", "pressure"};
  if (g_cfg.get< tag::steady >()) nodefieldnames.push_back( "mach" );

  using tk::operator/=;
  auto r = m_u.extract(0);
  auto u = m_u.extract(1);  u /= r;
  auto v = m_u.extract(2);  v /= r;
  auto w = m_u.extract(3);  w /= r;
  auto e = m_u.extract(4);  e /= r;
  std::vector< tk::real > pr( m_u.nunk() ), ma;
  if (g_cfg.get< tag::steady >()) ma.resize( m_u.nunk() );
  for (std::size_t i=0; i<pr.size(); ++i) {
    auto vv = u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
    pr[i] = eos::pressure( r[i]*(e[i] - 0.5*vv) );
    if (g_cfg.get< tag::steady >()) {
      ma[i] = std::sqrt(vv) / eos::soundspeed( r[i], pr[i] );
    }
  }

  std::vector< std::vector< tk::real > > nodefields{
    std::move(r), std::move(u), std::move(v), std::move(w), std::move(e),
    std::move(pr) };
  if (g_cfg.get< tag::steady >()) nodefields.push_back( std::move(ma) );

  for (std::size_t c=0; c<ncomp-5; ++c) {
    nodefieldnames.push_back( "c" + std::to_string(c) );
    nodefields.push_back( m_u.extract(5+c) );
  }

  // query function to evaluate analytic solution (if defined)
  auto sol = problems::SOL();

  if (sol) {
    const auto& coord = d->Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    auto an = m_u;
    std::vector< tk::real > ap( m_u.nunk() );
    for (std::size_t i=0; i<an.nunk(); ++i) {
      auto s = sol( x[i], y[i], z[i], d->T() );
      s[1] /= s[0];
      s[2] /= s[0];
      s[3] /= s[0];
      s[4] /= s[0];
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
      s[4] -= 0.5*(s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
      ap[i] = eos::pressure( s[0]*s[4] );
    }
    for (std::size_t c=0; c<5; ++c) {
      nodefieldnames.push_back( nodefieldnames[c] + "_analytic" );
      nodefields.push_back( an.extract(c) );
    }
    nodefieldnames.push_back( nodefieldnames[5] + "_analytic" );
    nodefields.push_back( std::move(ap) );
    for (std::size_t c=0; c<ncomp-5; ++c) {
      nodefieldnames.push_back( nodefieldnames[6+c] + "_analytic" );
      nodefields.push_back( an.extract(5+c) );
    }
  }

  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Surface output

  std::vector< std::string > nodesurfnames;
  std::vector< std::vector< tk::real > > nodesurfs;

  const auto& f = g_cfg.get< tag::fieldout >();

  if (!f.empty()) {
    nodesurfnames.push_back( "density" );
    nodesurfnames.push_back( "velocityx" );
    nodesurfnames.push_back( "velocityy" );
    nodesurfnames.push_back( "velocityz" );
    nodesurfnames.push_back( "energy" );
    nodesurfnames.push_back( "pressure" );

    for (std::size_t c=0; c<ncomp-5; ++c) {
      nodesurfnames.push_back( "c" + std::to_string(c) );
    }
    if (g_cfg.get< tag::steady >()) {
      nodesurfnames.push_back( "mach" );
    }

    auto bnode = tk::bfacenodes( m_bface, m_triinpoel );
    std::set< int > outsets( begin(f), end(f) );
    for (auto sideset : outsets) {
      auto b = bnode.find(sideset);
      if (b == end(bnode)) continue;
      const auto& nodes = b->second;
      auto i = nodesurfs.size();
      auto ns = ncomp + 1;
      if (g_cfg.get< tag::steady >()) ++ns;
      nodesurfs.insert( end(nodesurfs), ns,
                        std::vector< tk::real >( nodes.size() ) );
      std::size_t j = 0;
      for (auto n : nodes) {
        const auto s = m_u[n];
        std::size_t p = 0;
        nodesurfs[i+(p++)][j] = s[0];
        nodesurfs[i+(p++)][j] = s[1]/s[0];
        nodesurfs[i+(p++)][j] = s[2]/s[0];
        nodesurfs[i+(p++)][j] = s[3]/s[0];
        nodesurfs[i+(p++)][j] = s[4]/s[0];
        auto vv = (s[1]*s[1] + s[2]*s[2] + s[3]*s[3])/s[0]/s[0];
        auto ei = s[4]/s[0] - 0.5*vv;
        auto sp = eos::pressure( s[0]*ei );
        nodesurfs[i+(p++)][j] = sp;
        for (std::size_t c=0; c<ncomp-5; ++c) nodesurfs[i+(p++)+c][j] = s[5+c];
        if (g_cfg.get< tag::steady >()) {
          nodesurfs[i+(p++)][j] = std::sqrt(vv) / eos::soundspeed( s[0], sp );
        }
        ++j;
      }
    }
  }

  // Send mesh and fields data (solution dump) for output to file
  d->write( d->Inpoel(), d->Coord(), m_bface, tk::remap(m_bnode,d->Lid()),
            m_triinpoel, {}, nodefieldnames, {}, nodesurfnames,
            {}, nodefields, {}, nodesurfs, cb );
}

void
KozCG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    auto ncomp = m_u.nprop();
    const auto& inpoel = d->Inpoel();
    std::vector< std::vector< tk::real > > hist( d->Hist().size() );
    std::size_t j = 0;
    for (const auto& p : d->Hist()) {
      auto e = p.get< tag::elem >();        // host element id
      const auto& n = p.get< tag::fn >();   // shapefunctions evaluated at point
      hist[j].resize( ncomp+1, 0.0 );
      for (std::size_t i=0; i<4; ++i) {
        const auto u = m_u[ inpoel[e*4+i] ];
        hist[j][0] += n[i] * u[0];
        hist[j][1] += n[i] * u[1]/u[0];
        hist[j][2] += n[i] * u[2]/u[0];
        hist[j][3] += n[i] * u[3]/u[0];
        hist[j][4] += n[i] * u[4]/u[0];
        auto ei = u[4]/u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0]/u[0];
        hist[j][5] += n[i] * eos::pressure( u[0]*ei );
        for (std::size_t c=5; c<ncomp; ++c) hist[j][c+1] += n[i] * u[c];
      }
      ++j;
    }
    d->history( std::move(hist) );
  }

  // Field data
  if (d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished) {
    writeFields( CkCallback(CkIndex_KozCG::integrals(), thisProxy[thisIndex]) );
  } else {
    integrals();
  }
}

void
KozCG::integrals()
// *****************************************************************************
// Compute integral quantities for output
// *****************************************************************************
{
  auto d = Disc();

  if (d->integiter() or d->integtime() or d->integrange()) {

    using namespace integrals;
    std::vector< std::map< int, tk::real > > ints( NUMINT );

    // Prepend integral vector with metadata on the current time step:
    // current iteration count, current physical time, time step size
    ints[ ITER ][ 0 ] = static_cast< tk::real >( d->It() );
    ints[ TIME ][ 0 ] = d->T();
    ints[ DT ][ 0 ] = d->Dt();

    // Compute integrals requested for surfaces requested
    const auto& reqv = g_cfg.get< tag::integout_integrals >();
    std::unordered_set< std::string > req( begin(reqv), end(reqv) );
    if (req.count("mass_flow_rate")) {
      for (const auto& [s,sint] : m_surfint) {
        auto& mfr = ints[ MASS_FLOW_RATE ][ s ];
        const auto& nodes = sint.first;
        const auto& ndA = sint.second;
        auto n = ndA.data();
        for (auto p : nodes) {
          mfr += n[0]*m_u(p,1) + n[1]*m_u(p,2) + n[2]*m_u(p,3);
          n += 3;
        }
      }
    }

    auto stream = serialize( d->MeshId(), ints );
    d->contribute( stream.first, stream.second.get(), IntegralsMerger,
      CkCallback(CkIndex_Transporter::integrals(nullptr), d->Tr()) );

  } else {

    step();

  }
}

void
KozCG::evalLB( int nrestart )
// *****************************************************************************
// Evaluate whether to do load balancing
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  auto d = Disc();

  // Detect if just returned from a checkpoint and if so, zero timers and
  // finished flag
  if (d->restarted( nrestart )) m_finished = 0;

  // Load balancing if user frequency is reached or after the second time-step
  if (d->lb()) {

    AtSync();
    if (g_cfg.get< tag::nonblocking >()) dt();

  } else {

    dt();

  }
}

void
KozCG::evalRestart()
// *****************************************************************************
// Evaluate whether to save checkpoint/restart
// *****************************************************************************
{
  auto d = Disc();

  const auto rsfreq = g_cfg.get< tag::rsfreq >();
  const auto benchmark = g_cfg.get< tag::benchmark >();

  if ( !benchmark && (d->It()) % rsfreq == 0 ) {

    std::vector< std::size_t > meshdata{ /* finished = */ 0, d->MeshId() };
    contribute( meshdata, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,checkpoint), d->Tr()) );

  } else {

    evalLB( /* nrestart = */ -1 );

  }
}

void
KozCG::step()
// *****************************************************************************
// Evaluate whether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  if (thisIndex == 0) d->status();

  if (not m_finished) {

    evalRestart();

  } else {

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/kozcg.def.h"
