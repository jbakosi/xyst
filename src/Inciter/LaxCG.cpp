// *****************************************************************************
/*!
  \file      src/Inciter/LaxCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LaxCG: Time-derivative preconditioning for all Ma
*/
// *****************************************************************************

#include "XystBuildConfig.hpp"
#include "LaxCG.hpp"
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
#include "Lax.hpp"
#include "Dt.hpp"
#include "Problems.hpp"
#include "EOS.hpp"
#include "BC.hpp"

namespace inciter {

extern ctr::Config g_cfg;

static CkReduction::reducerType IntegralsMerger;

//! Runge-Kutta coefficients
static const std::array< tk::real, 3 > rkcoef{{ 1.0/3.0, 1.0/2.0, 1.0 }};

} // inciter::

using inciter::g_cfg;
using inciter::LaxCG;

LaxCG::LaxCG( const CProxy_Discretization& disc,
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
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_u( Disc()->Gid().size(), g_cfg.get< tag::problem_ncomp >() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_grad( m_u.nunk(), m_u.nprop()*3 ),
  m_stage( 0 ),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_cfg.get< tag::t0 >() ),
  m_finished( 0 )
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

  // Compute total box IC volume
  d->boxvol();

  // Activate SDAG wait for initially computing integrals
  thisProxy[ thisIndex ].wait4int();
}

void
LaxCG::setupBC()
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
        n.insert( m_triinpoel[f*3+0] );
        n.insert( m_triinpoel[f*3+1] );
        n.insert( m_triinpoel[f*3+2] );
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
    if (n != end(dir))
      for (auto p : n->second) {
        auto& m = dirbcset[p];
        if (m.empty()) m.resize( ncomp, 0 );
        for (std::size_t c=0; c<ncomp; ++c)
          if (!m[c]) m[c] = mask[c+1];  // overwrite mask if 0 -> 1
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
          n.insert( m_triinpoel[f*3+0] );
          n.insert( m_triinpoel[f*3+1] );
          n.insert( m_triinpoel[f*3+2] );
        }
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
        n.insert( m_triinpoel[f*3+0] );
        n.insert( m_triinpoel[f*3+1] );
        n.insert( m_triinpoel[f*3+2] );
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
        n.insert( m_triinpoel[f*3+0] );
        n.insert( m_triinpoel[f*3+1] );
        n.insert( m_triinpoel[f*3+2] );
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
}

void
LaxCG::feop()
// *****************************************************************************
// Start (re-)computing finite element domain and boundary operators
// *****************************************************************************
{
  auto d = Disc();

  // Prepare boundary conditions data structures
  setupBC();

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
LaxCG::bndint()
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
  tk::destroy( m_bndedgeint );

  for (const auto& [ setid, faceids ] : m_bface) { // for all side sets
    for (auto f : faceids) { // for all side set triangles
      const auto N = m_triinpoel.data() + f*3;
      const std::array< tk::real, 3 >
        ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
        ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
      auto n = tk::cross( ba, ca );
      auto A2 = tk::length( n );
      n[0] /= A2;
      n[1] /= A2;
      n[2] /= A2;
      const tk::real centroid[3] = {
        (x[N[0]] + x[N[1]] + x[N[2]]) / 3.0,
        (y[N[0]] + y[N[1]] + y[N[2]]) / 3.0,
        (z[N[0]] + z[N[1]] + z[N[2]]) / 3.0 };
      for (const auto& [i,j] : tk::lpoet) {
        auto p = N[i];
        auto q = N[j];
        tk::real r = invdistsq( centroid, p );
        auto& v = m_bnorm[setid];      // associate side set id
        auto& bpn = v[gid[p]];         // associate global node id of bnd pnt
        bpn[0] += r * n[0];            // inv.dist.sq-weighted normal
        bpn[1] += r * n[1];
        bpn[2] += r * n[2];
        bpn[3] += r;                   // inv.dist.sq of node from centroid
        auto& b = m_bndpoinint[gid[p]];// assoc global id of bnd point
        b[0] += n[0] * A2 / 6.0;        // bnd-point integral
        b[1] += n[1] * A2 / 6.0;
        b[2] += n[2] * A2 / 6.0;
        tk::UnsMesh::Edge ed{ gid[p], gid[q] };
        tk::real sig = 1.0;
        if (ed[0] > ed[1]) {
          std::swap( ed[0], ed[1] );
          sig = -1.0;
        }
        auto& e = m_bndedgeint[ ed ];
        e[0] += sig * n[0] * A2 / 24.0; // bnd-edge integral
        e[1] += sig * n[1] * A2 / 24.0;
        e[2] += sig * n[2] * A2 / 24.0;
      }
    }
  }
}

void
LaxCG::domint()
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
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::cross( ca, da );
    grad[2] = tk::cross( da, ba );
    grad[3] = tk::cross( ba, ca );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    for (const auto& [p,q] : tk::lpoed) {
      tk::UnsMesh::Edge ed{ gid[N[p]], gid[N[q]] };
      tk::real sig = 1.0;
      if (ed[0] > ed[1]) {
        std::swap( ed[0], ed[1] );
        sig = -1.0;
      }
      auto& n = m_domedgeint[ ed ];
      n[0] += sig * (grad[p][0] - grad[q][0]) / 48.0;
      n[1] += sig * (grad[p][1] - grad[q][1]) / 48.0;
      n[2] += sig * (grad[p][2] - grad[q][2]) / 48.0;
    }
  }
}

void
LaxCG::comnorm( const decltype(m_bnorm)& inbnd )
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
LaxCG::registerReducers()
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
LaxCG::ResumeFromSync()
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
LaxCG::setup( tk::real v )
// *****************************************************************************
// Start setup for solution
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  Disc()->Boxvol() = v;

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
LaxCG::start()
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
LaxCG::bnorm()
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
LaxCG::streamable()
// *****************************************************************************
// Convert integrals into streamable data structures
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();

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
  tk::destroy( m_bndpoinint );

  // Generate boundary superedges
  bndsuped();
  tk::destroy( m_bndedgeint );

  // Generate domain superedges
  domsuped();
  tk::destroy( m_domedgeint );

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
LaxCG::bndsuped()
// *****************************************************************************
// Generate superedge-groups for boundary-edge loops
//! \see See Lohner, Sec. 15.1.6.2, An Introduction to Applied CFD Techniques,
//!      Wiley, 2008.
// *****************************************************************************
{
  #ifndef NDEBUG
  auto nbedge = m_bndedgeint.size();
  #endif

  const auto& lid = Disc()->Lid();
  const auto& gid = Disc()->Gid();

  tk::destroy( m_bsupedge[0] );
  tk::destroy( m_bsupedge[1] );

  tk::destroy( m_bsupint[0] );
  tk::destroy( m_bsupint[1] );

  for (const auto& [setid, tri] : m_bface) {
    for (auto e : tri) {
      std::size_t N[3] = { m_triinpoel[e*3+0], m_triinpoel[e*3+1],
                           m_triinpoel[e*3+2] };
      int f = 0;
      tk::real sig[3];
      decltype(m_bndedgeint)::const_iterator b[3];
      for (const auto& [p,q] : tk::lpoet) {
        tk::UnsMesh::Edge ed{ gid[N[p]], gid[N[q]] };
        sig[f] = ed[0] < ed[1] ? 1.0 : -1.0;
        b[f] = m_bndedgeint.find( ed );
        if (b[f] == end(m_bndedgeint)) break; else ++f;
      }
      if (f == 3) {
        m_bsupedge[0].push_back( N[0] );
        m_bsupedge[0].push_back( N[1] );
        m_bsupedge[0].push_back( N[2] );
        m_bsupedge[0].push_back( m_symbcnodeset.count(N[0]) );
        m_bsupedge[0].push_back( m_symbcnodeset.count(N[1]) );
        m_bsupedge[0].push_back( m_symbcnodeset.count(N[2]) );
        for (int ed=0; ed<3; ++ed) {
          m_bsupint[0].push_back( sig[ed] * b[ed]->second[0] );
          m_bsupint[0].push_back( sig[ed] * b[ed]->second[1] );
          m_bsupint[0].push_back( sig[ed] * b[ed]->second[2] );
          m_bndedgeint.erase( b[ed] );
        }
      }
    }
  }

  m_bsupedge[1].resize( m_bndedgeint.size()*4 );
  m_bsupint[1].resize( m_bndedgeint.size()*3 );
  std::size_t k = 0;
  for (const auto& [ed,b] : m_bndedgeint) {
    auto p = tk::cref_find( lid, ed[0] );
    auto q = tk::cref_find( lid, ed[1] );
    auto e = m_bsupedge[1].data() + k*4;
    e[0] = p;
    e[1] = q;
    e[2] = m_symbcnodeset.count(p);
    e[3] = m_symbcnodeset.count(q);
    auto i = m_bsupint[1].data() + k*3;
    i[0] = b[0];
    i[1] = b[1];
    i[2] = b[2];
    ++k;
  }

  //std::cout << std::setprecision(2)
  //          << "superedges: ntri:" << m_bsupedge[0].size()/6
  //          << "(nedge:" << m_bsupedge[0].size()/3 << ","
  //          << 100.0 * static_cast< tk::real >( m_bsupedge[0].size()/3 ) /
  //                     static_cast< tk::real >( nbedge )
  //          << "%) + nedge:"
  //          << m_bsupedge[1].size()/4 << "("
  //          << 100.0 * static_cast< tk::real >( m_bsupedge[1].size()/4 ) /
  //                     static_cast< tk::real >( nbedge )
  //          << "%) = " << m_bsupedge[0].size()/2 + m_bsupedge[1].size()/4
  //          << " of "<< nbedge << " total boundary edges\n";

  Assert( m_bsupedge[0].size()/2 + m_bsupedge[1].size()/4 == nbedge,
          "Not all boundary edges accounted for in superedge groups" );
}

void
LaxCG::domsuped()
// *****************************************************************************
// Generate superedge-groups for domain-edge loops
//! \see See Lohner, Sec. 15.1.6.2, An Introduction to Applied CFD Techniques,
//!      Wiley, 2008.
// *****************************************************************************
{
  Assert( !m_domedgeint.empty(), "No domain edges to group" );

  #ifndef NDEBUG
  auto nedge = m_domedgeint.size();
  #endif

  const auto& inpoel = Disc()->Inpoel();
  const auto& lid = Disc()->Lid();
  const auto& gid = Disc()->Gid();

  tk::destroy( m_dsupedge[0] );
  tk::destroy( m_dsupedge[1] );
  tk::destroy( m_dsupedge[2] );

  tk::destroy( m_dsupint[0] );
  tk::destroy( m_dsupint[1] );
  tk::destroy( m_dsupint[2] );

  tk::UnsMesh::FaceSet untri;
  for (std::size_t e=0; e<inpoel.size()/4; e++) {
    std::size_t N[4] = {
      inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    for (const auto& [a,b,c] : tk::lpofa) untri.insert( { N[a], N[b], N[c] } );
  }

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    std::size_t N[4] = {
      inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    int f = 0;
    tk::real sig[6];
    decltype(m_domedgeint)::const_iterator d[6];
    for (const auto& [p,q] : tk::lpoed) {
      tk::UnsMesh::Edge ed{ gid[N[p]], gid[N[q]] };
      sig[f] = ed[0] < ed[1] ? 1.0 : -1.0;
      d[f] = m_domedgeint.find( ed );
      if (d[f] == end(m_domedgeint)) break; else ++f;
    }
    if (f == 6) {
      m_dsupedge[0].push_back( N[0] );
      m_dsupedge[0].push_back( N[1] );
      m_dsupedge[0].push_back( N[2] );
      m_dsupedge[0].push_back( N[3] );
      for (const auto& [a,b,c] : tk::lpofa) untri.erase( { N[a], N[b], N[c] } );
      for (int ed=0; ed<6; ++ed) {
        m_dsupint[0].push_back( sig[ed] * d[ed]->second[0] );
        m_dsupint[0].push_back( sig[ed] * d[ed]->second[1] );
        m_dsupint[0].push_back( sig[ed] * d[ed]->second[2] );
        m_domedgeint.erase( d[ed] );
      }
    }
  }

  for (const auto& N : untri) {
    int f = 0;
    tk::real sig[3];
    decltype(m_domedgeint)::const_iterator d[3];
    for (const auto& [p,q] : tk::lpoet) {
      tk::UnsMesh::Edge ed{ gid[N[p]], gid[N[q]] };
      sig[f] = ed[0] < ed[1] ? 1.0 : -1.0;
      d[f] = m_domedgeint.find( ed );
      if (d[f] == end(m_domedgeint)) break; else ++f;
    }
    if (f == 3) {
      m_dsupedge[1].push_back( N[0] );
      m_dsupedge[1].push_back( N[1] );
      m_dsupedge[1].push_back( N[2] );
      for (int ed=0; ed<3; ++ed) {
        m_dsupint[1].push_back( sig[ed] * d[ed]->second[0] );
        m_dsupint[1].push_back( sig[ed] * d[ed]->second[1] );
        m_dsupint[1].push_back( sig[ed] * d[ed]->second[2] );
        m_domedgeint.erase( d[ed] );
      }
    }
  }

  m_dsupedge[2].resize( m_domedgeint.size()*2 );
  m_dsupint[2].resize( m_domedgeint.size()*3 );
  std::size_t k = 0;
  for (const auto& [ed,d] : m_domedgeint) {
    auto e = m_dsupedge[2].data() + k*2;
    e[0] = tk::cref_find( lid, ed[0] );
    e[1] = tk::cref_find( lid, ed[1] );
    auto i = m_dsupint[2].data() + k*3;
    i[0] = d[0];
    i[1] = d[1];
    i[2] = d[2];
    ++k;
  }

  //std::cout << std::setprecision(2)
  //          << "superedges: ntet:" << m_dsupedge[0].size()/4 << "(nedge:"
  //          << m_dsupedge[0].size()/4*6 << ","
  //          << 100.0 * static_cast< tk::real >( m_dsupedge[0].size()/4*6 ) /
  //                     static_cast< tk::real >( nedge )
  //          << "%) + ntri:" << m_dsupedge[1].size()/3
  //          << "(nedge:" << m_dsupedge[1].size() << ","
  //          << 100.0 * static_cast< tk::real >( m_dsupedge[1].size() ) /
  //                     static_cast< tk::real >( nedge )
  //          << "%) + nedge:"
  //          << m_dsupedge[2].size()/2 << "("
  //          << 100.0 * static_cast< tk::real >( m_dsupedge[2].size()/2 ) /
  //                     static_cast< tk::real >( nedge )
  //          << "%) = " << m_dsupedge[0].size()/4*6 + m_dsupedge[1].size() +
  //             m_dsupedge[2].size()/2 << " of "<< nedge << " total edges\n";

  Assert( m_dsupedge[0].size()/4*6 + m_dsupedge[1].size() +
          m_dsupedge[2].size()/2 == nedge,
          "Not all edges accounted for in superedge groups" );
}

void
// cppcheck-suppress unusedFunction
LaxCG::merge()
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
  BC( d->T() );

  if (d->Initial()) {
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_LaxCG::start(), thisProxy[thisIndex]) );
  } else {
    feop_complete();
  }
}

void
LaxCG::BC( tk::real t )
// *****************************************************************************
// Apply boundary conditions
//! \param[in] t Physical time
// *****************************************************************************
{
  auto d = Disc();

  // Apply Dirichlet BCs
  physics::dirbc( m_u, t, d->Coord(), d->BoxNodes(), m_dirbcmasks );

  // Apply symmetry BCs
  physics::symbc( m_u, m_symbcnodes, m_symbcnorms );

  // Apply farfield BCs
  physics::farbc( m_u, m_farbcnodes, m_farbcnorms );

  // Apply pressure BCs
  physics::prebc( m_u, m_prebcnodes, m_prebcvals );
}

void
LaxCG::dt()
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

    if (g_cfg.get< tag::steady >()) {
      physics::dt( d->Vol(), m_u, m_dtp );
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );
    } else {
      mindt = physics::dt( d->Vol(), m_u );
    }

  }

  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  // Contribute to minimum dt across all chares and advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(LaxCG,advance), thisProxy) );
}

void
LaxCG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
// *****************************************************************************
{
  // Set new time step size
  if (m_stage == 0) Disc()->setdt( newdt );

  grad();
}

void
LaxCG::grad()
// *****************************************************************************
// Compute gradients for next time step
// *****************************************************************************
{
  auto d = Disc();

  lax::grad( m_bpoin, m_bpint, m_dsupedge, m_dsupint, m_bsupedge, m_bsupint,
             m_u, m_grad );

  // Send gradient contributions to neighbor chares
  if (d->NodeCommMap().empty()) {
    comgrad_complete();
  } else {
    const auto& lid = d->Lid();
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::unordered_map< std::size_t, std::vector< tk::real > > exp;
      for (auto g : n) exp[g] = m_grad[ tk::cref_find(lid,g) ];
      thisProxy[c].comgrad( exp );
    }
  }
  owngrad_complete();
}

void
LaxCG::comgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to node gradients on chare-boundaries
//! \param[in] ingrad Partial contributions to chare-boundary nodes. Key: 
//!   global mesh node IDs, value: contributions for all scalar components.
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
LaxCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();
  const auto steady = g_cfg.get< tag::steady >();

  // Combine own and communicated contributions to gradients
  for (const auto& [g,r] : m_gradc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_grad(i,c) += r[c];
  }
  tk::destroy(m_gradc);

  // divide weak result in gradients by nodal volume
  const auto& vol = d->Vol();
  for (std::size_t p=0; p<m_grad.nunk(); ++p)
    for (std::size_t c=0; c<m_grad.nprop(); ++c)
      m_grad(p,c) /= vol[p];

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];

  if (steady) {
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  }

  lax::rhs( m_dsupedge, m_dsupint, m_bsupedge, m_bsupint, m_bpoin, m_bpint,
            m_bpsym, d->Coord(), m_grad, m_u, d->V(), d->T(), m_tp, m_rhs );

  if (steady) {
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];
  }

  // Communicate rhs to other chares on chare-boundary
  if (d->NodeCommMap().empty()) {
    comrhs_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::unordered_map< std::size_t, std::vector< tk::real > > exp;
      for (auto g : n) exp[g] = m_rhs[ tk::cref_find(lid,g) ];
      thisProxy[c].comrhs( exp );
    }
  }
  ownrhs_complete();
}

void
LaxCG::comrhs(
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
// cppcheck-suppress unusedFunction
LaxCG::solve()
// *****************************************************************************
//  Advance systems of equations
// *****************************************************************************
{
  auto d = Disc();
  const auto lid = d->Lid();
  const auto steady = g_cfg.get< tag::steady >();

  // Combine own and communicated contributions to rhs
  for (const auto& [g,r] : m_rhsc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_rhs(i,c) += r[c];
  }
  tk::destroy(m_rhsc);

  // Update state at time n
  if (m_stage == 0) m_un = m_u;

  // Advance solution
  auto dt = d->Dt();
  const auto& vol = d->Vol();
  for (std::size_t i=0; i<m_u.nunk(); ++i) {
    if (steady) dt = m_dtp[i];
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(i,c) = m_un(i,c) - rkcoef[m_stage] * dt * m_rhs(i,c) / vol[i];
    }
  }

  // Configure and apply scalar source to solution (if defined)
  auto src = problems::PHYS_SRC();
  if (src) src( d->Coord(), d->T(), m_u );

  // Enforce boundary conditions
  BC( d->T() + rkcoef[m_stage] * d->Dt() );

  if (m_stage < 2) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4grad();
    thisProxy[ thisIndex ].wait4rhs();

    // start next time step stage
    stage();

  } else {

    // Activate SDAG waits for finishing this time step stage
    thisProxy[ thisIndex ].wait4stage();
    // Compute diagnostics, e.g., residuals
    auto diag = m_diag.compute( *d, m_u, m_un, g_cfg.get< tag::diag_iter >() );
    // Increase number of iterations and physical time
    d->next();
    // Advance physical time for local time stepping
    if (steady) {
      using tk::operator+=;
      m_tp += m_dtp;
    }
    // Evaluate residuals
    if (!diag) evalres( std::vector< tk::real >( m_u.nprop(), 1.0 ) );

  }
}

void
LaxCG::evalres( const std::vector< tk::real >& l2res )
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
LaxCG::refine()
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
LaxCG::resizePostAMR(
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
  m_un.rm( removedNodes );
  m_rhs.rm( removedNodes );
  m_grad.rm( removedNodes );

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
LaxCG::writeFields( CkCallback cb )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] cb Function to continue with after the write
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) { cb.send(); return; }

  auto d = Disc();
  auto ncomp = m_u.nprop();
  auto steady = g_cfg.get< tag::steady >();

  // Field output

  std::vector< std::string > nodefieldnames
    {"density", "velocityx", "velocityy", "velocityz", "energy", "pressure"};
  if (steady) nodefieldnames.push_back( "mach" );

  using tk::operator/=;
  auto r = m_u.extract(0);
  auto u = m_u.extract(1);  u /= r;
  auto v = m_u.extract(2);  v /= r;
  auto w = m_u.extract(3);  w /= r;
  auto e = m_u.extract(4);  e /= r;
  std::vector< tk::real > pr( m_u.nunk() ), ma;
  if (steady) ma.resize( m_u.nunk() );
  for (std::size_t i=0; i<pr.size(); ++i) {
    auto vv = u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
    pr[i] = eos::pressure( r[i]*(e[i] - 0.5*vv) );
    if (steady) ma[i] = std::sqrt(vv) / eos::soundspeed( r[i], pr[i] );
  }

  std::vector< std::vector< tk::real > > nodefields{
    std::move(r), std::move(u), std::move(v), std::move(w), std::move(e),
    std::move(pr) };
  if (steady) nodefields.push_back( std::move(ma) );

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

    if (steady) nodesurfnames.push_back( "mach" );

    auto bnode = tk::bfacenodes( m_bface, m_triinpoel );
    std::set< int > outsets( begin(f), end(f) );
    for (auto sideset : outsets) {
      auto b = bnode.find(sideset);
      if (b == end(bnode)) continue;
      const auto& nodes = b->second;
      auto i = nodesurfs.size();
      auto ns = ncomp + 1;
      if (steady) ++ns;
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
        if (steady) {
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
LaxCG::out()
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
    writeFields( CkCallback(CkIndex_LaxCG::integrals(), thisProxy[thisIndex]) );
  } else {
    integrals();
  }
}

void
LaxCG::integrals()
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
    // Compute mass flow rate for surfaces requested
    for (const auto& [s,sint] : m_surfint) {
      // cppcheck-suppress unreadVariable
      auto& mfr = ints[ MASS_FLOW_RATE ][ s ];
      const auto& nodes = sint.first;
      const auto& ndA = sint.second;
      for (std::size_t i=0; i<nodes.size(); ++i) {
        auto p = nodes[i];
        mfr += ndA[i*3+0] * m_u(p,1)
             + ndA[i*3+1] * m_u(p,2)
             + ndA[i*3+2] * m_u(p,3);
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
LaxCG::stage()
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
LaxCG::evalLB( int nrestart )
// *****************************************************************************
// Evaluate whether to do load balancing
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  auto d = Disc();

  // Detect if just returned from a checkpoint and if so, zero timers and
  // finished flag
  if (d->restarted( nrestart )) m_finished = 0;

  const auto lbfreq = g_cfg.get< tag::lbfreq >();
  const auto nonblocking = g_cfg.get< tag::nonblocking >();

  // Load balancing if user frequency is reached or after the second time-step
  if ( (d->It()) % lbfreq == 0 || d->It() == 2 ) {

    AtSync();
    if (nonblocking) dt();

  } else {

    dt();

  }
}

void
LaxCG::evalRestart()
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
LaxCG::step()
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

#include "NoWarning/laxcg.def.h"
