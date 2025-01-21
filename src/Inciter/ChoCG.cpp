// *****************************************************************************
/*!
  \file      src/Inciter/ChoCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     ChoCG: Projection-based solver for incompressible flow
*/
// *****************************************************************************

#include "XystBuildConfig.hpp"
#include "ChoCG.hpp"
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
#include "Chorin.hpp"
#include "Problems.hpp"
#include "EOS.hpp"
#include "BC.hpp"
#include "Print.hpp"

namespace inciter {

extern ctr::Config g_cfg;

static CkReduction::reducerType IntegralsMerger;

//! Runge-Kutta coefficients
//! Runge-Kutta coefficients
static const std::array< std::vector< tk::real >, 4 > rkcoef{{
  { 1.0 },
  { 1.0/2.0, 1.0 },
  { 1.0/3.0, 1.0/2.0, 1.0 },
  { 1.0/4.0, 1.0/3.0, 1.0/2.0, 1.0 }
}};

} // inciter::

using inciter::g_cfg;
using inciter::ChoCG;

ChoCG::ChoCG( const CProxy_Discretization& disc,
              const tk::CProxy_ConjugateGradients& cgpre,
              const tk::CProxy_ConjugateGradients& cgmom,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_cgpre( cgpre ),
  m_cgmom( cgmom ),
  m_nrhs( 0 ),
  m_nnorm( 0 ),
  m_naec( 0 ),
  m_nalw( 0 ),
  m_nlim( 0 ),
  m_nsgrad( 0 ),
  m_npgrad( 0 ),
  m_nvgrad( 0 ),
  m_nflux( 0 ),
  m_ndiv( 0 ),
  m_nbpint( 0 ),
  m_np( 0 ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_u( Disc()->Gid().size(), g_cfg.get< tag::problem_ncomp >() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_pr( m_u.nunk(), 0.0 ),
  m_p( m_u.nunk(), m_u.nprop()*2 ),
  m_q( m_u.nunk(), m_u.nprop()*2 ),
  m_a( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_sgrad( m_u.nunk(), 3UL ),
  m_pgrad( m_u.nunk(), 3UL ),
  m_vgrad( m_u.nunk(), m_u.nprop()*3 ),
  m_flux( m_u.nunk(), 3UL ),
  m_div( m_u.nunk() ),
  m_stage( 0 ),
  m_finished( 0 ),
  m_rkcoef( rkcoef[ g_cfg.get< tag::rk >() - 1 ] )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] cgpre ConjugateGradients Charm++ proxy for pressure solve
//! \param[in] cgmom ConjugateGradients Charm++ proxy for momentum solve
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
  // Recompute points surrounding points
  psup = tk::genPsup( d->Inpoel(), 4, tk::genEsup( d->Inpoel(), 4 ) );

  // Compute total box IC volume
  d->boxvol();

  // Setup LHS matrix for pressure solve
  m_cgpre[ thisIndex ].insert( prelhs( psup ),
                               d->Gid(),
                               d->Lid(),
                               d->NodeCommMap() );

  // Setup empty LHS matrix for momentum solve if needed
  if (g_cfg.get< tag::theta >() > std::numeric_limits< tk::real >::epsilon()) {
    if (g_cfg.get< tag::fct >()) Throw( "Implicit FCT not implemented" );
    m_cgmom[ thisIndex ].insert( momlhs( psup ),
                                 d->Gid(),
                                 d->Lid(),
                                 d->NodeCommMap() );
  }

  // Activate SDAG waits for setup
  thisProxy[ thisIndex ].wait4int();
}

std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
ChoCG::prelhs( const std::pair< std::vector< std::size_t >,
                                std::vector< std::size_t > >& psup )
// *****************************************************************************
//  Setup lhs matrix for pressure solve
//! \param[in] psup Points surrounding points
//! \return { A, x, b } in linear system A * x = b to solve for pressure
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];

  // Matrix with compressed sparse row storage
  tk::CSR A( /*DOF=*/ 1, psup );

  // Fill matrix with Laplacian
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) * 6.0;
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::cross( ca, da );
    grad[2] = tk::cross( da, ba );
    grad[3] = tk::cross( ba, ca );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    for (std::size_t a=0; a<4; ++a)
      for (std::size_t b=0; b<4; ++b)
        A(N[a],N[b]) -= tk::dot( grad[a], grad[b] ) / J;
  }

  auto nunk = X.size();
  std::vector< tk::real > x( nunk, 0.0 ), b( nunk, 0.0 );

  return { std::move(A), std::move(x), std::move(b) };
}

std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
ChoCG::momlhs( const std::pair< std::vector< std::size_t >,
                                std::vector< std::size_t > >& psup )
// *****************************************************************************
//  Setup empty lhs matrix for momentum solve
//! \param[in] psup Points surrounding points
//! \return { A, x, b } in linear system A * x = b to solve for momentum
// *****************************************************************************
{
  auto ncomp = m_u.nprop();

  // Matrix with compressed sparse row storage
  tk::CSR A( /*DOF=*/ ncomp, psup );

  auto nunk = (psup.second.size() - 1) * ncomp;
  std::vector< tk::real > x( nunk, 0.0 ), b( nunk, 0.0 );

  return { std::move(A), std::move(x), std::move(b) };
}

void
ChoCG::setupDirBC( const std::vector< std::vector< int > >& cfgmask,
                   const std::vector< std::vector< double > >& cfgval,
                   std::size_t ncomp,
                   std::vector< std::size_t >& mask,
                   std::vector< double >& val )
// *****************************************************************************
//  Prepare Dirichlet boundary condition data structures
//! \param[in] cfgmask Boundary condition mask config data to use
//! \param[in] cfgval Boundary condition values config data to use
//! \param[in] ncomp Number of scalar component BCs expected per mesh node
//! \param[in,out] mask Mesh nodes and their Dirichlet BC masks
//! \param[in,out] val Mesh nodes and their Dirichlet BC values
// *****************************************************************************
{
  // Query Dirichlet BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > dir;
  for (const auto& s : cfgmask) {
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
  for (const auto& s : cfgmask) {
    auto k = m_bnode.find(s[0]);
    if (k != end(m_bnode)) {
      auto& n = dir[ k->first ];
      for (auto g : k->second) {
        n.insert( tk::cref_find(lid,g) );
      }
    }
  }

  // Associate sidesets to Dirichlet BC values if configured by user
  std::unordered_map< int, std::vector< double > > dirval;
  for (const auto& s : cfgval) {
    auto k = dir.find( static_cast<int>(s[0]) );
    if (k != end(dir)) {
      auto& v = dirval[ k->first ];
      v.resize( s.size()-1 );
      for (std::size_t i=1; i<s.size(); ++i) v[i-1] = s[i];
    }
  }

  // Collect unique set of nodes + Dirichlet BC components mask and value
  auto nmask = ncomp + 1;
  std::unordered_map< std::size_t,
                      std::pair< std::vector< int >,
                                 std::vector< double > > > dirbcset;
  for (const auto& vec : cfgmask) {
    ErrChk( vec.size() == nmask, "Incorrect Dirichlet BC mask ncomp" );
    auto n = dir.find( vec[0] );
    if (n != end(dir)) {
      std::vector< double > v( ncomp, 0.0 );
      auto m = dirval.find( vec[0] );
      if (m != end(dirval)) {
        ErrChk( m->second.size() == ncomp, "Incorrect Dirichlet BC val ncomp" );
        v = m->second;
      }
      for (auto p : n->second) {
        auto& mv = dirbcset[p]; // mask & value
        mv.second = v;
        auto& mval = mv.first;
        if (mval.empty()) mval.resize( ncomp, 0 );
        for (std::size_t c=0; c<ncomp; ++c)
          if (!mval[c]) mval[c] = vec[c+1];  // overwrite mask if 0 -> 1
      }
    }
  }

  // Compile streamable list of nodes + Dirichlet BC components mask and values
  tk::destroy( mask );
  for (const auto& [p,mv] : dirbcset) {
    mask.push_back( p );
    mask.insert( end(mask), begin(mv.first), end(mv.first) );
    val.push_back( static_cast< double >( p ) );
    val.insert( end(val), begin(mv.second), end(mv.second) );
  }

  ErrChk( mask.size() % nmask == 0, "Dirichlet BC mask incomplete" );
  ErrChk( val.size() % nmask == 0, "Dirichlet BC val incomplete" );
  ErrChk( mask.size() == val.size(), "Dirichlet BC mask & val size mismatch" );
}

void
ChoCG::feop()
// *****************************************************************************
// Start (re-)computing finite element domain and boundary operators
// *****************************************************************************
{
  auto d = Disc();

  // Prepare Dirichlet boundary conditions data structures
  setupDirBC( g_cfg.get< tag::bc_dir >(), g_cfg.get< tag::bc_dirval >(),
              m_u.nprop(), m_dirbcmask, m_dirbcval );
  setupDirBC( g_cfg.get< tag::pre_bc_dir >(), g_cfg.get< tag::pre_bc_dirval >(),
              1, m_dirbcmaskp, m_dirbcvalp );

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
ChoCG::bndint()
// *****************************************************************************
//  Compute local contributions to boundary normals and integrals
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
        tk::real r = invdistsq( centroid, p );
        auto& v = m_bnorm[setid];      // associate side set id
        auto& bpn = v[gid[p]];         // associate global node id of bnd pnt
        bpn[0] += r * n[0];            // inv.dist.sq-weighted normal
        bpn[1] += r * n[1];
        bpn[2] += r * n[2];
        bpn[3] += r;                   // inv.dist.sq of node from centroid
        auto& b = m_bndpoinint[gid[p]];// assoc global id of bnd point
        b[0] += n[0] * A2 / 6.0;       // bnd-point integral
        b[1] += n[1] * A2 / 6.0;
        b[2] += n[2] * A2 / 6.0;
      }
    }
  }
}

void
ChoCG::domint()
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
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );
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
      n[3] += J / 120.0;
      n[4] += tk::dot( grad[p], grad[q] ) / J / 6.0;
    }
  }
}

void
ChoCG::comnorm( const decltype(m_bnorm)& inbnd )
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
ChoCG::registerReducers()
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
ChoCG::ResumeFromSync()
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
ChoCG::setup( tk::real v )
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
ChoCG::bnorm()
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
ChoCG::streamable()
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
        auto t = m_triinpoel.data() + f*3;
        n.push_back( t[0] );                    // nodes on side set
        n.push_back( t[1] );
        n.push_back( t[2] );
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
    auto a = ndA.data();
    for (auto p : nodes) {
      const auto& b = tk::cref_find( m_bndpoinint, gid[p] );
      a[0] = b[0];      // store ni * dA
      a[1] = b[1];
      a[2] = b[2];
      a += 3;
    }
  }
  tk::destroy( m_bndpoinint );

  // Generate domain superedges
  domsuped();

  //  Prepare symmetry boundary condition data structures

  // Query symmetry BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > sym;
  for (auto s : g_cfg.get< tag::bc_sym >()) {
    auto k = m_bface.find(s);
    if (k != end(m_bface)) {
      auto& n = sym[ k->first ];
      for (auto f : k->second) {
        const auto& t = m_triinpoel.data() + f*3;
        n.insert( t[0] );
        n.insert( t[1] );
        n.insert( t[2] );
      }
    }
  }

  // Generate unique set of symmetry BC nodes of all symmetryc BC side sets
  std::set< std::size_t > symbcnodeset;
  for (const auto& [s,n] : sym) symbcnodeset.insert( begin(n), end(n) );

  // Generate symmetry BC data as streamable data structures
  tk::destroy( m_symbcnodes );
  tk::destroy( m_symbcnorms );
  for (auto p : symbcnodeset) {
    for (const auto& s : g_cfg.get< tag::bc_sym >()) {
      auto m = m_bnorm.find( s );
      if (m != end(m_bnorm)) {
        auto r = m->second.find( p );
        if (r != end(m->second)) {
          m_symbcnodes.push_back( p );
          m_symbcnorms.push_back( r->second[0] );
          m_symbcnorms.push_back( r->second[1] );
          m_symbcnorms.push_back( r->second[2] );
        }
      }
    }
  }
  tk::destroy( m_bnorm );

  //  Prepare noslip boundary condition data structures

  // Query noslip BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > noslip;
  for (auto s : g_cfg.get< tag::bc_noslip >()) {
    auto k = m_bface.find(s);
    if (k != end(m_bface)) {
      auto& n = noslip[ k->first ];
      for (auto f : k->second) {
        const auto& t = m_triinpoel.data() + f*3;
        n.insert( t[0] );
        n.insert( t[1] );
        n.insert( t[2] );
      }
    }
  }

  // Generate unique set of noslip BC nodes of all noslip BC side sets
  std::set< std::size_t > noslipbcnodeset;
  for (const auto& [s,n] : noslip) noslipbcnodeset.insert( begin(n), end(n) );

  // Generate noslip BC data as streamable data structures
  tk::destroy( m_noslipbcnodes );
  m_noslipbcnodes.insert( m_noslipbcnodes.end(),
                          begin(noslipbcnodeset), end(noslipbcnodeset) );
}

void
ChoCG::domsuped()
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
        const auto& ded = d[ed]->second;
        m_dsupint[0].push_back( sig[ed] * ded[0] );
        m_dsupint[0].push_back( sig[ed] * ded[1] );
        m_dsupint[0].push_back( sig[ed] * ded[2] );
        m_dsupint[0].push_back( ded[3] );
        m_dsupint[0].push_back( ded[4] );
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
        const auto& ded = d[ed]->second;
        m_dsupint[1].push_back( sig[ed] * ded[0] );
        m_dsupint[1].push_back( sig[ed] * ded[1] );
        m_dsupint[1].push_back( sig[ed] * ded[2] );
        m_dsupint[1].push_back( ded[3] );
        m_dsupint[1].push_back( ded[4] );
        m_domedgeint.erase( d[ed] );
      }
    }
  }

  m_dsupedge[2].resize( m_domedgeint.size()*2 );
  m_dsupint[2].resize( m_domedgeint.size()*5 );
  std::size_t k = 0;
  for (const auto& [ed,d] : m_domedgeint) {
    auto e = m_dsupedge[2].data() + k*2;
    e[0] = tk::cref_find( lid, ed[0] );
    e[1] = tk::cref_find( lid, ed[1] );
    auto i = m_dsupint[2].data() + k*5;
    i[0] = d[0];
    i[1] = d[1];
    i[2] = d[2];
    i[3] = d[3];
    i[4] = d[4];
    ++k;
  }

  if (g_cfg.get< tag::fct >()) {
    const auto ncomp = m_u.nprop();
    m_dsuplim[0].resize( m_dsupedge[0].size() * 6 * ncomp );
    m_dsuplim[1].resize( m_dsupedge[1].size() * 3 * ncomp );
    m_dsuplim[2].resize( m_dsupedge[2].size() * ncomp );
  }

  tk::destroy( m_domedgeint );

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
ChoCG::merge()
// *****************************************************************************
// Combine own and communicated portions of the integrals
// *****************************************************************************
{
  // Combine own and communicated contributions to boundary point normals
  bnorm();

  // Convert integrals into streamable data structures
  streamable();

  // Enforce boundary conditions on initial conditions
  BC( m_u, Disc()->T() );

  // Start measuring initial div-free time
  m_timer.emplace_back();

  // Compute initial momentum flux
  thisProxy[ thisIndex ].wait4div();
  thisProxy[ thisIndex ].wait4sgrad();
  div( m_u );
}

void
ChoCG::fingrad( tk::Fields& grad,
  std::unordered_map< std::size_t, std::vector< tk::real > >& gradc )
// *****************************************************************************
//  Finalize computing gradient
//! \param[in,out] grad Gradient to finalize
//! \param[in,out] gradc Gradient communication buffer to finalize
// *****************************************************************************
{
  auto d = Disc();
  const auto lid = d->Lid();

  // Combine own and communicated contributions
  for (const auto& [g,r] : gradc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) grad(i,c) += r[c];
  }
  tk::destroy(gradc);

  // Divide weak result by nodal volume
  const auto& vol = d->Vol();
  for (std::size_t p=0; p<grad.nunk(); ++p) {
    for (std::size_t c=0; c<grad.nprop(); ++c) {
      grad(p,c) /= vol[p];
    }
  }
}

void
ChoCG::div( const tk::Fields& u )
// *****************************************************************************
//  Start computing divergence
// \para[in] u Vector field whose divergence to compute
// *****************************************************************************
{
  auto d = Disc();
  const auto lid = d->Lid();

  // Finalize momentum flux communications if needed
  if (m_np == 1) {
    fingrad( m_flux, m_fluxc );
    physics::symbc( m_flux, m_symbcnodes, m_symbcnorms, /*pos=*/0 );
  }

  // Compute divergence
  std::fill( begin(m_div), end(m_div), 0.0 );
  chorin::div( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel,
               d->Dt(), m_pr, m_pgrad, u, m_div, m_np>1 );

  // Communicate velocity divergence to other chares on chare-boundary
  if (d->NodeCommMap().empty()) {
    comdiv_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_divc) exp;
      for (auto g : n) exp[g] = m_div[ tk::cref_find(lid,g) ];
      thisProxy[c].comdiv( exp );
    }
  }
  owndiv_complete();
}

void
ChoCG::comdiv( const std::unordered_map< std::size_t, tk::real >& indiv )
// *****************************************************************************
//  Receive contributions to velocity divergence on chare-boundaries
//! \param[in] indiv Partial contributions of velocity divergence to
//!   chare-boundary nodes. Key: global mesh node IDs, value: contribution.
//! \details This function receives contributions to m_div, which stores the
//!   velocity divergence at mesh nodes. While m_div stores own contributions,
//!   m_divc collects the neighbor chare contributions during communication.
//!   This way work on m_div and m_divc is overlapped.
// *****************************************************************************
{
  for (const auto& [g,r] : indiv) m_divc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_ndiv == Disc()->NodeCommMap().size()) {
    m_ndiv = 0;
    comdiv_complete();
  }
}

void
ChoCG::velgrad()
// *****************************************************************************
//  Start computing velocity gradient
// *****************************************************************************
{
  auto d = Disc();

  // Compute momentum flux
  m_vgrad.fill( 0.0 );
  chorin::vgrad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_u, m_vgrad );

  // Communicate velocity divergence to other chares on chare-boundary
  const auto& lid = d->Lid();
  if (d->NodeCommMap().empty()) {
    comvgrad_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_vgradc) exp;
      for (auto g : n) exp[g] = m_vgrad[ tk::cref_find(lid,g) ];
      thisProxy[c].comvgrad( exp );
    }
  }
  ownvgrad_complete();
}

void
ChoCG::comvgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to velocity gradients on chare-boundaries
//! \param[in] ingrad Partial contributions of momentum flux to
//!   chare-boundary nodes. Key: global mesh node IDs, values: contributions.
//! \details This function receives contributions to m_vgrad, which stores the
//!   velocity gradients at mesh nodes. While m_vgrad stores own contributions,
//!   m_vgradc collects the neighbor chare contributions during communication.
//!   This way work on m_vgrad and m_vgradc is overlapped.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : ingrad) m_vgradc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nvgrad == Disc()->NodeCommMap().size()) {
    m_nvgrad = 0;
    comvgrad_complete();
  }
}

void
ChoCG::flux()
// *****************************************************************************
//  Start computing momentum flux
// *****************************************************************************
{
  auto d = Disc();

  // Finalize computing velocity gradients
  fingrad( m_vgrad, m_vgradc );

  // Compute momentum flux
  m_flux.fill( 0.0 );
  chorin::flux( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_u, m_vgrad,
                m_flux );

  // Communicate velocity divergence to other chares on chare-boundary
  const auto& lid = d->Lid();
  if (d->NodeCommMap().empty()) {
    comflux_complete();
  } else {
    for (const auto& [c,n] : d->NodeCommMap()) {
      decltype(m_fluxc) exp;
      for (auto g : n) exp[g] = m_flux[ tk::cref_find(lid,g) ];
      thisProxy[c].comflux( exp );
    }
  }
  ownflux_complete();
}

void
ChoCG::comflux(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& influx )
// *****************************************************************************
//  Receive contributions to momentum flux on chare-boundaries
//! \param[in] influx Partial contributions of momentum flux to
//!   chare-boundary nodes. Key: global mesh node IDs, values: contributions.
//! \details This function receives contributions to m_flux, which stores the
//!   momentum flux at mesh nodes. While m_flux stores own contributions,
//!   m_fluxc collects the neighbor chare contributions during communication.
//!   This way work on m_flux and m_fluxc is overlapped.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : influx) m_fluxc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nflux == Disc()->NodeCommMap().size()) {
    m_nflux = 0;
    comflux_complete();
  }
}

void
ChoCG::pinit()
// *****************************************************************************
//  Initialize Poisson solve
// *****************************************************************************
{
  auto d = Disc();
  const auto lid = d->Lid();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Combine own and communicated contributions to velocity divergence
  for (const auto& [g,r] : m_divc) m_div[ tk::cref_find(lid,g) ] += r;
  tk::destroy(m_divc);

  // Divide Poisson rhs by dt if solving for time-increment
  if (m_np > 1) for (auto& div : m_div) div /= d->Dt();

  // Configure Poisson BCs
  std::unordered_map< std::size_t,
    std::vector< std::pair< int, tk::real > > > dirbc;
  std::vector< tk::real > neubc;

  if (m_np < 3) {       // only during startup and first time step
    // Configure Dirichlet BCs
    if (!g_cfg.get< tag::pre_bc_dir >().empty()) {
      auto ic = problems::PRESSURE_IC();
      std::size_t nmask = 1 + 1;
      Assert( m_dirbcmaskp.size() % nmask == 0, "Size mismatch" );
      for (std::size_t i=0; i<m_dirbcmaskp.size()/nmask; ++i) {
        auto p = m_dirbcmaskp[i*nmask+0];     // local node id
        auto mask = m_dirbcmaskp[i*nmask+1];
        if (mask == 1) {                                  // mask == 1: IC value
          auto val = m_np>1 ? 0.0 : ic( x[p], y[p], z[p] );
          dirbc[p] = {{ { 1, val } }};
        } else if (mask == 2 && !m_dirbcvalp.empty()) {   // mask == 2: BC value
          auto val = m_np>1 ? 0.0 : m_dirbcvalp[i*nmask+1];
          dirbc[p] = {{ { 1, val } }};
        }
      }
    }

    // Configure Neumann BCs
    auto pg = problems::PRESSURE_GRAD();
    if (pg) {
      // Collect Neumann BC elements
      std::vector< std::uint8_t > besym( m_triinpoel.size(), 0 );
      for (auto s : g_cfg.get< tag::pre_bc_sym >()) {
        auto k = m_bface.find(s);
        if (k != end(m_bface)) for (auto f : k->second) besym[f] = 1;
      }
      // Setup Neumann BCs
      neubc.resize( x.size(), 0.0 );
      for (std::size_t e=0; e<m_triinpoel.size()/3; ++e) {
        if (besym[e]) {
          const auto N = m_triinpoel.data() + e*3;
          tk::real n[3];
          tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                        x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                        n[0], n[1], n[2] );
          auto g = pg( x[N[0]], y[N[0]], z[N[0]] );
          neubc[ N[0] ] -= n[0]*g[0] + n[1]*g[1] + n[2]*g[2];
          g = pg( x[N[1]], y[N[1]], z[N[1]] );
          neubc[ N[1] ] -= n[0]*g[0] + n[1]*g[1] + n[2]*g[2];
          g = pg( x[N[2]], y[N[2]], z[N[2]] );
          neubc[ N[2] ] -= n[0]*g[0] + n[1]*g[1] + n[2]*g[2];
        }
      }
    }

    // Set hydrostat
    auto h = g_cfg.get< tag::pre_hydrostat >();
    if (h != std::numeric_limits< uint64_t >::max()) {
      auto pi = lid.find( h );
      if (pi != end(lid)) {
        auto p = pi->second;
        auto ic = problems::PRESSURE_IC();
        auto val = m_np>1 ? 0.0 : ic( x[p], y[p], z[p] );
        auto& b = dirbc[p];
        if (b.empty()) b = {{ { 1, val }} };
      }
    }

    // Configure right hand side
    auto pr = problems::PRESSURE_RHS();
    if (pr) {
      const auto& vol = d->Vol();
      for (std::size_t i=0; i<x.size(); ++i) {
        m_div[i] = pr( x[i], y[i], z[i] ) * vol[i];
      }
    }
  }

  // Initialize Poisson solve setting ICs and BCs
  m_cgpre[ thisIndex ].ckLocal()->
    init( {}, m_div, neubc, dirbc, m_np<3, g_cfg.get< tag::pre_pc >(),
          CkCallback( CkIndex_ChoCG::psolve(), thisProxy[thisIndex] ) );
}

void
ChoCG::psolve()
// *****************************************************************************
//  Solve Poisson equation
// *****************************************************************************
{
  auto iter = g_cfg.get< tag::pre_iter >();
  auto tol = g_cfg.get< tag::pre_tol >();
  auto verbose = g_cfg.get< tag::pre_verbose >();

  auto c = m_np != 1 ?
           CkCallback( CkIndex_ChoCG::sgrad(), thisProxy[thisIndex] ) :
           CkCallback( CkIndex_ChoCG::psolved(), thisProxy[thisIndex] );

  m_cgpre[ thisIndex ].ckLocal()->solve( iter, tol, thisIndex, verbose, c );
}

void
ChoCG::sgrad()
// *****************************************************************************
// Compute recent conjugate gradients solution gradient
// *****************************************************************************
{
  auto d = Disc();

  auto sol = m_cgpre[ thisIndex ].ckLocal()->solution();
  m_sgrad.fill( 0.0 );
  chorin::grad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, sol, m_sgrad );

  // Send gradient contributions to neighbor chares
  if (d->NodeCommMap().empty()) {
    comsgrad_complete();
  } else {
    const auto& lid = d->Lid();
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::unordered_map< std::size_t, std::vector< tk::real > > exp;
      for (auto g : n) exp[g] = m_sgrad[ tk::cref_find(lid,g) ];
      thisProxy[c].comsgrad( exp );
    }
  }
  ownsgrad_complete();
}

void
ChoCG::comsgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to conjugrate gradients solution gradient
//! \param[in] ingrad Partial contributions to chare-boundary nodes. Key: 
//!   global mesh node IDs, value: contributions for all scalar components.
//! \details This function receives contributions to m_sgrad, which stores the
//!   gradients at mesh nodes. While m_sgrad stores own contributions, m_sgradc
//!   collects the neighbor chare contributions during communication. This way
//!   work on m_sgrad and m_sgradc is overlapped.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : ingrad) m_sgradc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nsgrad == Disc()->NodeCommMap().size()) {
    m_nsgrad = 0;
    comsgrad_complete();
  }
}

void
ChoCG::psolved()
// *****************************************************************************
// Continue setup after Poisson solve and gradient computation
// *****************************************************************************
{
  auto d = Disc();

  if (thisIndex == 0) d->pit( m_cgpre[ thisIndex ].ckLocal()->it() );

  if (m_np != 1) {
    // Finalize gradient communications
    fingrad( m_sgrad, m_sgradc );
    // Project velocity to divergence-free subspace
    auto dt = m_np > 1 ? d->Dt() : 1.0;
    for (std::size_t i=0; i<m_u.nunk(); ++i) {
      m_u(i,0) -= dt * m_sgrad(i,0);
      m_u(i,1) -= dt * m_sgrad(i,1);
      m_u(i,2) -= dt * m_sgrad(i,2);
    }
    // Enforce boundary conditions
    BC( m_u, d->T() + d->Dt() );
  }

  if (d->Initial()) {

    if (g_cfg.get< tag::nstep >() == 1) {  // test first Poisson solve only

      m_pr = m_cgpre[ thisIndex ].ckLocal()->solution();
      thisProxy[ thisIndex ].wait4step();
      writeFields( CkCallback(CkIndex_ChoCG::diag(), thisProxy[thisIndex]) );

    } else {

      if (++m_np < 2) {
        // Compute momentum flux for initial pressure-Poisson rhs
        thisProxy[ thisIndex ].wait4vgrad();
        thisProxy[ thisIndex ].wait4flux();
        thisProxy[ thisIndex ].wait4div();
        velgrad();
      } else {
        if (thisIndex == 0) {
          tk::Print() << "Initial div-free time: " << m_timer[0].dsec()
                      << " sec\n";
        }
        // Assign initial pressure and compute its gradient
        m_pr = m_cgpre[ thisIndex ].ckLocal()->solution();
        pgrad();
      }

    }

  } else {

    // Update pressure and compute its gradient
    using tk::operator+=;
    m_pr += m_cgpre[ thisIndex ].ckLocal()->solution();
    pgrad();

  }
}

void
ChoCG::pgrad()
// *****************************************************************************
//  Compute pressure gradient
// *****************************************************************************
{
  auto d = Disc();

  thisProxy[ thisIndex ].wait4pgrad();

  m_pgrad.fill( 0.0 );
  chorin::grad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_pr, m_pgrad );

  // Send gradient contributions to neighbor chares
  if (d->NodeCommMap().empty()) {
    compgrad_complete();
  } else {
    const auto& lid = d->Lid();
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::unordered_map< std::size_t, std::vector< tk::real > > exp;
      for (auto g : n) exp[g] = m_pgrad[ tk::cref_find(lid,g) ];
      thisProxy[c].compgrad( exp );
    }
  }
  ownpgrad_complete();
}

void
ChoCG::compgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to pressure gradient  on chare-boundaries
//! \param[in] ingrad Partial contributions to chare-boundary nodes. Key:
//!   global mesh node IDs, value: contributions for all scalar components.
//! \details This function receives contributions to m_pgrad, which stores the
//!   gradients at mesh nodes. While m_pgrad stores own contributions, m_pgradc
//!   collects the neighbor chare contributions during communication. This way
//!   work on m_pgrad and m_pgradc is overlapped.
// *****************************************************************************
{
  using tk::operator+=;
  for (const auto& [g,r] : ingrad) m_pgradc[g] += r;

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_npgrad == Disc()->NodeCommMap().size()) {
    m_npgrad = 0;
    compgrad_complete();
  }
}

void
ChoCG::finpgrad()
// *****************************************************************************
//  Compute pressure gradient
// *****************************************************************************
{
  auto d = Disc();

  // Finalize pressure gradient communications
  fingrad( m_pgrad, m_pgradc );

  if (d->Initial()) {
    writeFields( CkCallback(CkIndex_ChoCG::start(), thisProxy[thisIndex]) );
  } else {
    diag();
  }
}

void
ChoCG::start()
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
// cppcheck-suppress unusedFunction
ChoCG::aec()
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

  // tetrahedron superedges
  for (std::size_t e=0; e<m_dsupedge[0].size()/4; ++e) {
    const auto N = m_dsupedge[0].data() + e*4;
    const auto D = m_dsupint[0].data();
    std::size_t i = 0;
    for (const auto& [p,q] : tk::lpoed) {
      auto dif = D[(e*6+i)*5+3];
      for (std::size_t c=0; c<ncomp; ++c) {
        auto aec = -dif * ctau * (m_u(N[p],c) - m_u(N[q],c));
        auto a = c*2;
        auto b = a+1;
        if (aec > 0.0) std::swap(a,b);
        m_p(N[p],a) -= aec;
        m_p(N[q],b) += aec;
      }
      ++i;
    }
  }

  // triangle superedges
  for (std::size_t e=0; e<m_dsupedge[1].size()/3; ++e) {
    const auto N = m_dsupedge[1].data() + e*3;
    const auto D = m_dsupint[1].data();
    std::size_t i = 0;
    for (const auto& [p,q] : tk::lpoet) {
      auto dif = D[(e*3+i)*5+3];
      for (std::size_t c=0; c<ncomp; ++c) {
        auto aec = -dif * ctau * (m_u(N[p],c) - m_u(N[q],c));
        auto a = c*2;
        auto b = a+1;
        if (aec > 0.0) std::swap(a,b);
        m_p(N[p],a) -= aec;
        m_p(N[q],b) += aec;
      }
      ++i;
    }
  }

  // edges
  for (std::size_t e=0; e<m_dsupedge[2].size()/2; ++e) {
    const auto N = m_dsupedge[2].data() + e*2;
    const auto dif = m_dsupint[2][e*5+3];
    for (std::size_t c=0; c<ncomp; ++c) {
      auto aec = -dif * ctau * (m_u(N[0],c) - m_u(N[1],c));
      auto a = c*2;
      auto b = a+1;
      if (aec > 0.0) std::swap(a,b);
      m_p(N[0],a) -= aec;
      m_p(N[1],b) += aec;
    }
  }

  // Apply symmetry BCs on AEC
  for (std::size_t i=0; i<m_symbcnodes.size(); ++i) {
    auto p = m_symbcnodes[i];
    auto n = m_symbcnorms.data() + i*3;
    auto rvnp = m_p(p,0)*n[0] + m_p(p,2)*n[1] + m_p(p,4)*n[2];
    auto rvnn = m_p(p,1)*n[0] + m_p(p,3)*n[1] + m_p(p,5)*n[2];
    m_p(p,0) -= rvnp * n[0];
    m_p(p,1) -= rvnn * n[0];
    m_p(p,2) -= rvnp * n[1];
    m_p(p,3) -= rvnn * n[1];
    m_p(p,4) -= rvnp * n[2];
    m_p(p,5) -= rvnn * n[2];
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
ChoCG::comaec( const std::unordered_map< std::size_t,
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
ChoCG::alw()
// *****************************************************************************
// Compute allowed limits, Q+/-
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& lid = d->Lid();
  const auto& vol = d->Vol();

  // Combine own and communicated contributions to antidiffusive contributions
  // and low-order solution
  for (const auto& [g,p] : m_pc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<p.size(); ++c) m_p(i,c) += p[c];
  }
  tk::destroy(m_pc);

  // Finish computing antidiffusive contributions and low-order solution
  auto dt = m_rkcoef[m_stage] * d->Dt();
  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      m_p(i,a) /= vol[i];
      m_p(i,b) /= vol[i];
      // low-order solution
      m_rhs(i,c) = m_un(i,c) - dt*m_rhs(i,c)/vol[i] - m_p(i,a) - m_p(i,b);
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

  // tetrahedron superedges
  for (std::size_t e=0; e<m_dsupedge[0].size()/4; ++e) {
    const auto N = m_dsupedge[0].data() + e*4;
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      for (const auto& [p,q] : tk::lpoed) {
        tk::real alwp, alwn;
        if (g_cfg.get< tag::fctclip >()) {
          alwp = max( m_rhs(N[p],c), m_rhs(N[q],c) );
          alwn = min( m_rhs(N[p],c), m_rhs(N[q],c) );
        } else {
          alwp = max( max(m_rhs(N[p],c), m_u(N[p],c)),
                      max(m_rhs(N[q],c), m_u(N[q],c)) );
          alwn = min( min(m_rhs(N[p],c), m_u(N[p],c)),
                      min(m_rhs(N[q],c), m_u(N[q],c)) );
        }
        m_q(N[p],a) = max(m_q(N[p],a), alwp);
        m_q(N[p],b) = min(m_q(N[p],b), alwn);
        m_q(N[q],a) = max(m_q(N[q],a), alwp);
        m_q(N[q],b) = min(m_q(N[q],b), alwn);
      }
    }
  }

  // triangle superedges
  for (std::size_t e=0; e<m_dsupedge[1].size()/3; ++e) {
    const auto N = m_dsupedge[1].data() + e*3;
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      for (const auto& [p,q] : tk::lpoet) {
        tk::real alwp, alwn;
        if (g_cfg.get< tag::fctclip >()) {
          alwp = max( m_rhs(N[p],c), m_rhs(N[q],c) );
          alwn = min( m_rhs(N[p],c), m_rhs(N[q],c) );
        } else {
          alwp = max( max(m_rhs(N[p],c), m_u(N[p],c)),
                      max(m_rhs(N[q],c), m_u(N[q],c)) );
          alwn = min( min(m_rhs(N[p],c), m_u(N[p],c)),
                      min(m_rhs(N[q],c), m_u(N[q],c)) );
        }
        m_q(N[p],a) = max(m_q(N[p],a), alwp);
        m_q(N[p],b) = min(m_q(N[p],b), alwn);
        m_q(N[q],a) = max(m_q(N[q],a), alwp);
        m_q(N[q],b) = min(m_q(N[q],b), alwn);
      }
    }
  }

  // edges
  for (std::size_t e=0; e<m_dsupedge[2].size()/2; ++e) {
    const auto N = m_dsupedge[2].data() + e*2;
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      tk::real alwp, alwn;
      if (g_cfg.get< tag::fctclip >()) {
        alwp = max( m_rhs(N[0],c), m_rhs(N[1],c) );
        alwn = min( m_rhs(N[0],c), m_rhs(N[1],c) );
      } else {
        alwp = max( max(m_rhs(N[0],c), m_u(N[0],c)),
                    max(m_rhs(N[1],c), m_u(N[1],c)) );
        alwn = min( min(m_rhs(N[0],c), m_u(N[0],c)),
                    min(m_rhs(N[1],c), m_u(N[1],c)) );
      }
      m_q(N[0],a) = max(m_q(N[0],a), alwp);
      m_q(N[0],b) = min(m_q(N[0],b), alwn);
      m_q(N[1],a) = max(m_q(N[1],a), alwp);
      m_q(N[1],b) = min(m_q(N[1],b), alwn);
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
ChoCG::comalw( const std::unordered_map< std::size_t,
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
      auto a = c*2;
      auto b = a+1;
      q[a] = std::max( q[a], alw[a] );
      q[b] = std::min( q[b], alw[b] );
    }
  }

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nalw == Disc()->NodeCommMap().size()) {
    m_nalw = 0;
    comalw_complete();
  }
}

void
ChoCG::lim()
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
      auto a = c*2;
      auto b = a+1;
      m_q(i,a) = max( m_q(i,a), alw[a] );
      m_q(i,b) = min( m_q(i,b), alw[b] );
    }
  }
  tk::destroy(m_qc);

  // Finish computing allowed limits
  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      m_q(i,a) -= m_rhs(i,c);
      m_q(i,b) -= m_rhs(i,c);
    }
  }

  // Limit coefficients, C

  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      auto a = c*2;
      auto b = a+1;
      auto eps = std::numeric_limits< tk::real >::epsilon();
      m_q(i,a) = m_p(i,a) <  eps ? 0.0 : min(1.0, m_q(i,a)/m_p(i,a));
      m_q(i,b) = m_p(i,b) > -eps ? 0.0 : min(1.0, m_q(i,b)/m_p(i,b));
    }
  }

  // Limited antidiffusive contributions

  auto ctau = g_cfg.get< tag::fctdif >();
  m_a.fill( 0.0 );

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

  // tetrahedron superedges
  for (std::size_t e=0; e<m_dsupedge[0].size()/4; ++e) {
    const auto N = m_dsupedge[0].data() + e*4;
    const auto D = m_dsupint[0].data();
    auto C = m_dsuplim[0].data();
    std::size_t i = 0;
    for (const auto& [p,q] : tk::lpoed) {
      auto dif = D[(e*6+i)*5+3];
      auto coef = C + (e*6+i)*ncomp;
      tk::real aec[ncomp];
      for (std::size_t c=0; c<ncomp; ++c) {
        aec[c] = -dif * ctau * (m_u(N[p],c) - m_u(N[q],c));
        auto a = c*2;
        auto b = a+1;
        coef[c] = min( aec[c] < 0.0 ? m_q(N[p],a) : m_q(N[p],b),
                       aec[c] > 0.0 ? m_q(N[q],a) : m_q(N[q],b) );
      }
      tk::real cs = 1.0;
      for (auto c : fctsys) cs = min( cs, coef[c] );
      for (auto c : fctsys) coef[c] = cs;
      for (std::size_t c=0; c<ncomp; ++c) {
        aec[c] *= coef[c];
        m_a(N[p],c) -= aec[c];
        m_a(N[q],c) += aec[c];
      }
      ++i;
    }
  }

  // triangle superedges
  for (std::size_t e=0; e<m_dsupedge[1].size()/3; ++e) {
    const auto N = m_dsupedge[1].data() + e*3;
    const auto D = m_dsupint[1].data();
    auto C = m_dsuplim[0].data();
    std::size_t i = 0;
    for (const auto& [p,q] : tk::lpoet) {
      auto dif = D[(e*3+i)*5+3];
      auto coef = C + (e*3+i)*ncomp;
      tk::real aec[ncomp];
      for (std::size_t c=0; c<ncomp; ++c) {
        aec[c] = -dif * ctau * (m_u(N[p],c) - m_u(N[q],c));
        auto a = c*2;
        auto b = a+1;
        coef[c] = min( aec[c] < 0.0 ? m_q(N[p],a) : m_q(N[p],b),
                       aec[c] > 0.0 ? m_q(N[q],a) : m_q(N[q],b) );
      }
      tk::real cs = 1.0;
      for (auto c : fctsys) cs = min( cs, coef[c] );
      for (auto c : fctsys) coef[c] = cs;
      for (std::size_t c=0; c<ncomp; ++c) {
        aec[c] *= coef[c];
        m_a(N[p],c) -= aec[c];
        m_a(N[q],c) += aec[c];
      }
      ++i;
    }
  }

  // edges
  for (std::size_t e=0; e<m_dsupedge[2].size()/2; ++e) {
    const auto N = m_dsupedge[2].data() + e*2;
    const auto dif = m_dsupint[2][e*5+3];
    auto coef = m_dsuplim[2].data() + e*ncomp;
    tk::real aec[ncomp];
    for (std::size_t c=0; c<ncomp; ++c) {
      aec[c] = -dif * ctau * (m_u(N[0],c) - m_u(N[1],c));
      auto a = c*2;
      auto b = a+1;
      coef[c] = min( aec[c] < 0.0 ? m_q(N[0],a) : m_q(N[0],b),
                     aec[c] > 0.0 ? m_q(N[1],a) : m_q(N[1],b) );
    }
    tk::real cs = 1.0;
    for (auto c : fctsys) cs = min( cs, coef[c] );
    for (auto c : fctsys) coef[c] = cs;
    for (std::size_t c=0; c<ncomp; ++c) {
      aec[c] *= coef[c];
      m_a(N[0],c) -= aec[c];
      m_a(N[1],c) += aec[c];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif

  // Communicate limited antidiffusive contributions
  if (d->NodeCommMap().empty()){
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
ChoCG::comlim( const std::unordered_map< std::size_t,
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
ChoCG::BC( tk::Fields& u, tk::real t )
// *****************************************************************************
// Apply boundary conditions
//! \param[in,out] u Solution to apply BCs to
//! \param[in] t Physical time
// *****************************************************************************
{
  auto d = Disc();

  physics::dirbc( u, t, d->Coord(), d->BoxNodes(), m_dirbcmask, m_dirbcval );
  physics::symbc( u, m_symbcnodes, m_symbcnorms, /*pos=*/0 );
  physics::noslipbc( u, m_noslipbcnodes, /*pos=*/0 );
}

void
ChoCG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  auto d = Disc();
  const auto& vol = d->Vol();

  tk::real mindt = std::numeric_limits< tk::real >::max();
  auto const_dt = g_cfg.get< tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  // use constant dt if configured
  if (std::abs(const_dt) > eps) {

    // cppcheck-suppress redundantInitialization
    mindt = const_dt;

  } else {

    auto cfl = g_cfg.get< tag::cfl >();
    auto mu = g_cfg.get< tag::mat_dyn_viscosity >();
    auto large = std::numeric_limits< tk::real >::max();

    for (std::size_t i=0; i<m_u.nunk(); ++i) {
      auto u = m_u(i,0);
      auto v = m_u(i,1);
      auto w = m_u(i,2);
      auto vel = std::sqrt( u*u + v*v + w*w );
      auto L = std::cbrt( vol[i] );
      auto euler_dt = L / std::max( vel, 1.0e-8 );
      mindt = std::min( mindt, euler_dt );
      auto visc_dt = mu > eps ? L * L / mu : large;
      mindt = std::min( mindt, visc_dt );
    }
    mindt *= cfl;

  }

  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4rhs();
  thisProxy[ thisIndex ].wait4aec();
  thisProxy[ thisIndex ].wait4alw();
  thisProxy[ thisIndex ].wait4sol();
  thisProxy[ thisIndex ].wait4div();
  thisProxy[ thisIndex ].wait4sgrad();
  thisProxy[ thisIndex ].wait4step();

  // Contribute to minimum dt across all chares and advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(ChoCG,advance), thisProxy) );
}

void
ChoCG::advance( tk::real newdt )
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

  // Compute lhs and rhs of transport equations
  lhs();
  rhs();
}

void
ChoCG::lhs()
// *****************************************************************************
// Fill lhs matrix of transport equations
// *****************************************************************************
{
  auto theta = g_cfg.get< tag::theta >();
  if (theta < std::numeric_limits< tk::real >::epsilon()) return;

  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];
  const auto ncomp = m_u.nprop();
  const auto mu = g_cfg.get< tag::mat_dyn_viscosity >();

  auto dt = d->Dt();
  auto& A = Lhs();
  A.zero();

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::cross( ca, da );
    grad[2] = tk::cross( da, ba );
    grad[3] = tk::cross( ba, ca );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    for (std::size_t a=0; a<4; ++a) {
      for (std::size_t b=0; b<4; ++b) {
        auto v = J/dt/120.0 * ((a == b) ? 2.0 : 1.0);
        v += theta * mu * tk::dot(grad[a],grad[b]) / J / 6.0;
        for (std::size_t c=0; c<ncomp; ++c) A(N[a],N[b],c) -= v;
      }
    }
  }
}

void
ChoCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

  // Compute own portion of right-hand side for all equations
  auto dt = m_rkcoef[m_stage] * d->Dt();
  chorin::rhs( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel,
               dt, m_pr, m_u, m_vgrad, m_pgrad, m_rhs );

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
ChoCG::comrhs(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& inrhs )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] inrhs Partial contributions of RHS to chare-boundary nodes. Key:
//!   global mesh node IDs, value: contributions for all scalar components.
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped.
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
ChoCG::fct()
// *****************************************************************************
// Continue with flux-corrected transport if enabled
// *****************************************************************************
{
  // Combine own and communicated contributions to rhs
  const auto lid = Disc()->Lid();
  for (const auto& [g,r] : m_rhsc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_rhs(i,c) += r[c];
  }
  tk::destroy(m_rhsc);

  auto eps = std::numeric_limits< tk::real >::epsilon();
  if (g_cfg.get< tag::theta >() < eps and g_cfg.get< tag::fct >())
    aec();
  else
    solve();
}

void
// cppcheck-suppress unusedFunction
ChoCG::solve()
// *****************************************************************************
//  Advance systems of equations
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& vol = d->Vol();

  // Combine own and communicated contributions to limited antidiffusive
  // contributions
  for (const auto& [g,a] : m_ac) {
    auto i = tk::cref_find( d->Lid(), g );
    for (std::size_t c=0; c<a.size(); ++c) m_a(i,c) += a[c];
  }
  tk::destroy(m_ac);

  if (m_stage == 0) m_un = m_u;

  auto eps = std::numeric_limits<tk::real>::epsilon();
  if (g_cfg.get< tag::theta >() < eps || m_stage+1 < m_rkcoef.size()) {

    // Apply rhs in explicit solve
    if (g_cfg.get< tag::fct >()) {
      for (std::size_t i=0; i<npoin; ++i)
        for (std::size_t c=0; c<ncomp; ++c)
          m_a(i,c) = m_rhs(i,c) + m_a(i,c)/vol[i];
    }
    else {
      auto dt = m_rkcoef[m_stage] * d->Dt();
      for (std::size_t i=0; i<npoin; ++i)
        for (std::size_t c=0; c<ncomp; ++c)
          m_a(i,c) = m_un(i,c) - dt*m_rhs(i,c)/vol[i];
    }

    // Continue to advective-diffusive prediction
    pred();

  } else {

    // Configure momentum BCs
    std::unordered_map< std::size_t,
      std::vector< std::pair< int, tk::real > > > dirbc;

    if (m_np < 3) {       // only during startup and first time step
      // Configure Dirichlet BCs
      std::size_t nmask = ncomp + 1;
      Assert( m_dirbcmask.size() % nmask == 0, "Size mismatch" );
      for (std::size_t i=0; i<m_dirbcmask.size()/nmask; ++i) {
        auto p = m_dirbcmask[i*nmask+0];     // local node id
        auto& bc = dirbc[p];
        bc.resize( ncomp );
        for (std::size_t c=0; c<ncomp; ++c) {
          bc[c] = { m_dirbcmask[i*nmask+1+c], 0.0 };
        }
      }
      for (auto p : m_noslipbcnodes) {
        auto& bc = dirbc[p];
        bc.resize( ncomp );
        for (std::size_t c=0; c<ncomp; ++c) {
          bc[c] = { 1, 0.0 };
        }
      }
    }

    // Initialize semi-implicit momentum/transport solve
    m_cgmom[ thisIndex ].ckLocal()->
      init( {}, m_rhs.vec(), {}, dirbc, m_np<3,  g_cfg.get< tag::mom_pc >(),
            CkCallback( CkIndex_ChoCG::msolve(), thisProxy[thisIndex] ) );

  }
}

void
ChoCG::msolve()
// *****************************************************************************
//  Solve for momentum/transport system of equations
// *****************************************************************************
{
  auto iter = g_cfg.get< tag::mom_iter >();
  auto tol = g_cfg.get< tag::mom_tol >();
  auto verbose = g_cfg.get< tag::mom_verbose >();

  m_cgmom[ thisIndex ].ckLocal()->solve( iter, tol, thisIndex, verbose,
    CkCallback( CkIndex_ChoCG::msolved(), thisProxy[thisIndex] ) );
}

void
ChoCG::msolved()
// *****************************************************************************
// Continue after momentum/transport solve in semi-implcit solve
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();

  if (thisIndex == 0) d->mit( m_cgmom[ thisIndex ].ckLocal()->it() );

  // Update momentum/transport solution in semi-implicit solve
  auto& du = m_cgmom[ thisIndex ].ckLocal()->solution();
  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      m_a(i,c) = m_un(i,c) + du[i*ncomp+c];
    }
  }

  // Continue to advective-diffusive prediction
  pred();
}

void
ChoCG::pred()
// *****************************************************************************
//  Compute advective-diffusive prediction of momentum/transport
// *****************************************************************************
{
  auto d = Disc();

  // Configure and apply scalar source to solution (if defined)
  auto src = problems::PHYS_SRC();
  if (src) src( d->Coord(), d->T(), m_a );

  // Enforce boundary conditions
  BC( m_a, d->T() + m_rkcoef[m_stage] * d->Dt() );

  // Update momentum/transport solution
  m_u = m_a;
  m_a.fill( 0.0 );

  // Compute velocity gradients if needed
  if (g_cfg.get< tag::flux >() == "damp4") {
    thisProxy[ thisIndex ].wait4vgrad();
    velgrad();
  } else {
    corr();
  }
}

void
ChoCG::corr()
// *****************************************************************************
//  Compute pressure correction
// *****************************************************************************
{
  // Finalize computing velocity gradients
  if (g_cfg.get< tag::flux >() == "damp4") fingrad( m_vgrad, m_vgradc );

  if (++m_stage < m_rkcoef.size()) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4rhs();
    thisProxy[ thisIndex ].wait4aec();
    thisProxy[ thisIndex ].wait4alw();
    thisProxy[ thisIndex ].wait4sol();
    // Continue to next time stage of momentum/transport prediction
    rhs();

  } else {

    // Reset Runge-Kutta stage counter
    m_stage = 0;
    // Continue to pressure correction and projection
    div( m_u );

  }
}

void
ChoCG::diag()
// *****************************************************************************
//  Compute diagnostics
// *****************************************************************************
{
  auto d = Disc();

  // Increase number of iterations and physical time
  d->next();

  // Compute diagnostics, e.g., residuals
  auto diag_iter = g_cfg.get< tag::diag_iter >();
  const auto& dp = m_cgpre[ thisIndex ].ckLocal()->solution();
  auto diagnostics = m_diag.precompute( *d, m_u, m_un, m_pr, dp, diag_iter );

  // Evaluate residuals
  if (!diagnostics) evalres( std::vector< tk::real >( m_u.nprop(), 1.0 ) );
}

void
ChoCG::evalres( const std::vector< tk::real >& )
// *****************************************************************************
//  Evaluate residuals
// *****************************************************************************
{
  refine();
}

void
ChoCG::refine()
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
ChoCG::resizePostAMR(
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
  //m_pr.rm( removedNodes );
  m_rhs.rm( removedNodes );

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  m_u.resize( npoin );
  m_pr.resize( npoin );
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
ChoCG::writeFields( CkCallback cb )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] cb Function to continue with after the write
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) { cb.send(); return; }

  auto d = Disc();

  // Field output

  std::vector< std::string > nodefieldnames{
    "velocityx", "velocityy", "velocityz", "divergence", "pressure" };

  std::vector< std::vector< tk::real > > nodefields;

  nodefields.push_back( m_u.extract(0) );
  nodefields.push_back( m_u.extract(1) );
  nodefields.push_back( m_u.extract(2) );

  // Divide weak result by nodal volume
  const auto& vol = d->Vol();
  for (std::size_t i=0; i<m_div.size(); ++i) m_div[i] /= vol[i];
  nodefields.push_back( m_div );

  nodefields.push_back( m_pr ) ;

  //nodefieldnames.push_back( "dp/dx" );
  //nodefieldnames.push_back( "dp/dy" );
  //nodefieldnames.push_back( "dp/dz" );
  //nodefields.push_back( m_pgrad.extract(0) );
  //nodefields.push_back( m_pgrad.extract(1) );
  //nodefields.push_back( m_pgrad.extract(2) );

  //nodefieldnames.push_back( "fx" );
  //nodefieldnames.push_back( "fy" );
  //nodefieldnames.push_back( "fz" );
  //nodefields.push_back( m_flux.extract(0) );
  //nodefields.push_back( m_flux.extract(1) );
  //nodefields.push_back( m_flux.extract(2) );

  //nodefieldnames.push_back( "du/dx" );
  //nodefieldnames.push_back( "du/dy" );
  //nodefieldnames.push_back( "du/dz" );
  //nodefieldnames.push_back( "dv/dx" );
  //nodefieldnames.push_back( "dv/dy" );
  //nodefieldnames.push_back( "dv/dz" );
  //nodefieldnames.push_back( "dw/dx" );
  //nodefieldnames.push_back( "dw/dy" );
  //nodefieldnames.push_back( "dw/dz" );
  //nodefields.push_back( m_vgrad.extract(0) );
  //nodefields.push_back( m_vgrad.extract(1) );
  //nodefields.push_back( m_vgrad.extract(2) );
  //nodefields.push_back( m_vgrad.extract(3) );
  //nodefields.push_back( m_vgrad.extract(4) );
  //nodefields.push_back( m_vgrad.extract(5) );
  //nodefields.push_back( m_vgrad.extract(6) );
  //nodefields.push_back( m_vgrad.extract(7) );
  //nodefields.push_back( m_vgrad.extract(8) );

  auto ncomp = m_u.nprop();
  for (std::size_t c=0; c<ncomp-3; ++c) {
    nodefieldnames.push_back( "c" + std::to_string(c) );
    nodefields.push_back( m_u.extract(3+c) );
  }

  // query function to evaluate analytic solution (if defined)
  auto sol = problems::SOL();

  if (sol) {
    const auto& coord = d->Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    auto an = m_u;
    for (std::size_t i=0; i<an.nunk(); ++i) {
      auto s = sol( x[i], y[i], z[i], d->T() );
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
    }
    nodefieldnames.push_back( "velocity_analyticx" );
    nodefields.push_back( an.extract(0) );
    nodefieldnames.push_back( "velocity_analyticy" );
    nodefields.push_back( an.extract(1) );
    nodefieldnames.push_back( "velocity_analyticz" );
    nodefields.push_back( an.extract(2) );
    for (std::size_t c=0; c<ncomp-3; ++c) {
      nodefieldnames.push_back( nodefieldnames[5+c] + "_analytic" );
      nodefields.push_back( an.extract(3+c) );
    }
  }

  // also output analytic pressure (if defined)
  auto pressure_sol = problems::PRESSURE_SOL();
  if (pressure_sol) {
    const auto& coord = d->Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    auto ap = m_pr;
    for (std::size_t i=0; i<ap.size(); ++i) {
      ap[i] = pressure_sol( x[i], y[i], z[i] );
    }
    nodefieldnames.push_back( "analytic" );
    nodefields.push_back( ap );
  }

  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Surface output

  std::vector< std::string > nodesurfnames;
  std::vector< std::vector< tk::real > > nodesurfs;

  const auto& f = g_cfg.get< tag::fieldout >();

  if (!f.empty()) {
    std::size_t nc = 5;
    nodesurfnames.push_back( "velocityx" );
    nodesurfnames.push_back( "velocityy" );
    nodesurfnames.push_back( "velocityz" );
    nodesurfnames.push_back( "divergence" );
    nodesurfnames.push_back( "pressure" );

    auto bnode = tk::bfacenodes( m_bface, m_triinpoel );
    std::set< int > outsets( begin(f), end(f) );
    for (auto sideset : outsets) {
      auto b = bnode.find(sideset);
      if (b == end(bnode)) continue;
      const auto& nodes = b->second;
      auto i = nodesurfs.size();
      nodesurfs.insert( end(nodesurfs), nc,
                        std::vector< tk::real >( nodes.size() ) );
      std::size_t j = 0;
      for (auto n : nodes) {
        const auto s = m_u[n];
        std::size_t p = 0;
        nodesurfs[i+(p++)][j] = s[0];
        nodesurfs[i+(p++)][j] = s[1];
        nodesurfs[i+(p++)][j] = s[2];
        nodesurfs[i+(p++)][j] = m_div[n];
        nodesurfs[i+(p++)][j] = m_pr[n];
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
ChoCG::out()
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
    writeFields( CkCallback(CkIndex_ChoCG::integrals(), thisProxy[thisIndex]) );
  } else {
    integrals();
  }
}

void
ChoCG::integrals()
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
          mfr += n[0]*m_u(p,0) + n[1]*m_u(p,1) + n[2]*m_u(p,2);
          n += 3;
        }
      }
    }
    if (req.count("force")) {
      auto mu = g_cfg.get< tag::mat_dyn_viscosity >();
      for (const auto& [s,sint] : m_surfint) {
        auto& fx = ints[ FORCE_X ][ s ];
        auto& fy = ints[ FORCE_Y ][ s ];
        auto& fz = ints[ FORCE_Z ][ s ];
        const auto& nodes = sint.first;
        const auto& ndA = sint.second;
        auto n = ndA.data();
        for (auto p : nodes) {
          // pressure force
          fx -= n[0]*m_pr[p];
          fy -= n[1]*m_pr[p];
          fz -= n[2]*m_pr[p];
          // viscous force
          fx += mu*(m_vgrad(p,0)*n[0] + m_vgrad(p,1)*n[1] + m_vgrad(p,2)*n[2]);
          fy += mu*(m_vgrad(p,3)*n[0] + m_vgrad(p,4)*n[1] + m_vgrad(p,5)*n[2]);
          fz += mu*(m_vgrad(p,6)*n[0] + m_vgrad(p,7)*n[1] + m_vgrad(p,8)*n[2]);
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
ChoCG::evalLB( int nrestart )
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
ChoCG::evalRestart()
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
ChoCG::step()
// *****************************************************************************
// Evaluate whether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  if(thisIndex == 0) d->status();

  if (not m_finished) {

    evalRestart();

  } else {

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/chocg.def.h"
