// *****************************************************************************
/*!
  \file      src/Inciter/LohCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LohCG: Artificial compressibility solver for incompressible flow
*/
// *****************************************************************************

#include "XystBuildConfig.hpp"
#include "LohCG.hpp"
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
#include "Lohner.hpp"
#include "Problems.hpp"
#include "EOS.hpp"
#include "BC.hpp"
#include "Print.hpp"

namespace inciter {

extern ctr::Config g_cfg;

static CkReduction::reducerType IntegralsMerger;

//! Runge-Kutta coefficients
static const std::array< std::vector< tk::real >, 4 > rkcoef{{
  { 1.0 },
  { 1.0/2.0, 1.0 },
  { 1.0/3.0, 1.0/2.0, 1.0 },
  { 1.0/4.0, 1.0/3.0, 1.0/2.0, 1.0 }
}};

} // inciter::

using inciter::g_cfg;
using inciter::LohCG;

LohCG::LohCG( const CProxy_Discretization& disc,
              const tk::CProxy_ConjugateGradients& cgpre,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel )
try :
  m_disc( disc ),
  m_cgpre( cgpre ),
  m_nrhs( 0 ),
  m_nnorm( 0 ),
  m_ngrad( 0 ),
  m_nsgrad( 0 ),
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
  m_grad( m_u.nunk(), ngradcomp() ),
  m_vgrad( m_u.nunk(), 9UL ),
  m_flux( m_u.nunk(), 3UL ),
  m_div( m_u.nunk() ),
  m_sgrad( m_u.nunk(), 3UL ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_stage( 0 ),
  m_finished( 0 ),
  m_rkcoef( rkcoef[ g_cfg.get< tag::rk >() - 1 ] )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] cgpre ConjugateGradients Charm++ proxy for initial pressure solve
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

  // Activate SDAG waits for setup
  thisProxy[ thisIndex ].wait4int();

} // Catch std::exception
  catch (std::exception& se) {
    // (re-)throw tk::Excpetion
    Throw( std::string("RUNTIME ERROR in CSR constructor: ") + se.what() );
  }


std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
LohCG::prelhs( const std::pair< std::vector< std::size_t >,
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

std::size_t
LohCG::ngradcomp() const
// *****************************************************************************
//  Compute number of scalar components for gradients
//! \return Number scalar components required for gradients
// *****************************************************************************
{
  std::size_t n = 0;
  const auto& req = g_cfg.get< tag::integout, tag::integrals >();

  if (g_cfg.get< tag::flux >() == "damp4" or
      std::find( begin(req), end(req), "force") != end(req))
  {
    n += m_u.nprop() * 3;     // (p,u,v,w,c0,...) x 3
  }

  return n;
}

void
LohCG::setupDirBC( const std::vector< std::vector< int > >& cfgmask,
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
        if (mv.second.empty()) mv.second = v;
        auto& mval = mv.first;
        if (mval.empty()) mval.resize( ncomp, 0 );
        for (std::size_t c=0; c<ncomp; ++c) {
          if (!mval[c]) mval[c] = vec[c+1];  // overwrite mask if 0 -> nonzero
          //mval[c] = std::max( mval[c], vec[c+1] );
          //std::cout << p << ": " << mval[c] << '\n';
        }
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
LohCG::feop()
// *****************************************************************************
// Start (re-)computing finite element domain and boundary operators
// *****************************************************************************
{
  auto d = Disc();

  bool multi = g_cfg.get< tag::input >().size() > 1;
  auto meshid = d->MeshId();

  // Prepare Dirichlet boundary conditions data structures
  const auto& bc_dir =
    multi ? g_cfg.get< tag::bc_dir_ >()[ meshid ]
          : g_cfg.get< tag::bc_dir >();
  const auto& bc_dirval =
    multi ? g_cfg.get< tag::bc_dirval_ >()[ meshid ]
          : g_cfg.get< tag::bc_dirval >();

  setupDirBC( bc_dir, bc_dirval, m_u.nprop(), m_dirbcmask, m_dirbcval );

  // Prepare pressure Dirichlet boundary conditions data structures
  const auto& tp =
    multi ? g_cfg.get< tag::pressure_ >()[ meshid ]
          : g_cfg.get< tag::pressure >();
  setupDirBC( tp.get< tag::bc_dir >(), tp.get< tag::bc_dirval >(),
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
LohCG::bndint()
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
LohCG::domint()
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
      n[3] += tk::dot( grad[p], grad[q] ) / J / 6.0;
    }
  }
}

void
LohCG::comnorm( const decltype(m_bnorm)& inbnd )
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
LohCG::registerReducers()
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
LohCG::ResumeFromSync()
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
LohCG::setup( tk::real v )
// *****************************************************************************
// Start setup for solution
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions
  problems::initialize( d->Coord(), m_u, d->T(), d->MeshId(), d->BoxNodes() );

  // Query time history field output labels from all PDEs integrated
  if (!g_cfg.get< tag::histout, tag::points >().empty()) {
    std::vector< std::string > var
      {"density", "xvelocity", "yvelocity", "zvelocity", "energy", "pressure"};
    auto ncomp = m_u.nprop();
    for (std::size_t c=5; c<ncomp; ++c)
      var.push_back( "c" + std::to_string(c-5) );
    d->histheader( std::move(var) );
  }

  // Setup hole data structures (if coupled)
  d->hole( m_bface, m_triinpoel,
    CkCallback( CkIndex_LohCG::transferFL(), thisProxy[thisIndex] ) );

  // Compute finite element operators
  feop();
}

void
LohCG::transferFL()
// *****************************************************************************
// Initiate transfer of transfer flags (if coupled)
// *****************************************************************************
{
  auto d = Disc();

  // Prepare integrid-boundary data structures
  d->intergrid( m_bnode );

  // Find mesh nodes within holes
  d->holefind();

  // Initiate transfer of transfer flags
  auto c = CkCallback(CkIndex_LohCG::transfer_complete(), thisProxy[thisIndex]);
  d->transfer( m_u, c, /* trflag = */ true );
}

void
LohCG::bnorm()
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
LohCG::prep_surfint()
// *****************************************************************************
// Prepare surface integral data strurctures
// *****************************************************************************
{
  // Query surface integral output nodes
  std::unordered_map< int, std::vector< std::size_t > > surfintnodes;
  const auto& is = g_cfg.get< tag::integout, tag::sidesets >();
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
}

void
LohCG::prep_symbc()
// *****************************************************************************
//  Prepare symmetry boundary condition data structures
// *****************************************************************************
{
  bool multi = g_cfg.get< tag::input >().size() > 1;
  auto meshid = Disc()->MeshId();

  const auto& bc_sym =
    multi ? g_cfg.get< tag::bc_sym_ >()[ meshid ]
          : g_cfg.get< tag::bc_sym >();

  // Query symmetry BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > sym;
  for (auto s : bc_sym) {
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
    for (const auto& s : bc_sym) {
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
}

void
LohCG::prep_noslipbc()
// *****************************************************************************
//  Prepare no-slip boundary condition data structures
// *****************************************************************************
{
  bool multi = g_cfg.get< tag::input >().size() > 1;
  auto meshid = Disc()->MeshId();

  const auto& bc_noslip =
    multi ? g_cfg.get< tag::bc_noslip_ >()[ meshid ]
          : g_cfg.get< tag::bc_noslip >();

  // Query noslip BC nodes associated to side sets
  std::unordered_map< int, std::unordered_set< std::size_t > > noslip;
  for (auto s : bc_noslip) {
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
LohCG::streamable()
// *****************************************************************************
// Convert integrals into streamable data structures
// *****************************************************************************
{
  // Prepare surface integral data strurctures
  prep_surfint();

  // Generate domain superedges
  domsuped();

  // Prepare symmetry boundary condition data structures
  prep_symbc();

  //  Prepare no-slip boundary condition data structures
  prep_noslipbc();
}

void
LohCG::domsuped()
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
        m_domedgeint.erase( d[ed] );
      }
    }
  }

  m_dsupedge[2].resize( m_domedgeint.size()*2 );
  m_dsupint[2].resize( m_domedgeint.size()*4 );
  std::size_t k = 0;
  for (const auto& [ed,d] : m_domedgeint) {
    auto e = m_dsupedge[2].data() + k*2;
    e[0] = tk::cref_find( lid, ed[0] );
    e[1] = tk::cref_find( lid, ed[1] );
    auto i = m_dsupint[2].data() + k*4;
    i[0] = d[0];
    i[1] = d[1];
    i[2] = d[2];
    i[3] = d[3];
    ++k;
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
LohCG::holeset()
// *****************************************************************************
// Set solution in holes (if coupled)
// *****************************************************************************
{
  bool multi = g_cfg.get< tag::input >().size() > 1;
  if (!multi) return;

  auto d = Disc();
  if (d->MeshId()) return;

  const auto& flag = d->TransferFlag();
  Assert( flag.size() == m_u.nunk(), "Size mismatch" );
  //auto eps = std::numeric_limits< tk::real >::epsilon();

  for (std::size_t i=0; i<m_u.nunk(); ++i) {
    if (flag[i] > 0.0) {  // zero solution in holes
      for (std::size_t c=0; c<m_u.nprop(); ++c) m_u(i,c) = 0.0;
    }
  }
}

void
// cppcheck-suppress unusedFunction
LohCG::merge()
// *****************************************************************************
// Combine own and communicated portions of the integrals
// *****************************************************************************
{
//std::cout << Disc()->MeshId() << ": " << thisIndex << ":\n";
  // Combine own and communicated contributions to boundary point normals
  bnorm();

  // Convert integrals into streamable data structures
  streamable();

  // Enforce boundary conditions on initial conditions
  auto d = Disc();
  auto t = d->T() + d->Dt();
  physics::dirbc( d->MeshId(), m_u, t, d->Coord(), d->BoxNodes(), m_dirbcmask,
                  m_dirbcval );
  physics::symbc( m_u, m_symbcnodes, m_symbcnorms, /*pos=*/1 );
  physics::noslipbc( m_u, m_noslipbcnodes, /*pos=*/1 );

  // Set solution in holes (if coupled) after initial conditions transfer
  holeset();

  // Start measuring initial div-free time
  m_timer.emplace_back();

  // Compute initial momentum flux
  thisProxy[ thisIndex ].wait4div();
  thisProxy[ thisIndex ].wait4sgrad();
  div( m_u, /*pos=*/1 );
}

void
LohCG::fingrad( tk::Fields& grad,
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
LohCG::div( const tk::Fields& u, std::size_t pos )
// *****************************************************************************
//  Start computing divergence
//! \param[in] u Vector field whose divergence to compute
//! \param[in] pos Position at which the three vector components are in u
// *****************************************************************************
{
  auto d = Disc();
  const auto lid = d->Lid();

  // Finalize momentum flux communications if needed
  if (m_np == 1) {
    fingrad( m_flux, m_fluxc );
    physics::symbc( m_flux, m_symbcnodes, m_symbcnorms, /*pos=*/0 );
  }

  // Compute velocity divergence
  std::fill( begin(m_div), end(m_div), 0.0 );
  lohner::div( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, u, m_div, pos );

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
LohCG::comdiv( const std::unordered_map< std::size_t, tk::real >& indiv )
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
LohCG::velgrad()
// *****************************************************************************
//  Start computing velocity gradient
// *****************************************************************************
{
  auto d = Disc();

  // Compute momentum flux
  m_vgrad.fill( 0.0 );
  lohner::vgrad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_u, m_vgrad );

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
LohCG::comvgrad(
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
LohCG::flux()
// *****************************************************************************
//  Start computing momentum flux
// *****************************************************************************
{
  auto d = Disc();

  // Finalize computing velocity gradients
  fingrad( m_vgrad, m_vgradc );

  // Compute momentum flux
  m_flux.fill( 0.0 );
  lohner::flux( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_u, m_vgrad,
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
LohCG::comflux(
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
LohCG::pinit()
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
  const auto& tp = g_cfg.get< tag::pressure >();

  // Combine own and communicated contributions to velocity divergence
  for (const auto& [g,r] : m_divc) m_div[ tk::cref_find(lid,g) ] += r;
  tk::destroy(m_divc);

  // Configure Dirichlet BCs
  std::unordered_map< std::size_t,
    std::vector< std::pair< int, tk::real > > > dirbc;
  if (!tp.get< tag::bc_dir >().empty()) {
    auto ic = problems::PRESSURE_IC();
    std::size_t nmask = 1 + 1;
    Assert( m_dirbcmaskp.size() % nmask == 0, "Size mismatch" );
    for (std::size_t i=0; i<m_dirbcmaskp.size()/nmask; ++i) {
      auto p = m_dirbcmaskp[i*nmask+0];     // local node id
      auto mask = m_dirbcmaskp[i*nmask+1];
      if (mask == 1) {                                  // mask == 1: IC value
        auto val = ic( x[p], y[p], z[p], /*meshid=*/0 );
        dirbc[p] = {{ { 1, val } }};
      } else if (mask == 2 && !m_dirbcvalp.empty()) {   // mask == 2: BC value
        auto val = m_dirbcvalp[i*nmask+1];
        dirbc[p] = {{ { 1, val } }};
      }
    }
  }

  // Configure Neumann BCs
  std::vector< tk::real > neubc;
  auto pg = problems::PRESSURE_GRAD();
  if (pg) {
    // Collect Neumann BC elements
    std::vector< std::uint8_t > besym( m_triinpoel.size(), 0 );
    for (auto s : tp.get< tag::bc_sym >()) {
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
  auto h = tp.get< tag::hydrostat >();
  if (h != std::numeric_limits< uint64_t >::max()) {
    auto pi = lid.find( h );
    if (pi != end(lid)) {
      auto p = pi->second;
      auto ic = problems::PRESSURE_IC();
      auto val = m_np>1 ? 0.0 : ic( x[p], y[p], z[p], /*meshid=*/0 );
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

  // Initialize Poisson solve
  const auto& pc = tp.get< tag::pc >();
  m_cgpre[ thisIndex ].ckLocal()->init( {}, m_div, neubc, dirbc, true, pc,
    CkCallback( CkIndex_LohCG::psolve(), thisProxy[thisIndex] ) );
}

void
LohCG::psolve()
// *****************************************************************************
//  Solve Poisson equation
// *****************************************************************************
{
  const auto& tp = g_cfg.get< tag::pressure >();
  auto iter = tp.get< tag::iter >();
  auto tol = tp.get< tag::tol >();
  auto verbose = tp.get< tag::verbose >();

  auto c = m_np != 1 ?
           CkCallback( CkIndex_LohCG::sgrad(), thisProxy[thisIndex] ) :
           CkCallback( CkIndex_LohCG::psolved(), thisProxy[thisIndex] );

  m_cgpre[ thisIndex ].ckLocal()->solve( iter, tol, thisIndex, verbose, c );
}

void
LohCG::sgrad()
// *****************************************************************************
// Compute recent conjugate gradients solution gradient
// *****************************************************************************
{
  auto d = Disc();

  auto sol = m_cgpre[ thisIndex ].ckLocal()->solution();
  m_sgrad.fill( 0.0 );
  lohner::grad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, sol, m_sgrad );

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
LohCG::comsgrad(
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
LohCG::psolved()
// *****************************************************************************
// Continue setup after Poisson solve and gradient computation
// *****************************************************************************
{
  auto d = Disc();

  if (thisIndex == 0) d->pit( m_cgpre[ thisIndex ].ckLocal()->it() );

  const auto& problem = inciter::g_cfg.get< tag::problem >();
  bool test_overset = problem.find("overset") != std::string::npos;

  if (m_np != 1 and !test_overset and !d->MeshId()) {
    // Finalize gradient communications
    fingrad( m_sgrad, m_sgradc );
    // Project velocity to divergence-free subspace
    for (std::size_t i=0; i<m_u.nunk(); ++i) {
      m_u(i,1) -= m_sgrad(i,0);
      m_u(i,2) -= m_sgrad(i,1);
      m_u(i,3) -= m_sgrad(i,2);
    }
    // Enforce boundary conditions
    auto t = d->T() + d->Dt();
    physics::dirbc( d->MeshId(), m_u, t, d->Coord(), d->BoxNodes(), m_dirbcmask,
                    m_dirbcval);
    physics::symbc( m_u, m_symbcnodes, m_symbcnorms, /*pos=*/1 );
    physics::noslipbc( m_u, m_noslipbcnodes, /*pos=*/1 );
    // Set solution in holes (if coupled)
    holeset();
  }

  // Initiate transfer of initial conditions (if coupled)
  auto c = CkCallback( CkIndex_LohCG::transferIC(), thisProxy[thisIndex] );
  d->transfer( m_u, c, /* trflag = */ false );
}

void
LohCG::transferIC()
// *****************************************************************************
// Initiate transfer of initial conditions (if coupled)
// *****************************************************************************
{
  auto d = Disc();

  if (g_cfg.get< tag::nstep >() == 1) {  // test first Poisson solve only

    auto p = m_cgpre[ thisIndex ].ckLocal()->solution();
    for (std::size_t i=0; i<m_u.nunk(); ++i) m_u(i,0) = p[i];
    thisProxy[ thisIndex ].wait4step();
    writeFields( CkCallback(CkIndex_LohCG::diag(), thisProxy[thisIndex]) );

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
      // Assign initial pressure and start timestepping
      if (!d->MeshId()) {
        auto p = m_cgpre[ thisIndex ].ckLocal()->solution();
        for (std::size_t i=0; i<m_u.nunk(); ++i) m_u(i,0) = p[i];
        holeset();      // set solution in holes (if coupled)
      }
      writeFields( CkCallback(CkIndex_LohCG::start(), thisProxy[thisIndex]) );
    }

  }
}

void
LohCG::diag()
// *****************************************************************************
//  Compute diagnostics
// *****************************************************************************
{
  auto d = Disc();

  // Increase number of iterations and physical time
  d->next();

  // Compute diagnostics, e.g., residuals
  auto diag_iter = g_cfg.get< tag::diag_iter >();
  auto diagnostics = m_diag.accompute( *d, m_u, m_un, diag_iter );

  // Evaluate residuals
  if (!diagnostics) evalres( std::vector< tk::real >( m_u.nprop(), 1.0 ) );
}

void
LohCG::start()
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
LohCG::dt()
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
    auto large = std::numeric_limits< tk::real >::max();
    auto c = g_cfg.get< tag::soundspeed >();
    auto mu = g_cfg.get< tag::mat_dyn_viscosity >();
    auto dif = g_cfg.get< tag::mat_dyn_diffusivity >();
    dif = std::max( mu, dif );

    for (std::size_t i=0; i<m_u.nunk(); ++i) {
      auto u = m_u(i,1);
      auto v = m_u(i,2);
      auto w = m_u(i,3);
      auto vel = std::sqrt( u*u + v*v + w*w );
      auto L = std::cbrt( vol[i] );
      auto euler_dt = L / std::max( vel+c, 1.0e-8 );
      mindt = std::min( mindt, euler_dt );
      auto visc_dt = dif > eps ? L * L / dif : large;
      mindt = std::min( mindt, visc_dt );
    }
    mindt *= cfl;

  }

  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4step();

  // Contribute to minimum dt across all chares (and meshes if coupled)
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(Transporter,transfer_dt), d->Tr()) );
}

void
LohCG::advance( tk::real newdt )
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

  // Start next time step stage
  stage();
}

void
LohCG::stage()
// *****************************************************************************
// Start next time step stage
// *****************************************************************************
{
  // Activate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  // Compute gradients
  grad();
}

void
LohCG::grad()
// *****************************************************************************
//  Compute gradients
// *****************************************************************************
{
  auto d = Disc();

  if (m_grad.nprop()) {
    m_grad.fill( 0.0 );
    lohner::grad( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, m_u, m_grad );
  }

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
LohCG::comgrad(
  const std::unordered_map< std::size_t, std::vector< tk::real > >& ingrad )
// *****************************************************************************
//  Receive contributions to gradients on chare-boundaries
//! \param[in] ingrad Partial contributions to chare-boundary nodes. Key:
//!   global mesh node IDs, value: contributions for all scalar components.
//! \details This function receives contributions to m_grad, which stores the
//!   gradients at mesh nodes. While m_grad stores own contributions, m_gradc
//!   collects the neighbor chare contributions during communication. This way
//!   work on m_grad and m_gradc is overlapped.
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
LohCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

  // Combine own and communicated contributions to gradients
  if (m_grad.nprop()) fingrad( m_grad, m_gradc );

  // Compute own portion of right-hand side for all equations
  lohner::rhs( m_dsupedge, m_dsupint, d->Coord(), m_triinpoel, d->V(), d->T(),
               m_u, m_grad, m_rhs );

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
LohCG::comrhs(
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
// cppcheck-suppress unusedFunction
LohCG::solve()
// *****************************************************************************
//  Advance systems of equations
// *****************************************************************************
{
  auto d = Disc();
  const auto npoin = m_u.nunk();
  const auto ncomp = m_u.nprop();
  const auto& vol = d->Vol();

  // Combine own and communicated contributions to rhs
  const auto lid = d->Lid();
  for (const auto& [g,r] : m_rhsc) {
    auto i = tk::cref_find( lid, g );
    for (std::size_t c=0; c<r.size(); ++c) m_rhs(i,c) += r[c];
  }
  tk::destroy(m_rhsc);

  if (m_stage == 0) m_un = m_u;

  // Advance system
  auto dt = m_rkcoef[m_stage] * d->Dt();
  for (std::size_t i=0; i<npoin; ++i) {
    for (std::size_t c=0; c<ncomp; ++c) {
      m_u(i,c) = m_un(i,c) - dt*m_rhs(i,c)/vol[i];
    }
  }

  // Configure and apply scalar source to solution (if defined)
  auto src = problems::PHYS_SRC();
  if (src) src( d->Coord(), d->T(), m_u );

  // Enforce boundary conditions
  auto t = d->T() + m_rkcoef[m_stage] * d->Dt();
  auto meshid = d->MeshId();
  physics::dirbc( meshid, m_u, t, d->Coord(), d->BoxNodes(), m_dirbcmask,
                  m_dirbcval );
  physics::dirbcp( meshid, m_u, d->Coord(), m_dirbcmaskp, m_dirbcvalp );
  physics::symbc( m_u, m_symbcnodes, m_symbcnorms, /*pos=*/1 );
  physics::noslipbc( m_u, m_noslipbcnodes, /*pos=*/1 );

  // Initiate transfer of updated solution (if coupled)
  auto c = CkCallback( CkIndex_LohCG::solved(), thisProxy[thisIndex] );
  d->transfer( m_u, c, /* trflag = */ false );
}

void
LohCG::solved()
// *****************************************************************************
//  Solution has been updated
// *****************************************************************************
{
  auto d = Disc();
  auto t = d->T() + m_rkcoef[m_stage] * d->Dt();
  auto meshid = d->MeshId();
  physics::dirbc( meshid, m_u, t, d->Coord(), d->BoxNodes(), m_dirbcmask,
                  m_dirbcval );
  physics::dirbcp( meshid, m_u, d->Coord(), m_dirbcmaskp, m_dirbcvalp );
  physics::symbc( m_u, m_symbcnodes, m_symbcnorms, /*pos=*/1 );
  physics::noslipbc( m_u, m_noslipbcnodes, /*pos=*/1 );

  // Set solution in holes (if coupled)
  holeset();

  if (++m_stage < m_rkcoef.size()) {

    // Start next time step stage
    stage();

  } else {

    // Reset Runge-Kutta stage counter
    m_stage = 0;
    // Compute diagnostics, e.g., residuals
    diag();

  }
}

void
LohCG::evalres( const std::vector< tk::real >& )
// *****************************************************************************
//  Evaluate residuals
// *****************************************************************************
{
  refine();
}

void
LohCG::refine()
// *****************************************************************************
// Optionally refine/derefine mesh
// *****************************************************************************
{
  auto d = Disc();

  // See if this is the last time step
  if (d->finished()) m_finished = 1;

  const auto& ht = g_cfg.get< tag::href >();
  auto dtref = ht.get< tag::dt >();
  auto dtfreq = ht.get< tag::dtfreq >();

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
LohCG::resizePostAMR(
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
LohCG::writeFields( CkCallback cb )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] cb Function to continue with after the write
// *****************************************************************************
{
  if (g_cfg.get< tag::benchmark >()) { cb.send(); return; }

  auto d = Disc();

  // Field output

  std::vector< std::string > nodefieldnames{
    "pressure", "velocityx", "velocityy", "velocityz" };

  std::vector< std::vector< tk::real > > nodefields;

  nodefields.push_back( m_u.extract(0) );
  nodefields.push_back( m_u.extract(1) );
  nodefields.push_back( m_u.extract(2) );
  nodefields.push_back( m_u.extract(3) );

  auto ncomp = m_u.nprop();
  for (std::size_t c=0; c<ncomp-4; ++c) {
    nodefieldnames.push_back( "c" + std::to_string(c) );
    nodefields.push_back( m_u.extract(4+c) );
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
      auto s = sol( x[i], y[i], z[i], d->T(), /*meshid=*/0 );
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
    }
    nodefieldnames.push_back( "pressure_analytic" );
    nodefields.push_back( an.extract(0) );
    nodefieldnames.push_back( "velocity_analyticx" );
    nodefields.push_back( an.extract(1) );
    nodefieldnames.push_back( "velocity_analyticy" );
    nodefields.push_back( an.extract(2) );
    nodefieldnames.push_back( "velocity_analyticz" );
    nodefields.push_back( an.extract(3) );
    for (std::size_t c=0; c<ncomp-4; ++c) {
      nodefieldnames.push_back( nodefieldnames[4+c] + "_analytic" );
      nodefields.push_back( an.extract(4+c) );
    }
  }

  bool multi = g_cfg.get< tag::input >().size() > 1;
  if (multi) {
    nodefieldnames.push_back( "flag" );
    nodefields.push_back( d->TransferFlag() );
  }

  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Surface output

  std::vector< std::string > nodesurfnames;
  std::vector< std::vector< tk::real > > nodesurfs;

  const auto& ft = multi ? g_cfg.get< tag::fieldout_ >()[ d->MeshId() ]
                         : g_cfg.get< tag::fieldout >();
  const auto& f = ft.get< tag::sidesets >();

  if (!f.empty()) {
    auto ns = ncomp + 1;
    nodesurfnames.push_back( "pressure" );
    nodesurfnames.push_back( "velocityx" );
    nodesurfnames.push_back( "velocityy" );
    nodesurfnames.push_back( "velocityz" );
    if (multi) {
      ++ns;
      nodesurfnames.push_back( "flag" );
    }

    auto bnode = tk::bfacenodes( m_bface, m_triinpoel );
    std::set< int > outsets( begin(f), end(f) );
    for (auto sideset : outsets) {
      auto b = bnode.find(sideset);
      if (b == end(bnode)) continue;
      const auto& nodes = b->second;
      auto i = nodesurfs.size();
      nodesurfs.insert( end(nodesurfs), ns,
                        std::vector< tk::real >( nodes.size() ) );
      std::size_t j = 0;
      for (auto n : nodes) {
        const auto s = m_u[n];
        std::size_t p = 0;
        nodesurfs[i+(p++)][j] = s[0];
        nodesurfs[i+(p++)][j] = s[1];
        nodesurfs[i+(p++)][j] = s[2];
        nodesurfs[i+(p++)][j] = s[3];
        if (multi) nodesurfs[i+(p++)][j] = d->TransferFlag()[n];
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
LohCG::out()
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
    writeFields( CkCallback(CkIndex_LohCG::integrals(), thisProxy[thisIndex]) );
  } else {
    integrals();
  }
}

void
LohCG::integrals()
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
    const auto& reqv = g_cfg.get< tag::integout, tag::integrals >();
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
          fx -= n[0]*m_u(p,0);
          fy -= n[1]*m_u(p,0);
          fz -= n[2]*m_u(p,0);
          // viscous force
          fx += mu * (m_grad(p,3)*n[0] + m_grad(p,4)*n[1] + m_grad(p,5)*n[2]);
          fy += mu * (m_grad(p,6)*n[0] + m_grad(p,7)*n[1] + m_grad(p,8)*n[2]);
          fz += mu * (m_grad(p,9)*n[0] + m_grad(p,10)*n[1] + m_grad(p,11)*n[2]);
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
LohCG::evalLB( int nrestart )
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
LohCG::evalRestart()
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
LohCG::step()
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

#include "NoWarning/lohcg.def.h"
