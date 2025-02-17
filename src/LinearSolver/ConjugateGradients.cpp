// *****************************************************************************
/*!
  \file      src/LinearSolver/ConjugateGradients.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients.
  \details   Charm++ chare array for asynchronous distributed
             conjugate gradients linear solver for the symmetric linear system
             A * x = b, where A is in compressed sparse row storage, containing
             only the upper triangular part of A.
  \see Y. Saad, Iterative Methods for Sparse Linear Systems: 2nd Edition,
    ISBN 9780898718003, 2003.

    Algorithm 'Preconditioned Conjugate Gradient':
    ```
    Compute r0 := b - A x0, z0 := M^{-1} r0, p0 := z0
    For j=0,1,..., until convergence, do
      alpha_j := (r_j,z_j) / (A p_j,p_j)
      x_{j+1} := x_j + alpha_j p_j
      r_{j+1} := r_j - alpha_j A p_j
      z_{j+1} := M^{-1} r_{j+1}
      beta_j  := (r_{j+1},z_{j+1}) / (r_j,z_j)
      p_{j+1} := z_{j+1} + beta_j p_j
    end
    ```
*/
// *****************************************************************************

#include <numeric>

#include "Exception.hpp"
#include "ConjugateGradients.hpp"
#include "Vector.hpp"
#include "ContainerUtil.hpp"
#include "Reorder.hpp"
#include "Print.hpp"

using tk::ConjugateGradients;

ConjugateGradients::ConjugateGradients(
  const CSR& A,
  const std::vector< tk::real >& x,
  const std::vector< tk::real >& b,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::unordered_map< int,
          std::unordered_set< std::size_t > >& nodecommmap ) :
  m_A( A ),
  m_An( A ),
  m_x( x ),
  m_b( b ),
  m_pc( "none" ),
  m_gid( gid ),
  m_lid( lid ),
  m_nodeCommMap( nodecommmap ),
  m_r( m_A.rsize(), 0.0 ),
  m_z( m_A.rsize(), 0.0 ),
  m_d( m_A.rsize(), 0.0 ),
  m_nr( 0 ),
  m_na( 0 ),
  m_nb( 0 ),
  m_p( m_A.rsize(), 0.0 ),
  m_q( m_A.rsize(), 0.0 ),
  m_nq( 0 ),
  m_nd( 0 ),
  m_initres(),
  m_solved(),
  m_normb( 0.0 ),
  m_it( 0 ),
  m_maxit( 0 ),
  m_finished( false ),
  m_verbose( 0 ),
  m_rho( 0.0 ),
  m_rho0( 0.0 ),
  m_alpha( 0.0 ),
  m_converged( false ),
  m_nx( 0 ),
  m_apply( false )
// *****************************************************************************
//  Constructor
//! \param[in] A Left hand side matrix of the linear system to solve in Ax=b
//! \param[in] x Solution (initial guess) of the linear system to solve in Ax=b
//! \param[in] b Right hand side of the linear system to solve in Ax=b
//! \param[in] gid Global node ids
//! \param[in] lid Local node ids associated to global ones
//! \param[in] nodecommmap Global mesh node IDs shared with other chares
//!   associated to their chare IDs
// *****************************************************************************
{
  // Fill in gid and lid for serial solve
  if (gid.empty() || lid.empty() || nodecommmap.empty()) {
    m_gid.resize( m_A.rsize()/m_A.Ncomp() );
    std::iota( begin(m_gid), end(m_gid), 0 );
    for (auto g : m_gid) m_lid[g] = g;
  }

  Assert( m_A.rsize() == m_gid.size()*A.Ncomp(), "Size mismatch" );
  Assert( m_x.size() == m_gid.size()*A.Ncomp(), "Size mismatch" );
  Assert( m_b.size() == m_gid.size()*A.Ncomp(), "Size mismatch" );
}

void
ConjugateGradients::setup( CkCallback c )
// *****************************************************************************
//  Setup solver
//! \param[in] c Call to continue with
//! \details This function initiates computing the residual (r=b-A*x), its dot
//!   product, and the rhs norm.
// *****************************************************************************
{
  m_initres = c;
  m_converged = false;
  m_finished = false;

  // initiate computing A * x (for the initial residual)
  thisProxy[ thisIndex ].wait4res();
  residual();
  pc();

  // initiate computing norm of right hand side
  dot( m_b, m_b,
       CkCallback( CkReductionTarget(ConjugateGradients,normb), thisProxy ) );
}

void
ConjugateGradients::dot( const std::vector< tk::real >& a,
                         const std::vector< tk::real >& b,
                         CkCallback c )
// *****************************************************************************
//  Initiate computation of dot product of two vectors
//! \param[in] a 1st vector of dot product
//! \param[in] b 2nd vector of dot product
//! \param[in] c Callback to target with the final result
// *****************************************************************************
{
  Assert( a.size() == b.size(), "Size mismatch" );

  tk::real D = 0.0;
  auto ncomp = m_A.Ncomp();
  for (std::size_t i=0; i<a.size()/ncomp; ++i) {
    auto incomp = i*ncomp;
    if (not slave(m_nodeCommMap,m_gid[i],thisIndex))
      for (std::size_t d=0; d<ncomp; ++d)
        D += a[incomp+d] * b[incomp+d];
  }

  contribute( sizeof(tk::real), &D, CkReduction::sum_double, c );
}

void
ConjugateGradients::normb( tk::real n )
// *****************************************************************************
// Compute the norm of the right hand side
//! \param[in] n Norm of right hand side (aggregated across all chares)
// *****************************************************************************
{
  m_normb = std::sqrt(n);
  normb_complete();
}

void
ConjugateGradients::residual()
// *****************************************************************************
//  Initiate A * x for computing the initial residual, r = b - A * x
// *****************************************************************************
{
  // Compute own contribution to r = A * x
  m_A.mult( m_x, m_r );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comres_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > rc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nr( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nr[d] = m_r[ i*ncomp+d ];
        rc[j++] = std::move(nr);
      }
      thisProxy[c].comres( std::vector<std::size_t>(begin(n),end(n)), rc );
    }
  }

  ownres_complete();
}

void
ConjugateGradients::comres( const std::vector< std::size_t >& gid,
                            const std::vector< std::vector< tk::real > >& rc )
// *****************************************************************************
//  Receive contributions to A * x on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] rc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( rc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_rc[ gid[i] ] += rc[i];

  if (++m_nr == m_nodeCommMap.size()) {
    m_nr = 0;
    comres_complete();
  }
}

void
ConjugateGradients::pc()
// *****************************************************************************
//  Setup preconditioner
// *****************************************************************************
{
  auto ncomp = m_A.Ncomp();

  if (m_pc == "none") {

    for (std::size_t i=0; i<m_q.size()/ncomp; ++i) {
      auto c = tk::count( m_nodeCommMap, m_gid[i] );
      for (std::size_t d=0; d<ncomp; ++d) {
        m_q[i*ncomp+d] = 1.0 / c;
      }
    }

  } else if (m_pc == "jacobi") {

    // Extract Jacobi preconditioner from matrix
    for (std::size_t i=0; i<m_q.size()/ncomp; ++i) {
      for (std::size_t d=0; d<ncomp; ++d) {
        m_q[i*ncomp+d] = m_A(i,i,d);
      }
    }

  }

  // Send partial preconditioner on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comd_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > qc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nq( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nq[d] = m_q[ i*ncomp+d ];
        qc[j++] = std::move(nq);
      }
      thisProxy[c].comd( std::vector<std::size_t>(begin(n),end(n)), qc );
    }
  }

  ownd_complete();
}

void
ConjugateGradients::comd( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& qc )
// *****************************************************************************
//  Receive contributions to preconditioner on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] qc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( qc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_qc[ gid[i] ] += qc[i];

  if (++m_nd == m_nodeCommMap.size()) {
    m_nd = 0;
    comd_complete();
  }
}

void
ConjugateGradients::initres()
// *****************************************************************************
// Finish computing the initial residual, r = b - A * x
// *****************************************************************************
{
  auto ncomp = m_A.Ncomp();

  // Combine own and communicated contributions to r = A * x
  for (const auto& [gid,r] : m_rc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c) m_r[i*ncomp+c] += r[c];
  }
  tk::destroy( m_rc );

  // Finish computing the initial residual, r = b - A * x
  for (auto& r : m_r) r *= -1.0;
  m_r += m_b;

  // Initialize p
  m_p = m_r;

  // Combine own and communicated contributions to q = A * p
  for (const auto& [gid,q] : m_qc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c)
      m_q[i*ncomp+c] += q[c];
  }
  tk::destroy( m_qc );

  // Set preconditioner
  m_d = m_q;

  // initialize preconditioned solve: compute z = M^{-1} * r
  for (std::size_t i=0; i<m_z.size(); ++i) m_z[i] = m_r[i] / m_d[i];

  // initiate computing the dot product of the initial residual, rho = (r,z)
  dot( m_r, m_z,
       CkCallback( CkReductionTarget(ConjugateGradients,rho), thisProxy ) );
}

void
ConjugateGradients::rho( tk::real r )
// *****************************************************************************
// Compute rho = (r,z)
//! \param[in] r Dot product, rho = (r,z) (aggregated across all chares)
// *****************************************************************************
{
  // store dot product of residual
  m_rho = r;

  // send back rhs norm to caller
  m_initres.send( CkDataMsg::buildNew( sizeof(tk::real), &m_normb ) );
}

void
ConjugateGradients::init(
  const std::vector< tk::real >& x,
  const std::vector< tk::real >& b,
  const std::vector< tk::real >& neubc,
  const std::unordered_map< std::size_t,
          std::vector< std::pair< int, tk::real > > >& dirbc,
  bool apply,
  const std::string& pc,
  CkCallback cb )
// *****************************************************************************
//  Initialize linear solve: set initial guess and boundary conditions
//! \param[in] x Initial guess
//! \param[in] b Right hand side vector
//! \param[in] neubc Right hand side vector with Neumann BCs partially applied
//! \param[in] dirbc Local node ids and associated Dirichlet BCs
//! \param[in] apply True to apply boundary conditions
//! \param[in] pc Preconditioner to use
//! \param[in] cb Call to continue with when initialized and ready for a solve
//! \details This function allows setting the initial guess and boundary
//!   conditions, followed by computing the initial residual and the rhs norm.
//!   It also performs communication of BC data.
// *****************************************************************************
{
  // Configure preconditioner
  m_pc = pc;

  // Optionally set initial guess
  if (!x.empty()) m_x = x; //else std::fill( begin(m_x), end(m_x), 0.0 );

  // Optionally set rhs, assumed complete in parallel
  if (!b.empty()) m_b = b;

  // Save if applied BCs
  m_apply = apply;
  if (not m_apply) {

    // Recompute initial residual (r=b-A*x), its dot product, and the rhs norm
    setup( cb );

  }
  else {

    // Save matrix state
    m_An = m_A;

    // Set Neumann BCs, partial in parallel, communication below
    if (!neubc.empty()) m_q = neubc; else std::fill(begin(m_q), end(m_q), 0.0);

    // Store incoming Dirichlet BCs, partial in parallel, communication below
    tk::destroy(m_dirbc);
    for (auto&& [i,bcval] : dirbc) m_dirbc[i] = std::move(bcval);

    // Get ready to communicate boundary conditions. This is necessary because
    // there can be nodes a chare contributes to but does not apply BCs on.
    // This happens if a node is in the node communication map but not on the
    // list of incoming BCs on this chare. To have all chares share the same
    // view on all BC nodes, we send the global node ids together with the
    // Dirichlet BCs at which BCs are set to those fellow chares that also
    // contribute to those BC nodes. Only after this communication step will we
    // apply the BCs, which then will correctly setup the BC rows of the matrix
    // and on the rhs that (1) may contain a nonzero value, (2) may have partial
    // contributions due to Neumann BCs, and (3) will be modified to apply
    // Dirichlet BCs.
    thisProxy[ thisIndex ].wait4bc();
    thisProxy[ thisIndex ].wait4r();

    // Send boundary conditions to those who contribute to those rows
    if (m_nodeCommMap.empty()) {
      combc_complete();
    } else {
      auto ncomp = m_A.Ncomp();
      for (const auto& [c,n] : m_nodeCommMap) {
        decltype(m_dirbc) dbc;
        std::vector< std::vector< tk::real > > qc( n.size() );
        std::size_t k = 0;
        for (auto g : n) {
          auto i = tk::cref_find( m_lid, g );
          auto j = m_dirbc.find(i);
          if (j != end(m_dirbc)) dbc[g] = j->second;
          std::vector< tk::real > nq( ncomp );
          for (std::size_t d=0; d<ncomp; ++d) nq[d] = m_q[ i*ncomp+d ];
          qc[k++] = std::move(nq);
        }
        std::vector< std::size_t > gid( begin(n), end(n) );
        thisProxy[c].combc( dbc, gid, qc );
      }
    }

    ownbc_complete( cb );

  }
}

void
ConjugateGradients::combc(
  const std::map< std::size_t, std::vector< std::pair< int, tk::real > > >& dbc,
  const std::vector< std::size_t >& gid,
  const std::vector< std::vector< tk::real > >& qc )
// *****************************************************************************
//  Receive contributions to boundary conditions and rhs on chare-boundaries
//! \param[in] dbc Contributions to Dirichlet boundary conditions
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] qc Contributions to Neumann boundary conditions
// *****************************************************************************
{
  for (const auto& [g,dirbc] : dbc) m_dirbcc[g] = dirbc;
  for (std::size_t i=0; i<gid.size(); ++i) m_qc[ gid[i] ] += qc[i];

  if (++m_nb == m_nodeCommMap.size()) {
    m_nb = 0;
    combc_complete();
  }
}

void
ConjugateGradients::apply( CkCallback cb )
// *****************************************************************************
//  Apply boundary conditions
//! \param[in] cb Call to continue with after applying the BCs is complete
// *****************************************************************************
{
  auto ncomp = m_A.Ncomp();

  // Merge own and received contributions to Dirichlet boundary conditions
  for (const auto& [g,dirbc] : m_dirbcc) {
    m_dirbc[ tk::cref_find(m_lid,g) ] = dirbc;
  }
  tk::destroy( m_dirbcc );

  // Merge own and received contributions to Neumann boundary conditions
  for (const auto& [gid,q] : m_qc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c) m_q[i*ncomp+c] += q[c];
  }
  tk::destroy( m_qc );

  // Apply Neumann BCs on rhs
  m_b += m_q;

  // Apply Dirichlet BCs on matrix and rhs (with decreasing local id)
  std::fill( begin(m_r), end(m_r), 0.0 );
  for (auto bi = m_dirbc.rbegin(); bi != m_dirbc.rend(); ++bi) {
    auto i = bi->first;
    const auto& dirbc = bi->second;
    for (std::size_t j=0; j<ncomp; ++j) {
      if (dirbc[j].first) {
        m_A.dirichlet( i, dirbc[j].second, m_r, m_gid, m_nodeCommMap, j );
      }
    }
  }

  // Send rhs with Dirichlet BCs partially applied to fellow chares
  if (m_nodeCommMap.empty()) {
    comr_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > rc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nr( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nr[d] = m_r[ i*ncomp+d ];
        rc[j++] = std::move(nr);
      }
      thisProxy[c].comr( std::vector<std::size_t>(begin(n),end(n)), rc );
    }
  }

  ownr_complete( cb );
}

void
ConjugateGradients::comr( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& rc )
// *****************************************************************************
//  Receive contributions to rhs with Dirichlet BCs applied on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] rc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( rc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_rc[ gid[i] ] += rc[i];

  if (++m_na == m_nodeCommMap.size()) {
    m_na = 0;
    comr_complete();
  }
}

void
ConjugateGradients::r( CkCallback cb )
// *****************************************************************************
//  Finish computing rhs with BCs applied
//! \param[in] cb Call to continue with after applying the BCs is complete
// *****************************************************************************
{
  auto ncomp = m_A.Ncomp();

  // Combine own and communicated contributions to Dirichlet BCs applied to rhs
  for (const auto& [gid,r] : m_rc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c) m_r[i*ncomp+c] += r[c];
  }
  tk::destroy( m_rc );

  // Subtract matrix columns * BC values from rhs at Dirichlet BC nodes
  m_b -= m_r;

  // Apply Dirichlet BCs on rhs
  for (const auto& [i,dirbc] : m_dirbc) {
    for (std::size_t j=0; j<ncomp; ++j) {
      if (dirbc[j].first) {
        m_b[i*ncomp+j] = dirbc[j].second;
      }
    }
  }

  // Recompute initial residual (r=b-A*x), its dot product, and the rhs norm
  setup( cb );
}

void
ConjugateGradients::solve( std::size_t maxit,
                           tk::real tol,
                           int pe,
                           uint64_t verbose,
                           CkCallback c )
// *****************************************************************************
//  Solve linear system
//! \param[in] maxit Max iteration count
//! \param[in] tol Stop tolerance
//! \param[in] pe Processing element
//! \param[in] verbose Verbose output
//! \param[in] c Call to continue with after solve is complete
// *****************************************************************************
{
  m_maxit = maxit;
  m_tol = tol;
  m_solved = c;
  m_it = 0;
  m_verbose = pe == 0 ? verbose : 0;

  if (m_verbose > 1) tk::Print() << "Xyst> Conjugate gradients start\n";

  if (m_converged) x(); else next();
}

void
ConjugateGradients::next()
// *****************************************************************************
//  Start next linear solver iteration
// *****************************************************************************
{
  if (m_it == 0) m_alpha = 0.0; else m_alpha = m_rho/m_rho0;
  m_rho0 = m_rho;

  // compute p = z + alpha * p
  for (std::size_t i=0; i<m_p.size(); ++i) m_p[i] = m_z[i] + m_alpha * m_p[i];

  // initiate computing q = A * p
  thisProxy[ thisIndex ].wait4q();
  qAp();
}


void
ConjugateGradients::qAp()
// *****************************************************************************
//  Initiate computing q = A * p
// *****************************************************************************
{
  // Compute own contribution to q = A * p
  m_A.mult( m_p, m_q );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comq_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > qc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nq( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nq[d] = m_q[ i*ncomp+d ];
        qc[j++] = std::move(nq);
      }
      thisProxy[c].comq( std::vector<std::size_t>(begin(n),end(n)), qc );
    }
  }

  ownq_complete();
}

void
ConjugateGradients::comq( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& qc )
// *****************************************************************************
//  Receive contributions to q = A * p on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] qc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( qc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_qc[ gid[i] ] += qc[i];

  if (++m_nq == m_nodeCommMap.size()) {
    m_nq = 0;
    comq_complete();
  }
}

void
ConjugateGradients::q()
// *****************************************************************************
// Finish computing q = A * p
// *****************************************************************************
{
  // Combine own and communicated contributions to q = A * p
  auto ncomp = m_A.Ncomp();
  for (const auto& [gid,q] : m_qc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c)
      m_q[i*ncomp+c] += q[c];
  }
  tk::destroy( m_qc );

  // initiate computing (p,q)
  dot( m_p, m_q,
       CkCallback( CkReductionTarget(ConjugateGradients,pq), thisProxy ) );
}

void
ConjugateGradients::pq( tk::real d )
// *****************************************************************************
// Compute the dot product (p,q)
//! \param[in] d Dot product of (p,q) (aggregated across all chares)
// *****************************************************************************
{
  // If (p,q)=0, then p and q are orthogonal and the system either has a trivial
  // solution, x=x0, or the BCs are incomplete or wrong, in either case the
  // solve cannot continue.
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::abs(d) < eps) {
    if (m_verbose > 1) {
      tk::Print() << "Xyst> Conjugate gradients it " << m_it << ", (p,q) = 0\n";
    }
    m_finished = true;
    m_alpha = 0.0;
  } else {
    m_alpha = m_rho / d;
  }

  // compute r = r - alpha * q
  for (std::size_t i=0; i<m_r.size(); ++i) m_r[i] -= m_alpha * m_q[i];

  // apply preconditioner: compute z = M^{-1} * r
  for (std::size_t i=0; i<m_z.size(); ++i) m_z[i] = m_r[i] / m_d[i];

  // initiate computing (r,z)
  dot( m_r, m_z,
       CkCallback( CkReductionTarget(ConjugateGradients,rz), thisProxy ) );
  // initiate computing norm of residual: (r,r)
  dot( m_r, m_r,
       CkCallback( CkReductionTarget(ConjugateGradients,normres), thisProxy ) );
}

void
ConjugateGradients::normres( tk::real r )
// *****************************************************************************
// Compute norm of residual: (r,r)
//! \param[in] r Dot product, (r,r) (aggregated across all chares)
// *****************************************************************************
{
  m_normr = r;
  normr_complete();
}

void
ConjugateGradients::rz( tk::real rz )
// *****************************************************************************
// Compute (r,z)
//! \param[in] rz Dot product, (r,z) (aggregated across all chares)
// *****************************************************************************
{
  m_rho = rz;

  // Advance solution: x = x + alpha * p
  for (std::size_t i=0; i<m_x.size(); ++i) m_x[i] += m_alpha * m_p[i];

  // Communicate solution
  thisProxy[ thisIndex ].wait4x();

  // Send solution on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comx_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > xc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nx( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nx[d] = m_x[ i*ncomp+d ];
        xc[j++] = std::move(nx);
      }
      thisProxy[c].comx( std::vector<std::size_t>(begin(n),end(n)), xc );
    }
  }

  ownx_complete();
}

void
ConjugateGradients::comx( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& xc )
// *****************************************************************************
//  Receive contributions to final solution on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] xc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( xc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_xc[ gid[i] ] += xc[i];

  if (++m_nx == m_nodeCommMap.size()) {
    m_nx = 0;
    comx_complete();
  }
}

void
ConjugateGradients::x()
// *****************************************************************************
//  Assemble solution on chare boundaries and decide what's next
// *****************************************************************************
{
  // Assemble solution on chare boundaries by averaging
  auto ncomp = m_A.Ncomp();
  for (const auto& [g,x] : m_xc) {
    auto i = tk::cref_find(m_lid,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] += x[d];
    auto c = tk::count(m_nodeCommMap,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] /= c;
  }
  tk::destroy( m_xc );

  ++m_it;
  auto normb = m_normb > 1.0e-14 ? m_normb : 1.0;
  auto normr = std::sqrt( m_normr );

  if ( m_finished || normr < m_tol*normb || m_it >= m_maxit ) {

    m_converged = normr > m_tol*normb ? false : true;

    if (m_verbose) {
      std::stringstream c;
      if (m_converged) {
        c << " < " << m_tol << ", converged";
      } else {
        c << " > " << m_tol << ", not converged, ||r||/||b|| = "
          << normr << '/' << normb;
      }
      tk::Print() << "Xyst> Conjugate gradients it " << m_it
                  << ", norm = " << normr/normb << c.str() << '\n';
    }

    // Restore matrix to state before init()
    if (m_apply) m_A = m_An;

    m_solved.send( CkDataMsg::buildNew( sizeof(tk::real), &normr ) );

  } else {

    if (m_verbose > 1) {
      tk::Print() << "Xyst> Conjugate gradients it "
                  << m_it << ", norm = " << normr/normb << '\n';
    }

    next();

  }
}

#include "NoWarning/conjugategradients.def.h"
