// *****************************************************************************
/*!
  \file      src/LinearSolver/ConjugateGradients.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients
  \details   Charm++ chare array for asynchronous distributed
    conjugate gradients linear solver.

    There are a potentially large number of ConjugateGradients Charm++ chares.
    Each ConjugateGradient chare gets a chunk of the full load, due to partiting
    the mesh, on which the solve is performed.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#pragma once

#include "Types.hpp"
#include "CSR.hpp"

#include "NoWarning/conjugategradients.decl.h"

namespace tk {

//! \brief ConjugateGradients Charm++ chare array used to perform a distributed
//!   linear solve with the conjugate gradients algorithm
class ConjugateGradients : public CBase_ConjugateGradients {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    ConjugateGradients_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #endif

    //! Constructor
    explicit ConjugateGradients(
      const CSR& A,
      const std::vector< tk::real >& x,
      const std::vector< tk::real >& b,
      const std::vector< std::size_t >& gid = {},
      const std::unordered_map< std::size_t, std::size_t >& lid = {},
      const std::unordered_map< int,
              std::unordered_set< std::size_t > >& nodecommmap = {} );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif

    //! Constructor taking a tuple of {A,x,b} by rvalue reference
    explicit ConjugateGradients(
      std::tuple< tk::CSR,
                  std::vector< tk::real >,
                  std::vector< tk::real > >&& system,
      const std::vector< std::size_t >& gid,
      const std::unordered_map< std::size_t, std::size_t >& lid,
      const std::unordered_map< int,
              std::unordered_set< std::size_t > >& nodecommmap ) :
      ConjugateGradients( std::move(std::get<0>(system)),
                          std::move(std::get<1>(system)),
                          std::move(std::get<2>(system)),
                          gid, lid, nodecommmap ) {}

    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Migrate constructor
    explicit ConjugateGradients( CkMigrateMessage* m )
     : CBase_ConjugateGradients( m ) {}

    //! Solve linear system
    void solve( std::size_t maxit,
                tk::real tol,
                int pe,
                uint64_t verbose,
                CkCallback c );

    //! Initialize linear solve: set initial guess and boundary conditions
    void init( const std::vector< tk::real >& x,
               const std::vector< tk::real >& b,
               const std::vector< tk::real >& neubc,
               const std::unordered_map< std::size_t,
                       std::vector< std::pair< int, tk::real > > >& dirbc,
               const std::string& pc,
               CkCallback cb );

    //! Setup solver
    void setup( CkCallback c );

    //! Compute the norm of the right hand side
    void normb( tk::real n );

    //! Compute rho = (r,z)
    void rho( tk::real r );

    //! Receive contributions to r = b - A * x on chare-boundaries
    void comres( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& rc );

    //! Receive contributions to boundary conditions and rhs on chare-boundaries
    void combc( const std::map< std::size_t,
                        std::vector< std::pair< int, tk::real > > >& dbc,
                const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& qc );

    //! \brief Receive contributions to rhs with Dirichlet BCs applied on
    //!        chare-boundaries
    void comr( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& rc );

    //! Receive contributions to preconditioner chare-boundaries
    void comd( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& qc );

    //! Receive contributions to q = A * p on chare-boundaries
    void comq( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& qc );

    //! Receive contributions to final solution on chare-boundaries
    void comx( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& xc );

    //! Compute the dot product (p,q)
    void pq( tk::real d );

    //! Compute the norm of the residual (r,r)
    void normres( tk::real r );

    //! Compute the dot product (r,z)
    void rz( tk::real rz );

    //! Access solution
    const std::vector< tk::real >& solution() const { return m_x; }

    //! Return convergence flag
    bool converged() const { return m_converged; }

    //! Return number of iterations taken
    std::size_t it() const { return m_it; }

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_A;
      p | m_An;
      p | m_x;
      p | m_b;
      p | m_pc;
      p | m_gid;
      p | m_lid;
      p | m_nodeCommMap;
      p | m_r;
      p | m_z;
      p | m_d;
      p | m_rc;
      p | m_nr;
      p | m_na;
      p | m_dirbc;
      p | m_dirbcc;
      p | m_nb;
      p | m_p;
      p | m_q;
      p | m_qc;
      p | m_nq;
      p | m_nd;
      p | m_initres;
      p | m_solved;
      p | m_normb;
      p | m_it;
      p | m_maxit;
      p | m_finished;
      p | m_verbose;
      p | m_tol;
      p | m_rho;
      p | m_rho0;
      p | m_alpha;
      p | m_converged;
      p | m_xc;
      p | m_nx;
      p | m_normr;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c ConjugateGradients object reference
    friend void operator|( PUP::er& p, ConjugateGradients& c ) { c.pup(p); }
    ///@}

  private:
    //! Sparse matrix
    CSR m_A;
    //! Sparse matrix before boundary conditions
    CSR m_An;
    //! Solution/unknown
    std::vector< tk::real > m_x;
    //! Right hand side
    std::vector< tk::real > m_b;
    //! Preconditioner to use
    std::string m_pc;
    //! Global node IDs
    std::vector< std::size_t > m_gid;
    //! Local node IDs associated to global ones
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Global mesh node IDs shared with other chares associated to chare IDs
    std::unordered_map< int, std::unordered_set< std::size_t > > m_nodeCommMap;
    //! Auxiliary vector for CG solve
    std::vector< tk::real > m_r;
    //! Receive buffer for communication of r = b - A * x
    std::unordered_map< std::size_t, std::vector< tk::real > > m_rc;
    //! Auxiliary vector for preconditioned CG solve
    std::vector< tk::real > m_z;
    //! Jacobi preconditioner
    std::vector< tk::real > m_d;
    //! Counter for assembling m_r
    std::size_t m_nr;
    //! Counter for assembling m_r (rhs with BCs applied)
    std::size_t m_na;
    //! Dirichlet boundary conditions
    std::map< std::size_t, std::vector< std::pair<int,tk::real> > > m_dirbc;
    //! Dirichlet boundary conditions communication buffer
    std::map< std::size_t, std::vector< std::pair<int,tk::real> > > m_dirbcc;
    //! Counter for assembling boundary conditions
    std::size_t m_nb;
    //! Auxiliary vector for CG solve
    std::vector< tk::real > m_p;
    //! Auxiliary vector for CG solve
    std::vector< tk::real > m_q;
    //! Receive buffer for communication of q = A * p
    std::unordered_map< std::size_t, std::vector< tk::real > > m_qc;
    //! Counter for assembling m_q
    std::size_t m_nq;
    //! Counter for assembling the preconditioner
    std::size_t m_nd;
    //! Charm++ callback to continue with when the setup is complete
    CkCallback m_initres;
    //! Charm++ callback to continue with when the solve is complete
    CkCallback m_solved;
    //! L2 norm of the right hand side
    tk::real m_normb;
    //! Iteration count
    std::size_t m_it;
    //! Max iteration count
    std::size_t m_maxit;
    //! True if finished
    bool m_finished;
    //! Verbose output
    uint64_t m_verbose;
    //! Stop tolerance
    tk::real m_tol;
    //! Helper scalar for CG algorithm
    tk::real m_rho;
    //! Helper scalar for CG algorithm
    tk::real m_rho0;
    //! Helper scalar for CG algorithm
    tk::real m_alpha;
    //! Convergence flag: true if linear smoother converged to tolerance
    bool m_converged;
    //! Receive buffer for solution
    std::unordered_map< std::size_t, std::vector< tk::real > > m_xc;
    //! Counter for assembling the solution on chare boundaries
    std::size_t m_nx;
    //! Norm of the residual
    tk::real m_normr;

    //! Initiate computationa of dot product of two vectors
    void dot( const std::vector< tk::real >& a,
              const std::vector< tk::real >& b,
              CkCallback c );

    //! Initiate A * x for computing the residual, r = b - A * x
    void residual();

    //! Finish computing the initial residual, r = b - A * x
    void initres();

    //! Setup preconditioner
    void pc();

    //! Apply boundary conditions
    void apply( CkCallback cb );

    //! Finish computing rhs with applied BCs
    void r( CkCallback cb );

    //! Initiate computing q = A * p
    void qAp();

    //! Finish computing q = A * p
    void q();

    //! Start next linear solver iteration
    void next();

    //! Assemble solution on chare boundaries and decide what's next
    void x();
};

} // tk::
