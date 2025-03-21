// *****************************************************************************
/*!
  \file      src/Inciter/LohCG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LohCG: Artificial compressibility solver for incompressible flow
*/
// *****************************************************************************

#pragma once

#include <vector>
#include <map>

#include "Types.hpp"
#include "Fields.hpp"
#include "Table.hpp"
#include "DerivedData.hpp"
#include "NodeDiagnostics.hpp"
#include "PUPUtil.hpp"

#include "NoWarning/lohcg.decl.h"

namespace inciter {

//! LohCG Charm++ chare array used to advance PDEs in time with LohCG
class LohCG : public CBase_LohCG {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    LohCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit LohCG( const CProxy_Discretization& disc,
                    const tk::CProxy_ConjugateGradients& cgpre,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit LohCG( CkMigrateMessage* m ) : CBase_LohCG( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Return from migration
    void ResumeFromSync() override;

    //! Start setup for solution
    void setup( tk::real v );

    //! Initialize Poisson solve
    void pinit();

    //! Solve Poisson equation
    void psolve();

    //! Continue after Poisson solve
    void psolved();

    // Start time stepping
    void start();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Evaluate diagnostics
    void diag();

    //! Start (re-)computing domain and boundary integrals
    void feop();

    //! Receive contributions to boundary point normals on chare-boundaries
    void comnorm( const std::unordered_map< int,
      std::unordered_map< std::size_t, std::array<tk::real,4> > >& inbnd );

    //! Receive contributions to velocity gradients
    void comvgrad( const std::unordered_map< std::size_t,
                           std::vector< tk::real > >& ingrad );

    //! Receive contributions to momentum flux on chare-boundaries
    void comflux( const std::unordered_map< std::size_t,
                          std::vector< tk::real > >& influx );

    //! Receive contributions to conjugate gradients solution gradient
    void comsgrad( const std::unordered_map< std::size_t,
                            std::vector< tk::real > >& ingrad );

    //! Receive contributions to gradient
    void comgrad( const std::unordered_map< std::size_t,
                          std::vector< tk::real > >& ingrad );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inrhs );

    //! Receive contributions to velocity divergence on chare-boundaries
    void comdiv( const std::unordered_map< std::size_t, tk::real >& indiv );

    //! Solution has been updated
    void solved();

    //! Evaluate residuals
    void evalres( const std::vector< tk::real >& l2res );

    //! Receive new mesh from Refiner
    void resizePostAMR(
      const std::vector< std::size_t >& ginpoel,
      const tk::UnsMesh::Chunk& chunk,
      const tk::UnsMesh::Coords& coord,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
      const std::unordered_map< std::size_t, std::size_t >& addedTets,
      const std::set< std::size_t >& removedNodes,
      const std::unordered_map< int, std::unordered_set< std::size_t > >&
        nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& bnode,
      const std::vector< std::size_t >& triinpoel );

    //! Const-ref access to current solution
    //! \return Const-ref to current solution
    const tk::Fields& solution() const { return m_u; }

    //! Compute integral quantities for output
    void integrals();

    //! Compute recent conjugate gradients solution gradient
    void sgrad();

    //! Evaluate whether to continue with next time step
    void step();

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_cgpre;
      p | m_nrhs;
      p | m_nnorm;
      p | m_ngrad;
      p | m_nsgrad;
      p | m_nvgrad;
      p | m_nflux;
      p | m_ndiv;
      p | m_nbpint;
      p | m_np;
      p | m_bnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_u;
      p | m_un;
      p | m_grad;
      p | m_gradc;
      p | m_vgrad;
      p | m_vgradc;
      p | m_flux;
      p | m_fluxc;
      p | m_div;
      p | m_divc;
      p | m_sgrad;
      p | m_sgradc;
      p | m_rhs;
      p | m_rhsc;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormc;
      p | m_bndpoinint;
      p | m_domedgeint;
      p | m_bpint;
      p | m_dsupedge;
      p | m_dsupint;
      p | m_dirbcmask;
      p | m_dirbcval;
      p | m_dirbcmaskp;
      p | m_dirbcvalp;
      p | m_symbcnodes;
      p | m_symbcnorms;
      p | m_noslipbcnodes;
      p | m_surfint;
      p | m_stage;
      p | m_finished;
      p | m_rkcoef;
      p | m_timer;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i LohCG object reference
    friend void operator|( PUP::er& p, LohCG& i ) { i.pup(p); }
    //@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Conjugate Gradients Charm++ proxy for pressure solve
    tk::CProxy_ConjugateGradients m_cgpre;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for receiving boundary point normals
    std::size_t m_nnorm;
    //! Counter for receiving gradient
    std::size_t m_ngrad;
    //! Counter for receiving conjugrate gradient solution gradient
    std::size_t m_nsgrad;
    //! Counter for receiving velocity gradient
    std::size_t m_nvgrad;
    //! Counter for receiving momentum flux
    std::size_t m_nflux;
    //! Counter for receiving boundary velocity divergences
    std::size_t m_ndiv;
    //! Counter for receiving boundary point integrals
    std::size_t m_nbpint;
    //! Count number of Poisson solves during setup
    std::size_t m_np;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary triangle face connecitivity where BCs are set by user
    std::vector< std::size_t > m_triinpoel;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Unknown/solution vector at mesh nodes at previous time step
    tk::Fields m_un;
    //! Gradient in mesh nodes
    tk::Fields m_grad;
    //! Gradient receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_gradc;
    //! Velocity gradient in mesh nodes
    tk::Fields m_vgrad;
    //! Velocity gradient receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_vgradc;
    //! Momentum flux in mesh nodes
    tk::Fields m_flux;
    //! Momentum flux receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_fluxc;
    //! Velocity divergence
    std::vector< tk::real > m_div;
    //! Receive buffer for communication of the velocity divergence
    //! \details Key: global node id, value: velocity divergence
    std::unordered_map< std::size_t, tk::real > m_divc;
    //! Conjugate gradient solution gradient in mesh nodes
    tk::Fields m_sgrad;
    //! Conjugate gradient solution gradient receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_sgradc;
    //! Right-hand side vector (for the high order system)
    tk::Fields m_rhs;
    //! Receive buffer for communication of the right hand side
    //! \details Key: global node id, value: rhs for all scalar components per
    //!   node.
    std::unordered_map< std::size_t, std::vector< tk::real > > m_rhsc;
    //! Diagnostics object
    NodeDiagnostics m_diag;
    //! Boundary point normals
    //! \details Outer key: side set id. Inner key: global node id of boundary
    //!   point, value: weighted normal vector, inverse distance square.
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnorm;
    //! Boundary point normals receive buffer
    //! \details Outer key: side set id. Inner key: global node id of boundary
    //!   point, value: weighted normals and inverse distance square.
    decltype(m_bnorm) m_bnormc;
    //! Boundary point integrals
    //! \details Key: global node id of boundary point, value: boundary point
    //!   integral contributions.
    std::unordered_map< std::size_t, std::array<tk::real,3> > m_bndpoinint;
    //! Domain edge integrals
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 4 >,
      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_domedgeint;
    //! Streamable boundary point integrals
    std::vector< tk::real > m_bpint;
    //! Superedge (tet, face, edge) end points with local ids for domain edges
    std::array< std::vector< std::size_t >, 3 > m_dsupedge;
    //! Superedge (tet, face, edge) domain edge integrals
    std::array< std::vector< tk::real >, 3 > m_dsupint;
    //! Nodes and their Dirichlet BC masks
    std::vector< std::size_t > m_dirbcmask;
    //! Nodes and their Dirichlet BC values
    std::vector< double > m_dirbcval;
    //! Nodes and their pressure Dirichlet BC masks
    std::vector< std::size_t > m_dirbcmaskp;
    //! Nodes and their pressure Dirichlet BC values
    std::vector< double > m_dirbcvalp;
    //! Streamable nodes at which symmetry BCs are set
    std::vector< std::size_t > m_symbcnodes;
    //! Streamable normals at nodes at which symmetry BCs are set
    std::vector< tk::real > m_symbcnorms;
    //! Streamable nodes at which noslip BCs are set
    std::vector< std::size_t > m_noslipbcnodes;
    //! Streamable surface integral nodes and normals * dA on surfaces
    std::map< int, std::pair< std::vector< std::size_t >,
                              std::vector< tk::real > > > m_surfint;
    //! Runge-Kutta stage counter
    std::size_t m_stage;
    //! True in the last time step
    int m_finished;
    //! Runge-Kutta coefficients
    std::vector< tk::real > m_rkcoef;
    //! Timer
    std::vector< tk::Timer > m_timer;

    //! Compute number of scalar components for gradients
    std::size_t ngradcomp() const;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

   //! Prepare Dirichlet boundary condition data structures
   void setupDirBC( const std::vector< std::vector< int > >& cfgmask,
                    const std::vector< std::vector< double > >& cfgval,
                    std::size_t ncomp,
                    std::vector< std::size_t >& mask,
                    std::vector< double >& val );

    //! Start computing velocity divergence
    void div( const tk::Fields& u, std::size_t pos = 0 );

    //! Start computing velocity gradient
    void velgrad();

    //! Start computing momentum flux
    void flux();

    //! Finalize computing gradient
    void fingrad( tk::Fields& grad,
      std::unordered_map< std::size_t, std::vector< tk::real > >& gradc );

    //! Compute local contributions to domain edge integrals
    void domint();

    //! Setup lhs matrix for pressure solve
    std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
    prelhs( const std::pair< std::vector< std::size_t >,
                             std::vector< std::size_t > >& psup );

    //! Compute chare-boundary edges
    void bndEdges();

    //! Compute local contributions to boundary normals and integrals
    void bndint();

    //! Combine own and communicated portions of the boundary point normals
    void bnorm();

    //! Prepare surface integral data strurctures
    void prep_surfint();

    //! Prepare symmetry boundary condition data structures
    void prep_symbc();

    //! Prepare no-slip boundary condition data structures
    void prep_noslipbc();

    //! Prepare integrid-boundary data structures (if coupled)
    void prep_intergrid();

    //! Convert integrals into streamable data structures
    void streamable();

    //! Generate superedge-groups for domain-edge loops
    void domsuped();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback cb );

    //! Combine own and communicated portions of the integrals
    void merge();

    //! Compute gradients
    void grad();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Advance systems of equations
    void solve();

    //! Start next time step stage
    void stage();

    //! Optionally refine/derefine mesh
    void refine();

    //! Compute time step size
    void dt();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Apply scalar source to solution
    void src();
};

} // inciter::
