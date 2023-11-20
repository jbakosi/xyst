// *****************************************************************************
/*!
  \file      src/Inciter/KozCG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     KozCG: Taylor-Galerkin, FCT, element-based continuous Galerkin
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

#include "NoWarning/kozcg.decl.h"

namespace inciter {

//! KozCG Charm++ chare array used to advance PDEs in time with KozCG
class KozCG : public CBase_KozCG {

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
    KozCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit KozCG( const CProxy_Discretization& disc,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit KozCG( CkMigrateMessage* m ) : CBase_KozCG( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Return from migration
    void ResumeFromSync() override;

    //! Start setup for solution
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v );

    // Start time stepping
    void start();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Start (re-)computing domain and boundary integrals
    void feop();

    //! Receive contributions to boundary point normals on chare-boundaries
    void comnorm( const std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& inbnd );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inrhs );

    //! Receive antidiffusive and low-order contributions on chare-boundaries
    void comaec( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inaec );


    //! Receive allowed limits contributions on chare-boundaries
    void comalw( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inalw );

    //! Receive limited antidiffusive contributions on chare-boundaries
    void comlim( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inlim );

    //! Optionally refine/derefine mesh
    void refine( const std::vector< tk::real >& l2res );

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

    //! Evaluate whether to continue with next time step
    void step();

    //! Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Continue to next time step
    void next();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_nrhs;
      p | m_nnorm;
      p | m_naec;
      p | m_nalw;
      p | m_nlim;
      p | m_bnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_u;
      p | m_p;
      p | m_pc;
      p | m_q;
      p | m_qc;
      p | m_a;
      p | m_ac;
      // do not pup these, will recompute after migration anyway
      if (p.isUnpacking()) {
        m_rhs.resize( m_u.nunk(), m_u.nprop() );
      }
      p | m_rhsc;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormc;
      p | m_bndpoinint;
      p | m_dirbcmasks;
      p | m_prebcnodes;
      p | m_prebcvals;
      p | m_symbcnodeset;
      p | m_symbcnodes;
      p | m_symbcnorms;
      p | m_farbcnodeset;
      p | m_farbcnodes;
      p | m_farbcnorms;
      p | m_surfint;
      p | m_dtp;
      p | m_tp;
      p | m_finished;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i KozCG object reference
    friend void operator|( PUP::er& p, KozCG& i ) { i.pup(p); }
    //@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for receiving boundary point normals
    std::size_t m_nnorm;
    //! Counter for receiving antidiffusive contributions
    std::size_t m_naec;
    //! Counter for receiving allowed limits
    std::size_t m_nalw;
    //! Counter for receiving limited antidiffusive contributions
    std::size_t m_nlim;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary triangle face connecitivity where BCs are set by user
    std::vector< std::size_t > m_triinpoel;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Max/min antidiffusive edge contributions at mesh nodes
    tk::Fields m_p;
    //! Receive buffer for max/min antidiffusive edge contributions
    //! \details Key: global node id, value: max/min antidiff edge contributions
    //!   in nodes.
    std::unordered_map< std::size_t, std::vector< tk::real > > m_pc;
    //! Max/min allowed limits at mesh nodes
    tk::Fields m_q;
    //! Receive buffer for max/min allowed limits
    //! \details Key: global node id, value: max/min allowed limits in nodes.
    std::unordered_map< std::size_t, std::vector< tk::real > > m_qc;
    //! Limited antidiffusive contributions at mesh nodes
    tk::Fields m_a;
    //! Receive buffer for limited antidiffusive contributions
    //! \details Key: global node id, value: limited antidiffusive contributions
    //!     in nodes.
    std::unordered_map< std::size_t, std::vector< tk::real > > m_ac;
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
    //!   point, value: weighted normals, inverse distance square, nodal area.
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > > m_bnorm;
    //! Boundary point normals receive buffer
    //! \details Outer key: side set id. Inner key: global node id of boundary
    //!   point, value: weighted normals and inverse distance square.
    decltype(m_bnorm) m_bnormc;
    //! Boundary point integrals
    //! \details Key: global node id of boundary point, value: boundary point
    //!   integral contributions.
    std::unordered_map< std::size_t, std::array< tk::real, 3 > > m_bndpoinint;
    //! Streamable boundary point symmetry BC flags
    std::vector< std::uint8_t > m_bpsym;
    //! Nodes and their Dirichlet BC masks
    std::vector< std::size_t > m_dirbcmasks;
    //! Nodes at pressure BCs
    std::vector< std::size_t > m_prebcnodes;
    //! Density and pressure values at pressure BCs
    std::vector< tk::real > m_prebcvals;
    //! Unique set of ordered nodes at which symmetry BCs are set
    std::set< std::size_t > m_symbcnodeset;
    //! Streamable nodes at which symmetry BCs are set
    std::vector< std::size_t > m_symbcnodes;
    //! Streamable normals at nodes at which symmetry BCs are set
    std::vector< tk::real > m_symbcnorms;
    //! Unique set of ordered nodes at which farfield BCs are set
    std::set< std::size_t > m_farbcnodeset;
    //! Streamable nodes at which farfield BCs are set
    std::vector< std::size_t > m_farbcnodes;
    //! Streamable normals at nodes at which farfield BCs are set
    std::vector< tk::real > m_farbcnorms;
    //! Streamable surface integral nodes and normals * dA on surfaces
    std::map< int, std::pair< std::vector< std::size_t >,
                              std::vector< tk::real > > > m_surfint;
    //! Time step size for each mesh node
    std::vector< tk::real > m_dtp;
    //! Physical time for each mesh node
    std::vector< tk::real > m_tp;
    //! True in the last time step
    int m_finished;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Prepare boundary condition data structures
    void setupBC();

    //! Compute chare-boundary edges
    void bndEdges();

    //! Compute boundary point normals
    void bndint();

    //! Combine own and communicated portions of the boundary point normals
    void bnorm();

    //! Convert integrals into streamable data structures
    void streamable();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback cb );

    //! Combine own and communicated portions of the integrals
    void merge();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Compute antidiffusive contributions: P+/-,  low-order solution: ul
    void aec();

    //! Compute allowed limits, Q+/-
    void alw();

    //! Compute limit coefficients
    void lim();

    //! Advance systems of equations
    void solve();

    //! Compute time step size
    void dt();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Apply boundary conditions
    void BC( tk::Fields& u, tk::real t );

    //! Apply scalar source to solution
    void src();
};

} // inciter::
