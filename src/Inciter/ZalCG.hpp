// *****************************************************************************
/*!
  \file      src/Inciter/ZalCG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     ZalCG: Taylor-Galerkin, FCT, edge-based continuous Galerkin
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
#include "InciterConfig.hpp"

#include "NoWarning/zalcg.decl.h"

namespace inciter {

extern ctr::Config g_cfg;

//! ZalCG Charm++ chare array used to advance PDEs in time with ZalCG
class ZalCG : public CBase_ZalCG {

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
    ZalCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit ZalCG( const CProxy_Discretization& disc,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit ZalCG( CkMigrateMessage* m ) : CBase_ZalCG( m ) {}
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

    //! Receive contributions to node gradients on chare-boundaries
    void comgrad( const std::unordered_map< std::size_t,
                          std::vector< tk::real > >& ingrad );

    //! Receive contributions to stabilization contributions on chare-boundaries
    void comstab( const std::unordered_map< std::size_t, tk::real >& instab );

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

    //! Evaluate residuals
    void evalres( const std::vector< tk::real >& l2res );

    //! Receive activation request
    void comrea( int reactivate );

    //! Receive activation status
    void comact( int ch, int deactivated );

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

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Continue to next time step
    void next();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_ngrad;
      p | m_nstab;
      p | m_nrhs;
      p | m_nnorm;
      p | m_naec;
      p | m_nalw;
      p | m_nlim;
      p | m_nrea;
      p | m_nact;
      p | m_todeactivate;
      p | m_toreactivate;
      p | m_deactivated;
      p | m_inactive;
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
        m_grad.resize( g_cfg.get< tag::stab4 >() ? m_u.nunk() : 0,
                       3 + m_u.nprop() );
        m_stab.resize( m_grad.nunk(), 1UL );
      }
      p | m_rhsc;
      p | m_gradc;
      p | m_stabc;
      p | m_vol;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormc;
      p | m_bndpoinint;
      p | m_bndedgeint;
      p | m_domedgeint;
      p | m_bpoin;
      p | m_bpint;
      p | m_bsupedge;
      p | m_bsupint;
      p | m_dsupedge;
      p | m_dsupint;
      p | m_chbndedge;
      p | m_besym;
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
    //! \param[in,out] i ZalCG object reference
    friend void operator|( PUP::er& p, ZalCG& i ) { i.pup(p); }
    //@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for receiving gradients
    std::size_t m_ngrad;
    //! Counter for receiving stabilization coefficients
    std::size_t m_nstab;
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
    //! Counter for receiving reactivation requests
    std::size_t m_nrea;
    //! Counter for receiving activation status communications
    std::size_t m_nact;
    //! Flag: 1 if chare desires to deactivate
    int m_todeactivate;
    //! Flag: 1 if chare desires to reactivate
    int m_toreactivate;
    //! Flag: 1 if chare is deactivated, 0 if active
    int m_deactivated;
    //! Deactived chares this chare communicates with
    std::unordered_set< int > m_inactive;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Chare-boundary triangle face connecitivity
    std::vector< std::size_t > m_triinpoel;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Max/min antidiffusive edge contributions at mesh nodes
    tk::Fields m_p;
    //! Receive buffer for max/min antidiffusive edge contributions
    //! \details Key: global node id, value: max/min antidiff edge contributions
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
    //! Boundary edge integrals
    //! \details Key: boundary edge-end points with global node ids, value:
    //!   boundary edge integral contributions.
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_bndedgeint;
    //! Domain edge integrals
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 4 >,
      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_domedgeint;
    //! Streamable boundary point local ids
    std::vector< std::size_t > m_bpoin;
    //! Streamable boundary point integrals
    std::vector< tk::real > m_bpint;
    //! Superedge (face, edge) end points with local ids for boundary edges
    std::array< std::vector< std::size_t >, 2 > m_bsupedge;
    //! Superedge (tet, face, edge) boundary edge integrals
    std::array< std::vector< tk::real >, 2 > m_bsupint;
    //! Superedge (tet, face, edge) end points with local ids for domain edges
    std::array< std::vector< std::size_t >, 3 > m_dsupedge;
    //! Superedge (tet, face, edge) domain edge integrals
    std::array< std::vector< tk::real >, 3 > m_dsupint;
    //! Chare-boundary edge end-points with integrals
    //! \details Outer key: neighbor chare id, value: contents of domain-edge
    //!   integral associated to the edge with local node ids
    std::unordered_map< int, decltype(m_domedgeint) > m_chbndedge;
    //! Streamable boundary point symmetry BC flags
    std::vector< std::uint8_t > m_besym;
    //! Gradients in mesh nodes
    tk::Fields m_grad;
    //! Gradients receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_gradc;
    //! Stabilization coefficients in mesh nodes
    tk::Fields m_stab;
    //! Stabilization coefficients receive buffer
    std::unordered_map< std::size_t, tk::real > m_stabc;
    //! Nodal volumes dynamically adjusted for deactivated chares
    std::vector< tk::real > m_vol;
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

    //! Compute local contributions to domain edge integrals
    void domint();

    //! Compute chare-boundary edges
    void bndEdges();

    //! Compute boundary point normals
    void bndint();

    //! Combine own and communicated portions of the boundary point normals
    void bnorm();

    //! Convert integrals into streamable data structures
    void streamable();

    //! Generate superedge-groups for boundary-edge integrals
    void bndsuped();

    //! Generate superedge-groups for domain-edge integral
    void domsuped();

    //! Generate edges along chare boundary
    void chbnded();

    //! Apply diffusion on active hull
    void huldif();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback cb );

    //! Combine own and communicated portions of the integrals
    void merge();

    //! Compute next time step
    void compute();

    //! Compute gradients for next time step
    void grad();

    //! Compute stabilization coefficients for next time step
    void stab();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Continue with flux-corrected transport if enabled
    void fct();

    //! Compute antidiffusive contributions: P+/-
    void aec();

    //! Compute allowed limits, Q+/-
    void alw();

    //! Compute limit coefficients
    void lim();

    //! Advance systems of equations
    void solve();

    //! Adjust node volumes along inactive neighbor chares
    void deavol();

    //! Decide if edge is active
    int active( std::size_t p,
                std::size_t q,
                tk::real tol,
                const std::vector< uint64_t >& sys );

    //! Decide whether to deactivate this chare
    int dea( const std::vector< uint64_t >& sys );

    //! Decide whether to teactivate a neighbor chare
    std::unordered_map< int, int > rea( const std::vector< uint64_t >& sys );

    //! Deactivate regions
    void deactivate();

    //! Refine/derefine mesh
    void refine();

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
