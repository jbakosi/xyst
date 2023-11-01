// *****************************************************************************
/*!
  \file      src/Inciter/RieCG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     RieCG: Riemann, MUSCL, Runge-Kutta, edge-based continuous Galerkin
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

#include "NoWarning/riecg.decl.h"

namespace inciter {

//! RieCG Charm++ chare array used to advance PDEs in time with RieCG+RK
class RieCG : public CBase_RieCG {

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
    RieCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit RieCG( const CProxy_Discretization& disc,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit RieCG( CkMigrateMessage* m ) : CBase_RieCG( m ) {}
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
      std::unordered_map< std::size_t, std::array<tk::real,4> > >& inbnd );

    //! Receive contributions to node gradients on chare-boundaries
    void comgrad( const std::unordered_map< std::size_t,
                          std::vector< tk::real > >& ingrad );

    //! Receive contributions to right-hand side vector on chare-boundaries
    void comrhs( const std::unordered_map< std::size_t,
                         std::vector< tk::real > >& inrhs );

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

    // Evaluate whether to do load balancing
    void evalLB( int nrestart );

    //! Evaluate whether to continue with next time step stage
    void stage();

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
      p | m_nbpint;
      p | m_nbeint;
      p | m_ndeint;
      p | m_ngrad;
      p | m_bnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_bpoinid;
      p | m_bpoinin;
      p | m_u;
      // do not pup these, will recompute after migration anyway
      if (p.isUnpacking()) {
        m_un.resize( m_u.nunk(), m_u.nprop() );
        m_rhs.resize( m_u.nunk(), m_u.nprop() );
        m_grad.resize( m_u.nunk(), m_u.nprop()*3 );
      }
      p | m_rhsc;
      p | m_gradc;
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
      p | m_bpsym;
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
      p | m_stage;
      p | m_dtp;
      p | m_tp;
      p | m_finished;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i RieCG object reference
    friend void operator|( PUP::er& p, RieCG& i ) { i.pup(p); }
    //@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Counter for right-hand side vector nodes updated
    std::size_t m_nrhs;
    //! Counter for receiving boundary point normals
    std::size_t m_nnorm;
    //! Counter for receiving boundary point integrals
    std::size_t m_nbpint;
    //! Counter for receiving boundary edge integrals
    std::size_t m_nbeint;
    //! Counter for receiving domain edge integrals
    std::size_t m_ndeint;
    //! Counter for receiving gradients
    std::size_t m_ngrad;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary triangle face connecitivity where BCs are set by user
    std::vector< std::size_t > m_triinpoel;
    //! Streamable boundary point local ids
    std::vector< std::size_t > m_bpoinid;
    //! Streamable boundary point integrals
    std::vector< tk::real > m_bpoinin;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Unknown/solution vector at mesh nodes at previous time
    tk::Fields m_un;
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
      std::unordered_map< std::size_t, std::array<tk::real,4> > > m_bnorm;
    //! Boundary point normals receive buffer
    //! \details Outer key: side set id. Inner key: global node id of boundary
    //!   point, value: weighted normals and inverse distance square.
    decltype(m_bnorm) m_bnormc;
    //! Boundary point integrals
    //! \details Key: global node id of boundary point, value: boundary point
    //!   integral contributions.
    std::unordered_map< std::size_t, std::array<tk::real,3> > m_bndpoinint;
    //! Boundary edge integrals
    //! \details Key: boundary edge-end points with global node ids, value:
    //!   boundary edge integral contributions.
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_bndedgeint;
    //! Domain edge integrals
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
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
    //! Streamable boundary point symmetry BC flags
    std::vector< std::uint8_t > m_bpsym;
    //! Gradients in mesh nodes
    tk::Fields m_grad;
    //! Gradients receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_gradc;
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
    //! Runge-Kutta stage counter
    std::size_t m_stage;
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

    //! Generate superedge-groups for boundary-edge loops
    void bndsuped();

    //! Generate superedge-groups for domain-edge loops
    void domsuped();

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( CkCallback cb );

    //! Combine own and communicated portions of the integrals
    void merge();

    //! Compute righ-hand side vector of transport equations
    void rhs();

    //! Advance systems of equations
    void solve();

    //! Compute time step size
    void dt();

    //! Evaluate whether to save checkpoint/restart
    void evalRestart();

    //! Apply boundary conditions
    void BC( tk::real t );

    //! Compute gradients for next time step
    void grad();

    //! Apply scalar source to solution
    void src();
};

} // inciter::
