// *****************************************************************************
/*!
  \file      src/Inciter/AirCG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     AirCG for a PDE system with continuous Galerkin + ALE + RK
  \details   AirCG advances a system of partial differential equations (PDEs)
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme in the arbitrary Eulerian-Lagrangian
    reference frame.

    There are a potentially large number of AirCG Charm++ chares created by
    Transporter. Each AirCG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    ALE time-stepping is performed in an unsplit fashion, as opposed to
    Lagrange + remap. See also J. Waltz, N.R. Morgan, T.R. Canfield, M.R.J.
    Charest, L.D. Risinger, J.G. Wohlbier, A three-dimensional finite element
    arbitrary Lagrangian–Eulerian method for shock hydrodynamics on unstructured
    grids, Computers & Fluids, 92: 172-187, 2014.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
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
#include "Inciter/InputDeck/InputDeck.hpp"

#include "NoWarning/aircg.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! AirCG Charm++ chare array used to advance PDEs in time with AirCG+RK
class AirCG : public CBase_AirCG {

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
    AirCG_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit AirCG( const CProxy_Discretization& disc,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit AirCG( CkMigrateMessage* msg ) : CBase_AirCG( msg ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Return from migration
    void ResumeFromSync() override;

    //! Setup node-neighborhood (no-op)
    void nodeNeighSetup() {}

    //! Start setup for solution
    void setup();

    //! Receive total box IC volume and set conditions in box
    void box( tk::real v );

    // Start time stepping
    void start();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Start (re-)computing domain and boundary integrals
    void integrals();

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
      const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
      const tk::NodeCommMap& nodeCommMap,
      const std::map< int, std::vector< std::size_t > >& bface,
      const std::map< int, std::vector< std::size_t > >& bnode,
      const std::vector< std::size_t >& triinpoel );

    //! Extract field output to file
    void extractFieldOutput(
      const std::vector< std::size_t >& /* ginpoel */,
      const tk::UnsMesh::Chunk& /*chunk*/,
      const tk::UnsMesh::Coords& /*coord*/,
      const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /* addedNodes */,
      const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
      const tk::NodeCommMap& /*nodeCommMap*/,
      const std::map< int, std::vector< std::size_t > >& /*bface*/,
      const std::map< int, std::vector< std::size_t > >& /* bnode */,
      const std::vector< std::size_t >& /*triinpoel*/,
      CkCallback /*c*/ ) {}

    //! Const-ref access to current solution
    //! \return Const-ref to current solution
    const tk::Fields& solution() const { return m_u; }

    //! Resizing data sutrctures after mesh refinement has been completed
    void resized();

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
      p | m_ncomp;
      p | m_bnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_bndel;
      p | m_bpoinid;
      p | m_bpoinin;
      p | m_u;
      p | m_un;
      p | m_rhs;
      p | m_rhsc;
      p | m_diag;
      p | m_bnorm;
      p | m_bnormc;
      p | m_bndpoinint;
      p | m_bndpoinintc;
      p | m_bndedgeint;
      p | m_bndedgeintc;
      p | m_domedgeint;
      p | m_domedgeintc;
      p | m_bpoin;
      p | m_bpint;
      p | m_bedge;
      p | m_beint;
      p | m_dedge;
      p | m_deint;
      p | m_bpsym;
      p | m_besym;
      p | m_grad;
      p | m_gradc;
      p | m_dirbcnodes;
      p | m_symbcnodeset;
      p | m_symbcnodes;
      p | m_symbcnorms;
      p | m_farbcnodeset;
      p | m_farbcnodes;
      p | m_farbcnorms;
      p | m_stage;
      p | m_boxnodes;
      p | m_edgenode;
      p | m_dtp;
      p | m_tp;
      p | m_finished;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i AirCG object reference
    friend void operator|( PUP::er& p, AirCG& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

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
    //! Number of scalar components (flow:5 + transported scalars)
    std::size_t m_ncomp;
    //! Boundary node lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Boundary face lists mapped to side set ids used in the input file
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary triangle face connecitivity where BCs are set by user
    std::vector< std::size_t > m_triinpoel;
    //! Elements along mesh boundary
    std::vector< std::size_t > m_bndel;
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
    //!   point, value: weighted normals and inverse distance square.
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
    //! Boundary point integrals receive buffer
    //! \details Key: global node id of boundary point, value: boundary point
    //!   integral contributions.
    decltype(m_bndpoinint) m_bndpoinintc;
    //! Boundary edge integrals
    //! \details Key: boundary edge-end points with global node ids, value:
    //!   boundary edge integral contributions.
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_bndedgeint;
    //! Boundary edge integrals receive buffer
    //! \details Key: boundary edge-end points with global node ids, value:
    //!   boundary edge integral contributions.
    decltype(m_bndedgeint) m_bndedgeintc;
    //! Domain edge integrals
    std::unordered_map< tk::UnsMesh::Edge, std::array< tk::real, 3 >,
      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_domedgeint;
    //! Receive buffer for domain edge integrals along chare-boundary edges
    decltype(m_domedgeint) m_domedgeintc;
    //! Streamable boundary point local ids
    std::vector< std::size_t > m_bpoin;
    //! Streamable boundary point integrals
    std::vector< tk::real > m_bpint;
    //! Streamable boundary edges with local ids
    std::vector< std::size_t > m_bedge;
    //! Streamable boundary edge integrals
    std::vector< tk::real > m_beint;
    //! Streamable domain edge end points with local ids
    std::vector< std::size_t > m_dedge;
    //! Streamable domain edge integrals
    std::vector< tk::real > m_deint;
    //! Streamable boundary point symmetry BC flags
    std::vector< std::uint8_t > m_bpsym;
    //! Streamable boundary edge symmetry BC flags
    std::vector< std::uint8_t > m_besym;
    //! Gradients in mesh nodes
    tk::Fields m_grad;
    //! Gradients receive buffer
    std::unordered_map< std::size_t, std::vector< tk::real > > m_gradc;
    //! Streamable nodes at which Dirichlet BCs are set
    std::vector< std::size_t > m_dirbcnodes;
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
    //! Runge-Kutta stage counter
    std::size_t m_stage;
    //! Mesh node ids at which user-defined box ICs are defined (multiple boxes)
    std::vector< std::unordered_set< std::size_t > > m_boxnodes;
    //! Local node IDs of edges
    std::vector< std::size_t > m_edgenode;
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

    //! Compute local contributions to domain edge integrals
    void domint();

    //! Compute chare-boundary edges
    void bndEdges();

    //! Compute boundary point normals
    void bndint( const std::unordered_map< int,
                         std::unordered_set< std::size_t > >& bcnodes );

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
    void BC();

    //! Compute gradients for next time step
    void grad();
};

} // inciter::
