// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare partitioner nodegroup used to perform mesh
             partitioning
  \details   Charm++ chare partitioner nodegroup used to perform mesh read and
             partitioning, one worker per compute node.
*/
// *****************************************************************************
#ifndef Partitioner_h
#define Partitioner_h

#include <array>
#include <stddef.h>

#include "ContainerUtil.hpp"
#include "DerivedData.hpp"
#include "UnsMesh.hpp"
#include "Sorter.hpp"
#include "Refiner.hpp"
#include "Callback.hpp"

#include "NoWarning/partitioner.decl.h"

namespace inciter {

//! Partitioner Charm++ chare nodegroup class
//! \details Instantiations of Partitioner comprise a processor aware Charm++
//!   chare node group. When instantiated, a new object is created on each
//!   compute node and not more (as opposed to individual chares or chare array
//!   object elements). See also the Charm++ interface file partitioner.ci.
class Partitioner : public CBase_Partitioner {

  private:
    //! \brief Mesh data used for categorizing mesh chunks assigned to chares
    //!    after mesh partitioning and before mesh distribution across chares
    using MeshData =
      std::tuple<
        // Tetrahedron (domain element) connectivity
        std::vector< std::size_t >,
        // Boundary face connectivity for each side set
        std::unordered_map< int, std::vector< std::size_t > >,
        // Boundary node lists for each side set
        std::unordered_map< int, std::vector< std::size_t > > >;

  public:
    //! Constructor
    Partitioner( std::size_t meshid,
                 const std::string& filename,
                 const tk::PartitionerCallback& cbp,
                 const tk::RefinerCallback& cbr,
                 const tk::SorterCallback& cbs,
                 const CProxy_Transporter& host,
                 const CProxy_Refiner& refiner,
                 const CProxy_Sorter& sorter,
                 const tk::CProxy_MeshWriter& meshwriter,
                 const std::vector< CProxy_Discretization >& discretization,
                 const CProxy_RieCG& riecg,
                 const CProxy_LaxCG& laxcg,
                 const CProxy_ZalCG& zalcg,
                 const CProxy_KozCG& kozcg,
                 const CProxy_ChoCG& chocg,
                 const CProxy_LohCG& lohcg,
                 const tk::CProxy_ConjugateGradients& cgpre,
                 const tk::CProxy_ConjugateGradients& cgmom,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::map< int, std::vector< std::size_t > >& faces,
                 const std::map< int, std::vector< std::size_t > >& bnode );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Partitioner( CkMigrateMessage* m ) : CBase_Partitioner( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Aggregate mesh graph for mesh nodes owned
    void query( int fromnode,
                const std::unordered_map< std::size_t,
                        std::unordered_set< std::size_t > >& psup );
    //! Receive receipt of list of points surrounding points to query
    void recvquery();
    //! Respond to graph queries
    void response();
    //! Receive mesh graphs for our mesh chunk
    void psup( int fromnode, const std::unordered_map< std::size_t,
                                std::unordered_set< std::size_t > >& graph );
    //! Receive receipt of mesh nodes graphs
    void recvpsup();
    // Compute total load across distributed mesh
    void load();

    //! Reduction target to aggregate mesh partition assginments
    void parts( CkReductionMsg* msg );

    //! Partition the computational mesh into a number of chares
    void partition( int nchare );

    //! Receive mesh associated to chares we own after refinement
    void addMesh( int fromnode,
                  const std::unordered_map< int,
                    std::tuple<
                      std::vector< std::size_t >,
                      tk::UnsMesh::CoordMap,
                      std::unordered_map< int, std::vector< std::size_t > >,
                      std::unordered_map< int, std::vector< std::size_t > >
                    > >& chmesh );

    //! Acknowledge received mesh after initial mesh refinement
    void recvMesh();

    //! Optionally start refining the mesh
    void refine();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ nodegroup, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_meshid;
      p | m_npsup;
      p | m_cbp;
      p | m_cbr;
      p | m_cbs;
      p | m_host;
      p | m_refiner;
      p | m_sorter;
      p | m_meshwriter;
      p | m_discretization;
      p | m_riecg;
      p | m_laxcg;
      p | m_zalcg;
      p | m_kozcg;
      p | m_chocg;
      p | m_lohcg;
      p | m_cgpre;
      p | m_cgmom;
      p | m_ginpoel;
      p | m_graph;
      p | m_graphnode;
      p | m_coord;
      p | m_inpoel;
      p | m_lid;
      p | m_ndist;
      p | m_nchare;
      p | m_nface;
      p | m_linnodes;
      p | m_chinpoel;
      p | m_chcoordmap;
      p | m_chbface;
      p | m_chtriinpoel;
      p | m_chbnode;
      p | m_bface;
      p | m_triinpoel;
      p | m_bnode;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Partitioner object reference
    friend void operator|( PUP::er& p, Partitioner& i ) { i.pup(p); }
    //@}

  private:
    //! Mesh ID
    std::size_t m_meshid;
    //! Counter for number of nodes contributing to mesn graphs
    std::size_t m_npsup;
    //! Charm++ callbacks associated to compile-time tags for partitioner
    tk::PartitionerCallback m_cbp;
    //! Charm++ callbacks associated to compile-time tags for refiner
    tk::RefinerCallback m_cbr;
    //! Charm++ callbacks associated to compile-time tags for sorter
    tk::SorterCallback m_cbs;
    //! Host proxy
    CProxy_Transporter m_host;
    //! Mesh refiner proxy
    CProxy_Refiner m_refiner;
    //! Mesh sorter proxy
    CProxy_Sorter m_sorter;
    //! Mesh writer proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Discretization proxy for all meshes
    std::vector< CProxy_Discretization > m_discretization;
    //! Discretization scheme proxy
    CProxy_RieCG m_riecg;
    //! Discretization scheme proxy
    CProxy_LaxCG m_laxcg;
    //! Discretization scheme proxy
    CProxy_ZalCG m_zalcg;
    //! Discretization scheme proxy
    CProxy_KozCG m_kozcg;
    //! Discretization scheme proxy
    CProxy_ChoCG m_chocg;
    //! Discretization scheme proxy
    CProxy_LohCG m_lohcg;
    //! Conjugate Gradients Charm++ proxy for pressure solve
    tk::CProxy_ConjugateGradients m_cgpre;
    //! Conjugate Gradients Charm++ proxy for momentum solve
    tk::CProxy_ConjugateGradients m_cgmom;
    //! Element connectivity of this compute node's mesh chunk (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Aggregated mesh graph of owned nodes if graph-based partitioner is used
    std::unordered_map< std::size_t, std::vector< std::size_t > > m_graph;
    //! Graph->node map used to aggregate mesh graph
    //! \details Key: global mesh node id, value: 0: list of compute nodes
    //! contributing partial graph (points surround points) to the global mesh
    //! node id in key, 1: list of global mesh node ids surrounding node in key
    std::unordered_map< std::size_t,
      std::tuple< std::vector< int >,
                  std::unordered_set< std::size_t > > > m_graphnode;
    //! Coordinates of mesh nodes of this compute node's mesh chunk
    tk::UnsMesh::Coords m_coord;
    //! \brief Element connectivity with local node IDs of this compute node's
    //!   mesh chunk
    std::vector< std::size_t > m_inpoel;
    //! Global->local node IDs of elements of this compute node's mesh chunk
    //! \details Key: global node id, value: local node id
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Counter during mesh distribution
    std::size_t m_ndist;
    //! Total number of chares across all compute nodes
    int m_nchare;
    //! Counters (for each chare owned) for assigning face ids in parallel
    std::unordered_map< int, std::size_t > m_nface;
    //! \brief Map associating new node IDs (as in producing contiguous-row-id
    //!   linear system contributions) as map-values to old node IDs (as in
    //!   file) as map-keys
    std::unordered_map< std::size_t, std::size_t > m_linnodes;
    //! Mesh connectivity using global node IDs associated to chares owned
    std::unordered_map< int, std::vector< std::size_t > > m_chinpoel;
    //! Coordinates associated to global node IDs of our mesh chunk for chares
    std::unordered_map< int, tk::UnsMesh::CoordMap > m_chcoordmap;
    //! Side set id + boundary face id for each chare
    std::unordered_map< int,
      std::map< int, std::vector< std::size_t > > > m_chbface;
    //! Boundary face connectivity for each chare
    std::map< int, std::vector< std::size_t > > m_chtriinpoel;
    //! Side set id + boundary nodes for each chare
    std::unordered_map< int,
      std::map< int, std::vector< std::size_t > > > m_chbnode;
    //! Boundary face IDs associated associated to side set IDs
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;

    //!  Categorize mesh elements (given by their gobal node IDs) by target
    std::unordered_map< int, MeshData >
    categorize( const std::vector< std::size_t >& che ) const;

    //! Extract coordinates associated to global nodes of a mesh chunk
    tk::UnsMesh::CoordMap coordmap( const std::vector< std::size_t >& inpoel );

    //! Distribute mesh to target compute nodes after mesh partitioning
    void distribute( std::unordered_map< int, MeshData >&& mesh );

    //! Compute chare (partition) distribution across compute nodes
    std::array< int, 2 > distribution( int npart ) const;

    //! Return nodegroup id for chare id
    int node( int id ) const;

    //! Continue after partitioning finished
    void partitioned( std::vector< std::size_t >&& che );
};

} // inciter::

#endif // Partitioner_h
