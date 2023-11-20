// *****************************************************************************
/*!
  \file      src/Partition/Zoltan2Geom.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan2 library's geometric partitioners
*/
// *****************************************************************************

#include "NoWarning/Zoltan2_PartitioningProblem.hpp"

#include "ZoltanInterOp.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

//! GeometricMeshElemAdapter : Zoltan2::MeshAdapter
//! \details GeometricMeshElemAdapter specializes those virtual member functions
//!   of Zoltan2::MeshAdapter that are required for mesh-element-based
//!   geometric partitioning with Zoltan2
template< typename ZoltanTypes >
class GeometricMeshElemAdapter : public Zoltan2::MeshAdapter< ZoltanTypes > {

  private:
    using MeshEntityType = Zoltan2::MeshEntityType;
    using EntityTopologyType = Zoltan2::EntityTopologyType;

  public:
    using gno_t = typename Zoltan2::InputTraits< ZoltanTypes >::gno_t;
    using scalar_t = typename Zoltan2::InputTraits< ZoltanTypes >::scalar_t;
    using base_adapter_t = Zoltan2::MeshAdapter< ZoltanTypes >;

    //! Constructor
    //! \param[in] nelem Number of elements in mesh graph on this rank
    //! \param[in] centroid Mesh element coordinates (centroids)
    //! \param[in] elemid Mesh element global IDs
    GeometricMeshElemAdapter(
      std::size_t nelem,
      const std::array< std::vector< tk::real >, 3 >& centroid,
      const std::vector< long >& elemid )
    : m_nelem( nelem ),
      m_topology( EntityTopologyType::TETRAHEDRON ),
      m_centroid( centroid ),
      m_elemid( elemid )
    {}

    //! Returns the number of mesh entities on this rank
    //! \return Number of mesh elements on this rank
    // cppcheck-suppress unusedFunction
    std::size_t getLocalNumOf( MeshEntityType ) const override
    { return m_nelem; }

    //! Provide a pointer to this rank's identifiers
    //! \param[in,out] Ids Pointer to the list of global element Ids on this
    //!   rank
    // cppcheck-suppress unusedFunction
    void getIDsViewOf( MeshEntityType, const gno_t*& Ids) const override
    { Ids = m_elemid.data(); }

    //! Return dimensionality of the mesh
    //! \return Number of mesh dimension
    // cppcheck-suppress unusedFunction
    int getDimension() const override { return 3; }

    //! Provide a pointer to one dimension of mesh element coordinates
    //! \param[in] coords Pointer to a list of coordinate values for the
    //!   dimension
    //! \param[in,out] stride Describes the layout of the coordinate values in
    //!   the coords list. If stride is one, then the ith coordinate value is
    //!   coords[i], but if stride is two, then the ith coordinate value is
    //!   coords[2*i]
    //! \param dim Value from 0 to one less than getEntityCoordinateDimension()
    //!   specifying which dimension is being provided in the coords list
    // cppcheck-suppress unusedFunction
    void getCoordinatesViewOf( MeshEntityType,
                               const scalar_t*& coords,
                               int &stride,
                               int dim ) const override
    {
      coords = m_centroid[ static_cast<std::size_t>(dim) ].data();
      stride = 1;
    }

  private:
    //! Number of elements on this rank
    const std::size_t m_nelem;
    //! Mesh element topology types
    const EntityTopologyType m_topology;
    //! Mesh element coordinates (centroids)
    const std::array< std::vector< tk::real >, 3 >& m_centroid;
    //! Global mesh element ids
    const std::vector< long >& m_elemid;
};

static
std::array< std::vector< tk::real >, 3 >
centroids( const std::vector< std::size_t >& inpoel,
           const std::array< std::vector< tk::real >, 3 >& coord )
// *****************************************************************************
//  Compute element centroid coordinates
//! \param[in] inpoel Mesh connectivity with local ids
//! \param[in] coord Node coordinates
//! \return Centroids for all cells on this compute node
// *****************************************************************************
{
  Assert( tk::uniquecopy(inpoel).size() == coord[0].size(), "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Make room for element centroid coordinates
  std::array< std::vector< tk::real >, 3 > cent;
  auto& cx = cent[0];
  auto& cy = cent[1];
  auto& cz = cent[2];
  auto num = inpoel.size()/4;
  cx.resize( num );
  cy.resize( num );
  cz.resize( num );

  // Compute element centroids for mesh passed in
  for (std::size_t e=0; e<num; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    cx[e] = (x[A] + x[B] + x[C] + x[D]) / 4.0;
    cy[e] = (y[A] + y[B] + y[C] + y[D]) / 4.0;
    cz[e] = (z[A] + z[B] + z[C] + z[D]) / 4.0;
  }

  return cent;
}

std::vector< std::size_t >
geomPartMeshZ2( const char* alg,
                const std::vector< std::string >& zoltan_params,
                const std::vector< std::size_t >& inpoel,
                const std::array< std::vector< tk::real >, 3 >& coord,
                int npart )
// *****************************************************************************
//  Partition mesh using Zoltan2 with a geometric partitioner
//! \param[in] alg Partitioning algorithm to use
//! \param[in] zoltan_params Extra parameters pass to zoltan
//! \param[in] inpoel Mesh connectivity with local ids
//! \param[in] coord Node coordinates
//! \param[in] npart Number of desired partitions
//! \return Array of chare ownership IDs mapping elements to chares
//! \details This function uses Zoltan to partition the mesh in parallel.
//!   It assumes that the mesh is distributed among all the MPI ranks.
// *****************************************************************************
{
  // Set Zoltan parameters
  Teuchos::ParameterList params( "Zoltan parameters" );
  params.set( "algorithm", alg );
  params.set( "num_global_parts", npart );
  params.set( "objects_to_partition", "mesh_elements" );

  for (std::size_t i=0; i<zoltan_params.size()/2; ++i) {
    const auto p = zoltan_params.data() + i*2;
    params.set( p[0], p[1] );
  }

  // Define types for Zoltan2
  //  * 1st argument, 'scalar': the data type for element values, weights and
  //    coordinates
  //  * 2nd argument, 'lno' (local number): the integral data type used by
  //    the application and by Zoltan2 for local indices and local
  //    counts
  //  * 3rd argument 'gno' (global number): is the integral data type used by
  //    the application and Zoltan2 to represent global
  //    identifiers and global counts
  // See also
  // external/src/trilinos/packages/zoltan2/src/input/Zoltan2_InputTraits.hpp
  using ZoltanTypes = Zoltan2::BasicUserTypes< tk::real, long, long >;

  auto cent = centroids( inpoel, coord );
  std::vector< long > elemid( inpoel.size()/4 );
  std::iota( begin(elemid), end(elemid), 0 );

  // Create mesh adapter for Zoltan for mesh element partitioning
  using InciterZoltanAdapter = GeometricMeshElemAdapter< ZoltanTypes >;
  InciterZoltanAdapter adapter( elemid.size(), cent, elemid );

  // Create Zoltan2 partitioning problem using our mesh input adapter
  Zoltan2::PartitioningProblem< InciterZoltanAdapter >
    partitioner( &adapter, &params );

  // Perform partitioning using Zoltan
  partitioner.solve();

  // Copy over array of chare IDs corresponding to the ownership of elements
  // in our chunk of the mesh graph, i.e., the coloring or chare ids for the
  // mesh elements we operated on
  auto partlist = partitioner.getSolution().getPartListView();
  std::vector< std::size_t > chare( elemid.size() );
  for (std::size_t p=0; p<elemid.size(); ++p )
    chare[p] = static_cast< std::size_t >( partlist[p] );

  return chare;
}

} // inciter::
