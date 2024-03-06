// *****************************************************************************
/*!
  \file      src/Partition/ZoltanGeom.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library's geometric partitioners
*/
// *****************************************************************************

#include <numeric>

#include "Compiler.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsuggest-override"
  #pragma clang diagnostic ignored "-Wsuggest-destructor-override"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wdollar-in-identifier-extension"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-function-type"
#endif

#include "zoltan.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#include "ZoltanGeom.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

//! Zoltan mesh data structure
struct MESH_DATA {
  int numMyElems;
  ZOLTAN_ID_PTR myGlobalIDs;
  tk::real* x;
  tk::real* y;
  tk::real* z;
};

static int
get_number_of_objects( void* data, int* ierr ) {
  MESH_DATA* mesh = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;
  return mesh->numMyElems;
}

static void
get_object_list( void* data, int /* sizeGID */, int /* sizeLID */,
                 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                 int /* wgt_dim */, float* /*obj_wgts */, int* ierr )
{
  MESH_DATA* mesh = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;
  for (int i=0; i<mesh->numMyElems; ++i){
    globalID[i] = mesh->myGlobalIDs[i];
    localID[i] = static_cast< ZOLTAN_ID_TYPE >( i );
  }
}

static int
get_num_geometry( void* /*data */, int* ierr ) {
  *ierr = ZOLTAN_OK;
  return 3;
}

static void
get_geometry_list( void* data, int sizeGID, int sizeLID,
                   int num_obj, ZOLTAN_ID_PTR /* globalID */,
                   ZOLTAN_ID_PTR /* localID */,
                   int num_dim, double* geom_vec, int* ierr )
{
  MESH_DATA *mesh = static_cast< MESH_DATA* >( data );
  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)){
    *ierr = ZOLTAN_FATAL;
    return;
  }
  *ierr = ZOLTAN_OK;
  for (int i=0; i<num_obj; ++i) {
    geom_vec[3*i+0] = static_cast< tk::real >( mesh->x[i] );
    geom_vec[3*i+1] = static_cast< tk::real >( mesh->y[i] );
    geom_vec[3*i+2] = static_cast< tk::real >( mesh->z[i] );
  }
}

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
geomPartMesh( const char* alg,
              const std::vector< std::string >& zoltan_params,
              const std::vector< std::size_t >& inpoel,
              const std::array< std::vector< tk::real >, 3 >& coord,
              int npart )
// *****************************************************************************
//  Partition mesh using Zoltan with a geometric partitioner
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
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids,
                exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  std::vector< unsigned int > elemid( inpoel.size()/4 );
  std::iota( begin(elemid), end(elemid), 0 );

  MESH_DATA myMesh;
  myMesh.numMyElems = static_cast< int >( elemid.size() );
  myMesh.myGlobalIDs = const_cast< unsigned int* >( elemid.data() );

  auto cent = centroids( inpoel, coord );
  myMesh.x = const_cast< tk::real* >( cent[0].data() );
  myMesh.y = const_cast< tk::real* >( cent[1].data() );
  myMesh.z = const_cast< tk::real* >( cent[2].data() );

  Zoltan_Initialize( 0, nullptr, &ver );

  zz = Zoltan_Create( MPI_COMM_WORLD );

  Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0" );
  Zoltan_Set_Param( zz, "LB_METHOD", alg );
  Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );
  Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "RETURN_LISTS", "PART" );
  Zoltan_Set_Param( zz, "RCB_OUTPUT_LEVEL", "0" );
  Zoltan_Set_Param( zz, "AVERAGE_CUTS", "1" );
  Zoltan_Set_Param( zz, "NUM_GLOBAL_PARTS", std::to_string(npart).c_str() );

  for (std::size_t i=0; i<zoltan_params.size()/2; ++i) {
    const auto p = zoltan_params.data() + i*2;
    Zoltan_Set_Param( zz, p[0].c_str(), p[1].c_str() );
  }

  /* Query functions, to provide geometry to Zoltan */

  Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );
  Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );
  Zoltan_Set_Num_Geom_Fn( zz, get_num_geometry, &myMesh );
  Zoltan_Set_Geom_Multi_Fn( zz, get_geometry_list, &myMesh );

  Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */

  Assert( numExport == static_cast< int >( elemid.size() ), "Size mismatch" );

  // Copy over array of chare IDs corresponding to the ownership of elements in
  // our chunk of the mesh, i.e., the coloring or chare ids for the mesh
  // elements we operate on.
  std::vector< std::size_t > chare( elemid.size() );
  for (std::size_t p=0; p<static_cast<std::size_t>(numExport); ++p )
    chare[ exportLocalGids[p] ] = static_cast< std::size_t >( exportToPart[p] );

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  Zoltan_LB_Free_Part( &importGlobalGids, &importLocalGids,
                       &importProcs, &importToPart );
  Zoltan_LB_Free_Part( &exportGlobalGids, &exportLocalGids,
                       &exportProcs, &exportToPart );

  Zoltan_Destroy( &zz );

  return chare;
}

} // inciter::
