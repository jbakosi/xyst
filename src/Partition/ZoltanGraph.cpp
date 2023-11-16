// *****************************************************************************
/*!
  \file      src/Partition/ZoltanGraph.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library's graph partitioners.
*/
// *****************************************************************************

#include "Compiler.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
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

#include "ZoltanInterOp.hpp"
#include "ContainerUtil.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"

namespace zoltan {

//! Zoltan mesh data structure
struct MESH_DATA {
  int numMyVertices;            //!< number of vertices that I own initially
  ZOLTAN_ID_TYPE* vtxGID;       //!< globa ID of these vertices
  int numMyHEdges;              //!< number of my hyperedges
  int numAllNbors;              //!< number of vertices in my hyperedges
  ZOLTAN_ID_TYPE *edgeGID;      //!< global ID of each of my hyperedges
  int *nborIndex;               //!< index into nborGID array of edge ids
  ZOLTAN_ID_TYPE *nborGID;      //!< array of edge ids
};

static int
get_number_of_objects( void* data, int* ierr ) {
  MESH_DATA* mesh = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;
  return mesh->numMyVertices;
}

static void
get_object_list( void *data, int /* sizeGID */, int /* sizeLID */,
                 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                 int /* wgt_dim */, float* /*obj_wgts */, int* ierr )
{
  MESH_DATA* mesh = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;
  for (int i=0; i<mesh->numMyVertices; ++i){
    globalID[i] = mesh->vtxGID[i];
    localID[i] = static_cast< ZOLTAN_ID_TYPE >( i );
  }
}

static void
get_hypergraph_size( void* data, int* num_lists, int* num_nonzeros,
                     int* format, int* ierr )
{
  MESH_DATA* hg = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;
  *num_lists = hg->numMyHEdges;
  *num_nonzeros = hg->numAllNbors;
  *format = ZOLTAN_COMPRESSED_EDGE;
}

static void
get_hypergraph( void* data, int /* sizeGID */, int num_edges, int num_nonzeros,
                int format, ZOLTAN_ID_PTR edgeGID, int* vtxPtr,
                ZOLTAN_ID_PTR vtxGID, int* ierr )
{
  MESH_DATA* hg = static_cast< MESH_DATA* >( data );
  *ierr = ZOLTAN_OK;

  if ( (num_edges != hg->numMyHEdges) ||
       (num_nonzeros != hg->numAllNbors) ||
       (format != ZOLTAN_COMPRESSED_EDGE) )
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (int i=0; i<num_edges; ++i) {
    edgeGID[i] = hg->edgeGID[i];
    vtxPtr[i] = hg->nborIndex[i];
  }

  for (int i=0; i<num_nonzeros; ++i) vtxGID[i] = hg->nborGID[i];
}

static std::size_t
createHyperGraph( const std::vector< std::size_t >& inpoel,
                  const std::vector< std::size_t >& gid,
                  MESH_DATA& hg )
// *****************************************************************************
//  Create hypergraph data structure for Zoltan
//! \param[in] inpoel Mesh node connectivity with local ids
//! \param[in] gid Global node ids
//! \param[inout] hg Hypergraph data structure to fill
//! \return Number of hyperedges in graph (number of nodes in our mesh chunk)
// *****************************************************************************
{
  // Get number of points from graph
  const auto npoin = gid.size();

  // Create hypergraph data structure based on mesh graph
  hg.numMyVertices = static_cast< int >( npoin );
  hg.numMyHEdges = hg.numMyVertices;
  hg.vtxGID = static_cast< ZOLTAN_ID_PTR >(
                malloc(sizeof(ZOLTAN_ID_TYPE) * npoin) );
  hg.edgeGID = static_cast< ZOLTAN_ID_PTR >(
                 malloc(sizeof(ZOLTAN_ID_TYPE) * npoin) );
  hg.nborIndex = static_cast< int* >( malloc(sizeof(int) * (npoin+1)) );

  // Assign global point ids
  for (std::size_t i=0; i<npoin; ++i)
    hg.vtxGID[i] = static_cast< ZOLTAN_ID_TYPE >( gid[i] );

  // Generate points surrounding points of graph storing local node ids
  const auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  // Allocate data to store the hypergraph ids. The total number of vertices or
  // neighbors in all the hyperedges of the hypergraph, nhedge = all points
  // surrounding points + number of points, since psup does not store the
  // connection to the own point. In other words, here we need the number of
  // edges in the graph, independent of direction.
  auto nhedge = psup.first.size() - 1 + npoin;
  hg.numAllNbors = static_cast< int >( nhedge );
  hg.nborGID = static_cast< ZOLTAN_ID_PTR >(
                 malloc(sizeof(ZOLTAN_ID_TYPE) * nhedge) );

  // Fill up hypergraph edge ids and their indices
  hg.nborIndex[0] = 0;
  for (std::size_t p=0; p<npoin; ++p) {
    hg.edgeGID[p] = static_cast< ZOLTAN_ID_TYPE >( gid[p] );
    hg.nborGID[ hg.nborIndex[p] ] = static_cast< ZOLTAN_ID_TYPE >( gid[p] );
    int j = 0;
    for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i, ++j) {
      hg.nborGID[ hg.nborIndex[p] + 1 + j ] =
        static_cast< ZOLTAN_ID_TYPE >( gid[psup.first[i]] );
    }
    hg.nborIndex[p+1] = hg.nborIndex[p] + j;
  }

  return npoin;
}

static std::vector< std::size_t >
chElem( const std::vector< std::size_t >& chp,
        const std::vector< std::size_t >& inpoel )
// *****************************************************************************
//! Construct array of chare ownership IDs mapping mesh elements to chares
//! \param[in] chp Chares of points: array of chare ownership IDs mapping graph
//!   points to Charm++ chares. Size: number of points in chunk of mesh graph.
//! \param[in] inpoel Mesh tetrahedron element connectivity with local IDs
//! \return Vector of chare IDs mapped to mesh elements
//! \details This function constructs the vector che, 'chares of elements'. This
//!   is the equivalent of chp, 'chares of points', but stores the chare ids of
//!   mesh elements. The element ownership is computed based on the point
//!   ownership. If all points of an element are owned by a chare, then that
//!   chare owns the element. If not all points of an element are owned by the
//!   a chare, than the chare with lower chare id gets to own the element.
//! \note This function operates on all MPI ranks, working on different chunks
//!   of the mesh.
//! \author J. Bakosi
// ******************************************************************************
{
  std::vector< std::size_t > che( inpoel.size()/4 );

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const auto i = inpoel.data() + e*4;
    che[e] = std::min( chp[i[3]],
               std::min( chp[i[2]], std::min( chp[i[0]], chp[i[1]] ) ) );
  }

  return che;
}

std::vector< std::size_t >
graphPartMesh( const std::vector< std::size_t >& ginpoel,
               const std::vector< std::string >& zoltan_params,
               int npart )
// *****************************************************************************
//  Partition mesh using Zoltan with a geometric partitioner
//! \param[in] ginpoel Mesh connectivity with global ids
//! \param[in] zoltan_params Extra parameters pass to zoltan
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

  Zoltan_Initialize( 0, nullptr, &ver );

  zz = Zoltan_Create( MPI_COMM_WORLD );

  Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0" );
  Zoltan_Set_Param( zz, "PHG_OUTPUT_LEVEL", "0" );
  Zoltan_Set_Param( zz, "CHECK_HYPERGRAPH", "0" );
  Zoltan_Set_Param( zz, "LB_METHOD", "HYPERGRAPH" );
  Zoltan_Set_Param( zz, "HYPERGRAPH_PACKAGE", "PHG" );
  Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );
  Zoltan_Set_Param( zz, "PHG_MULTILEVEL", "0" );
  Zoltan_Set_Param( zz, "PHG_CUT_OBJECTIVE", "CONNECTIVITY" );
  Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "RETURN_LISTS", "PART" );
  Zoltan_Set_Param( zz, "RCB_OUTPUT_LEVEL", "0" );
  Zoltan_Set_Param( zz, "NUM_GLOBAL_PARTS", std::to_string(npart).c_str() );

  for (std::size_t i=0; i<zoltan_params.size()/2; ++i) {
    const auto p = zoltan_params.data() + i*2;
    Zoltan_Set_Param( zz, p[0].c_str(), p[1].c_str() );
  }

  // Generate element connectivity storing local node ids
  const auto& [ inpoel, gid, lid ] = tk::global2local( ginpoel );

  MESH_DATA myMesh;
  auto npoin = createHyperGraph( inpoel, gid, myMesh );

  // Set Zoltan query functions
  Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );
  Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );
  Zoltan_Set_HG_Size_CS_Fn( zz, get_hypergraph_size, &myMesh );
  Zoltan_Set_HG_CS_Fn( zz, get_hypergraph, &myMesh );

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

  Assert( numExport == static_cast< int >( npoin ), "Size mismatch" );

  if (myMesh.numMyVertices > 0) free( myMesh.vtxGID );
  if (myMesh.numMyHEdges > 0) free( myMesh.edgeGID );
  if (myMesh.numMyHEdges >= 0) free( myMesh.nborIndex );
  if (myMesh.numAllNbors > 0) free( myMesh.nborGID );

  // Copy over array of chare IDs corresponding to the ownership of vertices in
  // our chunk of the mesh, i.e., the coloring or chare ids for the mesh nodes
  // we operate on.
  std::vector< std::size_t > chp( npoin );
  for (std::size_t p=0; p<static_cast<std::size_t>(numExport); ++p )
    chp[ exportLocalGids[p] ] = static_cast< std::size_t >( exportToPart[p] );

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  Zoltan_LB_Free_Part( &importGlobalGids, &importLocalGids,
                       &importProcs, &importToPart );
  Zoltan_LB_Free_Part( &exportGlobalGids, &exportLocalGids,
                       &exportProcs, &exportToPart );

  Zoltan_Destroy( &zz );

  return chElem( chp, inpoel );
}

} // zoltan::
