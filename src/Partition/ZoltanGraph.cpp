// *****************************************************************************
/*!
  \file      src/Partition/ZoltanGraph.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
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

#include "ZoltanGraph.hpp"
#include "ContainerUtil.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"

namespace inciter {

//! Zoltan mesh data structure
struct MESH_DATA {
  int numMyVertices;            //!< number of vertices that I own initially
  ZOLTAN_ID_TYPE* vtxGID;       //!< global ID of these vertices
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

static void
createHyperGraph( const std::vector< std::size_t >& gid,
                  const std::unordered_map< std::size_t,
                                            std::vector< std::size_t > >& graph,
                  MESH_DATA& hg )
// *****************************************************************************
//  Create hypergraph data structure for Zoltan
//! \param[in] gid Global node ids
//! \param[in] graph Aggregated mesh graph point connectivity
//! \param[inout] hg Hypergraph data structure to fill
//! \return Number of hyperedges in graph (number of nodes in our mesh chunk)
// *****************************************************************************
{
  // Get number of points from graph
  const auto npoin = graph.size();

  // Create hypergraph data structure based on mesh graph
  hg.numMyVertices = static_cast< int >( npoin );
  hg.numMyHEdges = hg.numMyVertices;
  hg.vtxGID = static_cast< ZOLTAN_ID_PTR >(
                malloc(sizeof(ZOLTAN_ID_TYPE) * npoin) );
  hg.edgeGID = static_cast< ZOLTAN_ID_PTR >(
                 malloc(sizeof(ZOLTAN_ID_TYPE) * npoin) );
  hg.nborIndex = static_cast< int* >( malloc(sizeof(int) * (npoin+1)) );

  // generate linked vectors for points surrounding points and their indices
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > > psup;
  auto& psup1 = psup.first;
  auto& psup2 = psup.second;
  psup1.resize( 1, 0 );
  psup2.resize( 1, 0 );

  std::vector< std::size_t > gp;
  for (std::size_t p=0; p<gid.size(); ++p) {
    auto i = graph.find( gid[p] );
    if (i == end(graph)) continue;
    gp.push_back( gid[p] );
    const auto& n = i->second;
    psup2.push_back( psup2.back() + n.size() );
    psup1.insert( end(psup1), begin(n), end(n) );
  }

  Assert( gp.size() == graph.size(), "Size mismatch" );

  // Compute sum of number of vertices of hyperedges and allocate memory
  auto nedge = psup.first.size() - 1 + npoin;
  hg.numAllNbors = static_cast< int >( nedge );
  hg.nborGID = static_cast< ZOLTAN_ID_PTR >(
                 malloc(sizeof(ZOLTAN_ID_TYPE) * nedge) );

   // Fill up hypergraph edge ids and their indices
   hg.nborIndex[0] = 0;
   for (std::size_t p=0; p<npoin; ++p) {
     auto g = static_cast< ZOLTAN_ID_TYPE >( gp[p] );
     hg.vtxGID[p] = hg.edgeGID[p] = hg.nborGID[ hg.nborIndex[p] ] = g;
     int j = 1;
     for (auto i=psup2[p]+1; i<=psup2[p+1]; ++i, ++j) {
       hg.nborGID[ hg.nborIndex[p] + j ] =
         static_cast< ZOLTAN_ID_TYPE >( psup1[i] );
     }
     hg.nborIndex[p+1] = hg.nborIndex[p] + j;
   }
}

std::unordered_map< std::size_t, std::size_t >
graphPartMesh( const std::vector< std::size_t >& ginpoel,
               const std::unordered_map< std::size_t,
                                std::vector< std::size_t > >& graph,
               const std::vector< std::string >& zoltan_params,
               int npart )
// *****************************************************************************
//  Partition mesh using Zoltan with a geometric partitioner
//! \param[in] ginpoel Mesh connectivity with global ids
//! \param[in] graph Mesh graph point connectivity
//! \param[in] zoltan_params Extra parameters pass to zoltan
//! \param[in] npart Number of desired partitions
//! \return Array of chare ownership IDs mapping points owned to chares
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
  Zoltan_Set_Param( zz, "LB_METHOD", "HYPERGRAPH" );
  Zoltan_Set_Param( zz, "HYPERGRAPH_PACKAGE", "PHG" );
  Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );
  Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "EDGE_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "RETURN_LISTS", "PART" );
  Zoltan_Set_Param( zz, "NUM_GLOBAL_PARTS", std::to_string(npart).c_str() );

  for (std::size_t i=0; i<zoltan_params.size()/2; ++i) {
    const auto p = zoltan_params.data() + i*2;
    Zoltan_Set_Param( zz, p[0].c_str(), p[1].c_str() );
  }

  // Generate element connectivity storing local node ids
  const auto& [ inpoel, gid, lid ] = tk::global2local( ginpoel );

  MESH_DATA myMesh;
  createHyperGraph( gid, graph, myMesh );

  // Set Zoltan query functions
  Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );
  Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );
  Zoltan_Set_HG_Size_CS_Fn( zz, get_hypergraph_size, &myMesh );
  Zoltan_Set_HG_CS_Fn( zz, get_hypergraph, &myMesh );
  //Zoltan_Set_HG_Size_CS_Fn( zz, get_hypergraph_size, &myMesh );

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

  Assert( numExport == static_cast< int >( graph.size() ), "Size mismatch" );

  if (myMesh.numMyVertices > 0) free( myMesh.vtxGID );
  if (myMesh.numMyHEdges > 0) free( myMesh.edgeGID );
  if (myMesh.numMyHEdges >= 0) free( myMesh.nborIndex );
  if (myMesh.numAllNbors > 0) free( myMesh.nborGID );

  // Copy over array of chare IDs corresponding to the ownership of vertices in
  // our chunk of the mesh, i.e., the coloring or chare ids for the mesh nodes
  // we operate on.
  std::unordered_map< std::size_t, std::size_t > chp;
  for (std::size_t p=0; p<static_cast<std::size_t>(numExport); ++p ) {
    chp[ exportGlobalGids[p] ] = static_cast< std::size_t >( exportToPart[p] );
  }

  // Free the arrays allocated by Zoltan_LB_Partition
  Zoltan_LB_Free_Part( &importGlobalGids, &importLocalGids,
                       &importProcs, &importToPart );
  Zoltan_LB_Free_Part( &exportGlobalGids, &exportLocalGids,
                       &exportProcs, &exportToPart );
  // Fee the storage allocated for the Zoltan structure
  Zoltan_Destroy( &zz );

  return chp;
}

} // inciter::
