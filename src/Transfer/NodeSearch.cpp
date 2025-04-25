// *****************************************************************************
/*!
  \file      src/Transfer/NodeSearch.cpp
  \copyright 2020-2021 Charmworks, Inc.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions pertaining to node search between 3d meshes
*/
// *****************************************************************************

#include <cmath>

#include "NodeSearch.hpp"
#include "Transfer.hpp"
#include "DerivedData.hpp"
#include "InciterConfig.hpp"

#include "NoWarning/collidecharm.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

PUPbytes( Collision );

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#include "NoWarning/transfer.decl.h"

namespace transfer {
  extern CollideHandle g_collideHandle;
  extern CProxy_Transfer g_transferProxy;
}

namespace inciter {
  extern ctr::Config g_cfg;
}

using transfer::NodeSearch;

NodeSearch::NodeSearch( CkArrayID p, MeshData mesh, CkCallback cb )
  : m_firstchunk( mesh.firstchunk )
// *****************************************************************************
//  Constructor
//! \param[in] firstchunk Chunk ID used for the collision detection library
//! \param[in] cb Callback to inform application that the library is ready
// *****************************************************************************
{
//std::cout << "NodeSearch: " << mesh.nchare << '\n';
  mesh.proxy = thisProxy;
  CollideRegister( g_collideHandle, /*ignored*/0 );
  g_transferProxy.ckLocalBranch()->setMesh( p, mesh );
  contribute( cb );
}

void
NodeSearch::setSourceTets( const std::vector< std::size_t >& inpoel,
                           const std::array< std::vector< double >, 3 >& coord,
                           const tk::Fields& u,
                           const std::vector< double >& flag,
                           bool dir,
                           CkCallback cb )
// *****************************************************************************
//  Set the data for the source tetrahedrons to be collided
//! \param[in] inpoel Pointer to the connectivity data for the source mesh
//! \param[in] coord Pointer to the coordinate data for the source mesh
//! \param[in] u Pointer to the solution data for the source mesh
//! \param[in] flag Transfer flags
//! \param[in] dir Transfer direction: 0: bg to overset, 1: overset to bg
//! \param[in] cb Callback to call when src side of transfer is done
// *****************************************************************************
{
  m_inpoel = const_cast< std::vector< std::size_t >* >( &inpoel );
  m_coord = const_cast< std::array< std::vector< double >, 3 >* >( &coord );
  m_u = const_cast< tk::Fields* >( &u );
  m_flag = const_cast< std::vector< double >* >( &flag );
  m_dir = dir;
  m_done = cb;
  m_srcnotified = 0;

  // Send tetrahedron data to the collision detection library
  collideTets();
}

void
NodeSearch::setDestPoints( const std::array< std::vector< double >, 3 >& coord,
                           tk::Fields& u,
                           std::vector< double >& flag,
                           bool trflag,
                           bool dir,
                           CkCallback cb )
// *****************************************************************************
//  Set the data for the destination points to be collided
//! \param[in] coord Pointer to the coordinate data for the destination mesh
//! \param[in,out] u Pointer to the solution data for the destination mesh
//! \param[in,out] flag Transfer flags
//! \param[in] trflag Transfer flags if true
//! \param[in] dir Transfer direction: 0: bg to overset, 1: overset to bg
//! \param[in] cb Callback to call when dst side of transfer is done
// *****************************************************************************
{
  m_coord = const_cast< std::array< std::vector< double >, 3 >* >( &coord );
  m_u = static_cast< tk::Fields* >( &u );
  m_flag = const_cast< std::vector< double >* >( &flag );
  m_trflag = trflag;
  m_dir = dir;
  m_done = cb;

  // Initialize msg counters, callback, and background solution data
  m_numsent = 0;
  m_numreceived = 0;
  //background();

  // Send vertex data to the collision detection library
  collideVertices();
}

void
NodeSearch::background()
// *****************************************************************************
// Initialize dest mesh solution with background data
//! \details This is useful to see what points did not receive solution.
// *****************************************************************************
{
  tk::Fields& u = *m_u;
  for (std::size_t i = 0; i < u.nunk(); ++i) u(i,0) = -1.0;
}

void
NodeSearch::collideVertices()
// *****************************************************************************
// Pass vertex information to the collision detection library
// *****************************************************************************
{
  const std::array< std::vector< double >, 3 >& coord = *m_coord;
  auto nVertices = coord[0].size();
  std::size_t nBoxes = 0;
  std::vector< bbox3d > boxes( nVertices );
  std::vector< int > prio( nVertices );
  auto firstchunk = static_cast< int >( m_firstchunk );
  for (std::size_t i=0; i<nVertices; ++i) {
    boxes[nBoxes].empty();
    boxes[nBoxes].add(CkVector3d(coord[0][i], coord[1][i], coord[2][i]));
    prio[nBoxes] = firstchunk;
    ++nBoxes;
  }

  CollideBoxesPrio( g_collideHandle, firstchunk + thisIndex,
                    static_cast<int>(nBoxes), boxes.data(), prio.data() );
}

void
NodeSearch::collideTets() const
// *****************************************************************************
// Pass tet information to the collision detection library
// *****************************************************************************
{
  const std::vector< std::size_t >& inpoel = *m_inpoel;
  const std::array< std::vector< double >, 3 >& coord = *m_coord;
  auto nBoxes = inpoel.size() / 4;
  std::vector< bbox3d > boxes( nBoxes );
  std::vector< int > prio( nBoxes );
  auto firstchunk = static_cast< int >( m_firstchunk );
  for (std::size_t i=0; i<nBoxes; ++i) {
    boxes[i].empty();
    prio[i] = firstchunk;
    for (std::size_t j=0; j<4; ++j) {
      // Get index of the jth point of the ith tet
      auto p = inpoel[i*4+j];
      // Add that point to the tets bounding box
      boxes[i].add( CkVector3d( coord[0][p], coord[1][p], coord[2][p] ) );
    }
  }

  CollideBoxesPrio( g_collideHandle, firstchunk + thisIndex,
                    static_cast<int>(nBoxes), boxes.data(), prio.data() );
}

void
NodeSearch::processCollisions( const MeshData& src,
                               int nColl,
                               Collision* colls )
// *****************************************************************************
//  Process potential collisions by sending my points to the source mesh chares
//  that they potentially collide with.
//! \param[in] src Source mesh config data
//! \param[in] nColl Number of potential collisions to process
//! \param[in] colls List of potential collisions
// *****************************************************************************
{
  const std::array< std::vector< double >, 3 >& coord = *m_coord;
  int mychunk = m_firstchunk + thisIndex;

  std::vector< std::vector< PotentialCollision > >
    pColls( static_cast<std::size_t>(src.nchare) );

  // Separate potential collisions into lists based on the source mesh chare
  // that is involved in the potential collision
  for (int i=0; i<nColl; ++i) {
    int chareindex;
    PotentialCollision pColl;
    if (colls[i].A.chunk == mychunk) {
      chareindex = colls[i].B.chunk - src.firstchunk;
      pColl.dst = static_cast<std::size_t>(colls[i].A.number);
      pColl.src = static_cast<std::size_t>(colls[i].B.number);
    } else {
      chareindex = colls[i].A.chunk - src.firstchunk;
      pColl.dst = static_cast<std::size_t>(colls[i].B.number);
      pColl.src = static_cast<std::size_t>(colls[i].A.number);
    }

    #if defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wdeprecated-copy"
    #endif
    auto dst = pColl.dst;
    pColl.point = { coord[0][dst], coord[1][dst], coord[2][dst] };
    #if defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #endif

    pColls[ static_cast<std::size_t>(chareindex) ].push_back( pColl );
  }

  // Send out the lists of potential collisions to the source mesh chares
  for (int i=0; i<src.nchare; ++i) {
    auto I = static_cast< std::size_t >( i );
    m_numsent++;
    src.proxy[i].determineActualCollisions( thisProxy, thisIndex,
      static_cast<int>(pColls[I].size()), pColls[I].data() );
  }
}

void
NodeSearch::determineActualCollisions( CProxy_NodeSearch proxy,
                                       int index,
                                       int nColls,
                                       PotentialCollision* colls )
// *****************************************************************************
//  Identify actual collisions by calling intet on all possible collisions, and
//  interpolate solution values to send back to the destination mesh.
//! \param[in] proxy The proxy of the destination mesh chare array
//! \param[in] index The index in proxy to return the solution data to
//! \param[in] nColls Number of collisions to be checked
//! \param[in] colls List of potential collisions
// *****************************************************************************
{
  const std::vector< std::size_t >& inpoel = *m_inpoel;
  const tk::Fields& u = *m_u;
  const std::vector< double >& f = *m_flag;
  //CkPrintf( "Source chare %i received data for %i potential collisions\n",
  //          thisIndex, nColls);

  // Slightly shift dest mesh if symmetry is specified
  auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto& sym = inciter::g_cfg.get< tag::overset, tag::sym_ >()[ 1 ];
  std::array< double, 3 > shift{ 0.0, 0.0, sym=="z" ? eps : 0.0 };

  // Iterate over my potential collisions and call intet() to determine
  // if an actual collision occurred, and if so interpolate solution to dest
  int numInTet = 0;
  std::vector< SolutionData > exp;
  for (int i=0; i<nColls; ++i) {
    const auto& p = colls[i].point;
    std::vector< double > point{ p.x, p.y, p.z };
    std::array< double, 4 > N;
    if (tk::intet( *m_coord, *m_inpoel, point, colls[i].src, N, shift )) {
      ++numInTet;
      SolutionData data;
      data.dst = colls[i].dst;
      auto e = colls[i].src;
      const auto A = inpoel[e*4+0];
      const auto B = inpoel[e*4+1];
      const auto C = inpoel[e*4+2];
      const auto D = inpoel[e*4+3];
      data.sol.resize( u.nprop() + 1 );
      for (std::size_t c=0; c<u.nprop(); ++c) {
        data.sol[c] = N[0]*u(A,c) + N[1]*u(B,c) + N[2]*u(C,c) + N[3]*u(D,c);
      }
      data.sol.back() = N[0]*f[A] + N[1]*f[B] + N[2]*f[C] + N[3]*f[D];
      exp.push_back( data );
    }
  }

  //if (numInTet) {
  //  CkPrintf( "Source chare %i found %i/%i actual collisions\n",
  //            thisIndex, numInTet, nColls );
  //}

  // Send the solution data for the actual collisions back to the dest mesh
  proxy[index].transferSolution( exp );
  if (not m_srcnotified) {
    m_done.send();
    m_srcnotified = 1;
  }
}

void
NodeSearch::transferSolution( const std::vector< SolutionData >& data )
// *****************************************************************************
//  Transfer the interpolated solution data to destination mesh
//! \param[in] data List of solutions at nodes (multiple scalars transferred)
// *****************************************************************************
{
  tk::Fields& u = *m_u;
  std::vector< double >& f = *m_flag;
  //auto eps = std::numeric_limits< tk::real >::epsilon();

  for (std::size_t i=0; i<data.size(); ++i) {
    const auto& d = data[i];
    if (m_trflag) {     // transfer flags only
      if (m_dir) f[d.dst] = d.sol.back();       // overset to background only
    }
    else {              // transfer solution only
      if (m_dir) {      // overset to background
        if (f[d.dst] > 0.5) {
          for (std::size_t c=0; c<u.nprop(); ++c) u(d.dst,c) = d.sol[c];
        }
      }
      else {            // background to overset
        if (f[d.dst] > -0.5 and f[d.dst] < 0.5) {
          for (std::size_t c=0; c<u.nprop(); ++c) u(d.dst,c) = d.sol[c];
        }
      }
    }
  }

  // Inform the caller if we've received all solution data
  m_numreceived++;
  if (m_numreceived == m_numsent) {
    m_done.send();
  }
}

#include "NoWarning/nodesearch.def.h"
