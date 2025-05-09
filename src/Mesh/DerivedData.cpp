// *****************************************************************************
/*!
  \file      src/Mesh/DerivedData.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <cstddef>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <cfenv>

#include "Exception.hpp"
#include "DerivedData.hpp"
#include "ContainerUtil.hpp"
#include "Vector.hpp"

namespace tk {

std::size_t
npoin_in_graph( const std::vector< std::size_t >& inpoel )
// *****************************************************************************
// Compute number of points (nodes) in mesh from connectivity
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//! \return Number of mesh points (nodes)
// *****************************************************************************
{
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  return *minmax.second + 1;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< std::size_t >& inpoel, std::size_t nnpe )
// *****************************************************************************
//  Generate derived data structure, elements surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \return Linked lists storing elements surrounding points
//! \warning It is not okay to call this function with an empty container or a
//!   non-positive number of nodes per element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esup1_ and _esup2_, where _esup2_ holds the
//!   indices at which _esup1_ holds the element ids surrounding points. Looping
//!   over all elements surrounding all points can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i)
//!          use element id esup.first[i]
//!   \endcode
//!     To find out the number of points, _npoin_, the mesh connectivity,
//!     _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsup() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // allocate one of the linked lists storing elements surrounding points: esup2
  // fill with zeros
  std::vector< std::size_t > esup2( npoin+1, 0 );

  // element pass 1: count number of elements connected to each point
  for (auto n : inpoel) ++esup2[ n + 1 ];

  // storage/reshuffling pass 1: update storage counter and store
  // also find out the maximum size of esup1 (mesup)
  auto mesup = esup2[0]+1;
  for (std::size_t i=1; i<npoin+1; ++i) {
    esup2[i] += esup2[i-1];
    if (esup2[i]+1 > mesup) mesup = esup2[i]+1;
  }

  // now we know mesup, so allocate the other one of the linked lists storing
  // elements surrounding points: esup1
  std::vector< std::size_t > esup1( mesup );

  // store the elements in esup1
  std::size_t e = 0;
  for (auto n : inpoel) {
    auto j = esup2[n]+1;
    esup2[n] = j;
    esup1[j] = e/nnpe;
    ++e;
  }

  // storage/reshuffling pass 2
  for (auto i=npoin; i>0; --i) esup2[i] = esup2[i-1];
  esup2[0] = 0;

  // Return (move out) linked lists
  return std::make_pair( std::move(esup1), std::move(esup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genPsup( const std::vector< std::size_t >& inpoel,
         std::size_t nnpe,
         const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, points surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing points surrounding points
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _psup1_ and _psup2_, where _psup2_ holds the
//!   indices at which _psup1_ holds the point ids surrounding points. Looping
//!   over all points surrounding all points can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i)
//!          use point id psup.first[i]
//!   \endcode
//!    To find out the number of points, _npoin_, the mesh connectivity,
//!    _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//!   or the length-1 of the generated index list:
//!   \code{.cpp}
//!     auto npoin = psup.second.size()-1;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genPsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genPsup() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genPsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genPsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate both of the linked lists storing points surrounding points, we
  // only know the size of psup2, put in a single zero in psup1
  std::vector< std::size_t > psup2( npoin+1 ), psup1( 1, 0 );

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // fill both psup1 and psup2
  psup2[0] = 0;
  std::size_t j = 0;
  for (std::size_t p=0; p<npoin; ++p) {
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i ) {
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          ++j;
          psup1.push_back( q );
          lpoin[q] = p+1;
        }
      }
    }
    psup2[p+1] = j;
  }

  // sort point ids for each point in psup1
  for (std::size_t p=0; p<npoin; ++p)
    std::sort(
      std::next( begin(psup1), static_cast<std::ptrdiff_t>(psup2[p]+1) ),
      std::next( begin(psup1), static_cast<std::ptrdiff_t>(psup2[p+1]+1) ) );

  // Return (move out) linked lists
  return std::make_pair( std::move(psup1), std::move(psup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEdsup( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, edges surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing edges (point ids p < q) emanating from points
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _edsup1_ and _edsup2_, where _edsup2_ holds
//!   the indices at which _edsup1_ holds the edge-end point ids emanating from
//!   points for all points. The generated data structure, linked lists edsup1
//!   and edsup2, are very similar to psup1 and psup2, generated by genPsup(),
//!   except here only unique edges are stored, i.e., for edges with point ids
//!   p < q, only ids q are stored that are still associated to point p. Looping
//!   over all unique edges can then be accomplished by the following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=edsup.second[p]+1; i<=edsup.second[p+1]; ++i)
//!         use edge with point ids p < edsup.first[i]
//!   \endcode
//!   To find out the number of points, _npoin_, the mesh connectivity,
//!   _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genInpoed for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEdsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEdsup() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genEdsup() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEdsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEdsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  const auto& esup1 = esup.first;
  const auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // map to contain stars, a point associated to points connected with edges
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate edge connectivity and store as stars where center id < spike id
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) star[p].push_back(q);
          lpoin[q] = p+1;
        }
      }

  // linked lists (vectors) to store edges surrounding points and their indices
  std::vector< std::size_t > edsup1( 1, 0 ), edsup2( 1, 0 );

  // sort non-center points of each star and store nodes and indices in vectors
  for (auto& p : star) {
    std::sort( begin(p.second), end(p.second) );
    edsup2.push_back( edsup2.back() + p.second.size() );
    edsup1.insert( end(edsup1), begin(p.second), end(p.second) );
  }
  // fill up index array with the last index for points with no new edges
  for (std::size_t i=0; i<npoin-star.size(); ++i)
    edsup2.push_back( edsup2.back() );

  // Return (move out) linked lists
  return std::make_pair( std::move(edsup1), std::move(edsup2) );
}

std::vector< std::size_t >
genInpoed( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, edge connectivity
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linear vector storing edge connectivity (point ids p < q)
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linear vector and is very
//!   similar to the linked lists, _edsup1_ and _edsup2, generated by
//!   genEdsup(). The difference is that in the linear vector, inpoed, generated
//!   here, both edge point ids are stored as a pair, p < q, as opposed to the
//!   linked lists edsup1 and edsup2, in which edsup1 only stores the edge-end
//!   point ids (still associated to edge-start point ids when used together
//!   with edsup2). The rationale is that while inpoed is larger in memory, it
//!   allows direct access to edges (pair of point ids making up an edge),
//!   edsup1 and edsup2 are smaller in memory, still allow accessing the same
//!   data (edge point id pairs) but only in a linear fashion, not by direct
//!   access to particular edges. Accessing all unique edges using the edge
//!   connectivity data structure, inpoed, generated here can be accomplished by
//!   \code{.cpp}
//!     for (std::size_t e=0; e<inpoed.size()/2; ++e) {
//!       use point id p of edge e = inpoed[e*2];
//!       use point id q of edge e = inpoed[e*2+1];
//!     }
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genEdsup for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInpoed() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInpoed() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genInpoed() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genInpoed() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genInpoed() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  const auto& esup1 = esup.first;
  const auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // map to contain stars, a point associated to points connected with edges,
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate edge connectivity and store as stars where center id < spike id
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) star[p].push_back( q );
          lpoin[q] = p+1;
        }
      }

  // linear vector to store edge connectivity and their indices
  std::vector< std::size_t > inpoed;

  // sort non-center points of each star and store both start and end points of
  // each star in linear vector
  for (auto& p : star) {
    std::sort( begin(p.second), end(p.second) );
    for (auto e : p.second) {
      inpoed.push_back( p.first );
      inpoed.push_back( e );
    }
  }

  // Return (move out) linear vector
  return inpoed;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsupel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding points of elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing elements surrounding points of elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esupel1_ and _esupel2_, where _esupel2_
//!   holds the indices at which _esupel1_ holds the element ids surrounding
//!   points of elements. Looping over all elements surrounding the points of
//!   all elements can then be accomplished by the following loop:
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e)
//!       for (auto i=esupel.second[e]+1; i<=esupel.second[e+1]; ++i)
//!          use element id esupel.first[i]
//!   \endcode
//!   To find out the number of elements, _nelem_, the size of the mesh
//!   connectivity vector, _inpoel_, can be devided by the number of nodes per
//!   elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsupel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsupel() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEsupel() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsupel() with empty esup2" );

  const auto& esup1 = esup.first;
  const auto& esup2 = esup.second;

  // linked lists storing elements surrounding points of elements, put in a
  // single zero in both
  std::vector< std::size_t > esupel2( 1, 0 ), esupel1( 1, 0 );

  std::size_t e = 0;
  std::set< std::size_t > esuel;
  for (auto p : inpoel) {       // loop over all points of all elements
    // collect unique element ids of elements surrounding points of element
    for (auto i=esup2[p]+1; i<=esup2[p+1]; ++i) esuel.insert( esup1[i] );
    if (++e%nnpe == 0) {        // when finished checking all nodes of element
      // erase element whose surrounding elements are considered
      esuel.erase( e/nnpe-1 );
      // store unique element ids in esupel1
      esupel1.insert( end(esupel1), begin(esuel), end(esuel) );
      // store end-index for element used to address into esupel1
      esupel2.push_back( esupel2.back() + esuel.size() );
      esuel.clear();
    }
  }

  // Return (move out) linked lists
  return std::make_pair( std::move(esupel1), std::move(esupel2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsuel( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing elements surrounding elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esuel1_ and _esuel2_, where _esuel2_ holds
//!   the indices at which _esuel1_ holds the element ids surrounding elements.
//!   Looping over elements surrounding elements can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e)
//!       for (auto i=esuel.second[e]+1; i<=esuel.second[e+1]; ++i)
//!          use element id esuel.first[i]
//!   \endcode
//!   To find out the number of elements, _nelem_, the size of the mesh
//!   connectivity vector, _inpoel_, can be devided by the number of nodes per
//!   elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsuel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsuel() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by four" );
  Assert( !esup.first.empty(), "Attempt to call genEsuel() with empty esuel1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsuel() with empty esuel2" );

  const auto& esup1 = esup.first;
  const auto& esup2 = esup.second;

  auto nelem = inpoel.size()/nnpe;

  // lambda that returns 1 if elements hel and gel share a face
  auto adj = [ &inpoel, nnpe ]( std::size_t hel, std::size_t gel ) -> bool {
    std::size_t sp = 0;
    for (std::size_t h=0; h<nnpe; ++h)
      for (std::size_t g=0; g<nnpe; ++g)
        if (inpoel[hel*nnpe+h] == inpoel[gel*nnpe+g]) ++sp;
    if (sp == nnpe-1) return true; else return false;
  };

  // map to associate unique elements and their surrounding elements
  std::map< std::size_t, std::vector< std::size_t > > es;

  for (std::size_t e=0; e<nelem; ++e) {
    std::set< std::size_t > faces; // will collect elem ids of shared faces
    for (std::size_t n=0; n<nnpe; ++n) {
      auto i = inpoel[ e*nnpe+n ];
      for (auto j=esup2[i]+1; j<=esup2[i+1]; ++j)
        if (adj( e, esup1[j] )) faces.insert( esup1[j] );
    }
    // store element ids of shared faces
    for (auto j : faces) es[e].push_back(j);
  }

  // storing elements surrounding elements
  std::vector< std::size_t > esuel1( 1, 0 ), esuel2( 1, 0 );

  // store elements surrounding elements in linked lists
  for (const auto& e : es) {
    esuel2.push_back( esuel2.back() + e.second.size() );
    esuel1.insert( end(esuel1), begin(e.second), end(e.second) );
  }

  // Return (move out) linked lists
  return std::make_pair( std::move(esuel1), std::move(esuel2) );
}

std::vector< std::size_t >
genInedel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::vector< std::size_t >& inpoed )
// *****************************************************************************
//  Generate derived data structure, edges of elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] inpoed Edge connectivity as linear vector, see tk::genInpoed
//! \return Linear vector storing all edge ids * 2 of all elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or inpoed or a non-positive number of nodes per element; it will
//!   throw an exception.
//! \details The data generated here is stored in a linear vector with all
//!   edge ids (as defined by inpoed) of all elements. The edge ids stored in
//!   inedel can be directly used to index the vector inpoed. Because the
//!   derived data structure generated here, inedel, is intended to be used in
//!   conjunction with the linear vector inpoed and not with the linked lists
//!   edsup1 and edsup2, this function takes inpoed as an argument. Accessing
//!   the edges of element e using the edge of elements data structure, inedel,
//!   generated here can be accomplished by
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e) {
//!       for (std::size_t i=0; i<nepe; ++i) {
//!         use edge id inedel[e*nepe+i] of element e, or
//!         use point ids p < q of edge id inedel[e*nepe+i] of element e as
//!           p = inpoed[ inedel[e*nepe+i]*2 ]
//!           q = inpoed[ inedel[e*nepe+i]*2+1 ]
//!       }
//!     }
//!   \endcode
//!   where _nepe_ denotes the number of edges per elements: 3 for triangles, 6
//!   for tetrahedra. To find out the number of elements, _nelem_, the size of
//!   the mesh connectivity vector, _inpoel_, can be devided by the number of
//!   nodes per elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInedel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInedel() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genInedel() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !inpoed.empty(), "Attempt to call genInedel() with empty inpoed" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // First, generate index of star centers. This is necessary to avoid a
  // brute-force search for point ids of edges when searching for element edges.
  // Note that this is the same as edsup2, generated by genEdsup(). However,
  // because the derived data structure generated here, inedel, is intended to
  // be used in conjunction with the linear vector inpoed and not with the
  // linked lists edsup1 and edsup2, this function takes inpoed as an argument,
  // and so edsup2 is temporarily generated here to avoid a brute-force search.

  // map to contain stars, a point associated to points connected with edges
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate stars from inpoed; starting with zero, every even is a star
  // center, every odd is a spike
  for (std::size_t i=0; i<inpoed.size()/2; ++i)
    star[ inpoed[i*2] ].push_back( inpoed[i*2+1] );

  // store index of star centers in vector; assume non-center points of each
  // star have already been sorted
  std::vector< std::size_t > edsup2( 1, 0 );
  for (const auto& p : star) edsup2.push_back(edsup2.back() + p.second.size());
  // fill up index array with the last index for points with no new edges
  for (std::size_t i=0; i<npoin-star.size(); ++i)
    edsup2.push_back( edsup2.back() );
  star.clear();

  // Second, generate edges of elements

  auto nelem = inpoel.size()/nnpe;

  // map associating elem id with vector of edge ids
  std::map< std::size_t, std::vector< std::size_t > > edges;

  // generate map of elements associated to edge ids
  for (std::size_t e=0; e<nelem; ++e)
    for (std::size_t n=0; n<nnpe; ++n) {
      auto p = inpoel[e*nnpe+n];
      for (auto i=edsup2[p]+1; i<=edsup2[p+1]; ++i)
         for (std::size_t j=0; j<nnpe; ++j)
            if (inpoed[(i-1)*2+1] == inpoel[e*nnpe+j])
              edges[e].push_back( i-1 );
    }

  // linear vector to store the edge ids of all elements
  std::vector< std::size_t > inedel( sumvalsize(edges) );

  // store edge ids of elements in linear vector
  std::size_t j = 0;
  for (const auto& e : edges) for (auto p : e.second) inedel[ j++ ] = p;

  // Return (move out) vector
  return inedel;
}

std::unordered_map< UnsMesh::Edge, std::vector< std::size_t >,
                    UnsMesh::Hash<2>, UnsMesh::Eq<2> >
genEsued( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding edges
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Associative container storing elements surrounding edges (value),
//!    assigned to edge-end points (key)
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details Looping over elements surrounding all edges can be accomplished by
//!   the following loop:
//!   \code{.cpp}
//!    for (const auto& [edge,surr_elements] : esued) {
//!      use element edge-end-point ids edge[0] and edge[1]
//!      for (auto e : surr_elements) {
//!         use element id e
//!      }
//!    }
//!   \endcode
//!   esued.size() equals the number of edges.
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsued() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsued() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genEsued() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEsued() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEsued() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  const auto& esup1 = esup.first;
  const auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // lambda that returns true if element e contains edge (p < q)
  auto has = [ &inpoel, nnpe ]( std::size_t e, std::size_t p, std::size_t q ) {
    int sp = 0;
    for (std::size_t n=0; n<nnpe; ++n)
      if (inpoel[e*nnpe+n] == p || inpoel[e*nnpe+n] == q) ++sp;
    return sp == 2;
  };

  // map to associate edges to unique surrounding element ids
  std::unordered_map< UnsMesh::Edge, std::vector< std::size_t >,
                      UnsMesh::Hash<2>, UnsMesh::Eq<2> > esued;

  // generate edges and associated vector of unique surrounding element ids
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) {  // for edge given point ids p < q
            for (std::size_t j=esup2[p]+1; j<=esup2[p+1]; ++j ) {
              auto e = esup1[j];
              if (has(e,p,q)) esued[{p,q}].push_back(e);
            }
          }
          lpoin[q] = p+1;
        }
      }

  // sort element ids surrounding edges for each edge
  for (auto& p : esued) std::sort( begin(p.second), end(p.second) );

  // Return elements surrounding edges data structure
  return esued;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEdpas( int mvecl, std::size_t nnpe, std::size_t npoin,
          const std::vector< std::size_t >& inpoed )
// *****************************************************************************
//  Generate vector-groups for edges
//! \param[in] mvecl Max vector length to target
//! \param[in] nnpe Number of nodes per (super-)edge
//! \param[in] npoin Number mesh points
//! \param[in] inpoed Edge connectivity as linear vector, see tk::genInpoed for
//!            nnpe=2
//! \return Linked lists storing edge-groups so that any point of a group is
//!   accessed only once within a group.
//! \warning It is not okay to call this function with an empty container or a
//!   non-positive number of nodes per element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _edpas1_ and _edpas2_, where _edpas2_ holds
//!   the indices at which _edpas1_ holds the edge ids of a vector group.
//!   Looping over all groups can then be accomplished by the following loop:
//!   \code{.cpp}
//!     for (std::size_t w=0; w<edpas.second.size()-1; ++w)
//!       for (auto i=edpas.second[w]+1; i<=edpas.second[w+1]; ++i)
//!          use edge id edpas.first[i]
//!   \endcode
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  if (inpoed.empty()) return {};

  Assert( mvecl > 0, "Attempt to call genEdpas() with non-positive veclen" );
  Assert( nnpe > 0, "Attempt to call genEdpas() with non-positive nnpe" );
  Assert( npoin > 0, "Attempt to call genEdpas() with non-positive npoin" );
  Assert( inpoed.size()%nnpe == 0, "Size of inpoed must be divisible by nnpe" );

  auto nedge = inpoed.size() / nnpe;

  std::vector< std::size_t > ledge( nedge, 0 );
  std::vector< std::size_t > lpoin( npoin, 0 );

  std::pair< std::vector< std::size_t >, std::vector< std::size_t > > edpas;
  edpas.first.resize( nedge+1, 0 );
  edpas.second.push_back( 0 );

  std::unordered_set< std::size_t > unedge( nedge );
  for (std::size_t e=0; e<nedge; ++e) unedge.insert( e );

  std::size_t nenew = 0, ngrou = 0;

  while (nenew < nedge) {
    int nvecl = 0;
    ++ngrou;
    edpas.second.emplace_back();
    for (auto ie = begin(unedge); ie != end(unedge); ) {
      auto e = *ie;
      const auto N = inpoed.data() + e*nnpe;
      std::size_t nsw = 0;
      for (std::size_t i=0; i<nnpe; ++i) {
        if (lpoin[N[i]] == ngrou) break; else ++nsw;
      }
      if (nsw == nnpe) {
        for (std::size_t i=0; i<nnpe; ++i) lpoin[N[i]] = ngrou;
        ledge[e] = ngrou;
        ++nenew;
        ++nvecl;
        edpas.first[nenew] = e;
        edpas.second[ngrou] = nenew;
        ie = unedge.erase( ie );
      } else {
        ++ie;
      }
      if (nvecl == mvecl) break;
    }
  }

  //std::size_t ne = 0;
  //for (std::size_t i=0; i<edpas.second.size()-1; ++i) {
  //  std::cout << i+1 << ": " << edpas.second[i+1] - edpas.second[i] << '\n';
  //  ne += edpas.second[i+1] - edpas.second[i];
  //}
  //std::cout << "edges grouped: " << ne << " of " << nedge << '\n';
  //std::cout << "edge groups:";
  //for (std::size_t g=1; g<=edpas.second.size(); ++g) {
  //  std::cout << '\n';
  //  for (std::size_t e=0; e<nedge; ++e) {
  //    if (ledge[e] == g) {
  //      const auto N = inpoed.data() + e*nnpe;
  //      std::cout << e << ":\t";
  //      for (std::size_t i=0; i<nnpe-1; ++i) std::cout << N[i] << '-';
  //      std::cout << N[nnpe-1] << "\tgrp: " << ledge[e] << '\n';
  //    }
  //  }
  //}
  //std::cout << '\n';

  //std::cout << "\nnew access loop:\n";
  //for (std::size_t g=0; g<edpas.second.size()-1; ++g) {
  //  //#pragma omp simd
  //  for (auto w=edpas.second[g]+1; w<=edpas.second[g+1]; ++w) {
  //    auto e = edpas.first[w];
  //    const auto N = inpoed.data() + e*nnpe;
  //    for (std::size_t i=0; i<nnpe-1; ++i) std::cout << N[i] << '-';
  //    std::cout << N[nnpe-1] << ' ';
  //  }
  //  std::cout << '\n';
  //}
  //std::cout << '\n';

  //std::cout << "old access loop:\n";
  //for (std::size_t e=0; e<nedge; ++e) {
  //  const auto N = inpoed.data() + e*nnpe;
  //  for (std::size_t i=0; i<nnpe-1; ++i) std::cout << N[i] << '-';
  //  std::cout << N[nnpe-1] << ' ';
  //}
  //std::cout << "\n\n";

  // Return linked lists
  return edpas;
}

std::size_t
genNbfacTet( std::size_t tnbfac,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& triinpoel_complete,
             const std::map< int, std::vector< std::size_t > >& bface_complete,
             const std::unordered_map< std::size_t, std::size_t >& lid,
             std::vector< std::size_t >& triinpoel,
             std::map< int, std::vector< std::size_t > >& bface )
// *****************************************************************************
//  Generate number of boundary-faces and triangle boundary-face connectivity
//! \param[in] tnbfac Total number of boundary faces in the entire mesh.
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh.
//! \param[in] triinpoel_complete Interconnectivity of points and boundary-face
//!   in the entire mesh.
//! \param[in] bface_complete Map of boundary-face lists mapped to corresponding 
//!   side set ids for the entire mesh.
//! \param[in] lid Mapping between the node indices used in the smaller inpoel
//!   connectivity (a subset of the entire triinpoel_complete connectivity),
//!   e.g., after mesh partitioning.
//! \param[inout] triinpoel Interconnectivity of points and boundary-face in
//!   this mesh-partition.
//! \param[inout] bface Map of boundary-face lists mapped to corresponding 
//!   side set ids for this mesh-partition
//! \return Number of boundary-faces on this chare/mesh-partition.
//! \details This function takes a mesh by its domain-element
//!   (tetrahedron-connectivity) in inpoel and a boundary-face (triangle)
//!   connectivity in triinpoel_complete. Based on these two arrays, it
//!   searches for those faces of triinpoel_complete that are also in inpoel
//!   and as a result it generates (1) the number of boundary faces shared with
//!   the mesh in inpoel and (2) the intersection of the triangle element
//!   connectivity whose faces are shared with inpoel. An example use case is
//!   where triinpoel_complete contains the connectivity for the boundary of the
//!   full problem/mesh and inpoel contains the connectivity for only a chunk of
//!   an already partitioned mesh. This function then intersects
//!   triinpoel_complete with inpoel and returns only those faces that share
//!   nodes with inpoel.
//! \warning This is for Triangular face-elements only.
// *****************************************************************************
{
  // cppcheck-suppress unreadVariable
  std::size_t nbfac = 0, nnpf = 3;

  if (tnbfac > 0)
  {

  Assert( !inpoel.empty(),
          "Attempt to call genNbfacTet() on empty inpoel container" );
  Assert( !triinpoel_complete.empty(),
          "Attempt to call genNbfacTet() on empty triinpoel_complete container" );
  Assert( triinpoel_complete.size()/nnpf == tnbfac, 
          "Incorrect size of triinpoel in genNbfacTet()" );

  auto nptet = inpoel;
  auto nptri = triinpoel_complete;

  unique( nptet );
  unique( nptri );

  std::unordered_set< std::size_t > snptet;

  // getting the reduced inpoel as a set for quick searches
  snptet.insert( begin(nptet), end(nptet));

  // vector to store boundary-face-nodes in this chunk
  std::vector< std::size_t > nptri_chunk;

  // getting the nodes of the boundary-faces in this chunk
  for (auto i : nptri)
    if (snptet.find(i) != end(snptet))
      nptri_chunk.push_back(i);

  std::size_t tag, icoun;

  // matching nodes in nptri_chunk with nodes in inpoel and 
  // triinpoel_complete to get the number of faces in this chunk
  for (const auto& ss : bface_complete)
  {
    for (auto f : ss.second)
    {
      icoun = f*nnpf;
      tag = 0;
      for (std::size_t i=0; i<nnpf; ++i) {
        for (auto j : nptri_chunk) {
          // cppcheck-suppress useStlAlgorithm
          if (triinpoel_complete[icoun+i] == j) ++tag;
        }
      }
      if (tag == nnpf)
      // this is a boundary face
      {
        for (std::size_t i=0; i<nnpf; ++i)
        {
          auto ip = triinpoel_complete[icoun+i];

          // find local renumbered node-id to store in triinpoel
          triinpoel.push_back( cref_find(lid,ip) );
        }

        bface[ss.first].push_back(nbfac);
        ++nbfac;
      }
    }
  }

  }

  return nbfac;
}

std::vector< int >
genEsuelTet( const std::vector< std::size_t >& inpoel,
             const std::pair< std::vector< std::size_t >,
                              std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding elements
//  as a fixed length data structure as a full vector, including
//  boundary elements as -1.
//  \warning This is for Tetrahedra only.
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Vector storing elements surrounding elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second; it will throw an exception.
//! \details The data generated here is stored in a single vector, with length
//!   nfpe * nelem. Note however, that nelem is not explicitly provided, but
//!   calculated from inpoel. For boundary elements, at the boundary face, this
//!   esuelTet stores value -1 indicating that this is outside the domain. The
//!   convention for numbering the local face (triangle) connectivity is very
//!   important, e.g., in generating the inpofa array later. This node ordering
//!   convention is stored in tk::lpofa. Thus function is specific to
//!   tetrahedra, which is reflected in the fact that nnpe and nfpe are being
//!   set here in the function rather than being input arguments. To find out
//!   the number of elements, _nelem_, the size of the mesh connectivity vector,
//!   _inpoel_, can be devided by the number of nodes per elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsuelTet() on empty container" );
  Assert( !esup.first.empty(), "Attempt to call genEsuelTet() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsuelTet() with empty esup2" );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // set tetrahedron geometry
  // cppcheck-suppress unreadVariable
  std::size_t nnpe = 4, nfpe = 4, nnpf = 3;

  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by four" );

  // get nelem and npoin
  // cppcheck-suppress unreadVariable
  auto nelem = inpoel.size()/nnpe;
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  std::vector< int > esuelTet(nfpe*nelem, -1);
  std::vector< std::size_t > lhelp(nnpf,0),
                             lpoin(npoin,0);

  for (std::size_t e=0; e<nelem; ++e)
  {
    auto mark = nnpe*e;
    for (std::size_t fe=0; fe<nfpe; ++fe)
    {
      // array which stores points on this face
      lhelp[0] = inpoel[mark+lpofa[fe][0]];
      lhelp[1] = inpoel[mark+lpofa[fe][1]];
      lhelp[2] = inpoel[mark+lpofa[fe][2]];

      // mark in this array
      lpoin[lhelp[0]] = 1;
      lpoin[lhelp[1]] = 1;
      lpoin[lhelp[2]] = 1;

      // select a point on this face
      auto ipoin = lhelp[0];

      // loop over elements around this point
      for (std::size_t j=esup2[ipoin]+1; j<=esup2[ipoin+1]; ++j )
      {
        auto jelem = esup1[j];
        // if this jelem is not e itself then proceed
        if (jelem != e)
        {
          for (std::size_t fj=0; fj<nfpe; ++fj)
          {
            std::size_t icoun(0);
            for (std::size_t jnofa=0; jnofa<nnpf; ++jnofa)
            {
              auto markj = jelem*nnpe;
              auto jpoin = inpoel[markj+lpofa[fj][jnofa]];
              if (lpoin[jpoin] == 1) { ++icoun; }
            }
            //store esuel if
            if (icoun == nnpf)
            {
              auto markf = nfpe*e;
              esuelTet[markf+fe] = static_cast<int>(jelem);

              markf = nfpe*jelem;
              esuelTet[markf+fj] = static_cast<int>(e);
            }
          }
        }
      }
      // reset this array
      lpoin[lhelp[0]] = 0;
      lpoin[lhelp[1]] = 0;
      lpoin[lhelp[2]] = 0;
    }
  }

  return esuelTet;
}

std::size_t
genNipfac( std::size_t nfpe,
           std::size_t nbfac,
           const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate number of internal and physical-boundary faces
//! \param[in] nfpe Number of faces per element.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Total number of faces in the mesh
//! \details The unsigned integer here gives the number of internal and
//!    physical-boundary faces in the mesh. The data structure does not include
//!    faces that are on partition/chare-boundaries.
// *****************************************************************************
{
  Assert( !esuelTet.empty(), "Attempt to call genNipfac() with empty esuelTet" );
  Assert( esuelTet.size()%nfpe == 0,
                  "Size of esuelTet must be divisible by nfpe" );
  Assert( nfpe > 0, "Attempt to call genNipfac() with zero faces per element" );

  auto nelem = esuelTet.size()/nfpe;

  std::size_t nifac = 0;

  // loop through elements surrounding elements to find number of internal faces
  for (std::size_t e=0; e<nelem; ++e)
  {
    for (std::size_t ip=nfpe*e; ip<nfpe*(e+1); ++ip)
    {
      if (esuelTet[ip] != -1)
      {
        if ( e<static_cast< std::size_t >(esuelTet[ip]) )
        {
          ++nifac;
        }
      }
    }
  }

  return nifac + nbfac;
}

std::vector< int >
genEsuf( std::size_t nfpe,
         std::size_t nipfac,
         std::size_t nbfac,
         const std::vector< std::size_t >& belem,
         const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate derived data structure, elements surrounding faces
//! \param[in] nfpe  Number of faces per element.
//! \param[in] nipfac Number of internal and physical-boundary faces.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] belem Boundary element vector.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Elements surrounding faces.
//! \details The unsigned integer vector gives the IDs of the elements to the
//    left and the right of each face in the mesh. The convention followed 
//    throughout is : The left element always has an ID smaller than the ID of
//    the right element.
// *****************************************************************************
{
  Assert( esuelTet.size()%nfpe == 0, 
                  "Size of esuelTet must be divisible by nfpe" );
  Assert( nfpe > 0, "Attempt to call genEsuf() with zero faces per element" );

  auto nelem = esuelTet.size()/nfpe;

  std::vector< int > esuf(2*nipfac);

  // counters for number of internal and boundary faces
  std::size_t icoun(2*nbfac), bcoun(0);

  // loop to get face-element connectivity for internal faces
  for (std::size_t e=0; e<nelem; ++e) {
    for (std::size_t ip=nfpe*e; ip<nfpe*(e+1); ++ip) {
      auto jelem = esuelTet[ip];
      if (jelem != -1)
      {
        if ( e < static_cast< std::size_t >(jelem) )
        {
          esuf[icoun] = static_cast< int >(e);
          esuf[icoun+1] = static_cast< int >(jelem);
          icoun = icoun + 2;
        }
      }
    }
  }

  // loop to get face-element connectivity for physical-boundary faces
  bcoun = 0;
  for (auto ie : belem) {
    esuf[bcoun] = static_cast< int >(ie);
    esuf[bcoun+1] = -1;  // outside domain
    bcoun = bcoun + 2;
  }

  return esuf;
}

std::vector< std::size_t >
genInpofaTet( std::size_t nipfac,
              std::size_t nbfac,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::size_t >& triinpoel,
              const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate derived data structure, points on faces for tetrahedra only
//! \param[in] nipfac Number of internal and physical-boundary faces.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] inpoel Element-node connectivity.
//! \param[in] triinpoel Face-node connectivity.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Points surrounding faces. The unsigned integer vector gives the
//!   elements to the left and to the right of each face in the mesh.
// *****************************************************************************
{
  std::vector< std::size_t > inpofa;

  // set tetrahedron geometry
  // cppcheck-suppress unreadVariable
  std::size_t nnpe(4), nfpe(4), nnpf(3);

  Assert( esuelTet.size()%nfpe == 0,
                  "Size of esuelTet must be divisible by nfpe" );
  Assert( inpoel.size()%nnpe == 0,
                  "Size of inpoel must be divisible by nnpe" );

  inpofa.resize(nnpf*nipfac);

  // counters for number of internal and boundary faces
  std::size_t icoun(nnpf*nbfac);

  // loop over elems to get nodes on faces
  // this fills the interior face-node connectivity part
  for (std::size_t e=0; e<inpoel.size()/nnpe; ++e)
  {
    auto mark = nnpe*e;
    for (std::size_t f=0; f<nfpe ; ++f)
    {
      auto ip = nfpe*e + f;
      auto jelem = esuelTet[ip];
      if (jelem != -1)
      {
        if ( e < static_cast< std::size_t >(jelem) )
        {
          inpofa[icoun]   = inpoel[mark+lpofa[f][0]];
          inpofa[icoun+1] = inpoel[mark+lpofa[f][1]];
          inpofa[icoun+2] = inpoel[mark+lpofa[f][2]];
          icoun = icoun + nnpf;
        }
      }
    }
  }

  // this fills the boundary face-node connectivity part
  // consistent with triinpoel
  for (std::size_t f=0; f<nbfac; ++f)
  {
    icoun = nnpf * f;
    inpofa[icoun+0] = triinpoel[icoun+0];
    inpofa[icoun+1] = triinpoel[icoun+1];
    inpofa[icoun+2] = triinpoel[icoun+2];
  }

  return inpofa;
}
        
std::vector< std::size_t >
genBelemTet( std::size_t nbfac,
              const std::vector< std::size_t >& inpofa,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data, boundary elements
//! \param[in] nbfac Number of boundary faces.
//! \param[in] inpofa Face-node connectivity.
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Host elements or boundary elements. The unsigned integer vector
//!   gives the elements to the left of each boundary face in the mesh.
//! \details The data structure generated here contains an array of elements
//!   which share one or more of their faces with the physical boundary, i.e.,
//!   where exodus specifies a side-set for faces. Such elements are sometimes
//!   also called host or boundary elements.
// *****************************************************************************
{
  std::vector< std::size_t > belem(nbfac);

  if (nbfac > 0)
  {

  // set tetrahedron geometry
  // cppcheck-suppress unreadVariable
  std::size_t nnpf = 3, tag = 0;

  // loop over all the boundary faces
  for(std::size_t f=0; f<nbfac; ++f)
  {
    belem[f] = 0;

    // array storing the element-cluster around face
    std::vector< std::size_t > elemcluster;

    // loop over the nodes of this boundary face
    for(std::size_t lp=0; lp<nnpf; ++lp)
    {
      auto gp = inpofa[nnpf*f + lp];

      Assert( gp < esup.second.size(), "Indexing out of esup2" );
      // loop over elements surrounding this node
      for (auto i=esup.second[gp]+1; i<=esup.second[gp+1]; ++i)
      {
        // form element-cluster vector
        elemcluster.push_back(esup.first[i]);
      }
    }

    // loop over element cluster to find repeating elements
    for(std::size_t i=0; i<elemcluster.size(); ++i)
    {
      auto ge = elemcluster[i];
      tag = 1;
      for(std::size_t j=0; j<elemcluster.size(); ++j)
      {
        if ( i != j && elemcluster[j] == ge )
        {
          tag++;
        }
      }
      if (tag == nnpf)
      {
        // this is the required boundary element
        belem[f] = ge;
        break;
      }
    }
  }
  }

  return belem;
}
        
bool
leakyPartition( const std::vector< int >& esueltet,
                const std::vector< std::size_t >& inpoel,
                const std::array< std::vector< real >, 3 >& coord )
// *****************************************************************************
// Perform leak-test on mesh (partition)
//! \param[in] esueltet Elements surrounding elements for tetrahedra, see
//!   tk::genEsueltet()
//! \param[in] inpoel Element connectivity
//! \param[in] coord Node coordinates
//! \details This function computes a surface integral over the boundary of the
//!   incoming mesh (partition). A non-zero vector result indicates a leak, e.g.,
//!   a hole in the mesh (partition), which indicates an error either in the
//    mesh geometry, mesh partitioning, or in the data structures that represent
//    faces.
//! \return True if partition leaks.
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Storage for surface integral over our mesh partition
  std::array< real, 3 > s{{ 0.0, 0.0, 0.0}};

  for (std::size_t e=0; e<esueltet.size()/4; ++e) {   // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)     // for all tet faces
      if (esueltet[mark+f] == -1) {     // if face has no outside-neighbor tet
        // 3 local node IDs of face
        auto A = inpoel[ mark + lpofa[f][0] ];
        auto B = inpoel[ mark + lpofa[f][1] ];
        auto C = inpoel[ mark + lpofa[f][2] ];
        // Compute face area and normal
        real nx, ny, nz;
        auto a = normal( x[A],x[B],x[C], y[A],y[B],y[C], z[A],z[B],z[C],
                         nx, ny, nz );
        // Sum up face area * face unit-normal
        s[0] += a * nx;
        s[1] += a * ny;
        s[2] += a * nz;
      }
  }

  auto eps = 1.0e-9;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

bool
conforming( const std::vector< std::size_t >& inpoel,
            const std::array< std::vector< real >, 3 >& coord,
            bool cerr,
            const std::vector< std::size_t >& rid )
// *****************************************************************************
// Check if mesh (partition) is conforming
//! \param[in] inpoel Element connectivity
//! \param[in] coord Node coordinates
//! \param[in] cerr True if hanging-node edge data should be output to
//!   std::cerr (true by default)
//! \param[in] rid AMR Lib node id map
//!   std::cerr (true by default)
//! \return True if mesh (partition) has no hanging nodes and thus the mesh is
//!   conforming, false if non-conforming.
//! \details A conforming mesh by definition has no hanging nodes. A node is
//!   hanging if an edge of one element coincides with two (or more) edges (of
//!   two or more other elements). Thus, testing for conformity relies on
//!   checking the coordinates of all vertices: if any vertex coincides with
//!   that of a mid-point node of an edge, that is a hanging node. Note that
//!   this assumes that hanging nodes can only be at the mid-point of edges.
//!   This may happen after a mesh refinement step, due to a problem/bug,
//!   within the mesh refinement algorithm given by J. Waltz, Parallel adaptive
//!   refinement for unsteady flow calculations on 3D unstructured grids,
//!   International Journal for Numerical Methods in Fluids, 46: 37–57, 2004,
//!   which always adds/removes vertices at the mid-points of edges of a
//!   tetrahedron mesh within a single refinement step. Thus this algorithm is
//!   intended for this specific case, i.e., test for conformity after a
//!   single refinement step and not after multiple ones or for detecting
//!   hanging nodes in an arbitrary mesh.
//*****************************************************************************
{
  Assert( !inpoel.empty(),
          "Attempt to call conforming() with empty mesh connectivity" );
  Assert( inpoel.size() % 4 == 0,
          "Size of inpoel must be divisible by nnpe" );
  Assert( *std::min_element( begin(inpoel), end(inpoel) ) == 0,
          "Node ids should start from zero" );
  Assert( !coord[0].empty() && !coord[1].empty() && !coord[2].empty(),
          "Attempt to call conforming() with empty coordinates container" );

  using Coord = UnsMesh::Coord;
  using Edge = UnsMesh::Edge;
  using Tet = UnsMesh::Tet;

  // Compare operator to be used as less-than for std::array< real, 3 >,
  // implemented as a lexicographic ordering.
  struct CoordLess {
    const real eps = std::numeric_limits< real >::epsilon();
    bool operator() ( const Coord& lhs, const Coord& rhs ) const {
      if (lhs[0] < rhs[0])
        return true;
      else if (std::abs(lhs[0]-rhs[0]) < eps && lhs[1] < rhs[1])
        return true;
      else if (std::abs(lhs[0]-rhs[0]) < eps &&
               std::abs(lhs[1]-rhs[1]) < eps &&
               lhs[2] < rhs[2])
        return true;
      else
        return false;
    }
  };

  // Map associating data on potential hanging nodes. Key: coordinates of nodes
  // of edge-half points, value: tet id (local if in parallel), tet connectivity
  // (using local ids if in parallel), edge ids (local if in parallel).
  std::map< Coord,                      // edge-half node coordinates: x, y, z
            std::tuple< std::size_t,    // element id of edge-half node
                        Tet,            // element node ids of edge-half node
                        Edge >,         // edge containing half-node
            CoordLess > edgeNodes;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  fenv_t fe;
  feholdexcept( &fe );

  // Compute coordinates of nodes of mid-points of all edges
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    std::array<Edge,6> edge{{ {{A,B}}, {{B,C}}, {{A,C}},
                              {{A,D}}, {{B,D}}, {{C,D}} }};
    for (const auto& n : edge) {
      Coord en{{ (x[n[0]] + x[n[1]]) / 2.0,
                 (y[n[0]] + y[n[1]]) / 2.0,
                 (z[n[0]] + z[n[1]]) / 2.0 }};
      edgeNodes[ en ] = std::tuple<std::size_t,Tet,Edge>{ e, {{A,B,C,D}}, n };
    }
  }

  feclearexcept( FE_UNDERFLOW );
  feupdateenv( &fe );

  // Find hanging nodes. If the coordinates of an element vertex coincide with
  // that of a mid-point node of an edge, that is a hanging node. If we find one
  // such node we print out some info on it.
  auto ix = x.cbegin();
  auto iy = y.cbegin();
  auto iz = z.cbegin();

  bool hanging_node = false;

  while (ix != x.cend()) {
    Coord n{{ *ix, *iy, *iz }};
    auto i = edgeNodes.find( n );
    if (i != end(edgeNodes)) {
      const auto& hanging_node_coord = i->first;
      const auto& hanging_node_info = i->second;
      auto tet_id = std::get< 0 >( hanging_node_info );
      const auto& tet = std::get< 1 >( hanging_node_info );
      const auto& edge = std::get< 2 >( hanging_node_info );
      if (cerr) {
        std::cerr
          << "Mesh conformity test found hanging node with coordinates"" ("
          << hanging_node_coord[0] << ", "
          << hanging_node_coord[1] << ", "
          << hanging_node_coord[2] << ") of tetrahedron element "
          << tet_id << " with connectivity (" << tet[0] << ','
          << tet[1] << ',' << tet[2] << ',' << tet[3] << ") on edge ("
          << edge[0] << ',' << edge[1] << ")"
          << "AMR lib node ids for this edge: " << rid[edge[0]] << ','
          << rid[edge[1]] << std::endl;
      }
      hanging_node = true;
    }
    ++ix; ++iy; ++iz;
  }

  if (hanging_node) return false;

  return true;
}

bool
intet( const std::array< std::vector< real >, 3 >& coord,
       const std::vector< std::size_t >& inpoel,
       const std::vector< real >& p,
       std::size_t e,
       std::array< real, 4 >& N,
       const std::array< real, 3 >& eps )
// *****************************************************************************
//  Determine if a point is in a tetrahedron
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] p Point coordinates
//! \param[in] e Mesh cell index
//! \param[in,out] N Shapefunctions evaluated at the point
//! \param[in] eps Optionally shift point coordinates in x,y,z
//! \return True if ppoint is in mesh cell
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( p.size() == 3, "Size mismatch" );

  // Tetrahedron node indices
  const auto A = inpoel[e*4+0];
  const auto B = inpoel[e*4+1];
  const auto C = inpoel[e*4+2];
  const auto D = inpoel[e*4+3];

  // Tetrahedron node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Point coordinates
  const auto& xp = p[0] + eps[0];
  const auto& yp = p[1] + eps[1];
  const auto& zp = p[2] + eps[2];

  // Evaluate linear shapefunctions at point locations using Cramer's Rule
  //    | xp |   | x1 x2 x3 x4 |   | N1 |
  //    | yp | = | y1 y2 y3 y4 | • | N2 |
  //    | zp |   | z1 z2 z3 z4 |   | N3 |
  //    | 1  |   | 1  1  1  1  |   | N4 |

  real DetX = (y[B]*z[C] - y[C]*z[B] - y[B]*z[D] + y[D]*z[B] +
    y[C]*z[D] - y[D]*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[C]*y[A]*z[D] + x[C]*y[D]*z[A] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] - x[B]*y[C]*z[D] + x[B]*y[D]*z[C] +
    x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - x[D]*y[B]*z[C] + x[D]*y[C]*z[B];

  real DetX1 = (y[D]*z[C] - y[C]*z[D] + y[C]*zp - yp*z[C] -
    y[D]*zp + yp*z[D])*x[B] + x[C]*y[B]*z[D] - x[C]*y[D]*z[B] -
    x[D]*y[B]*z[C] + x[D]*y[C]*z[B] - x[C]*y[B]*zp + x[C]*yp*z[B] +
    xp*y[B]*z[C] - xp*y[C]*z[B] + x[D]*y[B]*zp - x[D]*yp*z[B] -
    xp*y[B]*z[D] + xp*y[D]*z[B] + x[C]*y[D]*zp - x[C]*yp*z[D] -
    x[D]*y[C]*zp + x[D]*yp*z[C] + xp*y[C]*z[D] - xp*y[D]*z[C];

  real DetX2 = (y[C]*z[D] - y[D]*z[C] - y[C]*zp + yp*z[C] +
    y[D]*zp - yp*z[D])*x[A] + x[C]*y[D]*z[A] - x[C]*y[A]*z[D] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] + x[C]*y[A]*zp - x[C]*yp*z[A] -
    xp*y[A]*z[C] + xp*y[C]*z[A] - x[D]*y[A]*zp + x[D]*yp*z[A] +
    xp*y[A]*z[D] - xp*y[D]*z[A] - x[C]*y[D]*zp + x[C]*yp*z[D] +
    x[D]*y[C]*zp - x[D]*yp*z[C] - xp*y[C]*z[D] + xp*y[D]*z[C];

  real DetX3 = (y[D]*z[B] - y[B]*z[D] + y[B]*zp - yp*z[B] -
    y[D]*zp + yp*z[D])*x[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[B]*y[A]*zp + x[B]*yp*z[A] +
    xp*y[A]*z[B] - xp*y[B]*z[A] + x[D]*y[A]*zp - x[D]*yp*z[A] -
    xp*y[A]*z[D] + xp*y[D]*z[A] + x[B]*y[D]*zp - x[B]*yp*z[D] -
    x[D]*y[B]*zp + x[D]*yp*z[B] + xp*y[B]*z[D] - xp*y[D]*z[B];

  real DetX4 = (y[B]*z[C] - y[C]*z[B] - y[B]*zp + yp*z[B] +
    y[C]*zp - yp*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*zp - x[B]*yp*z[A] -
    xp*y[A]*z[B] + xp*y[B]*z[A] - x[C]*y[A]*zp + x[C]*yp*z[A] +
    xp*y[A]*z[C] - xp*y[C]*z[A] - x[B]*y[C]*zp + x[B]*yp*z[C] +
    x[C]*y[B]*zp - x[C]*yp*z[B] - xp*y[B]*z[C] + xp*y[C]*z[B];

  // Shape functions evaluated at point
  N[0] = DetX1/DetX;
  N[1] = DetX2/DetX;
  N[2] = DetX3/DetX;
  N[3] = DetX4/DetX;

  // if min( N^i, 1-N^i ) > 0 for all i, point is in cell
  if ( std::min(N[0],1.0-N[0]) > 0.0 && std::min(N[1],1.0-N[1]) > 0.0 &&
       std::min(N[2],1.0-N[2]) > 0.0 && std::min(N[3],1.0-N[3]) > 0.0 )
  {
    return true;
  } else {
    return false;
  }
}

} // tk::
