// *****************************************************************************
/*!
  \file      tests/unit/Mesh/TestReorder.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Mesh/Reorder
  \details   Unit tests for Mesh/Reorder. All unit tests start from simple mesh
     connectivities defined in the code. The tetrahedron mesh in Gmsh ASCII
     format is as follows. Note that ids start from zero in the code, but from
     one in Gmsh.
     \code{.sh}
       $MeshFormat
       2.2 0 8
       $EndMeshFormat
       $Nodes
       14
       1 0 0 0
       2 1 0 0
       3 1 1 0
       4 0 1 0
       5 0 0 1
       6 1 0 1
       7 1 1 1
       8 0 1 1
       9 0.5 0.5 0
       10 0.5 0.5 1
       11 0.5 0 0.5
       12 1 0.5 0.5
       13 0.5 1 0.5
       14 0 0.5 0.5
       $EndNodes
       $Elements
       24
       1 4 1 0 12 14 9 11
       2 4 1 0 10 14 13 12
       3 4 1 0 14 13 12 9
       4 4 1 0 10 14 12 11
       5 4 1 0 1 14 5 11
       6 4 1 0 7 6 10 12
       7 4 1 0 14 8 5 10
       8 4 1 0 8 7 10 13
       9 4 1 0 7 13 3 12
       10 4 1 0 1 4 14 9
       11 4 1 0 13 4 3 9
       12 4 1 0 3 2 12 9
       13 4 1 0 4 8 14 13
       14 4 1 0 6 5 10 11
       15 4 1 0 1 2 9 11
       16 4 1 0 2 6 12 11
       17 4 1 0 6 10 12 11
       18 4 1 0 2 12 9 11
       19 4 1 0 5 14 10 11
       20 4 1 0 14 8 10 13
       21 4 1 0 13 3 12 9
       22 4 1 0 7 10 13 12
       23 4 1 0 14 4 13 9
       24 4 1 0 14 1 9 11
       $EndElements
     \endcode
     Here is the simple triangle mesh used below by the unit tests in Gmsh ASCII
     format. Note that ids start from zero in the code, but from one in Gmsh.
     \code{.sh}
       $MeshFormat
       2.2 0 8
       $EndMeshFormat
       $Nodes
       14
       1 0 0 0
       2 1 0 0
       3 1 1 0
       4 0 1 0
       5 0 0 1
       6 1 0 1
       7 1 1 1
       8 0 1 1
       9 0.5 0.5 0
       10 0.5 0.5 1
       11 0.5 0 0.5
       12 1 0.5 0.5
       13 0.5 1 0.5
       14 0 0.5 0.5
       $EndNodes
       $Elements
       24
       1 2 2 0 1 1 9 2
       2 2 2 0 1 1 4 9
       3 2 2 0 1 2 9 3
       4 2 2 0 1 3 9 4
       5 2 2 0 2 5 6 10
       6 2 2 0 2 5 10 8
       7 2 2 0 2 6 7 10
       8 2 2 0 2 7 8 10
       9 2 2 0 3 1 2 11
       10 2 2 0 3 1 11 5
       11 2 2 0 3 2 6 11
       12 2 2 0 3 5 11 6
       13 2 2 0 4 2 3 12
       14 2 2 0 4 2 12 6
       15 2 2 0 4 3 7 12
       16 2 2 0 4 6 12 7
       17 2 2 0 5 3 4 13
       18 2 2 0 5 3 13 7
       19 2 2 0 5 4 8 13
       20 2 2 0 5 7 13 8
       21 2 2 0 6 1 14 4
       22 2 2 0 6 1 5 14
       23 2 2 0 6 4 14 8
       24 2 2 0 6 5 8 14
       $EndElements
     \endcode
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct Reorder_common {

  // Mesh node coordinates
  std::array< std::vector< tk::real >, 3 > tetcoord {{
    {{ 0, 1, 1, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1, 0.5, 0 }},
    {{ 0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0, 0.5, 1, 0.5 }},
    {{ 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0.5, 0.5, 0.5, 0.5 }} }};

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > tetinpoel { 12, 14,  9, 11,
                                         10, 14, 13, 12,
                                         14, 13, 12,  9,
                                         10, 14, 12, 11,
                                         1,  14,  5, 11,
                                         7,   6, 10, 12,
                                         14,  8,  5, 10,
                                         8,   7, 10, 13,
                                         7,  13,  3, 12,
                                         1,   4, 14,  9,
                                         13,  4,  3,  9,
                                         3,   2, 12,  9,
                                         4,   8, 14, 13,
                                         6,   5, 10, 11,
                                         1,   2,  9, 11,
                                         2,   6, 12, 11,
                                         6,  10, 12, 11,
                                         2,  12,  9, 11,
                                         5,  14, 10, 11,
                                         14,  8, 10, 13,
                                         13,  3, 12,  9,
                                         7,  10, 13, 12,
                                         14,  4, 13,  9,
                                         14,  1,  9, 11 };
};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using Reorder_group = test_group< Reorder_common, MAX_TESTS_IN_GROUP >;
using Reorder_object = Reorder_group::object;

//! Define test group
static Reorder_group Reorder( "Mesh/Reorder" );

//! Test definitions for group

//! Attempt to shift empty container using shiftToZero
template<> template<>
void Reorder_object::test< 1 >() {
  set_test_name( "shiftToZero graceful with empty inpoel" );

  // Attempt to shift node IDs with empty connectivity. If some error happens or an
  // exception is throw that will go to the screen; no further tests are
  // necessary.
  std::vector< std::size_t > empty;
  tk::shiftToZero( empty );
}

//! Shift node ids to zero in line mesh
template<> template<>
void Reorder_object::test< 2 >() {
  set_test_name( "shiftToZero for lines" );

  // Mesh connectivity for simple line-only mesh
  std::vector< std::size_t > inpoel { 1, 2,
                                      2, 3,
                                      3, 4,
                                      4, 1,
                                      5, 6,
                                      6, 7,
                                      7, 8,
                                      8, 5,
                                      1, 5,
                                      2, 6,
                                      3, 7,
                                      4, 8 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0UL );
}

//! Shift node ids to zero in triangle mesh
template<> template<>
void Reorder_object::test< 3 >() {
  set_test_name( "shiftToZero for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0UL );
}

//! Shift node ids to zero in tetrahedron-only mesh
template<> template<>
void Reorder_object::test< 4 >() {
  set_test_name( "shiftToZero for tetrahedra" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0UL );
}

//! Renumber triangle mesh
template<> template<>
void Reorder_object::test< 5 >() {
  set_test_name( "renumber triangle mesh" );

  // Mesh connectivity for simple triangle mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Renumber triangle mesh
  const auto psup = tk::genPsup( inpoel, 3, tk::genEsup( inpoel, 3 ) );
  auto map = tk::renumber( psup );
  tk::remap( inpoel, map );

  // Test the result of reordering
  std::vector< std::size_t > correct_renumbered_inpoel{ 0, 4, 1,
                                                        0, 2, 4,
                                                        1, 4, 7,
                                                        7, 4, 2,
                                                        3, 8, 12,
                                                        3, 12, 10,
                                                        8, 13, 12,
                                                        13, 10, 12,
                                                        0, 1, 5,
                                                        0, 5, 3,
                                                        1, 8, 5,
                                                        3, 5, 8,
                                                        1, 7, 9,
                                                        1, 9, 8,
                                                        7, 13, 9,
                                                        8, 9, 13,
                                                        7, 2, 11,
                                                        7, 11, 13,
                                                        2, 10, 11,
                                                        13, 11, 10,
                                                        0, 6, 2,
                                                        0, 3, 6,
                                                        2, 6, 10,
                                                        3, 10, 6 };
  ensure( "reordered triangle mesh incorrect",
          inpoel == correct_renumbered_inpoel );
}

//! Renumber tetrahedron mesh
template<> template<>
void Reorder_object::test< 6 >() {
  set_test_name( "renumber tetrahedron mesh" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Renumber triangle mesh
  const auto psup = tk::genPsup( inpoel, 4, tk::genEsup( inpoel, 4 ) );
  auto map = tk::renumber( psup );
  tk::remap( inpoel, map );

  // Test the result of reordering
  std::vector< std::size_t > correct_renumbered_inpoel{ 9, 6, 4, 5,
                                                        12, 6, 11, 9,
                                                        6, 11, 9, 4,
                                                        12, 6, 9, 5,
                                                        0, 6, 3, 5,
                                                        13, 8, 12, 9,
                                                        6, 10, 3, 12,
                                                        10, 13, 12, 11,
                                                        13, 11, 7, 9,
                                                        0, 2, 6, 4,
                                                        11, 2, 7, 4,
                                                        7, 1, 9, 4,
                                                        2, 10, 6, 11,
                                                        8, 3, 12, 5,
                                                        0, 1, 4, 5,
                                                        1, 8, 9, 5,
                                                        8, 12, 9, 5,
                                                        1, 9, 4, 5,
                                                        3, 6, 12, 5,
                                                        6, 10, 12, 11,
                                                        11, 7, 9, 4,
                                                        13, 12, 11, 9,
                                                        6, 2, 11, 4,
                                                        6, 0, 4, 5 };

  ensure( "reordered tetrahedron mesh incorrect",
          inpoel == correct_renumbered_inpoel );
}

//! Test all positive Jacbians in tetrahedron mesh
template<> template<>
void Reorder_object::test< 7 >() {
  set_test_name( "all positive Jacobians" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  ensure( "Not all Jacobians are positive",
          tk::positiveJacobians( inpoel, tetcoord ) );
}

//! Test not all positive Jacbians in tetrahedron mesh
template<> template<>
void Reorder_object::test< 8 >() {
  set_test_name( "not all positive Jacobians" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Switch two vertices of a tet
  std::swap( inpoel[4*4+0], inpoel[4*4+1] );

  ensure( "All Jacobians are positive",
          !tk::positiveJacobians( inpoel, tetcoord ) );
}

//! \brief Test if positiveJacobians throws on inpoel non-divisible by the
//!   number of nodes per elements
template<> template<>
void Reorder_object::test< 9 >() {
  set_test_name( "positiveJacobians throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::positiveJacobians( inpoel, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on empty inpoel
template<> template<>
void Reorder_object::test< 10 >() {
  set_test_name( "positivaJacobians throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > empty;
    tk::positiveJacobians( empty, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on empty coord
template<> template<>
void Reorder_object::test< 11 >() {
  set_test_name( "positivaJacobians throws with empty coord" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    decltype(tetcoord) empty;
    tk::positiveJacobians( tetinpoel, empty );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on unique(inpoel).size != coord.size
template<> template<>
void Reorder_object::test< 12 >() {
  set_test_name( "posJacobians throws w inconsistent inp & crd" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    // Shift node IDs to start from zero
    auto inpoel = tetinpoel;
    tk::shiftToZero( inpoel );
    // Remove last coordinate from coord[0]
    auto coord = tetcoord;
    coord[0].pop_back();
    tk::positiveJacobians( tetinpoel, coord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws with non-zero-based inpoel
template<> template<>
void Reorder_object::test< 13 >() {
  set_test_name( "positivaJacobians throws with non-zero-based inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    // Do not shift node IDs to start from zero
    auto inpoel = tetinpoel;
    tk::positiveJacobians( inpoel, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! In-place remap vector of reals using a vector
template<> template<>
void Reorder_object::test< 14 >() {
  set_test_name( "in-place remap vector of reals" );

  // feed good data, test correct result
  std::vector< tk::real > a2{ 1.1, 2.2, 3.3, 4.4 };
  std::vector< std::size_t > r2{ 3, 1, 0, 2 };
  tk::remap( a2, r2 );

  ensure( "in-place remap of real vector incorrect",
          a2 == std::vector< tk::real >{ 3.3, 2.2, 4.4, 1.1 } );
  ensure( "map after in-place remap of ulong vector modified",
          r2 == std::vector< std::size_t >{ 3, 1, 0, 2 } );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // feed empty data
  std::vector< tk::real > ae{};
  std::vector< std::size_t > r1{ 3, 1, 0, 2 };
  tk::remap( ae, r1 );
  // data still empty
  ensure( "in-place remap of empty real vector modified", ae.empty() );

  // feed empty map
  std::vector< tk::real > a1{ 1.1, 2.2, 3.3 };
  std::vector< std::size_t > re;
  tk::remap( a1, re );
  // data intact
  ensure( "in-place remap of real vector modified if map is empty",
          a1 == std::vector< tk::real >{ 1.1, 2.2, 3.3 } );

  // feed unequal-size data and map
  try {
    std::vector< tk::real > a{ 1.1, 2.2, 3.3 };
    std::vector< std::size_t > r{ 3, 1, 0, 2 };
    tk::remap( a, r );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }

  // feed bad map
  try {
    std::vector< tk::real > a3{ 1.1, 2.2, 3.3 };
    std::vector< std::size_t > r{ 4, 1, 0, 2 }; // 4 would index out of a3
    tk::remap( a3, r );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! In-place remap vector of unsigned longs using a vector
template<> template<>
void Reorder_object::test< 15 >() {
  set_test_name( "in-place remap vector of unsigned longs using a vector" );

  // feed good data, test correct result
  std::vector< std::size_t > a{ 0, 1, 2, 3 };
  const std::vector< std::size_t > r{ 3, 1, 0, 2 };
  tk::remap( a, r );

  ensure( "map after remap of ulong vector modified",
          r == std::vector< std::size_t >{ 3, 1, 0, 2 } );
  ensure( "remap of ulong vector incorrect",
          a == std::vector< std::size_t >{ 3, 1, 0, 2 } );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // feed empty data
  std::vector< std::size_t > ae{};
  const std::vector< std::size_t > r1{3,1,0,2};
  tk::remap( ae, r1 );
  // source/destination data still empty
  ensure( "remap of empty ulong vector modified", ae.empty() );
  // source map intact
  ensure( "map not intact after remap of empty ulong vector",
          r1 == std::vector< std::size_t >{3,1,0,2} );

  // feed empty map
  std::vector< std::size_t > a1{ 1, 2, 3 };
  const std::vector< std::size_t > re;
  tk::remap( a1, re );
  // source/destination data intact
  ensure( "data modified after remap with empty map",
          a1 == std::vector< std::size_t >{1,2,3} );
  // source map intact
  ensure( "map not empty after remap with empty map", re.empty() );

  // feed bad map
  try {
    std::vector< std::size_t > a2{ 1, 2, 4 };    // 4 will index out of map
    const std::vector< std::size_t > rr{ 3, 1, 2, 0 };
    tk::remap( a2, rr );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Out-of-place remap vector of unsigned longs using a vector
template<> template<>
void Reorder_object::test< 16 >() {
  set_test_name( "remap vector of unsigned longs using a vector" );

  // feed good data, test correct result
  const std::vector< std::size_t > a{ 0, 1, 2, 3 };
  const std::vector< std::size_t > r{ 3, 1, 0, 2 };
  auto b = tk::remap( a, r );

  ensure( "src data after remap of ulong vector modified",
          a == std::vector< std::size_t >{ 0, 1, 2, 3 } );
  ensure( "map after remap of ulong vector modified",
          r == std::vector< std::size_t >{ 3, 1, 0, 2 } );
  ensure( "remap of ulong vector incorrect",
          b == std::vector< std::size_t >{ 3, 1, 0, 2 } );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // feed empty data
  const std::vector< std::size_t > ae{};
  const std::vector< std::size_t > r1{3,1,0,2};
  auto b1 = tk::remap( ae, r1 );
  // return empty vector
  ensure( "return of non-empty vector from remap of empty input", b1.empty() );
  // source data still empty
  ensure( "remap of empty ulong vector modified", ae.empty() );
  // source map intact
  ensure( "map not intact after remap of empty ulong vector",
          r1 == std::vector< std::size_t >{3,1,0,2} );

  // feed empty map
  const std::vector< std::size_t > a1{ 1, 2, 3 };
  const std::vector< std::size_t > re;
  auto b2 = tk::remap( a1, re );
  // return source vector
  ensure( "src and returned vector differs after remap with empty map",
          b2 == a1 );
  // source data intact
  ensure( "src data modified after remap with empty map",
          a1 == std::vector< std::size_t >{1,2,3} );
  // source map intact
  ensure( "map not empty after remap with empty map", re.empty() );

  // feed bad map
  try {
    const std::vector< std::size_t > a2{ 1, 2, 4 };    // 4 will index out of map
    const std::vector< std::size_t > rr{ 3, 1, 2, 0 };
    auto b3 = tk::remap( a2, rr );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Out-of-place remap vector of unsigned longs using a hash-map
template<> template<>
void Reorder_object::test< 17 >() {
  set_test_name( "remap vector of unsigned longs using a hash-map" );

  // feed good data, test correct result
  const std::vector< std::size_t > a{ 0, 1, 2, 3 };
  const std::unordered_map< std::size_t, std::size_t >
          r{{0,3},{1,1},{2,0},{3,2}};
  auto b = tk::remap( a, r );

  ensure( "src data after remap of ulong vector modified",
          a == std::vector< std::size_t >{ 0, 1, 2, 3 } );
  ensure( "map after remap of ulong vector modified",
          r == std::unordered_map< std::size_t, std::size_t >
               {{0,3},{1,1},{2,0},{3,2}} );
  ensure( "remap of ulong vector incorrect",
          b == std::vector< std::size_t >{ 3, 1, 0, 2 } );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  // feed empty data
  const std::vector< std::size_t > ae{};
  const std::unordered_map< std::size_t, std::size_t >
    r1{{0,3}, {1,1}, {2,0}, {3,2}};
  auto b1 = tk::remap( ae, r1 );
  // data still empty
  ensure( "remap of empty ulong vector modified", ae.empty() );

  // feed empty map
  try {
    const std::vector< std::size_t > a1{ 1, 2, 3 };
    const std::unordered_map< std::size_t, std::size_t > re{};
    auto b2 = tk::remap( a1, re );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }

  // feed bad map
  try {
    const std::vector< std::size_t > a1{ 1, 2, 4 }; // 4 is not in map keys
    const std::unordered_map< std::size_t, std::size_t >
            rr{{0,3},{1,1},{2,0},{3,2}};
    auto b3 = tk::remap( a1, rr );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Out-of-place remap map of vectors of unsigned longs using a hash-map
template<> template<>
void Reorder_object::test< 18 >() {
  set_test_name( "remap map of vectors of ulongs using a hash-map" );

  // feed good data, test correct result
  std::map< int, std::vector< std::size_t > > a{
    {1,{0,1,2,3}}, {32,{1,2,3,0}}, {42,{2,0,1}}, {12,{3,1,2,0}} };
  std::unordered_map< std::size_t, std::size_t > r{{0,3},{1,1},{2,0},{3,2}};
  auto b = tk::remap( a, r );

  ensure( "src data after remap of map of ulong vectors modified",
          a == std::map< int, std::vector< std::size_t > >{
             {1,{0,1,2,3}}, {32,{1,2,3,0}}, {42,{2,0,1}}, {12,{3,1,2,0}} } );
  ensure( "map after remap of map of ulong vectors modified",
          r == std::unordered_map< std::size_t, std::size_t >
               {{0,3},{1,1},{2,0},{3,2}} );
  ensure( "remap of map of ulong vectors incorrect",
          b == std::map< int, std::vector< std::size_t > >{
             {1,{3,1,0,2}}, {32,{1,0,2,3}}, {42,{0,3,1}}, {12,{2,1,0,3}} } );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  // feed empty data
  std::map< int, std::vector< std::size_t > > ae{};
  std::unordered_map< std::size_t, std::size_t > r1{{0,3}, {1,1}, {2,0}, {3,2}};
  auto b1 = tk::remap( ae, r1 );
  // data still empty
  ensure( "remap of empty ulong vector modified", ae.empty() );
  // return container of ids empty
  ensure( "remap of empty map non-empty", b1.empty() );
  // map intact
  ensure( "map modified after remap of empty map of vectors of ulongs",
          r1 == std::unordered_map< std::size_t, std::size_t >
                                  {{0,3}, {1,1}, {2,0}, {3,2}} );

  // feed empty map
  try {
    std::map< int, std::vector< std::size_t > > a1{{0,{1}}, {42,{2}}, {23,{3}}};
    std::unordered_map< std::size_t, std::size_t > re;
    auto b2 = tk::remap( a1, re );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }

  // feed bad map
  try {
    // id 4 is not in map keys
    std::map< int, std::vector< std::size_t > >
      a1{ {0,{1}}, {42,{2,1}}, {23,{4}} };
    std::unordered_map< std::size_t, std::size_t > rr{{0,3},{1,1},{2,0},{3,2}};
    auto b3 = tk::remap( a1, rr );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
