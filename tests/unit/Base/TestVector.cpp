// *****************************************************************************
/*!
  \file      tests/unit/Base/TestVector.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/Vector.hpp
  \details   Unit tests for Base/Vector.hpp
*/
// *****************************************************************************

#include <unistd.h>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Vector.hpp"
#include "Types.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct Vector_common {
  // cppcheck-suppress unusedStructMember
  double precision = 1.0e-15;    // required floating-point precision
};

//! Test group shortcuts
using Vector_group = test_group< Vector_common, MAX_TESTS_IN_GROUP >;
using Vector_object = Vector_group::object;

//! Define test group
static Vector_group Vector( "Base/Vector" );

//! Test definitions for group

//! Test cross product
template<> template<>
void Vector_object::test< 1 >() {
  set_test_name( "cross product" );

  std::array< tk::real, 3 > v1{{ 3.0, -3.0, 1.0 }},
                            v2{{ 4.0, 9.0, 2.0 }},
                            correct_result{{ -15.0, -2.0, 39.0 }};

  const auto result = tk::cross( v1, v2 );
  ensure_equals( "cross product incorrect",
                 result[0], correct_result[0], precision );
  ensure_equals( "cross product incorrect",
                 result[1], correct_result[1], precision );
  ensure_equals( "cross product incorrect",
                 result[2], correct_result[2], precision );
}

//! Test cross product divided by scalar
template<> template<>
void Vector_object::test< 2 >() {
  set_test_name( "cross product divided by scalar" );

  std::array< tk::real, 3 > v1{{ 3.0, -3.0, 1.0 }},
                            v2{{ 4.0, 9.0, 2.0 }},
                            correct_result{{ -7.5, -1.0, 19.5 }};

  const auto result = tk::crossdiv( v1, v2, 2.0 );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[0], correct_result[0], precision );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[1], correct_result[1], precision );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[2], correct_result[2], precision );
}

//! Test dot product
template<> template<>
void Vector_object::test< 3 >() {
  set_test_name( "dot product" );

  std::array< tk::real, 3 > v1{{ 1.0, 2.0, 3.0 }}, v2{{ 4.0, -5.0, 6.0 }};
  tk::real correct_result = 12.0;

  const auto result = tk::dot( v1, v2 );
  ensure_equals( "dot product incorrect", result, correct_result, precision );
}

//! Test triple product
template<> template<>
void Vector_object::test< 4 >() {
  set_test_name( "triple product" );

  std::array< tk::real, 3 > v1{{ -1.0, 3.0, 3.0 }},
                            v2{{ -2.0, 3.0, 1.0 }},
                            v3{{  0.0, 4.0, 0.0 }};
  tk::real correct_result = -20.0;

  const auto result = tk::triple( v1, v2, v3 );
  ensure_equals( "triple product incorrect", result, correct_result,
                 precision );
}

//! Test vector length
template<> template<>
void Vector_object::test< 5 >() {
  set_test_name( "length" );

  std::array< tk::real, 3 > v1{{ -1.0, 3.0, 3.0 }}, v2{{  0.0, 4.0, 0.0 }};
  tk::real correct_result_1 = 4.358898943540674;
  tk::real correct_result_2 = 4.0;

  const auto result_1 = tk::length( v1 );
  const auto result_2 = tk::length( v2 );
  ensure_equals( "length incorrect", result_1, correct_result_1, precision );
  ensure_equals( "length incorrect", result_2, correct_result_2, precision );
}

//! Test normalizing a vector
template<> template<>
void Vector_object::test< 6 >() {
  set_test_name( "unit" );

  std::array< tk::real, 3 > v1{{ -1.0, 3.0, 3.0 }}, v2{{  0.0, 4.0, 0.0 }};

  tk::unit( v1 );
  tk::unit( v2 );
  ensure_equals( "unit incorrect", tk::length(v1), 1.0, precision );
  ensure_equals( "unit incorrect", tk::length(v2), 1.0, precision );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
