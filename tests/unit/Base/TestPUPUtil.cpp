// *****************************************************************************
/*!
  \file      tests/unit/Base/TestPUPUtil.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/PUPUtil.hpp
  \details   Unit tests for Base/PUPUtil.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "PUPUtil.hpp"
#include "MigratedTypes.hpp"

#include "XystConfig.hpp"
#include "NoWarning/tutsuite.decl.h"
#include "migrated.decl.h"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct PUPUtil_common {};

//! Test group shortcuts
using PUPUtil_group = test_group< PUPUtil_common, MAX_TESTS_IN_GROUP >;
using PUPUtil_object = PUPUtil_group::object;

//! Define test group
static PUPUtil_group PUPUtil( "Base/PUPUtil" );

//! Test definitions for group

//! Charm++ chare to test Pack/Unpack utilities during network migration
class Migrated : public CBase_Migrated {

  public:

  //! Constructor taking (and migrating) a default strongly-typed enum
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Enum_default e ) : m_enum_default(e)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 1,
                         "Charm:migrate enum 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      using underlying_type =
        typename std::underlying_type< charm::Enum_default >::type;
      std::string expected = std::to_string(
       static_cast< underlying_type >( charm::Enum_default::F1 ) );
      std::string actual = std::to_string(
       static_cast< underlying_type >( m_enum_default ) );
      // Evaluate test
      ensure( "strongly-typed enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_default == charm::Enum_default::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a uint8_t strongly-typed enum
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Enum_uint8_t e ) : m_enum_uint8_t(e)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 2,
                         "Charm:migrate uint8_t enum 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      using underlying_type =
        typename std::underlying_type< charm::Enum_uint8_t >::type;
      std::string expected = std::to_string(
       static_cast< underlying_type >( charm::Enum_uint8_t::F1 ) );
      std::string actual = std::to_string(
       static_cast< underlying_type >( m_enum_uint8_t ) );
      // Evaluate test
      ensure( "strongly-typed enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_uint8_t == charm::Enum_uint8_t::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a C-style enum
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Enum_cstyle e ) : m_enum_cstyle(e)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 3,
                         "Charm:migrate C-style enum 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = std::to_string( charm::Enum_cstyle::F1 );
      std::string actual = std::to_string( m_enum_cstyle );
      // Evaluate test
      ensure( "C-style enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_cstyle == charm::Enum_cstyle::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::pair<int,double>
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Pair p ) : m_pair(p)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 4,
                         "Charm:migrate std::pair<int,double> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      double precision = 1.0e-15;    // required floating-point precision
      // Evaluate test
      ensure_equals( "std::pair (1st) different after migrated: "
                     "expected `2` actual `" + std::to_string(m_pair.first) +
                     "`", m_pair.first, 2.0, precision );
      ensure_equals( "std::pair (2nd) different after migrated: "
                     "expected `3.14` actual `" + std::to_string(m_pair.second)
                     + "`", m_pair.second, 3.14, precision );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::vector< std::string >
  // cppcheck-suppress uninitMemberVar
  // cppcheck-suppress passedByValue
  explicit Migrated( charm::Vector v ) : m_vector(v)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 5,
                         "Charm:migrate std::vector< std::string > 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ \"1\", \"blah\", \"boohoo\" ]";
      std::string actual = "[ " + m_vector[0] +
                           ", " + m_vector[1] +
                           ", " + m_vector[2] + " ]";
      // Evaluate test
      ensure( "std::vector different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_vector == charm::Vector{ "1", "blah", "boohoo" } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::tuple
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Tuple t ) : m_tuple(t)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 6,
                         "Charm:migrate std::tuple 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ 2, 6.42, {\"woodoo\", \"boohoo\"}, "
                             "32, {{14, \"blah\"}, {15, \"blahblah\"}} ]";
      ensure_equals( "std::tuple different after migrated: vector size is "
                     "different", std::get< 2 >( m_tuple ).size(), 2UL );
      ensure_equals( "std::tuple different after migrated: map size is "
                     "different", std::get< 4 >( m_tuple ).size(), 2UL );
      using underlying_type_def =
        typename std::underlying_type< charm::Enum_default >::type;
      using underlying_type_uin =
        typename std::underlying_type< charm::Enum_uint8_t >::type;
      const auto& m = std::get< 4 >( m_tuple );
      const auto& i1 = m.find( charm::Enum_uint8_t::F1 );
      std::string m1, m2;
      if ( i1 != end(m) )
        m1 = "{" + std::to_string(
                     static_cast< underlying_type_uin >(
                       charm::Enum_uint8_t::F1 ) )
                 + ", \"" + i1->second + "\"}";
      else
        fail("std::tuple different after migrated: cannot find 1st key in map");
      const auto& i2 = m.find( charm::Enum_uint8_t::F2 );
      if ( i2 != end(m) )
        m2 = "{" + std::to_string(
                     static_cast< underlying_type_uin >(
                       charm::Enum_uint8_t::F2 ) )
                 + ", \"" + i2->second + "\"}";
      else
        fail("std::tuple different after migrated: cannot find 2nd key in map");
      std::string actual = "[ " + std::to_string( std::get< 0 >( m_tuple ) ) +
                           ", " + std::to_string( std::get< 1 >( m_tuple ) ) +
                        ", {\"" + std::get< 2 >( m_tuple )[ 0 ] +
                       "\", \"" + std::get< 2 >( m_tuple )[ 1 ] +
                        "\"}, " + std::to_string(
                                    static_cast< underlying_type_def >(
                                      std::get< 3 >( m_tuple ) ) ) +
                         ", {" + m1 + ", " + m2 + "} ]";
      // Evaluate test
      ensure( "std::tuple different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_tuple ==
                  charm::Tuple{ 2, 6.42, {"woodoo", "boohoo"},
                                charm::Enum_default::F1,
                                std::map< charm::Enum_uint8_t, std::string >{
                                  { charm::Enum_uint8_t::F1, "blah" },
                                  { charm::Enum_uint8_t::F2, "blahblah" } } } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::array
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Array a ) : m_array(a)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 7,
                         "Charm:migrate std::array<int,2> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ 12, 2 ]";
      std::string actual = "[ " + std::to_string( m_array[0] ) + ", "
                                + std::to_string( m_array[1] ) + " ]";
      // Evaluate test
      ensure( "std::array<int,2> different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_array == charm::Array{ { 12, 2 } } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::unordered_map
  // cppcheck-suppress uninitMemberVar
  // cppcheck-suppress passedByValue
  explicit Migrated( charm::UnorderedMap m ) : m_unordered_map(m)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 8,
                         "Charm:migrate std::unordered_map<int,str> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = R"([ {11, "eleven"} {12, "twelve"} ])";
      std::string actual = "[ ";
      for (const auto& n : m_unordered_map) {
        actual += "{" + std::to_string(n.first) + ", " + n.second + "} ";
      }
      actual += "]";
      // Evaluate test
      ensure( "std::unordered_map<int,str> different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_unordered_map ==
                charm::UnorderedMap{ {11,"eleven"}, {12,"twelve"} } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::unordered_set
  // cppcheck-suppress uninitMemberVar
  // cppcheck-suppress passedByValue
  explicit Migrated( charm::UnorderedSet s ) : m_unordered_set(s)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 9,
                         "Charm:migrate std::unordered_set<int> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = R"([ {11} {12} ])";
      std::string actual = "[ ";
      for (auto n : m_unordered_set) actual += "{" + std::to_string(n) + "} ";
      actual += "]";
      // Evaluate test
      ensure( "std::unordered_set<int> different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_unordered_set ==
                charm::UnorderedSet{ 11, 12 } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::optional< std::string >
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::OptionalStr o ) : m_optional_str(o)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 10,
                         "Charm:migrate std::optional<str> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "blah";
      std::string actual = o ? *o : "";
      // Evaluate test
      ensure( "std::optional<str> different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_optional_str == charm::OptionalStr{ { "blah" } } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) an uninitialized std::optional< int >
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::OptionalInt o ) : m_optional_int(o)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 11,
                         "Charm:migrate std::optional<int> 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "std::nullopt";
      std::string actual = o ? std::to_string(*o) : "std::nullopt";
      // Evaluate test
      ensure( "std::optional<int> different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_optional_int == charm::OptionalInt{ std::nullopt } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a tk::tuple::tagged_tuple
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::TaggedTuple t ) : m_tagged_tuple(t)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 12,
                         "Charm:migrate tk::tuple::tagged_tuple 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = R"([ "Bob", 32, "bob@bob.bob" ])";
      std::string actual = "[ " +
        m_tagged_tuple.get< charm::tag::name >() + ", " +
        std::to_string( m_tagged_tuple.get< charm::tag::age >() ) +
        m_tagged_tuple.get< charm::tag::email >() + " ]";
      // Evaluate test
      ensure("tk::tuple::tagged_tuple different after migrated: "
             "expected `" + expected + "` actual `" + actual + "`",
             m_tagged_tuple == charm::TaggedTuple{{"Bob", 32, "bob@bob.bob"}});
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! \brief Constructor taking (and migrating) a std::variant<int,double>
  //!   holding int
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Variant v, int value ) : m_variant(v)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 13,
                         "Charm:migrate std::variant<int,double>(int) 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ " + std::to_string(value) + " ]";
      std::string actual = "[ " +
        std::to_string( std::get< int >( m_variant ) ) + " ]";
      // Evaluate test
      ensure_equals( "std::variant<int,double>(int) different after migrated:"
                     " expected `" + expected + "` actual `" + actual + "`",
                     std::get< int >( m_variant ), value );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! \brief Constructor taking (and migrating) a std::variant<int,double>
  //!   holding double
  // cppcheck-suppress uninitMemberVar
  explicit Migrated( charm::Variant v, double value ) : m_variant(v)
  {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 14,
                         "Charm:migrate std::variant<int,double>(double) 2",
                         tut::test_result::result_type::ok );

    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ " + std::to_string(value) + " ]";
      std::string actual = "[ " +
        std::to_string( std::get< double >( m_variant ) ) + " ]";
      // Evaluate test
      ensure_equals( "rk::Variant<int,double>(double) different after "
                     "migrated: expected `" + expected + "` actual `" + actual +
                     "`", std::get< double >( m_variant ), value,
                     1.0e-15 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluateTest(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  charm::Enum_default m_enum_default;
  charm::Enum_uint8_t m_enum_uint8_t;
  charm::Enum_cstyle m_enum_cstyle;
  charm::Pair m_pair;
  charm::Vector m_vector;
  charm::Tuple m_tuple;
  charm::Array m_array;
  charm::UnorderedMap m_unordered_map;
  charm::UnorderedSet m_unordered_set;
  charm::OptionalStr m_optional_str;
  charm::OptionalInt m_optional_int;
  charm::TaggedTuple m_tagged_tuple;
  charm::Variant m_variant;
};

//! Test Pack/Unpack of a default strongly-typed enum
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 1 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_default::F1 );
}

//! Test Pack/Unpack of a uint8_t strongly-typed enum
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 2 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate uint8_t enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_uint8_t::F1 );
}

//! Test Pack/Unpack of a C-style enum
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 3 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate C-style enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_cstyle::F1 );
}

//! Test Pack/Unpack of a std::pair< int, double >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 4 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::pair<int,double> 1" );

  CProxy_Migrated::ckNew( charm::Pair( 2, 3.14 ) );
}

//! Test Pack/Unpack of a std::vector< std::string >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 5 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::vector< std::string > 1" );

  CProxy_Migrated::ckNew( charm::Vector{"1", "blah", "boohoo"} );
}

//! \brief Test Pack/Unpack of a std::tuple< int, double,
//!    std::vector< std::string >, strongly-typed enum,
//!    std::map< strongly-typed(uint8_t) enum, std::string > >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 6 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::tuple 1" );

  CProxy_Migrated::ckNew(
    charm::Tuple{ 2, 6.42, {"woodoo", "boohoo"},
                  charm::Enum_default::F1,
                  std::map< charm::Enum_uint8_t, std::string >{
                    { charm::Enum_uint8_t::F1, "blah" },
                    { charm::Enum_uint8_t::F2, "blahblah" } } } );
}

//! Test Pack/Unpack of a std::array< int, 2 >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 7 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::array<int,2> 1" );

  CProxy_Migrated::ckNew( charm::Array{ { 12, 2 } } );
}

//! Test Pack/Unpack of a std::unordered_map< int, std::string >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 8 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::unordered_map<int,str> 1" );

  CProxy_Migrated::ckNew( charm::UnorderedMap{ {11,"eleven"}, {12,"twelve"} } );
}

//! Test Pack/Unpack of a std::unordered_set< int >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 9 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::unordered_set<int> 1" );

  CProxy_Migrated::ckNew( charm::UnorderedSet{ 11, 12 } );
}

//! Test Pack/Unpack of a std::optional< std::string >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 10 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::optional<str> 1" );

  CProxy_Migrated::ckNew( charm::OptionalStr{ {"blah"} } );
}

//! Test Pack/Unpack of an uninitialized std::optional< int >
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 11 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::optional<none> 1" );

  CProxy_Migrated::ckNew( charm::OptionalInt{ std::nullopt } );
}

//! Test Pack/Unpack of tk::tuple::tagged_tuple
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 12 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate tk::tuple::tagged_tuple 1" );

  CProxy_Migrated::ckNew( charm::TaggedTuple{{"Bob", 32, "bob@bob.bob"}} );
}

//! Test Pack/Unpack of std::variant<int,double> holding an int
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 13 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::variant<int,double>(int) 1" );

  std::variant< int, double > v{213};
  CProxy_Migrated::ckNew( charm::Variant(v), 213 );
}

//! Test Pack/Unpack of std::variant<int,double> holding a double
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void PUPUtil_object::test< 14 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate std::variant<int,double>(double) 1" );

  std::variant< int, double > v{3.14};
  CProxy_Migrated::ckNew( charm::Variant(v), 3.14 );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT

#include "NoWarning/migrated.def.h"
