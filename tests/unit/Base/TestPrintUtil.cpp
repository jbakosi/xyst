// *****************************************************************************
/*!
  \file      tests/unit/Base/TestPrintUtil.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/PrintUtil.hpp
  \details   Unit tests for Base/PrintUtil.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "PrintUtil.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct PrintUtil_common {};

//! Test group shortcuts
using PrintUtil_group = test_group< PrintUtil_common, MAX_TESTS_IN_GROUP >;
using PrintUtil_object = PrintUtil_group::object;

//! Define test group
static PrintUtil_group PrintUtil( "Base/PrintUtil" );

//! Test definitions for group

//! Test tk::operator<< used for writing enum class value to output stream
template<> template<>
void PrintUtil_object::test< 1 >() {
  set_test_name( "tk::operator<<( enum class )" );

  enum class Enum { FIRST=3, SECOND, THIRD };
  using tk::operator<<;
  std::stringstream ss;
  ss << Enum::FIRST;
  ensure_equals( "enum class first item output to stream", ss.str(), "3" );
  ss.str("");
  ss << Enum::SECOND;
  ensure_equals( "enum class second item output to stream", ss.str(), "4" );
  ss.str("");
  ss << Enum::THIRD;
  ensure_equals( "enum class third item output to stream", ss.str(), "5" );
}

//! Test tk::operator<< used for writing non-enum-class to output stream
template<> template<>
void PrintUtil_object::test< 2 >() {
  set_test_name( "tk::operator<<( non-enum-class )" );

  using tk::operator<<;
  std::stringstream ss;
  ss << "blah";
  ensure_equals( "non-enum-class item output to stream", ss.str(), "blah" );
}

//! Test tk::operator<< used for concatenating to std::basic_string for lvalues
template<> template<>
void PrintUtil_object::test< 3 >() {
  set_test_name( "tk::operator<<( std::basic_string& )" );

  using tk::operator<<;
  std::string s;
  s << 123;
  ensure_equals( "concatenate integer to std::string", s, "123" );
  s << 3.14;
  ensure_equals( "concatenate double to std::string", s, "1233.14" );
  s << " blah";
  ensure_equals( "concatenate const char* to std::string", s, "1233.14 blah" );
  s << std::string("string");
  ensure_equals( "concatenat std::string to std::string",
                 s, "1233.14 blahstring" );
  s = "";
  ensure_equals( "std::string empty after concatenations", s.size(), 0UL );
}

//! Test tk::operator<< used for concatenating to std::basic_string for rvalues
template<> template<>
void PrintUtil_object::test< 4 >() {
  set_test_name( "tk::operator<<( std::basic_string&& )" );

  using tk::operator<<;
  // lhs to operator<< are rvalues
  ensure_equals( "concatenate integer to std::string",
                 std::string( "start" ) << 123, "start123" );
  ensure_equals( "concatenate double to std::string",
                 std::string( "start" ) << 3.14, "start3.14" );
  ensure_equals( "concatenate const char* to std::string",
                 std::string( "start" ) << " blah", "start blah" );
  ensure_equals( "concatenate const char* to std::string",
                 std::string( "start" ) << std::string("blah"), "startblah" );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
