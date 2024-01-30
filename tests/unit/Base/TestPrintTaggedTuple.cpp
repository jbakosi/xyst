// *****************************************************************************
/*!
  \file      tests/unit/Base/TestPrintTaggedTuple.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for TaggedTuple printer
*/
// *****************************************************************************

#include <sstream>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "PrintTaggedTuple.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct PrintTaggedTuple_common {
  // Tags
  DEFTAG( name );
  DEFTAG( age );
  DEFTAG( email );
  DEFTAG( tag1 );
  DEFTAG( tag2 );
  DEFTAG( tag3 );

  using MemberList = brigand::list<
    name,  std::string,
    age,   int,
    email, std::string,
    tag1,  tk::TaggedTuple< brigand::list <
              tag2, std::string,
              tag3, std::string > > >;

  // Define a tagged tuple: odd template arguments are tags, even ones are types
  using record = tk::TaggedTuple< MemberList >;

  // Constructor
  PrintTaggedTuple_common() {
    tup.get< name >() = "Bob";
    tup.get< age >() = 32;
    tup.get< email >() = "bob@google.com";
    auto& t1 = tup.get< tag1 >();
    t1.get< tag2 >() = "string2";
    t1.get< tag3 >() = "string3";
  }

  record tup;
};

//! Test group shortcuts
using PrintTaggedTuple_group =
  test_group< PrintTaggedTuple_common, MAX_TESTS_IN_GROUP >;
using PrintTaggedTuple_object = PrintTaggedTuple_group::object;

//! Define test group
static PrintTaggedTuple_group PrintTaggedTuple( "Base/PrintTaggedTuple" );

//! Test definitions for group

//! Test operator<< of TaggedTuple
template<> template<>
void PrintTaggedTuple_object::test< 1 >() {
  set_test_name( "operator<<" );

  std::stringstream s;
  s << tup;
  ensure_equals( "operator<<(TaggedTuple)", s.str(), "name = Bob age = 32 "
    "email = bob@google.com tag1 = { tag2 = string2 tag3 = string3 } " );
}

//! Test print()
template<> template<>
void PrintTaggedTuple_object::test< 2 >() {
  set_test_name( "print()" );

  std::stringstream s;
  tk::print( s, tup );
  ensure_equals( "print()", s.str(), "name = Bob age = 32 "
    "email = bob@google.com tag1 = { tag2 = string2 tag3 = string3 } " );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
