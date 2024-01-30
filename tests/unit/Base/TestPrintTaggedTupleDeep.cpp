// *****************************************************************************
/*!
  \file      tests/unit/Base/TestPrintTaggedTupleDeep.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for TaggedTuple deep printer
*/
// *****************************************************************************

#include <sstream>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "PrintTaggedTupleDeep.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct PrintTaggedTupleDeep_common {
  // Tags
  DEFTAG( name );
  DEFTAG( age );
  DEFTAG( email );
  DEFTAG( tag1 );
  DEFTAG( tag2 );
  DEFTAG( tag3 );

  // Define a tk::TaggedTuple by inheriting from TaggedTuple
  struct Cmd : public tk::TaggedTuple< brigand::list<
                        name,  std::string,
                        age,   int,
                        email, std::string,
                        tag1,  tk::TaggedTuple< brigand::list <
                                  tag2, std::string,
                                  tag3, std::string > > > >
  {};

  // Constructor
  PrintTaggedTupleDeep_common() {
    cmd.get< name >() = "Bob";
    cmd.get< age >() = 32;
    cmd.get< email >() = "bob@google.com";
    auto& t1 = cmd.get< tag1 >();
    t1.get< tag2 >() = "string2";
    t1.get< tag3 >() = "string3";
  }

  Cmd cmd;
};

//! Test group shortcuts
using PrintTaggedTupleDeep_group =
  test_group< PrintTaggedTupleDeep_common, MAX_TESTS_IN_GROUP >;
using PrintTaggedTupleDeep_object = PrintTaggedTupleDeep_group::object;

//! Define test group
static PrintTaggedTupleDeep_group
  PrintTaggedTupleDeep( "Base/PrintTaggedTupleDeep" );

//! Test definitions for group

//! Test print() of TaggedTuple with depth/indentation
template<> template<>
void PrintTaggedTupleDeep_object::test< 1 >() {
  set_test_name( "print()" );

  std::stringstream s;
  tk::print( s, cmd );
  ensure_equals( "print()", s.str(),
R"({
  name = Bob
  age = 32
  email = bob@google.com
  tag1 = {
    tag2 = string2
    tag3 = string3
  }
}
)" );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
