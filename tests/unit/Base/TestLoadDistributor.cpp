// *****************************************************************************
/*!
  \file      tests/unit/Base/TestLoadDistributor.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/LoadDistributor
  \details   Unit tests for Base/LoadDistributor
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "LoadDistributor.hpp"
#include "Exception.hpp"

namespace unittest {

extern std::string g_executable;

} // unittest::

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct LoadDistributor_common {};

//! Test group shortcuts
using LoadDistributor_group =
  test_group< LoadDistributor_common, MAX_TESTS_IN_GROUP >;
using LoadDistributor_object = LoadDistributor_group::object;

//! Define test group
static LoadDistributor_group LoadDistributor( "Base/LoadDistributor" );

//! Test definitions for group

//! Test if linear distributor does not throw on bounded virtualization
template<> template<>
void LoadDistributor_object::test< 1 >() {
  set_test_name( "linear doesn't throw on bounded virt" );

  uint64_t chunksize, remainder;
  tk::linearLoadDistributor( 0.5, 1234, 2, chunksize, remainder );
}

//! Test if linear distributor throws on too low virtualization parameter
template<> template<>
void LoadDistributor_object::test< 2 >() {
  set_test_name( "linear throws on too low virt" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {

    uint64_t chunksize, remainder;
    tk::linearLoadDistributor( -0.5, 1234, 2, chunksize, remainder );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "must be between" ) !=
              std::string::npos );
  }
  #endif
}

//! Test if linear distributor throws on to high virtualization parameter
template<> template<>
void LoadDistributor_object::test< 3 >() {
  set_test_name( "linear throws on too high virt" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {

    uint64_t chunksize, remainder;
    tk::linearLoadDistributor( 1.5, 1234, 2, chunksize, remainder );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "must be between" ) !=
              std::string::npos );
  }
  #endif
}

//! Test if linear distributor returns number of chares less than or equal load
template<> template<>
void LoadDistributor_object::test< 4 >() {
  set_test_name( "linear returns sane nchare" );

  uint64_t nchare, chunksize, remainder;
  nchare = tk::linearLoadDistributor( 0.5, 1234, 2, chunksize, remainder );
  // nchare should never be larger than the load
  ensure( "number of chares too large", nchare < 1235 );
}


//! Test if linear distributor returns chunksize less than or equal the load
template<> template<>
void LoadDistributor_object::test< 5 >() {
  set_test_name( "linear returns sane chunksize" );

  uint64_t chunksize, remainder;
  tk::linearLoadDistributor( 0.0, 1234, 2, chunksize, remainder );
  // chunksize should never be larger than the load
  ensure( "chunksize too large", chunksize < 1235 );
}

//! Test if linear distributor returns remainder less than or equal to chunksize
template<> template<>
void LoadDistributor_object::test< 6 >() {
  set_test_name( "linear returns sane remainder" );

  uint64_t chunksize, remainder;
  tk::linearLoadDistributor( 1.0, 1234, 2, chunksize, remainder );
  // remainder should never be larger than chunksize
  ensure( "remainder too large",  remainder < chunksize );
}

//! Test if linear distributor throws on garbage number of processing elements
template<> template<>
void LoadDistributor_object::test< 7 >() {
  set_test_name( "linear throws on garbage npe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {

    uint64_t chunksize, remainder;
    tk::linearLoadDistributor( 0.5, 1234, 0, chunksize, remainder );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "must be larger" ) !=
              std::string::npos );
  }
  #endif
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
