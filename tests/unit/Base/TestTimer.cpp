// *****************************************************************************
/*!
  \file      tests/unit/Base/TestTimer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for tk::Timer
  \details   Unit tests for tk::Timer
*/
// *****************************************************************************

#include <unistd.h>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Timer.hpp"
#include "ContainerUtil.hpp"

#include "NoWarning/tutsuite.decl.h"

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

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
struct Timer_common {
  // cppcheck-suppress unusedStructMember
  double precision = 1.0e-3;    // required precision in seconds for timings
};

//! Test group shortcuts
using Timer_group = test_group< Timer_common, MAX_TESTS_IN_GROUP >;
using Timer_object = Timer_group::object;

//! Define test group
static Timer_group Timer( "Base/Timer" );

//! Test definitions for group

//! Test timing a 0.1s duration as float with given precision
template<> template<>
void Timer_object::test< 1 >() {
  double prec = 1.0e-1; // only for this single test (to pass on Mac OS)
  set_test_name( "measure 0.1s using dsec() with " + std::to_string(prec) +
                 "s prec" );

  tk::Timer timer;
  usleep( 100000 );    // in micro-seconds, sleep for 0.1 second
  // test if time measured with at least 1/10th of a millisecond prec
  ensure_equals( "time 0.1s elapsed as float", timer.dsec(), 0.1, prec );
}

//! Test timing a 1.0 duration as h:m:s with given precision
template<> template<>
void Timer_object::test< 2 >() {
  set_test_name( "measure 1.0s using hms() with " + std::to_string(precision) +
                 "s prec" );

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  const auto stamp = timer.hms();
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 1.0s elapsed as hrs",
                 static_cast<tk::real>(stamp.hrs.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as min",
                 static_cast<tk::real>(stamp.min.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as sec",
                 static_cast<tk::real>(stamp.sec.count()), 1.0, precision );
}

//! Test estimated time elapsed and to accomplishment triggered by term
template<> template<>
void Timer_object::test< 3 >() {
  set_test_name( "ETE and ETA triggered by terminate time" );

  // Setup a duration case
  tk::real term = 5.0;      // time at which to terminate time stepping
  tk::real time = 1.0;      // current time
  uint64_t nstep = 1000;    // max number of time steps to take (large)
  uint64_t it = 1;          // current iteration

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  tk::Timer::Watch ete, eta;
  timer.eta( term, time, nstep, it, 0.0, 0.0, 0.0, ete, eta );
  // test estimated time elapsed with given precision
  ensure_equals( "estimated time elapsed in hrs",
                 static_cast<tk::real>(ete.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in min",
                 static_cast<tk::real>(ete.min.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in sec",
                 static_cast<tk::real>(ete.sec.count()), 1.0, precision );
  // test estimated time to accomplishment with given precision
  ensure_equals( "estimated time to accomplishment in hrs",
                 static_cast<tk::real>(eta.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time to accomplishment in min",
                 static_cast<tk::real>(eta.min.count()), 0.0, precision );
  ensure_equals( "estimated time to accomplishment in sec",
                 static_cast<tk::real>(eta.sec.count()), 4.0, precision );
}

//! Test estimated time elapsed and to accomplishment triggered by nstep
template<> template<>
void Timer_object::test< 4 >() {
  set_test_name( "ETE and ETA triggered by max number of steps" );

  // Setup a duration case
  tk::real term = 500.0;    // time at which to terminate time stepping (large)
  tk::real time = 1.0;      // current time
  uint64_t nstep = 100;     // max number of time steps to take
  uint64_t it = 1;          // current iteration

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  tk::Timer::Watch ete, eta;
  timer.eta( term, time, nstep, it, 0.0, 0.0, 0.0, ete, eta );
  // test estimated time elapsed with given precision
  ensure_equals( "estimated time elapsed in hrs",
                 static_cast<tk::real>(ete.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in min",
                 static_cast<tk::real>(ete.min.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in sec",
                 static_cast<tk::real>(ete.sec.count()), 1.0, precision );
  // test estimated time to accomplishment with given precision
  ensure_equals( "estimated time to accomplishment in hrs",
                 static_cast<tk::real>(eta.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time to accomplishment in min",
                 static_cast<tk::real>(eta.min.count()), 1.0, precision );
  ensure_equals( "estimated time to accomplishment in sec",
                 static_cast<tk::real>(eta.sec.count()), 39.0, precision );
}

//! Test estimated time elapsed and to accomplishment triggered by residuals
template<> template<>
void Timer_object::test< 5 >() {
  set_test_name( "ETE and ETA triggered by target residual" );

  // Setup a duration case
  tk::real term = 500.0;    // time at which to terminate time stepping (unused)
  tk::real time = 1.0;      // current time (unused)
  uint64_t nstep = 100;     // max number of time steps to take (unused)
  uint64_t it = 1;          // current iteration
  tk::real res0 = 1.0e-3;   // previous residual
  tk::real res = 0.9997e-3; // current residual
  tk::real rest = 1.0e-8;   // target residual

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  tk::Timer::Watch ete, eta;
  timer.eta( term, time, nstep, it, res0, res, rest, ete, eta );
  // test estimated time elapsed with given precision
  ensure_equals( "estimated time elapsed in hrs",
                 static_cast<tk::real>(ete.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in min",
                 static_cast<tk::real>(ete.min.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in sec",
                 static_cast<tk::real>(ete.sec.count()), 1.0, precision );
  // test estimated time to accomplishment with given precision
  ensure_equals( "estimated time to accomplishment in hrs",
                 static_cast<tk::real>(eta.hrs.count()), 10.0, 1.0 );
  ensure_equals( "estimated time to accomplishment in min",
                 static_cast<tk::real>(eta.min.count()), 39.0, 100.0 );
  ensure_equals( "estimated time to accomplishment in sec",
                 static_cast<tk::real>(eta.sec.count()), 33.0, 100.0 );
}

//! Test converting a 1.0s duration timed as a float to Timer::Watch
template<> template<>
void Timer_object::test< 6 >() {
  set_test_name( "convert time stamp in float to Watch" );

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  // convert time stamp in float to Timer::Watch
  const auto w = tk::hms( timer.dsec() );
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in hrs",
                 static_cast<tk::real>(w.hrs.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in min",
                 static_cast<tk::real>(w.min.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in sec",
                 static_cast<tk::real>(w.sec.count()), 1.0, precision );
}

//! Charm chare having a tk::Timer object
class CharmTimer : public CBase_CharmTimer {
  public:
  explicit CharmTimer( const tk::Timer& timer ) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/Timer", 7,
                         "Charm:migrate tk::Timer 2",
                         tut::test_result::result_type::ok );

    // Evaluate test: The incoming timer's time point is queried here that
    // includes the time elapsed before the Charm++ chare has been created + the
    // migration time, so tested with a somewhat lose precision that includes
    // the approximate (guessed) migration time.
    try {
     ensure( "timer different after migrated: ", timer.dsec() > 1.0 );
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
};

//! Test Charm++ migration of a tk::Timer object across the network
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both triggers a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void Timer_object::test< 7 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate tk::Timer 1" );
  
  tk::Timer timer;
  usleep( 1000000 );     // in micro-seconds, sleep for 1.0 second
  CProxy_CharmTimer::ckNew( timer );    // fire up Charm++ chare
}

//! Test querying a timer from a std::map
template<> template<>
void Timer_object::test< 8 >() {
  set_test_name( "query timer from map" );

  std::map< std::string, tk::Timer > timer;
  timer[ "some timer" ];// start timing, assign to label
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  const auto t = tk::ref_find( timer, "some timer" );

  ensure_equals( "timer different", t.dsec(), 1.0, 0.1 );
}

//! Test that querying timer from a map throws with garbage key
template<> template<>
void Timer_object::test< 9 >() {
  set_test_name( "query throws with non-existent key" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::map< std::string, tk::Timer > timer;
    timer[ "some timer" ];// start timing, assign to label
    tk::cref_find( timer, std::string("some non-existent timer") );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
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

#include "NoWarning/charmtimer.def.h"
