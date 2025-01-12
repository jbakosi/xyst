// *****************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Template Unit Test suite class definition
  \details   Template Unit Test suite class definition. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
// *****************************************************************************

#include <iostream>
#include <utility>
#include <map>

#include "NoWarning/tut_runner.hpp"

#include "TUTSuite.hpp"
#include "TUTTest.hpp"
#include "Print.hpp"

#include "NoWarning/unittest.decl.h"

extern CProxy_Main mainProxy;

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

} // unittest::

using unittest::TUTSuite;

TUTSuite::TUTSuite( const std::string& grp ) :
  m_group2run( grp ),
  m_nrun( 0 ),
  m_ngroup( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 ),
  m_nskip( 0 ),
  m_nwarn( 0 ),
  m_nexcp( 0 ),
  m_nspaw( 0 )
// *****************************************************************************
// Constructor
//! \param[in] grp Test group to run
// *****************************************************************************
{
  const auto& groups = g_runner.get().list_groups();

  // If only select groups to be run, see if there is any that will run
  bool work = false;
  if (grp.empty() ||
      std::any_of( groups.cbegin(), groups.cend(), [&]( const std::string& g ) {
        return g.find(grp) != std::string::npos; } ))
  {
    work = true;
  }

  if (!work) {  // Quit if there is no work to be done
    tk::Print() << "\nNo test groups match '" + grp + "'.\n";
    mainProxy.finalize( false );
  }
}

void
TUTSuite::run()
// *****************************************************************************
//  Run all tests
// *****************************************************************************
{
  tk::Print().unithead( "Running unit tests", m_group2run );

  const auto& groups = g_runner.get().list_groups();

  // Fire up all tests in all groups using the Charm++ runtime system
  for (const auto& g : groups) {
    if (m_group2run.empty()) {  // consider all test groups
      spawngrp( g );
    } else if (g.find(m_group2run) != std::string::npos) {
      spawngrp( g );            // spawn groups configured
    }
  }
}

void
TUTSuite::spawngrp( const std::string& g )
// *****************************************************************************
//  Fire up all tests in a test group
//! \param[in] g Name of the test group
// *****************************************************************************
{
  ++m_ngroup;         // increase number of test groups to run

  // Add up number of additionally-spawned tests (this is so we know how many
  // to expect results from)
  const auto it = m_nspawned.find( g );
  if (it != m_nspawned.end()) m_nspaw += it->second;
  
  // Asynchronously fire up all tests in test group
  for (int t=1; t<=g_maxTestsInGroup; ++t) {
    if (m_fromPE0.count(g)) {
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t, 0 );
    } else {
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
    }
  }
}

void
TUTSuite::evaluateTest( std::vector< std::string >&& status )
// *****************************************************************************
// Evaluate a unit test
//! \param[in] status Vector strings containing the test results. See
//!   unittest::TUTTest constructor for the expected structure of status.
// *****************************************************************************
{
  // Increase number tests run (including dummies)
  ++m_nrun;

  // Evaluate test
  if (status[2] != "8") {             // only care about non-dummy tests
    ++m_ncomplete;                    // count number of tests completed
    if (status[2] == "3")             // count number of tests with a warning
      ++m_nwarn;
    else if (status[2] == "7")        // count number of skipped tests
      ++m_nskip;
    else if (status[2] == "2")        // count number of tests throwing
      ++m_nexcp;
    else if (status[2] != "0")        // count number of failed tests
      ++m_nfail;
  }

  // Echo one-liner info on result of test
  tk::Print().test( m_ncomplete, m_nfail, status );

  // Wait for all tests to finish, then quit
  if (m_nrun == m_ngroup*static_cast<std::size_t>(g_maxTestsInGroup) + m_nspaw)
  {
    auto pass = assess();
    mainProxy.finalize( pass );
  }
}

bool
TUTSuite::assess()
// *****************************************************************************
// Echo final assessment after the full unit test suite has finished
//! \return True of all tests passed, false if there was at least a failure or
//!   an exception
// *****************************************************************************
{
  tk::Print print;

  if (!m_nfail && !m_nwarn && !m_nskip && !m_nexcp) {

    print << "\nAll " + std::to_string(m_ncomplete) + " tests passed\n";

  } else {

    std::string skip, warn, fail, excp;

    if (m_nwarn) {
      warn = "finished with a warning: " + std::to_string(m_nwarn);
    }

    if (m_nskip) {
      skip = std::string(m_nwarn ? ", " : "") +
             "(fully or partially) skipped: " + std::to_string(m_nskip);
    }

    if (m_nexcp) {
      excp = std::string(m_nskip || m_nwarn ? ", " : "") +
             "threw exception: " + std::to_string(m_nexcp);
    }

    if (m_nfail) {
      fail = std::string(m_nexcp || m_nskip || m_nwarn ? ", " : "")
             + "failed: " + std::to_string(m_nfail);
    }

    print << "\nOf " + std::to_string(m_ncomplete) + " tests total: "
          << warn << skip << excp << fail << '\n';

  }

  return (m_nfail || m_nexcp) ? false : true;
}

#include "NoWarning/tutsuite.def.h"
