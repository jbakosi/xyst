// *****************************************************************************
/*!
  \file      src/Main/InciterDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter driver
  \details   Inciter driver.
*/
// *****************************************************************************

#include "InciterDriver.hpp"
#include "Inciter/InputDeck/Parser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"
#include "Print.hpp"

#include "NoWarning/transporter.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const ctr::CmdLine& cmdline )
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
// *****************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

  // Create pretty printer
  tk::Print print;

  print.item( "Non-blocking migration, -" + *kw::nonblocking::alias(),
               cmdline.get< tag::nonblocking >() ? "on" : "off" );
  print.item( "Benchmark mode, -" + *kw::benchmark::alias(),
               cmdline.get< tag::benchmark >() ? "on" : "off" );
  print.item( "On-screen feedback, -" + *kw::feedback::alias(),
               cmdline.get< tag::feedback >() ? "on" : "off" );
  print.item( "Load-balancing frequency, -" + *kw::lbfreq::alias(),
               std::to_string(cmdline.get< tag::lbfreq >()) );
  print.item( "Checkpoint/restart frequency, -" + *kw::rsfreq::alias(),
               std::to_string(cmdline.get< tag::rsfreq >()) );

  // Parse input deck into g_inputdeck
  print.item( "Control file", cmdline.get< tag::io, tag::control >() );
  g_inputdeck = g_inputdeck_defaults;   // overwrite with defaults if restarted
  InputDeckParser inputdeckParser( cmdline, g_inputdeck );

  // Output command line object to file
  auto logfilename = tk::inciter_executable() + "_input.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), "inputdeck", g_inputdeck );
}

void
InciterDriver::execute() const
// *****************************************************************************
//  Run inciter
// *****************************************************************************
{
  // Instantiate Transporter chare on PE 0 which drives time-integration
  CProxy_Transporter::ckNew( 0 );
}
