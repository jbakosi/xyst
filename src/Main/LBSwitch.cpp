// *****************************************************************************
/*!
  \file      src/Main/LBSwitch.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare group for switching on/off load balancing
  \details   Charm++ chare group for switching on/off load balancing.
*/
// *****************************************************************************

#include <iostream>

#include "LBSwitch.hpp"
#include "Print.hpp"

using tk::LBSwitch;

LBSwitch::LBSwitch()
// *****************************************************************************
//  Constructor: turn on automatic load balancing
// *****************************************************************************
{
  TurnManualLBOff();
  if (CkMyPe() == 0) Print() << "Xyst> Load balancing on\n";
}

void
// cppcheck-suppress unusedFunction
LBSwitch::off()
// *****************************************************************************
//  Turn off automatic load balancing
//! \details Since this is a [initproc] routine, the runtime system executes the
//!   routine exactly once on every PE early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  TurnManualLBOn();

  if (CkMyPe() == 0) Print() << "Xyst> Load balancing off\n";
}

#include "NoWarning/lbswitch.def.h"
