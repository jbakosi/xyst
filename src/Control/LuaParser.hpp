// *****************************************************************************
/*!
  \file      src/Control/LuaParser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Lua parser header
*/
// *****************************************************************************
#pragma once

#include <string>

namespace inciter {

// Parse lua input file
std::string
parseLua( const char* inputfile );

} // ::inciter
