// *****************************************************************************
/*!
  \file      src/NoWarning/exodusII.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Include exodusII.h with turning off specific compiler warnings
*/
// *****************************************************************************
#pragma once

#include "Compiler.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wdollar-in-identifier-extension"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
  #pragma GCC diagnostic ignored "-Wcast-function-type"
#endif

#include "exodus/include/exodusII.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif
