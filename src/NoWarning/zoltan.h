// *****************************************************************************
/*!
  \file      src/NoWarning/zoltan.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Include zoltan.h with turning off specific compiler warnings
*/
// *****************************************************************************
#pragma once

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-function-type"
#endif

#include <zoltan.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif
