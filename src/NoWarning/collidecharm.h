// *****************************************************************************
/*!
  \file      src/NoWarning/collidecharm.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Include collidecharm.h with turning off some compiler warnings.
*/
// *****************************************************************************
#pragma once

#include "Compiler.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdollar-in-identifier-extension"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wshadow"
  #pragma clang diagnostic ignored "-Wsuggest-override"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wcast-qual"
  #pragma clang diagnostic ignored "-Wimplicit-int-conversion"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wextra-semi-stmt"
  #pragma clang diagnostic ignored "-Wshadow-field-in-constructor"
  #pragma clang diagnostic ignored "-Wnon-virtual-dtor"
  #pragma clang diagnostic ignored "-Wcovered-switch-default"
  #pragma clang diagnostic ignored "-Wswitch-enum"
  #pragma clang diagnostic ignored "-Wzero-length-array"
  #pragma clang diagnostic ignored "-Wdouble-promotion"
  #pragma clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
  #pragma clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wfloat-conversion"
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
  #pragma clang diagnostic ignored "-Wshadow-field"
  #pragma clang diagnostic ignored "-Wsuggest-destructor-override"
  #pragma clang diagnostic ignored "-Wheader-hygiene"
//#elif defined(STRICT_GNUC)
//  #pragma GCC diagnostic push
//  #pragma GCC diagnostic ignored "-Wunused-parameter"
//  #pragma GCC diagnostic ignored "-Wshadow"
//  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "collidecharm.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
//#elif defined(STRICT_GNUC)
//  #pragma GCC diagnostic pop
#endif
