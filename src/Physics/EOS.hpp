// *****************************************************************************
/*!
  \file      src/Physics/EOS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state functions
*/
// *****************************************************************************
#pragma once

#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace eos {

using inciter::g_cfg;

//! Compute pressure
//! \param[in] re Specific internal energy times density
//! \return Pressure from the ideal gas equation of state
inline double
pressure( double re ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return re * (g-1.0);
}

//! Compute speed of sound from density and pressure
//! \param[in] r Density
//! \param[in] p Pressure
//! \return Speed of sound from the ideal gas equation of state
inline double
soundspeed( double r, double p ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return std::sqrt( g * p / r );
}

//! Compute specific total energy from density, momentum, and pressure
//! \param[in] r Density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] p Pressure
//! \return Specific total energy from the ideal gas equation of state
inline double
totalenergy( double r, double u, double v, double w, double p ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return p / (g-1.0) + 0.5 * r * (u*u + v*v + w*w);
}

} // eos::
