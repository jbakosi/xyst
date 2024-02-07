// *****************************************************************************
/*!
  \file      src/Physics/EOS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
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
//! \return Pressure computed from the ideal gas equation of state
inline double
pressure( double re ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return re * (g-1.0);
}

//! Calculate speed of sound from the material density and material pressure
//! \param[in] rho Material density
//! \param[in] pr Material pressure
//! \return Material speed of sound using the stiffened-gas EOS
inline double
soundspeed( double rho, double pr ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return std::sqrt( g * pr / rho );
}

//! \brief Calculate material specific total energy from the material density,
//!   momentum and material pressure
//! \param[in] rho Material density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] pr Material pressure
//! \return Material specific total energy using the stiffened-gas EOS
inline double
totalenergy( double rho, double u, double v, double w, double pr ) {
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  return pr / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w);
}

} // eos::
