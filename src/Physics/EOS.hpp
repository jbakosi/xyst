// *****************************************************************************
/*!
  \file      src/Physics/EOS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state functions
*/
// *****************************************************************************
#pragma once

#include "InciterInputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace eos {

using inciter::g_inputdeck;

//! Compute pressure
//! \param[in] r Density
//! \param[in] e Specific internal energy
//! \return Pressure computed from the ideal gas equation of state
inline tk::real
pressure( tk::real r, tk::real e ) {
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
  return r * e * (g-1.0);
}

//! Calculate speed of sound from the material density and material pressure
//! \param[in] rho Material density
//! \param[in] pr Material pressure
//! \return Material speed of sound using the stiffened-gas EOS
inline tk::real
soundspeed( tk::real rho, tk::real pr ) {
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
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
inline tk::real
totalenergy( tk::real rho, tk::real u, tk::real v, tk::real w, tk::real pr ) {
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
  return pr / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w);
}

} // eos::
