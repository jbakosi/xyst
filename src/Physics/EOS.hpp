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

#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace physics {

//! \brief Calculate density from the material pressure and temperature using
//!   the stiffened-gas equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the stiffened-gas EOS
inline tk::real
eos_density( tk::real pr, tk::real temp ) {
  using inciter::g_inputdeck;
  auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
  auto cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0][0];
  return pr / ((g-1.0) * cv * temp);
}

//! Compute pressure
//! \param[in] r Density
//! \param[in] e Specific internal energy
//! \return Pressure computed from the ideal gas equation of state
inline tk::real
eos_pressure( tk::real r, tk::real e ) {
  using inciter::g_inputdeck;
  auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
  return r * e * (g-1.0);
}

//! Compute pressure
//! \param[in] r Density
//! \param[in] re Specific total energy
//! \return Pressure computed from the ideal gas equation of state
inline tk::real
eos_pressure( tk::real r, tk::real u, tk::real v, tk::real w, tk::real re ) {
  using inciter::g_inputdeck;
  auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
  return (re - 0.5 * r * (u*u + v*v + w*w)) * (g-1.0);
}

//! Calculate speed of sound from the material density and material pressure
//! \param[in] rho Material density
//! \param[in] pr Material pressure
//! \return Material speed of sound using the stiffened-gas EOS
inline tk::real
eos_soundspeed( tk::real rho, tk::real pr ) {
  using inciter::g_inputdeck;
  auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
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
eos_totalenergy( tk::real rho, tk::real u, tk::real v, tk::real w, tk::real pr )
{
  using inciter::g_inputdeck;
  auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0][0];
  return pr / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w);
}

} // physics::
