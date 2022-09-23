// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/EOS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state class
  \details   This file defines functions for equations of state for the
    compressible flow equations.
*/
// *****************************************************************************
#ifndef EOS_h
#define EOS_h

#include "Data.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = kw::ncomp::info::expect::type;

//! \brief Calculate density from the material pressure and temperature using
//!   the stiffened-gas equation of state
//! \tparam Eq Equation type to operate on, e.g., tag::compflow
//! \param[in] system Equation system index
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the stiffened-gas EOS
template< class Eq >
tk::real eos_density( ncomp_t system, tk::real pr, tk::real temp ) {
  auto g = g_inputdeck.get< tag::param, Eq, tag::gamma >()[ system ][0];
  auto cv = g_inputdeck.get< tag::param, Eq, tag::cv >()[ system ][0];
  return pr / ((g-1.0) * cv * temp);
}

//! \brief Calculate pressure from the material density, momentum and total
//!   energy using the stiffened-gas equation of state
//! \tparam Eq Equation type to operate on, e.g., tag::compflow
//! \param[in] system Equation system index
//! \param[in] rho Material density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] rhoE Material total energy
//! \return Material pressure calculated using the stiffened-gas EOS
template< class Eq >
tk::real eos_pressure( ncomp_t system,
                       tk::real rho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real rhoE )
{
  auto g = g_inputdeck.get< tag::param, Eq, tag::gamma >()[ system ][0];
  return (rhoE - 0.5 * rho * (u*u + v*v + w*w)) * (g-1.0);
}

//! Calculate speed of sound from the material density and material pressure
//! \tparam Eq Equation type to operate on, e.g., tag::compflow
//! \param[in] system Equation system index
//! \param[in] rho Material density
//! \param[in] pr Material pressure
//! \return Material speed of sound using the stiffened-gas EOS
template< class Eq >
tk::real eos_soundspeed( ncomp_t system, tk::real rho, tk::real pr ) {
  auto g = g_inputdeck.get< tag::param, Eq, tag::gamma >()[ system ][0];
  return std::sqrt( g * pr / rho );
}

//! \brief Calculate material specific total energy from the material density,
//!   momentum and material pressure
//! \tparam Eq Equation type to operate on, e.g., tag::compflow
//! \param[in] system Equation system index
//! \param[in] rho Material density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] pr Material pressure
//! \return Material specific total energy using the stiffened-gas EOS
template< class Eq >
tk::real eos_totalenergy( ncomp_t system,
                          tk::real rho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real pr )
{
  auto g = g_inputdeck.get< tag::param, Eq, tag::gamma >()[ system ][0];
  return pr / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w);
}

} //inciter::

#endif // EOS_h
