// *****************************************************************************
/*!
  \file      src/Physics/Dt.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Compute dt for next time step for edge-based continuous Galerkin
*/
// *****************************************************************************

#include "Dt.hpp"
#include "EOS.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace physics {

using inciter::g_cfg;

tk::real
dt( const std::vector< tk::real >& vol, const tk::Fields& U )
// *****************************************************************************
//  Compute the minimum time step size
//! \param[in] vol Nodal volume (with contributions from other chares)
//! \param[in] U Solution vector at recent time step
//! \return Minimum time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();
  for (std::size_t p=0; p<U.nunk(); ++p) {
    auto r = U(p,0);
    auto u = U(p,1)/r;
    auto v = U(p,2)/r;
    auto w = U(p,3)/r;
    auto pr = eos::pressure( U(p,4) - 0.5*r*(u*u + v*v + w*w) );
    auto c = eos::soundspeed( r, std::max(pr,0.0) );
    auto L = std::cbrt( vol[p] );
    auto vel = std::sqrt( u*u + v*v + w*w );
    auto euler_dt = L / std::max( vel+c, 1.0e-8 );
    mindt = std::min( mindt, euler_dt );
  }

  mindt *= g_cfg.get< tag::cfl >();

  return mindt;
}

void
dt( const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp )
// *****************************************************************************
//  Compute a time step size for each mesh node (toward stationary state)
//! \param[in] vol Nodal volume (with contributions from other chares)
//! \param[in] U Solution vector at recent time step
//! \param[in,out] dtp Time step size for each mesh node
// *****************************************************************************
{
  auto cfl = g_cfg.get< tag::cfl >();

  for (std::size_t p=0; p<U.nunk(); ++p) {
    auto r = U(p,0);
    auto u = U(p,1)/r;
    auto v = U(p,2)/r;
    auto w = U(p,3)/r;
    auto pr = eos::pressure( U(p,4) - 0.5*r*(u*u + v*v + w*w) );
    auto c = eos::soundspeed( r, std::max(pr,0.0) );
    auto L = std::cbrt( vol[p] );
    auto vel = std::sqrt( u*u + v*v + w*w );
    dtp[p] = L / std::max( vel+c, 1.0e-8 ) * cfl;
  }
}

} // physics::
