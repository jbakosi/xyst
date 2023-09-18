// *****************************************************************************
/*!
  \file      src/Physics/Dt.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Compute dt for next time step for edge-based continuous Galerkin
*/
// *****************************************************************************

#include "Dt.hpp"
#include "EOS.hpp"
#include "Vector.hpp"
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
    auto r  = U(p,0,0);
    auto u  = U(p,1,0) / r;
    auto v  = U(p,2,0) / r;
    auto w  = U(p,3,0) / r;
    auto re = U(p,4,0) / r - 0.5*(u*u + v*v + w*w);
    auto vel = tk::length( u, v, w );
    auto pr = eos::pressure( r, re );
    auto c = eos::soundspeed( r, pr );
    auto L = std::cbrt( vol[p] );
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
    auto r  = U(p,0,0);
    auto u  = U(p,1,0) / r;
    auto v  = U(p,2,0) / r;
    auto w  = U(p,3,0) / r;
    auto re = U(p,4,0) / r - 0.5*(u*u + v*v + w*w);
    auto vel = tk::length( u, v, w );
    auto pr = eos::pressure( r, re );
    auto c = eos::soundspeed( r, pr < 0 ? 0 : pr );
    auto L = std::cbrt( vol[p] );
    dtp[p] = L / std::max( vel+c, 1.0e-8 ) * cfl;
  }
}

} // physics::
