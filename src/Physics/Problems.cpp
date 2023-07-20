// *****************************************************************************
/*!
  \file      src/Physics/Problems.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem-specific functions. Initial conditions, source terms.
*/
// *****************************************************************************

#include "Problems.hpp"
#include "EOS.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace problems {

using inciter::g_inputdeck;

// *****************************************************************************
//  Determine nodes that lie inside user-defined IC box(es)
//! \param[in] coord Mesh node coordinates
//! \return inbox List of nodes at which box user ICs are set for each IC box
// *****************************************************************************
std::vector< std::unordered_set< std::size_t > >
boxnodes( const std::array< std::vector< tk::real >, 3 >& coord ) {

  const auto& icbox = g_inputdeck.get< tag::ic >();

  if (icbox.empty()) return {};

  std::vector< std::unordered_set< std::size_t > > inbox;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  std::size_t bcnt = 0;
  for (const auto& b : icbox) {
    inbox.emplace_back();
    const auto& bx = b.get< tag::x >();
    const auto& by = b.get< tag::y >();
    const auto& bz = b.get< tag::z >();
    std::vector< tk::real > box{ bx[0], bx[1], by[0], by[1], bz[0], bz[1] };

    const auto eps = std::numeric_limits< tk::real >::epsilon();
    // Determine which nodes lie in the IC box
    if ( std::any_of( begin(box), end(box), [=](auto p)
                      { return abs(p) > eps; } ) )
    {
      std::array< tk::real, 3 > b_min{{box[0], box[2], box[4]}};
      std::array< tk::real, 3 > b_max{{box[1], box[3], box[5]}};
      for (std::size_t i=0; i<x.size(); ++i) {
        std::array< tk::real, 3 > node{{ x[i], y[i], z[i] }};
        if ( node[0]>b_min[0] && node[0]<b_max[0] &&
             node[1]>b_min[1] && node[1]<b_max[1] &&
             node[2]>b_min[2] && node[2]<b_max[2] )
        {
          inbox[bcnt].insert( i );
        }
      }
    }
    ++bcnt;
  }

  return inbox;
}

namespace userdef {

static std::vector< tk::real >
ic( tk::real, tk::real, tk::real, tk::real )
// *****************************************************************************
//! Set homogeneous initial conditions for a generic user-defined problem
//! \return Values of conserved variables
// *****************************************************************************
{
  auto ic_density = g_inputdeck.get< tag::ic_density >();
  const auto& ic_velocity = g_inputdeck.get< tag::ic_velocity >();
  ErrChk( ic_velocity.size() == 3, "ic_velocity must have 3 components" );

  std::vector< tk::real > u( 5, 0.0 );

  u[0] = ic_density;
  u[1] = u[0] * ic_velocity[0];
  u[2] = u[0] * ic_velocity[1];
  u[3] = u[0] * ic_velocity[2];

  auto ic_pressure = g_inputdeck.get< tag::ic_pressure >();
  auto ic_energy = g_inputdeck.get< tag::ic_energy >();
  auto ic_temperature = g_inputdeck.get< tag::ic_temperature >();

  auto largereal = std::numeric_limits< double >::max();

  if (std::abs(ic_pressure - largereal) > 1.0e-12) {

    u[4] = eos::totalenergy( u[0], u[1]/u[0], u[2]/u[0], u[3]/u[0],
                             ic_pressure );

  } else if (std::abs(ic_energy - largereal) > 1.0e-12) {

    u[4] = u[0] * ic_energy;

  } else if (std::abs(ic_temperature - largereal) > 1.0e-12) {

    auto cv = g_inputdeck.get< tag::mat_spec_heat_const_vol >();
    if (std::abs(cv - largereal) > 1.0e-12) {
      u[4] = u[0] * ic_temperature * cv;
    }

  } else {

    Throw( "IC background energy cannot be computed. Must specify "
           "one of background pressure, energy, or velocity." );

  }

  return u;
}

} // userdef::

namespace nonlinear_energy_growth {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Set initial conditions prescribing nonlinear energy growth
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using std::cos;

  // manufactured solution parameters
  auto ce = g_inputdeck.get< tag::problem_ce >();
  auto r0 = g_inputdeck.get< tag::problem_r0 >();
  auto a = g_inputdeck.get< tag::problem_alpha >();
  auto k = g_inputdeck.get< tag::problem_kappa >();
  const auto& b = g_inputdeck.get< tag::problem_beta >();

  auto ec = [ ce, t ]( tk::real kappa, tk::real h, tk::real p ) {
    return std::pow( -3.0*(ce + kappa*h*h*t), p );
  };

  auto hx = [ x, y, z, b ]() {
    return cos(b[0]*M_PI*x) * cos(b[1]*M_PI*y) * cos(b[2]*M_PI*z);
  };

  // density
  auto r = r0 + std::exp(-a*t) * (1.0 - x*x - y*y - z*z);
  // energy
  auto re = r * ec(k,hx(),-1.0/3.0);

  return { r, 0.0, 0.0, 0.0, re };
}

static std::vector< tk::real >
src( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Compute and return source term for nonlinear energy growth
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \return Source for flow variables + transported scalars
// *****************************************************************************
{
  using std::sin; using std::cos; using std::pow;

  // manufactured solution parameters
  auto a = g_inputdeck.get< tag::problem_alpha >();
  const auto& b = g_inputdeck.get< tag::problem_beta >();
  auto ce = g_inputdeck.get< tag::problem_ce >();
  auto kappa = g_inputdeck.get< tag::problem_kappa >();
  auto r0 = g_inputdeck.get< tag::problem_r0 >();
  // ratio of specific heats
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
  // spatial component of density field
  auto gx = 1.0 - x*x - y*y - z*z;
  // derivative of spatial component of density field
  std::array< tk::real, 3 > dg{ -2.0*x, -2.0*y, -2.0*z };
  // spatial component of energy field
  auto h = cos(b[0]*M_PI*x) * cos(b[1]*M_PI*y) * cos(b[2]*M_PI*z);
  // derivative of spatial component of energy field
  std::array< tk::real, 3 >
    dh{ -b[0]*M_PI*sin(b[0]*M_PI*x)*cos(b[1]*M_PI*y)*cos(b[2]*M_PI*z),
        -b[1]*M_PI*cos(b[0]*M_PI*x)*sin(b[1]*M_PI*y)*cos(b[2]*M_PI*z),
        -b[2]*M_PI*cos(b[0]*M_PI*x)*cos(b[1]*M_PI*y)*sin(b[2]*M_PI*z) };
  // temporal function f and its derivative
  auto ft = std::exp(-a*t);
  auto dfdt = -a*ft;
  // density and its derivatives
  auto rho = r0 + ft*gx;
  std::array< tk::real, 3 > drdx{ ft*dg[0], ft*dg[1], ft*dg[2] };
  auto drdt = gx*dfdt;
  // internal energy and its derivatives
  auto ie = pow( -3.0*(ce + kappa*h*h*t), -1.0/3.0 );
  std::array< tk::real, 3 > dedx{ 2.0 * pow(ie,4.0) * kappa * h * dh[0] * t,
                                  2.0 * pow(ie,4.0) * kappa * h * dh[1] * t,
                                  2.0 * pow(ie,4.0) * kappa * h * dh[2] * t };
  const auto dedt = kappa * h * h * pow(ie,4.0);

  std::vector< tk::real > s( 5, 0.0 );
  // density source
  s[0] = drdt;
  // momentum source
  s[1] = (g-1.0)*(rho*dedx[0] + ie*drdx[0]);
  s[2] = (g-1.0)*(rho*dedx[1] + ie*drdx[1]);
  s[3] = (g-1.0)*(rho*dedx[2] + ie*drdx[2]);
  // energy source
  s[4] = rho*dedt + ie*drdt;

  return s;
}

} // nonlinear_energy_growth::

namespace rayleigh_taylor {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Set initial conditions prescribing a Rayleigh-Taylor flow
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< tag::problem_alpha >();
  const auto& b = g_inputdeck.get< tag::problem_beta >();
  auto p0 = g_inputdeck.get< tag::problem_p0 >();
  auto r0 = g_inputdeck.get< tag::problem_r0 >();
  auto k = g_inputdeck.get< tag::problem_kappa >();

  // spatial component of density and pressure fields
  tk::real gx = b[0]*x*x + b[1]*y*y + b[2]*z*z;
  // density
  tk::real r = r0 - gx;
  // velocity
  tk::real ft = cos(k*M_PI*t);
  tk::real u = ft * z * sin(M_PI*x);
  tk::real v = ft * z * cos(M_PI*y);
  tk::real w = ft * ( -0.5*M_PI*z*z*(cos(M_PI*x) - sin(M_PI*y)) );
  // total specific energy
  tk::real rE = eos::totalenergy( r, u, v, w, p0 + a*gx );

  return { r, r*u, r*v, r*w, rE };
}

static std::vector< tk::real >
src( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Compute and return source term for a Rayleigh-Taylor flow
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \return Source for flow variables + transported scalars
// *****************************************************************************
{
  using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< tag::problem_alpha >();
  const auto& b = g_inputdeck.get< tag::problem_beta >();
  auto k = g_inputdeck.get< tag::problem_kappa >();
  auto p0 = g_inputdeck.get< tag::problem_p0 >();
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();

  // evaluate solution at x,y,z,t
  auto U = ic( x, y, z, t );

  // density, velocity, energy, pressure
  auto rho = U[0];
  auto u = U[1]/U[0];
  auto v = U[2]/U[0];
  auto w = U[3]/U[0];
  auto E = U[4]/U[0];
  auto p = p0 + a*(b[0]*x*x + b[1]*y*y + b[2]*z*z);

  // spatial gradients
  std::array< tk::real, 3 > drdx{{ -2.0*b[0]*x, -2.0*b[1]*y, -2.0*b[2]*z }};
  std::array< tk::real, 3 > dpdx{{ 2.0*a*b[0]*x, 2.0*a*b[1]*y, 2.0*a*b[2]*z }};
  tk::real ft = cos(k*M_PI*t);
  std::array< tk::real, 3 > dudx{{ ft*M_PI*z*cos(M_PI*x),
                                   0.0,
                                   ft*sin(M_PI*x) }};
  std::array< tk::real, 3 > dvdx{{ 0.0,
                                   -ft*M_PI*z*sin(M_PI*y),
                                   ft*cos(M_PI*y) }};
  std::array< tk::real, 3 > dwdx{{ ft*M_PI*0.5*M_PI*z*z*sin(M_PI*x),
                                   ft*M_PI*0.5*M_PI*z*z*cos(M_PI*y),
                                  -ft*M_PI*z*(cos(M_PI*x) - sin(M_PI*y)) }};
  std::array< tk::real, 3 > dedx{{
    dpdx[0]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[0]
    + u*dudx[0] + v*dvdx[0] + w*dwdx[0],
    dpdx[1]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[1]
    + u*dudx[1] + v*dvdx[1] + w*dwdx[1],
    dpdx[2]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[2]
    + u*dudx[2] + v*dvdx[2] + w*dwdx[2] }};

  // time derivatives
  auto dudt = -k*M_PI*sin(k*M_PI*t)*z*sin(M_PI*x);
  auto dvdt = -k*M_PI*sin(k*M_PI*t)*z*cos(M_PI*y);
  auto dwdt =  k*M_PI*sin(k*M_PI*t)/2*M_PI*z*z*(cos(M_PI*x) - sin(M_PI*y));
  auto dedt = u*dudt + v*dvdt + w*dwdt;

  std::vector< tk::real > s( 5, 0.0 );
  // density source
  s[0] = u*drdx[0] + v*drdx[1] + w*drdx[2];
  // momentum source
  s[1] = rho*dudt+u*s[0]+dpdx[0] + U[1]*dudx[0]+U[2]*dudx[1]+U[3]*dudx[2];
  s[2] = rho*dvdt+v*s[0]+dpdx[1] + U[1]*dvdx[0]+U[2]*dvdx[1]+U[3]*dvdx[2];
  s[3] = rho*dwdt+w*s[0]+dpdx[2] + U[1]*dwdx[0]+U[2]*dwdx[1]+U[3]*dwdx[2];
  // energy source
  s[4] = rho*dedt + E*s[0] + U[1]*dedx[0]+U[2]*dedx[1]+U[3]*dedx[2]
       + u*dpdx[0]+v*dpdx[1]+w*dpdx[2];

  return s;
}

} // rayleigh_taylor::

namespace sedov {
static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real )
// *****************************************************************************
//! Set initial conditions prescribing the Sedov blast wave
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using std::abs;

  // pressure
  auto eps = std::numeric_limits< tk::real >::epsilon();
  tk::real p;
  if (abs(x) < eps && abs(y) < eps && abs(z) < eps) {
    p = g_inputdeck.get< tag::problem_p0 >();
  } else {
    p = 0.67e-4;
  }

  // density
  tk::real r = 1.0;
  // velocity
  tk::real u = 0.0;
  tk::real v = 0.0;
  tk::real w = 0.0;
  // total specific energy
  tk::real rE = eos::totalenergy( r, u, v, w, p );

  return { r, r*u, r*v, r*w, rE };

}
} // sedov::

namespace sod {
static std::vector< tk::real >
ic( tk::real x, tk::real, tk::real, tk::real )
// *****************************************************************************
//! Set initial conditions prescribing the Sod shocktube
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  tk::real r, p, u, v, w, rE;

  if (x<0.5) {
    // density
    r = 1.0;
    // pressure
    p = 1.0;
  }
  else {
    // density
    r = 0.125;
    // pressure
    p = 0.1;
  }

  // velocity
  u = 0.0;
  v = 0.0;
  w = 0.0;

  // total specific energy
  rE = eos::totalenergy( r, u, v, w, p );

  return { r, r*u, r*v, r*w, rE };
}
} // sod::

namespace taylor_green {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real, tk::real )
// *****************************************************************************
//! Set initial conditions prescribing the Taylor-Green vortex
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  // density
  tk::real r = 1.0;
  // pressure
  tk::real p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
  // velocity
  tk::real u =  sin(M_PI*x) * cos(M_PI*y);
  tk::real v = -cos(M_PI*x) * sin(M_PI*y);
  tk::real w = 0.0;
  // total specific energy
  auto rE = eos::totalenergy( r, u, v, w, p );

  return { r, r*u, r*v, r*w, rE };
}

static std::vector< tk::real >
src( tk::real x, tk::real y, tk::real, tk::real )
// *****************************************************************************
//! Compute and return source term for a the Taylor-Green vortex
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \return Source for flow variables + transported scalars
// *****************************************************************************
{
  using std::cos;

  std::vector< tk::real > s( 5, 0.0 );
  s[4] = 3.0*M_PI/8.0*( cos(3.0*M_PI*x)*cos(M_PI*y)
                      - cos(3.0*M_PI*y)*cos(M_PI*x) );

  return s;
}

} // taylor_green::

namespace vortical_flow {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real )
// *****************************************************************************
//! Set initial conditions prescribing vortical flow
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  // manufactured solution parameters
  tk::real a = g_inputdeck.get< tag::problem_alpha >();
  tk::real k = g_inputdeck.get< tag::problem_kappa >();
  tk::real p0 = g_inputdeck.get< tag::problem_p0 >();
  // ratio of specific heats
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
  // velocity
  tk::real ru = a*x - k*y;
  tk::real rv = k*x + a*y;
  tk::real rw = -2.0*a*z;
  // total specific energy
  tk::real rE = (ru*ru + rv*rv + rw*rw)/2.0 + (p0 - 2.0*a*a*z*z) / (g - 1.0);

  return { 1.0, ru, rv, rw, rE };
}

static std::vector< tk::real >
src( tk::real x, tk::real y, tk::real z, tk::real )
// *****************************************************************************
//! Compute and return source term for vortical flow
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \return Source for flow variables + transported scalars
// *****************************************************************************
{
  // manufactured solution parameters
  auto a = g_inputdeck.get< tag::problem_alpha >();
  auto k = g_inputdeck.get< tag::problem_kappa >();
  // ratio of specific heats
  auto g = g_inputdeck.get< tag::mat_spec_heat_ratio >();
  // evaluate solution at x,y,z
  auto u = ic( x, y, z, 0.0 );

  std::vector< tk::real > s( 5, 0.0 );
  // momentum source
  s[1] = a*u[1]/u[0] - k*u[2]/u[0];
  s[2] = k*u[1]/u[0] + a*u[2]/u[0];
  // energy source
  s[4] = (s[1]*u[1] + s[2]*u[2])/u[0] + 8.0*a*a*a*z*z/(g-1.0);

  return s;
}

} // vortical_flow::

namespace slot_cyl {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real, tk::real t )
// *****************************************************************************
//! Set initial conditions prescribing slotted cylinder, cone, Gauss hump
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using std::sin; using std::cos; using std::sqrt;

  // manufactured solution parameters
  tk::real p0 = 1.0;

  std::vector< tk::real > u( 6, 0.0 );

  // prescribed velocity: rotate in x-y plane
  u[0] = 1.0;
  u[1] = u[0] * (0.5 - y);
  u[2] = u[0] * (x - 0.5);
  u[3] = 0.0;
  u[4] = eos::totalenergy( u[0], u[1]/u[0], u[2]/u[0], u[3]/u[0], p0 );

  const tk::real R0 = 0.15;

  // center of the cone
  tk::real x0 = 0.5;
  tk::real y0 = 0.25;
  tk::real r = sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
  tk::real kx = 0.5 + r*sin( t );
  tk::real ky = 0.5 - r*cos( t );

  // center of the hump
  x0 = 0.25;
  y0 = 0.5;
  r = sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
  tk::real hx = 0.5 + r*sin( t-M_PI/2.0 ),
           hy = 0.5 - r*cos( t-M_PI/2.0 );

  // center of the slotted cylinder
  x0 = 0.5;
  y0 = 0.75;
  r = sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
  tk::real cx = 0.5 + r*sin( t+M_PI ),
           cy = 0.5 - r*cos( t+M_PI );

  // end points of the cylinder slot
  tk::real i1x = 0.525, i1y = cy - r*cos( std::asin(0.025/r) ),
           i2x = 0.525, i2y = 0.8,
           i3x = 0.475, i3y = 0.8;

  // rotate end points of cylinder slot
  tk::real ri1x = 0.5 + cos(t)*(i1x-0.5) - sin(t)*(i1y-0.5),
           ri1y = 0.5 + sin(t)*(i1x-0.5) + cos(t)*(i1y-0.5),
           ri2x = 0.5 + cos(t)*(i2x-0.5) - sin(t)*(i2y-0.5),
           ri2y = 0.5 + sin(t)*(i2x-0.5) + cos(t)*(i2y-0.5),
           ri3x = 0.5 + cos(t)*(i3x-0.5) - sin(t)*(i3y-0.5),
           ri3y = 0.5 + sin(t)*(i3x-0.5) + cos(t)*(i3y-0.5);

  // direction of slot sides
  tk::real v1x = ri2x-ri1x, v1y = ri2y-ri1y,
           v2x = ri3x-ri2x, v2y = ri3y-ri2y;

  // lengths of direction of slot sides vectors
  tk::real v1 = sqrt(v1x*v1x + v1y*v1y),
           v2 = sqrt(v2x*v2x + v2y*v2y);

  // cone
  r = sqrt((x-kx)*(x-kx) + (y-ky)*(y-ky)) / R0;
  if (r<1.0) u[5] = 0.6*(1.0-r);

  // hump
  r = sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy)) / R0;
  if (r<1.0) u[5] = 0.2*(1.0+cos(M_PI*std::min(r,1.0)));

  // cylinder
  r = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) / R0;
  const std::array< tk::real, 2 > r1{{ v1x, v1y }},
                                  r2{{ x-ri1x, y-ri1y }};
  const auto d1 = (r1[0]*r2[1] - r2[0]*r1[1]) / v1;
  const std::array< tk::real, 2 > r3{{ v2x, v2y }},
                                  r4{{ x-ri2x, y-ri2y }};
  const auto d2 = (r3[0]*r4[1] - r4[0]*r3[1]) / v2;
  if (r<1.0 && (d1>0.05 || d1<0.0 || d2<0.0)) u[5] = 0.6;

  return u;
}

static std::vector< tk::real >
src( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Compute and return source term for slotted cylinder, cone, Gauss hump
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \return Source for flow variables + transported scalars
// *****************************************************************************
{
  // evaluate solution at x,y,z,t
  auto u = ic( x, y, z, t );

  std::vector< tk::real > s( 6, 0.0 );
  // momentum source
  s[1] = -u[2];
  s[2] =  u[1];

  return s;
}

} // slot_cyl::

namespace point_source {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Set initial conditions for point source problem
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \return Values of conserved variables
// *****************************************************************************
{
  auto u = userdef::ic( x, y, z, t );
  u.push_back( 0.0 );
  return u;
}

static void
src( const std::array< std::vector< tk::real >, 3 >& coord,
     tk::real t,
     tk::Fields& U )
// *****************************************************************************
//! Apply point-source directly to numerical solution
//! \param[in] coord Mesh node coordinates
//! \param[in] t Physical time
//! \param[in,out] U Solution vector at recent time step
//! \note This is different from other source terms, because this directly
//!   modifies the solution instead of applied as a source term mathematically.
//!   Hence the function signature is also different.
// *****************************************************************************
{
  if (U.nprop() == 5) return;

  const auto& source = g_inputdeck.get< tag::problem_source >();
  const auto& location = source.get< tag::location >();
  auto radius = source.get< tag::radius >();
  auto release_time = source.get< tag::release_time >();
  auto largereal = std::numeric_limits< double >::max();

  if (location.size() != 3 ||
      std::abs(radius - largereal) < 1.0e-12 ||
      std::abs(release_time - largereal) < 1.0e-12)
  {
    return;
  }

  auto sx = location[0];
  auto sy = location[1];
  auto sz = location[2];
  auto sr = radius;
  auto st = release_time;

  if (t < st) return;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t i=0; i<U.nunk(); ++i) {
    auto rx = sx - x[i];
    auto ry = sy - y[i];
    auto rz = sz - z[i];
    if (rx*rx + ry*ry + rz*rz < sr*sr) U(i,5,0) = 1.0;
  }
}

} // point_source::

std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
IC()
// *****************************************************************************
//  Query user config and assign function to set initial conditions
//! \return The function to call to set initial conditions
// *****************************************************************************
{
  const auto& problem = inciter::g_inputdeck.get< tag::problem >();

  std::function< std::vector< tk::real >
               ( tk::real, tk::real, tk::real, tk::real ) > ic;

  if (problem == "userdef")
    ic = userdef::ic;
  else if (problem == "nonlinear_energy_growth")
    ic = nonlinear_energy_growth::ic;
  else if (problem == "rayleigh_taylor")
    ic = rayleigh_taylor::ic;
  else if (problem == "sedov")
    ic = sedov::ic;
  else if (problem == "sod")
    ic = sod::ic;
  else if (problem == "taylor_green")
    ic = taylor_green::ic;
  else if (problem == "vortical_flow")
    ic = vortical_flow::ic;
  else if (problem == "slot_cyl")
    ic = slot_cyl::ic;
  else if (problem == "point_src")
    ic = point_source::ic;
  else
    Throw( "problem type ic not hooked up" );

  return ic;
}

std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
SOL()
// *****************************************************************************
//  Query user config and assign function to query analytic solutions
//! \return The function to call to query analytic solutions
// *****************************************************************************
{
  const auto& problem = inciter::g_inputdeck.get< tag::problem >();

  if (problem == "userdef" ||
      problem == "sod" ||
      problem == "sedov" ||
      problem == "point_src")
    return {};
  else
    return IC();
}

void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t,
            const std::vector< std::unordered_set< std::size_t > >& inbox )
// *****************************************************************************
//  Set inital conditions
//! \param[in] coord Mesh node coordinates
//! \param[in,out] U Array of unknowns
//! \param[in] t Physical time
//! \param[in] inbox Nodes at which box user ICs are set (for each box IC)
// *****************************************************************************
{
  Assert( coord[0].size() == U.nunk(), "Size mismatch" );

  const auto& icbox = g_inputdeck.get< tag::ic >();

  auto ic = IC();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  auto largereal = std::numeric_limits< double >::max();

  // Set initial conditions dependeing on problem configured
  for (std::size_t i=0; i<x.size(); ++i) {

    // Set background ICs
    auto s = ic( x[i], y[i], z[i], t );

    // Initialize user-defined ICs in boxes
    for (std::size_t j=0; j<icbox.size(); ++j) {
      const auto& b = icbox[j];
      if (inbox.size() > j && inbox[j].find(i) != inbox[j].end()) {
        auto boxr = b.get< tag::ic_density >();
        const auto& boxv = b.get< tag::ic_velocity >();
        auto boxp = b.get< tag::ic_pressure >();
        auto boxe = b.get< tag::ic_energy >();
        auto boxt = b.get< tag::ic_temperature >();

        tk::real r = 0.0, ru = 0.0, rv = 0.0, rw = 0.0, re = 0.0;
        if (std::abs(boxr - largereal) > 1.0e-12 && boxr > 0.0) {
          r = boxr;
        }
        if (std::abs(boxv[0] - largereal) > 1.0e-12 &&
            std::abs(boxv[1] - largereal) > 1.0e-12 &&
            std::abs(boxv[2] - largereal) > 1.0e-12)
        {
          ru = r * boxv[0];
          rv = r * boxv[1];
          rw = r * boxv[2];
        }
        if (std::abs(boxp- largereal) > 1.0e-12 && boxp> 0.0) {
          re = eos::totalenergy( r, ru/r, rv/r, rw/r, boxp);
        }
        if (std::abs(boxe - largereal) > 1.0e-12 && boxe > 0.0) {
          auto ux = ru/r, uy = rv/r, uz = rw/r;
          auto ke = 0.5*(ux*ux + uy*uy + uz*uz);
          re = r * (boxe + ke);
        }
        if (std::abs(boxt - largereal) > 1.0e-12 && boxt > 0.0) {
          auto cv = g_inputdeck.get< tag::mat_spec_heat_const_vol >();
          if (std::abs(cv - largereal) > 1.0e-12 && cv > 0.0) {
            re = r * boxt * cv;
          }
        }
        s[0] = r;
        s[1] = ru;
        s[2] = rv;
        s[3] = rw;
        s[4] = re;
      }
    }

    // Set values for ICs
    for (std::size_t c=0; c<s.size(); ++c) U(i,c,0) = s[c];

  }
}

std::function< std::vector< tk::real >
                 ( tk::real, tk::real, tk::real, tk::real ) >
SRC()
// *****************************************************************************
//  Query user config and assign function to add a source term
//! \return The function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
  const auto& problem = inciter::g_inputdeck.get< tag::problem >();

  std::function<
    std::vector< tk::real >( tk::real, tk::real, tk::real, tk::real ) > src;

  if (problem == "nonlinear_energy_growth")
    src = nonlinear_energy_growth::src;
  else if (problem == "rayleigh_taylor")
    src = rayleigh_taylor::src;
  else if (problem == "taylor_green")
    src = taylor_green::src;
  else if (problem == "vortical_flow")
    src = vortical_flow::src;
  else if (problem == "slot_cyl")
    src = slot_cyl::src;

  return src;
}

std::function< void( const std::array< std::vector< tk::real >, 3 >&,
                     tk::real,
                     tk::Fields& ) >
PHYS_SRC()
// *****************************************************************************
//  Query user config and assign function to apply source to numerical solution
//! \return The function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
  const auto& problem = inciter::g_inputdeck.get< tag::problem >();

  std::function< void( const std::array< std::vector< tk::real >, 3 >&,
                       tk::real,
                       tk::Fields& ) > src;

  if (problem == "point_src") {
    src = point_source::src;
  }

  return src;
}

} // problems::
