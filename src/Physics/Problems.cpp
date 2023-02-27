// *****************************************************************************
/*!
  \file      src/Physics/Problems.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \bried     Problem-specific functions. Initial conditions, source terms.
*/
// *****************************************************************************

#include "Problems.hpp"
#include "EOS.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace problems {

namespace userdef {

static std::vector< tk::real >
ic( tk::real, tk::real, tk::real, tk::real )
// *****************************************************************************
//! Set initial conditions for a generic user-defined problem
//! \return Values of conserved variables
// *****************************************************************************
{
  using inciter::g_inputdeck;

  const auto& ic = g_inputdeck.get< tag::param, tag::compflow, tag::ic >();

  const auto& bgrhoic = ic.get< tag::density >();
  const auto& bgvelic = ic.get< tag::velocity >();

  ErrChk( bgrhoic.size() == 1 && bgrhoic.back().size() == 1,
          "Need a background density as IC" );
  ErrChk( bgvelic.size() == 1 && bgvelic.back().size() == 3,
          "Need three velocity components as background IC" );

  std::vector< tk::real > u( 5, 0.0 );

  u[0] = bgrhoic[0][0];
  u[1] = u[0] * bgvelic[0][0];
  u[2] = u[0] * bgvelic[0][1];
  u[3] = u[0] * bgvelic[0][2];

  const auto& bgpreic = ic.get< tag::pressure >();
  const auto& bgenic = ic.get< tag::energy >();
  const auto& bgtempic = ic.get< tag::temperature >();

  if (bgpreic.size() > 0 && !bgpreic[0].empty()) {

    u[4] = eos::totalenergy( u[0], u[1]/u[0], u[2]/u[0], u[3]/u[0],
                             bgpreic[0][0] );

  } else if (bgenic.size() > 0 && !bgenic[0].empty()) {

    u[4] = u[0] * bgenic[0][0];

  } else if (bgtempic.size() > 0 && !bgtempic[0].empty()) {

    auto c_v = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0][0];
    u[4] = u[0] * bgtempic[0][0] * c_v;

  } else {

    Throw( "IC background energy cannot be computed. User must specify "
           "one of background pressure, energy, or velocity." );

  }

  ErrChk( u[0] > 0.0, "Density IC zero" );

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
//! \param[in] t Time where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow; using tag::param; using std::cos;

  // manufactured solution parameters
  auto ce = g_inputdeck.get< param, compflow, tag::ce >()[0];
  auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[0];
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  auto k = g_inputdeck.get< param, compflow, tag::kappa >()[0];
  auto bx = g_inputdeck.get< param, compflow, tag::betax >()[0];
  auto by = g_inputdeck.get< param, compflow, tag::betay >()[0];
  auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[0];

  auto ec = [ ce, t ]( tk::real kappa, tk::real h, tk::real p ) {
    return std::pow( -3.0*(ce + kappa*h*h*t), p );
  };

  auto hx = [ x, y, z, bx, by, bz ]() {
    return cos(bx*M_PI*x) * cos(by*M_PI*y) * cos(bz*M_PI*z);
  };

  // density
  auto r = r0 + std::exp(-a*t) * (1.0 - x*x - y*y - z*z);
  // energy
  auto re = r * ec(k,hx(),-1.0/3.0);

  return { r, 0.0, 0.0, 0.0, re };
}

static void
src( tk::real x, tk::real y, tk::real z, tk::real t,
     tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re,
     tk::real& )
// *****************************************************************************
//! Compute and return source term for nonlinear energy growth
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \param[in,out] r Density source
//! \param[in,out] ru X momentum source
//! \param[in,out] rv Y momentum source
//! \param[in,out] rw Z momentum source
//! \param[in,out] re Specific total energy source
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow;
  using std::sin; using std::cos; using std::pow;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  auto bx = g_inputdeck.get< param, compflow, tag::betax >()[0];
  auto by = g_inputdeck.get< param, compflow, tag::betay >()[0];
  auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[0];
  auto ce = g_inputdeck.get< param, compflow, tag::ce >()[0];
  auto kappa = g_inputdeck.get< param, compflow, tag::kappa >()[0];
  auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[0];
  // ratio of specific heats
  auto g = g_inputdeck.get< param, compflow, tag::gamma >()[0][0];
  // spatial component of density field
  auto gx = 1.0 - x*x - y*y - z*z;
  // derivative of spatial component of density field
  std::array< tk::real, 3 > dg{ -2.0*x, -2.0*y, -2.0*z };
  // spatial component of energy field
  auto h = cos(bx*M_PI*x) * cos(by*M_PI*y) * cos(bz*M_PI*z);
  // derivative of spatial component of energy field
  std::array< tk::real, 3 >
    dh{ -bx*M_PI*sin(bx*M_PI*x)*cos(by*M_PI*y)*cos(bz*M_PI*z),
        -by*M_PI*cos(bx*M_PI*x)*sin(by*M_PI*y)*cos(bz*M_PI*z),
        -bz*M_PI*cos(bx*M_PI*x)*cos(by*M_PI*y)*sin(bz*M_PI*z) };
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
  // density source
  r = drdt;
  // momentum source
  ru = (g-1.0)*(rho*dedx[0] + ie*drdx[0]);
  rv = (g-1.0)*(rho*dedx[1] + ie*drdx[1]);
  rw = (g-1.0)*(rho*dedx[2] + ie*drdx[2]);
  // energy source
  re = rho*dedt + ie*drdt;
}

} // nonlinear_energy_growth::

namespace rayleigh_taylor {

static std::vector< tk::real >
ic( tk::real x, tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//! Set initial conditions prescribing a Rayleigh-Taylor flow
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of conserved variables
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow; using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  auto bx = g_inputdeck.get< param, compflow, tag::betax >()[0];
  auto by = g_inputdeck.get< param, compflow, tag::betay >()[0];
  auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[0];
  auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[0];
  auto r0 = g_inputdeck.get< param, compflow, tag::r0 >()[0];
  auto k = g_inputdeck.get< param, compflow, tag::kappa >()[0];

  // spatial component of density and pressure fields
  tk::real gx = bx*x*x + by*y*y + bz*z*z;
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

static void
src( tk::real x, tk::real y, tk::real z, tk::real t,
     tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re,
     tk::real& )
// *****************************************************************************
//! Compute and return source term for a Rayleigh-Taylor flow
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \param[in,out] r Density source
//! \param[in,out] ru X momentum source
//! \param[in,out] rv Y momentum source
//! \param[in,out] rw Z momentum source
//! \param[in,out] re Specific total energy source
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow; using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  auto bx = g_inputdeck.get< param, compflow, tag::betax >()[0];
  auto by = g_inputdeck.get< param, compflow, tag::betay >()[0];
  auto bz = g_inputdeck.get< param, compflow, tag::betaz >()[0];
  auto k = g_inputdeck.get< param, compflow, tag::kappa >()[0];
  auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[0];
  // ratio of specific heats
  auto g = g_inputdeck.get< param, compflow, tag::gamma >()[0][0];

  // evaluate solution at x,y,z,t
  auto s = ic( x, y, z, t );

  // density, velocity, energy, pressure
  auto rho = s[0];
  auto u = s[1]/s[0];
  auto v = s[2]/s[0];
  auto w = s[3]/s[0];
  auto E = s[4]/s[0];
  auto p = p0 + a*(bx*x*x + by*y*y + bz*z*z);

  // spatial gradients
  std::array< tk::real, 3 > drdx{{ -2.0*bx*x, -2.0*by*y, -2.0*bz*z }};
  std::array< tk::real, 3 > dpdx{{ 2.0*a*bx*x, 2.0*a*by*y, 2.0*a*bz*z }};
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

  // density source
  r = u*drdx[0] + v*drdx[1] + w*drdx[2];
  // momentum source
  ru = rho*dudt+u*r+dpdx[0] + s[1]*dudx[0]+s[2]*dudx[1]+s[3]*dudx[2];
  rv = rho*dvdt+v*r+dpdx[1] + s[1]*dvdx[0]+s[2]*dvdx[1]+s[3]*dvdx[2];
  rw = rho*dwdt+w*r+dpdx[2] + s[1]*dwdx[0]+s[2]*dwdx[1]+s[3]*dwdx[2];
  // energy source
  re = rho*dedt + E*r + s[1]*dedx[0]+s[2]*dedx[1]+s[3]*dedx[2]
       + u*dpdx[0]+v*dpdx[1]+w*dpdx[2];
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
  using inciter::g_inputdeck;
  using std::abs;

  // pressure
  auto eps = std::numeric_limits< tk::real >::epsilon();
  tk::real p;
  if (abs(x) < eps && abs(y) < eps && abs(z) < eps)
    p = g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[ 0 ];
  else
    p = 0.67e-4;

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

static void
src( tk::real x, tk::real y, tk::real, tk::real,
     tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re,
     tk::real& )
// *****************************************************************************
//! Compute and return source term for a the Taylor-Green vortex
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in,out] r Density source
//! \param[in,out] ru X momentum source
//! \param[in,out] rv Y momentum source
//! \param[in,out] rw Z momentum source
//! \param[in,out] re Specific total energy source
// *****************************************************************************
{
  using std::cos;

  r = ru = rv = rw = 0.0;
  re = 3.0*M_PI/8.0*( cos(3.0*M_PI*x)*cos(M_PI*y)
                    - cos(3.0*M_PI*y)*cos(M_PI*x) );
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
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow;

  // manufactured solution parameters
  tk::real a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  tk::real b = g_inputdeck.get< param, compflow, tag::beta >()[0];
  tk::real p0 = g_inputdeck.get< param, compflow, tag::p0 >()[0];
  // ratio of specific heats
  tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[0][0];
  // velocity
  tk::real ru = a*x - b*y;
  tk::real rv = b*x + a*y;
  tk::real rw = -2.0*a*z;
  // total specific energy
  tk::real rE = (ru*ru + rv*rv + rw*rw)/2.0 + (p0 - 2.0*a*a*z*z) / (g - 1.0);

  return { 1.0, ru, rv, rw, rE };
}

static void
src( tk::real x, tk::real y, tk::real z, tk::real,
     tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re,
     tk::real& )
// *****************************************************************************
//! Compute and return source term for vortical flow
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in,out] r Density source
//! \param[in,out] ru X momentum source
//! \param[in,out] rv Y momentum source
//! \param[in,out] rw Z momentum source
//! \param[in,out] re Specific total energy source
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::param; using tag::compflow;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[0];
  auto b = g_inputdeck.get< param, compflow, tag::beta >()[0];
  // ratio of specific heats
  auto g = g_inputdeck.get< param, compflow, tag::gamma >()[0][0];
  // evaluate solution at x,y,z
  auto s = ic( x, y, z, 0.0 );

  // density source
  r = 0.0;
  // momentum source
  ru = a*s[1]/s[0] - b*s[2]/s[0];
  rv = b*s[1]/s[0] + a*s[2]/s[0];
  rw = 0.0;
  // energy source
  re = (ru*s[1] + rv*s[2])/s[0] + 8.0*a*a*a*z*z/(g-1.0);
}

} // vortical_flow::

namespace point_source {

static void
src( tk::real x, tk::real y, tk::real z, tk::real t,
     tk::real&, tk::real&, tk::real&, tk::real&, tk::real&, tk::real& sc )
// *****************************************************************************
//! Compute and return source term for vortical flow
//! \param[in] x X coordinate where to evaluate the source
//! \param[in] y Y coordinate where to evaluate the source
//! \param[in] z Z coordinate where to evaluate the source
//! \param[in] t Time where to evaluate the source
//! \param[in,out] sc Scalar source
// *****************************************************************************
{
  using inciter::g_inputdeck;
  using tag::transport;

std::cout << "!";
  const auto& ns = g_inputdeck.get< tag::param, transport, tag::source >();
  if (ns.empty()) return;

  const auto& nc = g_inputdeck.get< tag::component >().get< transport >();
  if (nc.empty()) return;

  const auto& src = ns[0];

  ErrChk( src.size() % 5 == 0, "Source must contain n x 5 values" );

  auto sx = src[0];
  auto sy = src[1];
  auto sz = src[2];
  auto sr = src[3];
  auto st = src[4];

  if (t < st) return;

  if (((sx-x)*(sx-x) + (sy-y)*(sy-y) + (sz-z)*(sz-z)) < sr*sr) sc = 1.0;
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
  using inciter::g_inputdeck;
  using ProblemType = inciter::ctr::ProblemType;

  auto problem = g_inputdeck.get< tag::problem >();

  std::function< std::vector< tk::real >
               ( tk::real, tk::real, tk::real, tk::real ) > ic;

  if (problem == ProblemType::USER_DEFINED)
    ic = userdef::ic;
  else if (problem == ProblemType::NONLINEAR_ENERGY_GROWTH)
    ic = nonlinear_energy_growth::ic;
  else if (problem == ProblemType::RAYLEIGH_TAYLOR)
    ic = rayleigh_taylor::ic;
  else if (problem == ProblemType::SEDOV)
    ic = sedov::ic;
  else if (problem == ProblemType::SOD)
    ic = sod::ic;
  else if (problem == ProblemType::TAYLOR_GREEN)
    ic = taylor_green::ic;
  else if (problem == ProblemType::VORTICAL_FLOW)
    ic = vortical_flow::ic;
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
  using inciter::g_inputdeck;
  using ProblemType = inciter::ctr::ProblemType;

  auto problem = g_inputdeck.get< tag::problem >();

  if (problem == ProblemType::USER_DEFINED ||
      problem == ProblemType::SOD)
    return {};
  else
    return IC();
}

void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t )
// *****************************************************************************
//  Initalize the compressible flow equations, prepare for time integration
//! \param[in] coord Mesh node coordinates
//! \param[in,out] U Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  Assert( coord[0].size() == U.nunk(), "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto ic = IC();

  // Set initial conditions dependeing on problem configured
  for (std::size_t i=0; i<x.size(); ++i) {
    auto s = ic( x[i], y[i], z[i], t );
    U(i,0,0) = s[0]; // rho
    U(i,1,0) = s[1]; // rho * u
    U(i,2,0) = s[2]; // rho * v
    U(i,3,0) = s[3]; // rho * w
    U(i,4,0) = s[4]; // rho * e, e: total = kinetic + internal
  }
}

std::function< void( tk::real, tk::real, tk::real, tk::real,
  tk::real&, tk::real&, tk::real&, tk::real&, tk::real&, tk::real& ) >
SRC()
// *****************************************************************************
//  Query user config and assign function to add a source term
//! \return The function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
  using inciter::g_inputdeck;
  auto problem = g_inputdeck.get< tag::problem >();

  std::function<
    void( tk::real, tk::real, tk::real, tk::real,
      tk::real&, tk::real&, tk::real&, tk::real&, tk::real&, tk::real& ) > src;

  using ProblemType = inciter::ctr::ProblemType;

  if (problem == ProblemType::NONLINEAR_ENERGY_GROWTH)
    src = nonlinear_energy_growth::src;
  else if (problem == ProblemType::RAYLEIGH_TAYLOR)
    src = rayleigh_taylor::src;
  else if (problem == ProblemType::TAYLOR_GREEN)
    src = taylor_green::src;
  else if (problem == ProblemType::VORTICAL_FLOW)
    src = vortical_flow::src;
  else if (problem == ProblemType::POINT_SOURCE)
    src = point_source::src;

  return src;
}

} // problems::
