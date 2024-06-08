// *****************************************************************************
/*!
  \file      src/Physics/Lax.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LaxCG: Time-derivative preconditioning for all Ma
  \see       Luo, Baum, Lohner, "Extension of Harten-Lax-van Leer Scheme for
             Flows at All Speeds", AIAA Journal, Vol. 43, No. 6, 2005
  \see       Weiss & Smith, "Preconditioning Applied to Variable and Constant
             Density Time-Accurate Flows on Unstructured Meshes", AIAA Journal,
             Vol. 33, No. 11, 1995, pp. 2050-2057.
*/
// *****************************************************************************

#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Lax.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace lax {

static const tk::real muscl_eps = 1.0e-9;
static const tk::real muscl_const = 1.0/3.0;

using inciter::g_cfg;

static void
muscl( std::size_t p,
       std::size_t q,
       const tk::UnsMesh::Coords& coord,
       const tk::Fields& G,
       tk::real& pL, tk::real& uL, tk::real& vL, tk::real& wL, tk::real& tL,
       tk::real& pR, tk::real& uR, tk::real& vR, tk::real& wR, tk::real& tR )
// *****************************************************************************
//! Compute MUSCL reconstruction in edge-end points for the flow variables
//! \param[in] p Left node id of edge-end
//! \param[in] q Right node id of edge-end
//! \param[in] coord Array of nodal coordinates
//! \param[in] G Gradient of all unknowns in mesh points
//! \param[in,out] pL Left pressure
//! \param[in,out] uL Left X velocity
//! \param[in,out] vL Left Y velocity
//! \param[in,out] wL Left Z velocity
//! \param[in,out] tL Left temperature
//! \param[in,out] pR Right pressure
//! \param[in,out] uR Right X velocity
//! \param[in,out] vR Right Y velocity
//! \param[in,out] wR Right Z velocity
//! \param[in,out] tR Right temperature
// *****************************************************************************
{
  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  tk::real vw[3] = { x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[5], delta2[5], delta3[5];
  tk::real ls[5] = { pL, uL, vL, wL, tL },
           rs[5] = { pR, uR, vR, wR, tR },
           url[5], urr[5];
  memcpy( url, ls, sizeof ls );
  memcpy( urr, rs, sizeof rs );

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<5; ++c) {

    auto g1 = G(p,c*3+0)*vw[0] + G(p,c*3+1)*vw[1] + G(p,c*3+2)*vw[2];
    auto g2 = G(q,c*3+0)*vw[0] + G(q,c*3+1)*vw[1] + G(q,c*3+2)*vw[2];

    delta2[c] = rs[c] - ls[c];
    delta1[c] = 2.0 * g1 - delta2[c];
    delta3[c] = 2.0 * g2 - delta2[c];

    // MUSCL extrapolation option 1:
    // ---------------------------------------------------------------------
    // See Waltz, J., Morgan, N. R., Canfield, T. R., Charest, M. R., Risinger,
    // L. D., & Wohlbier, J. G. (2014). A three-dimensional finite element
    // arbitrary Lagrangian–Eulerian method for shock hydrodynamics on
    // unstructured grids. Computers & Fluids, 92, 172-187.

    // van Leer limiter
    auto rcL = (delta2[c] + muscl_eps) / (delta1[c] + muscl_eps);
    auto rcR = (delta2[c] + muscl_eps) / (delta3[c] + muscl_eps);
    auto rLinv = (delta1[c] + muscl_eps) / (delta2[c] + muscl_eps);
    auto rRinv = (delta3[c] + muscl_eps) / (delta2[c] + muscl_eps);
    auto phiL = (std::abs(rcL) + rcL) / (std::abs(rcL) + 1.0);
    auto phiR = (std::abs(rcR) + rcR) / (std::abs(rcR) + 1.0);
    auto phi_L_inv = (std::abs(rLinv) + rLinv) / (std::abs(rLinv) + 1.0);
    auto phi_R_inv = (std::abs(rRinv) + rRinv) / (std::abs(rRinv) + 1.0);
    // update unknowns with reconstructed unknowns
    url[c] += 0.25*(delta1[c]*(1.0-muscl_const)*phiL +
                    delta2[c]*(1.0+muscl_const)*phi_L_inv);
    urr[c] -= 0.25*(delta3[c]*(1.0-muscl_const)*phiR +
                    delta2[c]*(1.0+muscl_const)*phi_R_inv);

    // MUSCL extrapolation option 2:
    // ---------------------------------------------------------------------
    // See Luo, H., Baum, J. D., & Lohner, R. (1994). Edge-based finite element
    // scheme for the Euler equations. AIAA journal, 32(6), 1183-1190.
    // van Leer, B. (1974). Towards the ultimate conservative difference
    // scheme. II. Monotonicity and conservation combined in a second-order
    // scheme. Journal of computational physics, 14(4), 361-370.
    // Derived from the flux limiter phi as: s = phi_inv - (1 - phi)

    // van Albada limiter
    //auto sL = std::max(0.0, (2.0*delta1[c]*delta2[c] + muscl_eps)
    //  /(delta1[c]*delta1[c] + delta2[c]*delta2[c] + muscl_eps));
    //auto sR = std::max(0.0, (2.0*delta3[c]*delta2[c] + muscl_eps)
    //  /(delta3[c]*delta3[c] + delta2[c]*delta2[c] + muscl_eps));
    //// update unknowns with reconstructed unknowns
    //url[c] += 0.25*sL*(delta1[c]*(1.0 - muscl_const*sL)
    //                 + delta2[c]*(1.0 + muscl_const*sL));
    //urr[c] -= 0.25*sR*(delta3[c]*(1.0 - muscl_const*sR)
    //                 + delta2[c]*(1.0 + muscl_const*sR));
  }

  // force first order if the reconstructions for pressure or temperature
  // would have allowed negative values
  if (ls[0] < delta1[0] || ls[4] < delta1[4]) memcpy( url, ls, sizeof ls );
  if (rs[0] < -delta3[0] || rs[4] < -delta3[4]) memcpy( urr, rs, sizeof rs );

  pL = url[0];
  uL = url[1];
  vL = url[2];
  wL = url[3];
  tL = url[4];

  pR = urr[0];
  uR = urr[1];
  vR = urr[2];
  wR = urr[3];
  tR = urr[4];
}

static void
muscl( std::size_t p, std::size_t q, const tk::UnsMesh::Coords& coord,
       const tk::Fields& G, tk::real uL[], tk::real uR[] )
// *****************************************************************************
//! Compute MUSCL reconstruction in edge-end points for transported scalars
//! \param[in] p Left node id of edge-end
//! \param[in] q Right node id of edge-end
//! \param[in] coord Array of nodal coordinates
//! \param[in] G Gradient of all unknowns in mesh points
//! \param[in,out] uL Primitive variables at left edge-end point
//! \param[in,out] uR Primitive variables at right edge-end point
// *****************************************************************************
{
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  auto ns = G.nprop() / 3 - 5;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  tk::real vw[3] = { x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[ns], delta2[ns], delta3[ns];

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<ns; ++c) {
    auto g = (5+c)*3;
    auto g1 = G(p,g+0)*vw[0] + G(p,g+1)*vw[1] + G(p,g+2)*vw[2];
    auto g2 = G(q,g+0)*vw[0] + G(q,g+1)*vw[1] + G(q,g+2)*vw[2];

    delta2[c] = uR[5+c] - uL[5+c];
    delta1[c] = 2.0 * g1 - delta2[c];
    delta3[c] = 2.0 * g2 - delta2[c];

    // van Leer limiter
    auto rL = (delta2[c] + muscl_eps) / (delta1[c] + muscl_eps);
    auto rR = (delta2[c] + muscl_eps) / (delta3[c] + muscl_eps);
    auto rLinv = (delta1[c] + muscl_eps) / (delta2[c] + muscl_eps);
    auto rRinv = (delta3[c] + muscl_eps) / (delta2[c] + muscl_eps);
    auto phiL = (std::abs(rL) + rL) / (std::abs(rL) + 1.0);
    auto phiR = (std::abs(rR) + rR) / (std::abs(rR) + 1.0);
    auto phi_L_inv = (std::abs(rLinv) + rLinv) / (std::abs(rLinv) + 1.0);
    auto phi_R_inv = (std::abs(rRinv) + rRinv) / (std::abs(rRinv) + 1.0);
    // update unknowns with reconstructed unknowns
    uL[5+c] += 0.25*(delta1[c]*(1.0-muscl_const)*phiL +
                     delta2[c]*(1.0+muscl_const)*phi_L_inv);
    uR[5+c] -= 0.25*(delta3[c]*(1.0-muscl_const)*phiR +
                     delta2[c]*(1.0+muscl_const)*phi_R_inv);
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
grad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const tk::Fields& U,
      tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of primitive variables in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Solution vector of primitive variables at recent time step
//! \param[in,out] G Nodal gradients
//! \return Gradients of primitive variables in all mesh points
// *****************************************************************************
{
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // cppcheck-suppress unreadVariable
  auto ncomp = U.nprop();

  Assert( G.nunk() == U.nunk(), "Size mismatch" );
  Assert( G.nprop() == ncomp*3, "Size mismatch" );
  G.fill( 0.0 );

  // domain integral

  // domain edge contributions: tetrahedron superedges
  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
    const auto N = dsupedge[0].data() + e*4;
    for (std::size_t c=0; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c), U(N[2],c), U(N[3],c) };
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[6];
        const auto d = dsupint[0].data();
        f[0] = d[(e*6+0)*3+j] * (u[1] + u[0]);
        f[1] = d[(e*6+1)*3+j] * (u[2] + u[1]);
        f[2] = d[(e*6+2)*3+j] * (u[0] + u[2]);
        f[3] = d[(e*6+3)*3+j] * (u[3] + u[0]);
        f[4] = d[(e*6+4)*3+j] * (u[3] + u[1]);
        f[5] = d[(e*6+5)*3+j] * (u[3] + u[2]);
        G(N[0],c*3+j) = G(N[0],c*3+j) - f[0] + f[2] - f[3];
        G(N[1],c*3+j) = G(N[1],c*3+j) + f[0] - f[1] - f[4];
        G(N[2],c*3+j) = G(N[2],c*3+j) + f[1] - f[2] - f[5];
        G(N[3],c*3+j) = G(N[3],c*3+j) + f[3] + f[4] + f[5];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    for (std::size_t c=0; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c), U(N[2],c) };
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[3];
        const auto d = dsupint[1].data();
        f[0] = d[(e*3+0)*3+j] * (u[1] + u[0]);
        f[1] = d[(e*3+1)*3+j] * (u[2] + u[1]);
        f[2] = d[(e*3+2)*3+j] * (u[0] + u[2]);
        G(N[0],c*3+j) = G(N[0],c*3+j) - f[0] + f[2];
        G(N[1],c*3+j) = G(N[1],c*3+j) + f[0] - f[1];
        G(N[2],c*3+j) = G(N[2],c*3+j) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*3;
    for (std::size_t c=0; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c) };
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = d[j] * (u[1] + u[0]);
        G(N[0],c*3+j) -= f;
        G(N[1],c*3+j) += f;
      }
    }
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto n = tk::cross( ba, ca );
    n[0] /= 12.0;
    n[1] /= 12.0;
    n[2] /= 12.0;
    for (std::size_t c=0; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c), U(N[2],c) };
      auto uab = (u[0] + u[1])/4.0;
      auto ubc = (u[1] + u[2])/4.0;
      auto uca = (u[2] + u[0])/4.0;
      tk::real g[] = { uab + uca + u[0],
                       uab + ubc + u[1],
                       ubc + uca + u[2] };
      for (std::size_t j=0; j<3; ++j) {
        G(N[0],c*3+j) += g[j] * n[j];
        G(N[1],c*3+j) += g[j] * n[j];
        G(N[2],c*3+j) += g[j] * n[j];
      }
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

tk::real
refvel( tk::real r, tk::real p, tk::real v )
// *****************************************************************************
//  Compute reference velocitity of the preconditioned system
//! \param[in] r Density
//! \param[in] p Pressure
//! \param[in] v Velocity magnitude
//! \return Reference velocity
// *****************************************************************************
{
  auto K = g_cfg.get< tag::turkel >();
  auto velinf = g_cfg.get< tag::velinf >();
  auto vinf = tk::length( velinf[0], velinf[1], velinf[2] );

  return std::min( eos::soundspeed( r, p ), std::max( v, K*vinf ) );
  //return eos::soundspeed( r, p );
}

static std::tuple< tk::real, tk::real >
sigvel( tk::real p, tk::real T, tk::real v, tk::real vn )
// *****************************************************************************
//  Compute signal velocitities of the preconditioned system
//! \param[in] p Pressure
//! \param[in] T Temperature
//! \param[in] v Velocity magnitude
//! \param[in] vn Face-normal velocity
//! \return v', c'
// *****************************************************************************
{
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  auto rgas = g_cfg.get< tag::mat_spec_gas_const >();
  auto cp = g*rgas/(g-1.0);

  auto r = p/T/rgas;
  auto rp = r/p;
  auto rt = -r/T;
  auto vr = refvel( r, p, v );
  auto vr2 = vr*vr;
  auto beta = rp + rt/r/cp;
  auto alpha = 0.5*(1.0 - beta*vr2);
  auto vpri = vn*(1.0 - alpha);
  auto cpri = std::sqrt( alpha*alpha*vn*vn + vr2 );

  return { vpri, cpri };
}

static void
rusanov( const tk::UnsMesh::Coords& coord,
         const tk::Fields& G,
         const tk::real dsupint[],
         std::size_t p,
         std::size_t q,
         const tk::real L[],
         const tk::real R[],
         tk::real f[] )
// *****************************************************************************
//! Compute advection fluxes on a single edge with Rusanov's flux
//! \param[in] coord Mesh node coordinates
//! \param[in] G Nodal gradients
//! \param[in] dsupint Domain superedge integral for this edge
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] L Left physics state variables
//! \param[in,out] R Rigth physics state variables
//! \param[in,out] f Flux computed
// *****************************************************************************
{
  auto ncomp = G.nprop() / 3;

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // will work on copies of physics variables
  tk::real l[ncomp], r[ncomp];
  memcpy( l, L, sizeof l );
  memcpy( r, R, sizeof r );

  // MUSCL reconstruction in edge-end points for flow variables
  muscl( p, q, coord, G, l[0], l[1], l[2], l[3], l[4],
                         r[0], r[1], r[2], r[3], r[4] );

  // dualface-normal velocities
  auto nx = dsupint[0];
  auto ny = dsupint[1];
  auto nz = dsupint[2];
  auto vnL = l[1]*nx + l[2]*ny + l[3]*nz;
  auto vnR = r[1]*nx + r[2]*ny + r[3]*nz;

  // pressure
  auto pL = l[0];
  auto pR = r[0];

  // preconditioned wave speed
  auto len = tk::length( nx, ny, nz );
  auto [vpL, cpL] = sigvel( l[0], l[4], tk::length(l[1],l[2],l[3]), vnL );
  auto [vpR, cpR] = sigvel( r[0], r[4], tk::length(r[1],r[2],r[3]), vnR );

  // convert to conserved variables
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  auto rgas = g_cfg.get< tag::mat_spec_gas_const >();
  l[0] = pL/l[4]/rgas;
  l[1] *= l[0];
  l[2] *= l[0];
  l[3] *= l[0];
  l[4] = pL/(g-1.0) + 0.5*(l[1]*l[1] + l[2]*l[2] + l[3]*l[3])/l[0];
  r[0] = pR/r[4]/rgas;
  r[1] *= r[0];
  r[2] *= r[0];
  r[3] *= r[0];
  r[4] = pR/(g-1.0) + 0.5*(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])/r[0];

  // dissipation: option 1
  //auto sL = std::abs(vpL) + cpL;
  //auto sR = std::abs(vpR) + cpR;
  //auto fw = std::max( sL, sR ) * len;

  // dissipation: option 2
  using std::abs;
  using std::max;
  auto sp = max(abs(vpL-cpL),max(abs(vpR-cpR),max(abs(vpL+cpL),abs(vpR+cpR))));
  auto sL = -sp;
  auto sR = +sp;
  auto fw = std::max( sL, sR ) * len;

  // flow fluxes
  f[0] = l[0]*vnL + r[0]*vnR + fw*(r[0] - l[0]);
  f[1] = l[1]*vnL + r[1]*vnR + (pL + pR)*nx + fw*(r[1] - l[1]);
  f[2] = l[2]*vnL + r[2]*vnR + (pL + pR)*ny + fw*(r[2] - l[2]);
  f[3] = l[3]*vnL + r[3]*vnR + (pL + pR)*nz + fw*(r[3] - l[3]);
  f[4] = (l[4] + pL)*vnL + (r[4] + pR)*vnR + fw*(r[4] - l[4]);

  // artificial viscosity
  const auto stab2 = g_cfg.get< tag::stab2 >();
  if (stab2) {
    auto stab2coef = g_cfg.get< tag::stab2coef >();
    auto fws = stab2coef * fw;
    f[0] -= fws*(l[0] - r[0]);
    f[1] -= fws*(l[1] - r[1]);
    f[2] -= fws*(l[2] - r[2]);
    f[3] -= fws*(l[3] - r[3]);
    f[4] -= fws*(l[4] - r[4]);
  }

  if (ncomp == 5) return;

  // MUSCL reconstruction in edge-end points for scalars
  muscl( p, q, coord, G, l, r );

  // scalar dissipation
  auto sw = std::max( std::abs(vnL), std::abs(vnR) );

  // scalar fluxes
  for (std::size_t c=5; c<ncomp; ++c) {
    f[c] = l[c]*vnL + r[c]*vnR + sw*(r[c] - l[c]);
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
hllc( const tk::UnsMesh::Coords& coord,
      const tk::Fields& G,
      const tk::real dsupint[],
      std::size_t p,
      std::size_t q,
      const tk::real L[],
      const tk::real R[],
      tk::real f[] )
// *****************************************************************************
//! Compute advection fluxes on a single edge with Harten-Lax-vanLeer-Contact
//! \param[in] coord Mesh node coordinates
//! \param[in] G Nodal gradients
//! \param[in] dsupint Domain superedge integral for this edge
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] L Left physics state variables
//! \param[in,out] R Rigth physics state variables
//! \param[in,out] f Flux computed
//! \see Toro, Riemann Solver and Numerical Methods for Fluid Dynamics. 3rd
//!      Edition, Springer, New York, 2009
//! \see Davis, Simplified Second-Order Godunov-Type Methods, SIAM J. Sci. Stat.
//!      Comput. 9:445-473, 1988
//! \see Roe, Approximate Riemann Solvers, Parameter Vectors, and Difference
//!      Schemes, J. Comput Phys., 43:357-372, 1981
// *****************************************************************************
{
  auto ncomp = G.nprop() / 3;

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // will work on copies of physics variables
  tk::real l[ncomp], r[ncomp];
  memcpy( l, L, sizeof l );
  memcpy( r, R, sizeof r );

  // MUSCL reconstruction in edge-end points for flow variables
  muscl( p, q, coord, G, l[0], l[1], l[2], l[3], l[4],
                         r[0], r[1], r[2], r[3], r[4] );

  // dualface-normal velocities
  auto nx = -dsupint[0];
  auto ny = -dsupint[1];
  auto nz = -dsupint[2];
  auto len = tk::length( nx, ny, nz );
  nx /= len;
  ny /= len;
  nz /= len;
  auto qL = l[1]*nx + l[2]*ny + l[3]*nz;
  auto qR = r[1]*nx + r[2]*ny + r[3]*nz;

  // pressure
  auto pL = l[0];
  auto pR = r[0];

  // preconditioned signal velocities
  auto [vpL, cpL] = sigvel( l[0], l[4], tk::length(l[1],l[2],l[3]), qL*len );
  auto [vpR, cpR] = sigvel( r[0], r[4], tk::length(r[1],r[2],r[3]), qR*len );

  // convert to conserved variables
  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  auto rgas = g_cfg.get< tag::mat_spec_gas_const >();
  l[0] = pL/l[4]/rgas;
  l[1] *= l[0];
  l[2] *= l[0];
  l[3] *= l[0];
  l[4] = pL/(g-1.0) + 0.5*(l[1]*l[1] + l[2]*l[2] + l[3]*l[3])/l[0];
  r[0] = pR/r[4]/rgas;
  r[1] *= r[0];
  r[2] *= r[0];
  r[3] *= r[0];
  r[4] = pR/(g-1.0) + 0.5*(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])/r[0];

  // preconditioned left and right wave speeds

  // option 1
  //auto sL = std::min( vpL - cpL, vpR - cpR );
  //auto sR = std::max( vpL + cpL, vpR + cpR );

  // option 2
  //auto sL = vpL - cpL;
  //auto sR = vpR + cpR;

  // option 3
  using std::abs;
  using std::max;
  auto sp = max(abs(vpL-cpL),max(abs(vpR-cpR),max(abs(vpL+cpL),abs(vpR+cpR))));
  auto sL = -sp;
  auto sR = +sp;

  //using std::sqrt;
  //auto hl = (l[4]/l[0] + pL) / l[0];
  //auto hr = (r[4]/r[0] + pR) / r[0];
  //auto srl = sqrt( l[0] );
  //auto srr = sqrt( r[0] );
  //auto hh = (srl*hl + srr*hr) / (srl + srr);
  //auto uh = (srl*vpL + srr*vpR) / (srl + srr);
  //auto ch = sqrt( (g-1.0)*(hh - 0.5*(uh*uh)) );
  //// option 4
  //auto sL = std::min( vpL - cpL, uh - ch );
  //auto sR = std::max( vpR + cpR, uh + ch );
  // option 5
  //auto sL = uh - ch;
  //auto sR = uh + ch;

  // no preconditioning
  //auto cL = eos::soundspeed(l[0],pL);
  //auto cR = eos::soundspeed(r[0],pR);
  //auto sL = fmin( qL - cL, qR - cR );
  //auto sR = fmax( qL + cL, qR + cR );

  // contact wave speed and pressure
  auto tL = sL - qL;
  auto tR = sR - qR;
  auto sM = (r[0]*qR*tR - l[0]*qL*tL + pL - pR) / (r[0]*tR - l[0]*tL);
  auto pS = pL - l[0]*tL*(qL - sM);

  // intermediate left-, and right-state conserved unknowns
  tk::real uL[ncomp], uR[ncomp];
  auto s = sL - sM;
  uL[0] = tL*l[0]/s;
  uL[1] = (tL*l[1] + (pS-pL)*nx)/s;
  uL[2] = (tL*l[2] + (pS-pL)*ny)/s;
  uL[3] = (tL*l[3] + (pS-pL)*nz)/s;
  uL[4] = (tL*l[4] - pL*qL + pS*sM)/s;
  s = sR - sM;
  uR[0] = tR*r[0]/s;
  uR[1] = (tR*r[1] + (pS-pR)*nx)/s;
  uR[2] = (tR*r[2] + (pS-pR)*ny)/s;
  uR[3] = (tR*r[3] + (pS-pR)*nz)/s;
  uR[4] = (tR*r[4] - pR*qR + pS*sM)/s;

  auto L2 = -2.0*len;
  nx *= L2;
  ny *= L2;
  nz *= L2;

  // flow fluxes
  if (sL > 0.0) {
    auto qL2 = qL * L2;
    f[0] = l[0]*qL2;
    f[1] = l[1]*qL2 + pL*nx;
    f[2] = l[2]*qL2 + pL*ny;
    f[3] = l[3]*qL2 + pL*nz;
    f[4] = (l[4] + pL)*qL2;
  }
  else if (sL <= 0.0 && sM > 0.0) {
    auto qL2 = qL * L2;
    auto sL2 = sL * L2;
    f[0] = l[0]*qL2 + sL2*(uL[0] - l[0]);
    f[1] = l[1]*qL2 + pL*nx + sL2*(uL[1] - l[1]);
    f[2] = l[2]*qL2 + pL*ny + sL2*(uL[2] - l[2]);
    f[3] = l[3]*qL2 + pL*nz + sL2*(uL[3] - l[3]);
    f[4] = (l[4] + pL)*qL2 + sL2*(uL[4] - l[4]);
  }
  else if (sM <= 0.0 && sR >= 0.0) {
    auto qR2 = qR * L2;
    auto sR2 = sR * L2;
    f[0] = r[0]*qR2 + sR2*(uR[0] - r[0]);
    f[1] = r[1]*qR2 + pR*nx + sR2*(uR[1] - r[1]);
    f[2] = r[2]*qR2 + pR*ny + sR2*(uR[2] - r[2]);
    f[3] = r[3]*qR2 + pR*nz + sR2*(uR[3] - r[3]);
    f[4] = (r[4] + pR)*qR2 + sR2*(uR[4] - r[4]);
  }
  else {
    auto qR2 = qR * L2;
    f[0] = r[0]*qR2;
    f[1] = r[1]*qR2 + pR*nx;
    f[2] = r[2]*qR2 + pR*ny;
    f[3] = r[3]*qR2 + pR*nz;
    f[4] = (r[4] + pR)*qR2;
  }

  // artificial viscosity
  //const auto stab2 = g_cfg.get< tag::stab2 >();
  //if (stab2) {
  //  auto stab2coef = g_cfg.get< tag::stab2coef >();
  //  auto sl = std::abs(vpL) + cpL;
  //  auto sr = std::abs(vpR) + cpR;
  //  auto fws = stab2coef * std::max(sl,sr) * len;
  //  f[0] -= fws*(l[0] - r[0]);
  //  f[1] -= fws*(l[1] - r[1]);
  //  f[2] -= fws*(l[2] - r[2]);
  //  f[3] -= fws*(l[3] - r[3]);
  //  f[4] -= fws*(l[4] - r[4]);
  //}

  if (ncomp == 5) return;

  // MUSCL reconstruction in edge-end points for scalars
  muscl( p, q, coord, G, l, r );

  // scalar fluxes
  //auto sw = std::max( std::abs(vpL), std::abs(vpR) ) * len;
  //for (std::size_t c=5; c<ncomp; ++c) {
  //  f[c] = (l[c]*qL + r[c]*qR)*len + sw*(r[c] - l[c]);
  //}

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static
std::function< void( const tk::UnsMesh::Coords&,
                     const tk::Fields&,
                     const tk::real[],
                     std::size_t,
                     std::size_t,
                     const tk::real[],
                     const tk::real[],
                     tk::real[] ) >
FLUX()
// *****************************************************************************
//! Configure flux function to use
//! \return Function to call to compute flux between two edge-end points
// *****************************************************************************
{
  const auto& flux = g_cfg.get< tag::flux >();

  if (flux == "rusanov")
    return rusanov;
  else if (flux == "hllc")
    return hllc;
  else
    Throw( "Flux not configured" );
}

static void
advdom( const tk::UnsMesh::Coords& coord,
        const std::array< std::vector< std::size_t >, 3 >& dsupedge,
        const std::array< std::vector< tk::real >, 3 >& dsupint,
        const tk::Fields& G,
        const tk::Fields& U,
        // cppcheck-suppress constParameter
        tk::Fields& R )
// *****************************************************************************
//! Compute domain integral for advection
//! \param[in] coord Mesh node coordinates
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] G Nodal gradients
//! \param[in] U Solution vector of primitive variables at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // number of transported scalars
  auto ncomp = U.nprop();

  // configure flux function
  auto flux = FLUX();

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // domain edge contributions: tetrahedron superedges
  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
    const auto N = dsupedge[0].data() + e*4;
    tk::real u[4][ncomp];
    for (std::size_t c=0; c<ncomp; ++c) {
      u[0][c] = U(N[0],c);
      u[1][c] = U(N[1],c);
      u[2][c] = U(N[2],c);
      u[3][c] = U(N[3],c);
    }
    // edge fluxes
    tk::real f[6][ncomp];
    const auto d = dsupint[0].data();
    flux( coord, G, d+(e*6+0)*3, N[0], N[1], u[0], u[1], f[0] );
    flux( coord, G, d+(e*6+1)*3, N[1], N[2], u[1], u[2], f[1] );
    flux( coord, G, d+(e*6+2)*3, N[2], N[0], u[2], u[0], f[2] );
    flux( coord, G, d+(e*6+3)*3, N[0], N[3], u[0], u[3], f[3] );
    flux( coord, G, d+(e*6+4)*3, N[1], N[3], u[1], u[3], f[4] );
    flux( coord, G, d+(e*6+5)*3, N[2], N[3], u[2], u[3], f[5] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c] - f[3][c];
      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c] - f[4][c];
      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c] - f[5][c];
      R(N[3],c) = R(N[3],c) + f[3][c] + f[4][c] + f[5][c];
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    tk::real u[3][ncomp];
    for (std::size_t c=0; c<ncomp; ++c) {
      u[0][c] = U(N[0],c);
      u[1][c] = U(N[1],c);
      u[2][c] = U(N[2],c);
    }
    // edge fluxes
    tk::real f[3][ncomp];
    const auto d = dsupint[1].data();
    flux( coord, G, d+(e*3+0)*3, N[0], N[1], u[0], u[1], f[0] );
    flux( coord, G, d+(e*3+1)*3, N[1], N[2], u[1], u[2], f[1] );
    flux( coord, G, d+(e*3+2)*3, N[2], N[0], u[2], u[0], f[2] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c];
      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c];
      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c];
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real u[2][ncomp];
    for (std::size_t c=0; c<ncomp; ++c) {
      u[0][c] = U(N[0],c);
      u[1][c] = U(N[1],c);
    }
    // edge fluxes
    tk::real f[ncomp];
    const auto d = dsupint[2].data();
    flux( coord, G, d+e*3, N[0], N[1], u[0], u[1], f );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) -= f[c];
      R(N[1],c) += f[c];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
advbnd( const std::vector< std::size_t >& triinpoel,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::vector< std::uint8_t >& besym,
        const tk::Fields& U,
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integral for advection
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] coord Mesh node coordinates
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector
// *****************************************************************************
{
  auto ncomp = U.nprop();

  auto g = g_cfg.get< tag::mat_spec_heat_ratio >();
  auto rgas = g_cfg.get< tag::mat_spec_gas_const >();

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;

    auto rA  = U(N[0],0)/U(N[0],4)/rgas;
    auto ruA = U(N[0],1) * rA;
    auto rvA = U(N[0],2) * rA;
    auto rwA = U(N[0],3) * rA;
    auto reA = U(N[0],0)/(g-1.0) + 0.5*(ruA*ruA + rvA*rvA + rwA*rwA)/rA;

    auto rB  = U(N[1],0)/U(N[1],4)/rgas;
    auto ruB = U(N[1],1) * rB;
    auto rvB = U(N[1],2) * rB;
    auto rwB = U(N[1],3) * rB;
    auto reB = U(N[1],0)/(g-1.0) + 0.5*(ruB*ruB + rvB*rvB + rwB*rwB)/rB;

    auto rC  = U(N[2],0)/U(N[2],4)/rgas;
    auto ruC = U(N[2],1) * rC;
    auto rvC = U(N[2],2) * rC;
    auto rwC = U(N[2],3) * rC;
    auto reC = U(N[2],0)/(g-1.0) + 0.5*(ruC*ruC + rvC*rvC + rwC*rwC)/rC;

    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto [nx,ny,nz] = tk::cross( ba, ca );
    nx /= 12.0;
    ny /= 12.0;
    nz /= 12.0;

    tk::real vn, f[ncomp][3];
    const auto sym = besym.data() + e*3;

    vn = sym[0] ? 0.0 : (nx*U(N[0],1) + ny*U(N[0],2) + nz*U(N[0],3));
    // flow
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + U(N[0],0)*nx;
    f[2][0] = rvA*vn + U(N[0],0)*ny;
    f[3][0] = rwA*vn + U(N[0],0)*nz;
    f[4][0] = (reA + U(N[0],0))*vn;
    // scalar
    for (std::size_t c=5; c<ncomp; ++c) f[c][0] = U(N[0],c)*vn;

    vn = sym[1] ? 0.0 : (nx*U(N[1],1) + ny*U(N[1],2) + nz*U(N[1],3));
    // flow
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + U(N[1],0)*nx;
    f[2][1] = rvB*vn + U(N[1],0)*ny;
    f[3][1] = rwB*vn + U(N[1],0)*nz;
    f[4][1] = (reB + U(N[1],0))*vn;
    // scalar
    for (std::size_t c=5; c<ncomp; ++c) f[c][1] = U(N[1],c)*vn;

    vn = sym[2] ? 0.0 : (nx*U(N[2],1) + ny*U(N[2],2) + nz*U(N[2],3));
    // flow
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + U(N[2],0)*nx;
    f[2][2] = rvC*vn + U(N[2],0)*ny;
    f[3][2] = rwC*vn + U(N[2],0)*nz;
    f[4][2] = (reC + U(N[2],0))*vn;
    // scalar
    for (std::size_t c=5; c<ncomp; ++c) f[c][2] = U(N[2],c)*vn;

    for (std::size_t c=0; c<ncomp; ++c) {
      auto fab = (f[c][0] + f[c][1])/4.0;
      auto fbc = (f[c][1] + f[c][2])/4.0;
      auto fca = (f[c][2] + f[c][0])/4.0;
      R(N[0],c) += fab + fca + f[c][0];
      R(N[1],c) += fab + fbc + f[c][1];
      R(N[2],c) += fbc + fca + f[c][2];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
src( const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R )
// *****************************************************************************
//  Compute source integral
//! \param[in] coord Mesh node coordinates
//! \param[in] v Nodal mesh volumes without contributions from other chares
//! \param[in] t Physical time
//! \param[in] tp Physical time for each mesh node
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  auto src = problems::SRC();
  if (!src) return;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t p=0; p<R.nunk(); ++p) {
    if (g_cfg.get< tag::steady >()) t = tp[p];
    auto s = src( x[p], y[p], z[p], t );
    for (std::size_t c=0; c<s.size(); ++c) R(p,c) -= s[c] * v[p];
  }
}

void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::uint8_t >& besym,
     const tk::Fields& G,
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] G Gradients in mesh nodes
//! \param[in] U Solution vector of primitive variables at recent time step
//! \param[in] v Nodal mesh volumes without contributions from other chares
//! \param[in] t Physical time
//! \param[in] tp Physical time for each mesh node
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  R.fill( 0.0 );
  advdom( coord, dsupedge, dsupint, G, U, R );
  advbnd( triinpoel, coord, besym, U, R );
  src( coord, v, t, tp, R );
}

} // lax::
