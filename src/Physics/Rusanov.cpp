// *****************************************************************************
/*!
  \file      src/Physics/Rusanov.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Rusanov, MUSCL, limiting for edge-based continuous Galerkin
*/
// *****************************************************************************

#include "Tags.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Rusanov.hpp"
#include "Problems.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace physics {

static const tk::real muscl_eps = 1.0e-9;
static const tk::real muscl_const = 1.0/3.0;

static void
muscl( std::size_t p,
       std::size_t q,
       const tk::UnsMesh::Coords& coord,
       const tk::Fields& G,
       tk::real& rL, tk::real& uL, tk::real& vL, tk::real& wL, tk::real& eL,
       tk::real& rR, tk::real& uR, tk::real& vR, tk::real& wR, tk::real& eR )
// *****************************************************************************
//! Compute MUSCL reconstruction in edge-end points for the flow variables
//! \param[in] p Left node id of edge-end
//! \param[in] q Right node id of edge-end
//! \param[in] coord Array of nodal coordinates
//! \param[in] G Gradient of all unknowns in mesh points
//! \param[in,out] rL Left density
//! \param[in,out] uL Left X velocity
//! \param[in,out] vL Left Y velocity
//! \param[in,out] wL Left Z velocity
//! \param[in,out] eL Left internal energy
//! \param[in,out] rR Right density
//! \param[in,out] uR Right X velocity
//! \param[in,out] vR Right Y velocity
//! \param[in,out] wR Right Z velocity
//! \param[in,out] eR Right internal energy
// *****************************************************************************
{
  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  std::array< tk::real, 3 > vw{ x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[5], delta2[5], delta3[5];
  std::array< tk::real, 5 > ls{ rL, uL, vL, wL, eL };
  std::array< tk::real, 5 > rs{ rR, uR, vR, wR, eR };
  auto url = ls;
  auto urr = rs;

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<5; ++c) {
    // gradients
    std::array< tk::real, 3 >
      g1{ G(p,c*3+0,0), G(p,c*3+1,0), G(p,c*3+2,0) },
      g2{ G(q,c*3+0,0), G(q,c*3+1,0), G(q,c*3+2,0) };

    delta2[c] = rs[c] - ls[c];
    delta1[c] = 2.0 * tk::dot(g1,vw) - delta2[c];
    delta3[c] = 2.0 * tk::dot(g2,vw) - delta2[c];

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

  // force first order if the reconstructions for density or internal energy
  // would have allowed negative values
  if (ls[0] < delta1[0] || ls[4] < delta1[4]) url = ls;
  if (rs[0] < -delta3[0] || rs[4] < -delta3[4]) urr = rs;

  rL = url[0];
  uL = url[1];
  vL = url[2];
  wL = url[3];
  eL = url[4];

  rR = urr[0];
  uR = urr[1];
  vR = urr[2];
  wR = urr[3];
  eR = urr[4];
}

static void
muscl( std::size_t p,
       std::size_t q,
       const tk::UnsMesh::Coords& coord,
       const tk::Fields& G,
       std::vector< tk::real >& uL,
       std::vector< tk::real >& uR )
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
  // number of transported scalars
  auto ns = G.nprop()/3 - 5;

  Assert( uL.size() == ns && uR.size() == ns, "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  std::array< tk::real, 3 > vw{ x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  std::vector< tk::real >
    delta1( ns, 0.0 ), delta2( ns, 0.0 ), delta3( ns, 0.0 );

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<ns; ++c) {
    // gradients
    std::array< tk::real, 3 >
      g1{ G(p,(5+c)*3+0,0), G(p,(5+c)*3+1,0), G(p,(5+c)*3+2,0) },
      g2{ G(q,(5+c)*3+0,0), G(q,(5+c)*3+1,0), G(q,(5+c)*3+2,0) };

    delta2[c] = uR[c] - uL[c];
    delta1[c] = 2.0 * tk::dot(g1,vw) - delta2[c];
    delta3[c] = 2.0 * tk::dot(g2,vw) - delta2[c];

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
    uL[c] += 0.25*(delta1[c]*(1.0-muscl_const)*phiL +
                   delta2[c]*(1.0+muscl_const)*phi_L_inv);
    uR[c] -= 0.25*(delta3[c]*(1.0-muscl_const)*phiR +
                   delta2[c]*(1.0+muscl_const)*phi_R_inv);
  }
}

void
grad( const std::vector< std::size_t >& dedge,
      const std::vector< tk::real >& deint,
      const std::vector< std::size_t >& bpoin,
      const std::vector< tk::real >& bpint,
      const std::vector< std::size_t >& bedge,
      const std::vector< tk::real >& beint,
      const tk::Fields& U,
      tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of primitive variables in all points
//! \param[in] U Solution vector at recent time step
//! \return Gradients of primitive variables in all mesh points
// *****************************************************************************
{
  auto ncomp = U.nprop();

  Assert( G.nunk() == U.nunk(), "Size mismatch" );
  Assert( G.nprop() == ncomp*3, "Size mismatch" );
  G.fill( 0.0 );

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // domain edge contributions
  for (std::size_t e=0; e<dedge.size()/2; ++e) {
    auto p = dedge[e*2+0];
    auto q = dedge[e*2+1];

    // primitive variables
    tk::real uL[ncomp], uR[ncomp];
    uL[0] = U(p,0,0);
    uL[1] = U(p,1,0) / uL[0];
    uL[2] = U(p,2,0) / uL[0];
    uL[3] = U(p,3,0) / uL[0];
    uL[4] = U(p,4,0) / uL[0] - 0.5*(uL[1]*uL[1] + uL[2]*uL[2] + uL[3]*uL[3]);
    uR[0] = U(q,0,0);
    uR[1] = U(q,1,0) / uR[0];
    uR[2] = U(q,2,0) / uR[0];
    uR[3] = U(q,3,0) / uR[0];
    uR[4] = U(q,4,0) / uR[0] - 0.5*(uR[1]*uR[1] + uR[2]*uR[2] + uR[3]*uR[3]);
    for (std::size_t c=5; c<ncomp; ++c) { uL[c] = U(p,c,0); uR[c] = U(q,c,0); }

    for (std::size_t c=0; c<ncomp; ++c) {
      auto f = uL[c] + uR[c];
      auto g = deint[e*3+0] * f;
      G(p,c*3+0,0) -= g;
      G(q,c*3+0,0) += g;
      g = deint[e*3+1] * f;
      G(p,c*3+1,0) -= g;
      G(q,c*3+1,0) += g;
      g = deint[e*3+2] * f;
      G(p,c*3+2,0) -= g;
      G(q,c*3+2,0) += g;
    }
  }

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto p = bpoin[b];

    // primitive variables
    tk::real u[ncomp];
    u[0] = U(p,0,0);
    u[1] = U(p,1,0) / u[0];
    u[2] = U(p,2,0) / u[0];
    u[3] = U(p,3,0) / u[0];
    u[4] = U(p,4,0) / u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
    for (std::size_t c=5; c<ncomp; ++c) u[c] = U(p,c,0);

    for (std::size_t c=0; c<ncomp; ++c) {
      G(p,c*3+0,0) += bpint[b*3+0] * u[c];
      G(p,c*3+1,0) += bpint[b*3+1] * u[c];
      G(p,c*3+2,0) += bpint[b*3+2] * u[c];
    }
  }

  // boundary edge contributions
  for (std::size_t e=0; e<bedge.size()/2; ++e) {
    auto p = bedge[e*2+0];
    auto q = bedge[e*2+1];

    // primitive variables
    tk::real uL[ncomp], uR[ncomp];
    uL[0] = U(p,0,0);
    uL[1] = U(p,1,0) / uL[0];
    uL[2] = U(p,2,0) / uL[0];
    uL[3] = U(p,3,0) / uL[0];
    uL[4] = U(p,4,0) / uL[0] - 0.5*(uL[1]*uL[1] + uL[2]*uL[2] + uL[3]*uL[3]);
    uR[0] = U(q,0,0);
    uR[1] = U(q,1,0) / uR[0];
    uR[2] = U(q,2,0) / uR[0];
    uR[3] = U(q,3,0) / uR[0];
    uR[4] = U(q,4,0) / uR[0] - 0.5*(uR[1]*uR[1] + uR[2]*uR[2] + uR[3]*uR[3]);
    for (std::size_t c=5; c<ncomp; ++c) { uL[c] = U(p,c,0); uR[c] = U(q,c,0); }

    for (std::size_t c=0; c<ncomp; ++c) {
      auto f = uL[c] + uR[c];
      auto g = beint[e*3+0] * f;
      G(p,c*3+0,0) -= g;
      G(q,c*3+0,0) += g;
      g = beint[e*3+1] * f;
      G(p,c*3+1,0) -= g;
      G(q,c*3+1,0) += g;
      g = beint[e*3+2] * f;
      G(p,c*3+2,0) -= g;
      G(q,c*3+2,0) += g;
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
advdom( const tk::UnsMesh::Coords& coord,
        const std::vector< std::size_t >& dedge,
        const std::vector< tk::real >& deint, 
        const tk::Fields& G,
        const tk::Fields& U,
        tk::Fields& R )
// *****************************************************************************
//! Compute domain-edge integral for advection
//! \param[in] coord Mesh node coordinates
//! \param[in] G Nodal gradients
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // number of transported scalars
  auto ns = U.nprop() - 5;

  // domain-edge integral: compute fluxes in edges
  for (std::size_t e=0; e<dedge.size()/2; ++e) {
    // edge-end points
    auto p = dedge[e*2+0];
    auto q = dedge[e*2+1];

    // primitive variables at edge-end points
    auto rL  = U(p,0,0);
    auto ruL = U(p,1,0) / rL;
    auto rvL = U(p,2,0) / rL;
    auto rwL = U(p,3,0) / rL;
    auto reL = U(p,4,0) / rL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL);
    auto rR  = U(q,0,0);
    auto ruR = U(q,1,0) / rR;
    auto rvR = U(q,2,0) / rR;
    auto rwR = U(q,3,0) / rR;
    auto reR = U(q,4,0) / rR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR);

    // MUSCL reconstruction in edge-end points for flow variables
    muscl( p, q, coord, G, rL, ruL, rvL, rwL, reL, rR, ruR, rvR, rwR, reR );

    // pressure
    auto pL = eos::pressure( rL, reL );
    auto pR = eos::pressure( rR, reR );

    // dualface-normal velocities
    auto nx = deint[e*3+0];
    auto ny = deint[e*3+1];
    auto nz = deint[e*3+2];
    auto vnL = ruL*nx + rvL*ny + rwL*nz;
    auto vnR = ruR*nx + rvR*ny + rwR*nz;

    // back to conserved variables
    reL = (reL + 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)) * rL;
    ruL *= rL;
    rvL *= rL;
    rwL *= rL;
    reR = (reR + 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)) * rR;
    ruR *= rR;
    rvR *= rR;
    rwR *= rR;

    // dissipation
    auto len = tk::length( nx, ny, nz );
    auto sl = std::abs(vnL) + eos::soundspeed(rL,pL)*len;
    auto sr = std::abs(vnR) + eos::soundspeed(rR,pR)*len;
    auto fw = std::max( sl, sr );

    // fluxes
    auto f = rL*vnL + rR*vnR + fw*(rR - rL);
    R(p,0,0) -= f;
    R(q,0,0) += f;
    f = ruL*vnL + ruR*vnR + (pL + pR)*nx + fw*(ruR - ruL);
    R(p,1,0) -= f;
    R(q,1,0) += f;
    f = rvL*vnL + rvR*vnR + (pL + pR)*ny + fw*(rvR - rvL);
    R(p,2,0) -= f;
    R(q,2,0) += f;
    f = rwL*vnL + rwR*vnR + (pL + pR)*nz + fw*(rwR - rwL);
    R(p,3,0) -= f;
    R(q,3,0) += f;
    f = (reL + pL)*vnL + (reR + pR)*vnR + fw*(reR - reL);
    R(p,4,0) -= f;
    R(q,4,0) += f;

    if (!ns) continue;

    // scalars at edge-end points
    std::vector< tk::real > uL( ns );
    std::vector< tk::real > uR( ns );
    for (std::size_t c=0; c<ns; ++c) {
      uL[c] = U(p,5+c,0);
      uR[c] = U(q,5+c,0);
    }

    // MUSCL reconstruction in edge-end points for scalars
    muscl( p, q, coord, G, uL, uR );

    // scalar dissipation
    auto sw = std::max( std::abs(vnL), std::abs(vnR) );

    // scalar fluxes
    for (std::size_t c=0; c<ns; ++c) {
      auto s = uL[c]*vnL + uR[c]*vnR + sw*(uR[c] - uL[c]);
      R(p,5+c,0) -= s;
      R(q,5+c,0) += s;
    }
  }
}

static void
advbnd( const tk::UnsMesh::Coords& coord,
        const std::vector< std::size_t >& bpoin,
        const std::vector< tk::real >& bpint,
        const std::vector< std::size_t >& bedge,
        const std::vector< tk::real >& beint,
        const std::vector< std::uint8_t >& bpsym,
        const std::vector< std::uint8_t >& besym,
        const tk::Fields& G,
        const tk::Fields& U,
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integrals for advection
//! \param[in] triinpoel Boundary triangle face connecitivity with local ids
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // number of transported scalars
  auto ns = U.nprop() - 5;

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto p = bpoin[b];

    // primitive variables at boundary point
    auto r = U(p,0,0);
    auto u = U(p,1,0) / r;
    auto v = U(p,2,0) / r;
    auto w = U(p,3,0) / r;
    auto pr = eos::pressure( r, u, v, w, U(p,4,0) );

    // boundary-normal velocity
    auto nx = bpint[b*3+0];
    auto ny = bpint[b*3+1];
    auto nz = bpint[b*3+2];
    auto vn = bpsym[b] ? 0.0 : (nx*u + ny*v + nz*w);

    // fluxes
    R(p,0,0) += U(p,0,0)*vn;
    R(p,1,0) += U(p,1,0)*vn + pr*nx;
    R(p,2,0) += U(p,2,0)*vn + pr*ny;
    R(p,3,0) += U(p,3,0)*vn + pr*nz;
    R(p,4,0) += (U(p,4,0) + pr)*vn;

    // scalar fluxes
    for (std::size_t c=5; c<U.nprop(); ++c) {
      R(p,c,0) += U(p,c,0)*vn;
    }
  }

  // boundary edge contributions
  for (std::size_t e=0; e<bedge.size()/2; ++e) {
    auto p = bedge[e*2+0];
    auto q = bedge[e*2+1];

    // primitive variables at boundary-edge end-points
    auto rL  = U(p,0,0);
    auto ruL = U(p,1,0) / rL;
    auto rvL = U(p,2,0) / rL;
    auto rwL = U(p,3,0) / rL;
    auto reL = U(p,4,0) / rL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL);
    auto rR  = U(q,0,0);
    auto ruR = U(q,1,0) / rR;
    auto rvR = U(q,2,0) / rR;
    auto rwR = U(q,3,0) / rR;
    auto reR = U(q,4,0) / rR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR);

    // MUSCL reconstruction in boundary-edge-end points for flow variables
    muscl( p, q, coord, G, rL, ruL, rvL, rwL, reL, rR, ruR, rvR, rwR, reR );

    // pressure
    auto pL = eos::pressure( rL, reL );
    auto pR = eos::pressure( rR, reR );

    // boundary-normal velocities in boundary-edge end-points
    auto nx = beint[e*3+0];
    auto ny = beint[e*3+1];
    auto nz = beint[e*3+2];
    auto vnL = besym[e*2+0] ? 0.0 : (nx*ruL + ny*rvL + nz*rwL);
    auto vnR = besym[e*2+1] ? 0.0 : (nx*ruR + ny*rvR + nz*rwR);

    // back to conserved variables
    reL = (reL + 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)) * rL;
    ruL *= rL;
    rvL *= rL;
    rwL *= rL;
    reR = (reR + 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)) * rR;
    ruR *= rR;
    rvR *= rR;
    rwR *= rR;

    // dissipation
    auto len = tk::length( nx, ny, nz );
    auto sl = std::abs(vnL) + eos::soundspeed(rL,pL)*len;
    auto sr = std::abs(vnR) + eos::soundspeed(rR,pR)*len;
    auto fw = std::max( sl, sr );

    // fluxes
    auto f = rL*vnL + rR*vnR + fw*(rR - rL);
    R(p,0,0) -= f;
    R(q,0,0) += f;
    f = ruL*vnL + ruR*vnR + (pL + pR)*nx + fw*(ruR - ruL);
    R(p,1,0) -= f;
    R(q,1,0) += f;
    f = rvL*vnL + rvR*vnR + (pL + pR)*ny + fw*(rvR - rvL);
    R(p,2,0) -= f;
    R(q,2,0) += f;
    f = rwL*vnL + rwR*vnR + (pL + pR)*nz + fw*(rwR - rwL);
    R(p,3,0) -= f;
    R(q,3,0) += f;
    f = (reL + pL)*vnL + (reR + pR)*vnR + fw*(reR - reL);
    R(p,4,0) -= f;
    R(q,4,0) += f;

    if (!ns) continue;

    // scalars at edge-end points
    std::vector< tk::real > uL( ns );
    std::vector< tk::real > uR( ns );
    for (std::size_t c=0; c<ns; ++c) {
      uL[c] = U(p,5+c,0);
      uR[c] = U(q,5+c,0);
    }

    // compute MUSCL reconstruction in boundary-edge-end points for scalars
    muscl( p, q, coord, G, uL, uR );

    // scalar dissipation
    auto sw = std::max( std::abs(vnL), std::abs(vnR) );

    // scalar fluxes
    for (std::size_t c=0; c<ns; ++c) {
      auto s = uL[c]*vnL + uR[c]*vnR + sw*(uR[c] - uL[c]);
      R(p,5+c,0) -= s;
      R(q,5+c,0) += s;
    }
  }
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
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  using inciter::g_inputdeck;

  auto src = problems::SRC();
  if (!src) return;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t p=0; p<R.nunk(); ++p) {
    std::array< tk::real, 6 > s;
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) t = tp[p];
    src( x[p], y[p], z[p], t, s[0], s[1], s[2], s[3], s[4], s[5] );
    for (std::size_t c=0; c<5; ++c) R(p,c,0) -= s[c] * v[p];
  }
}

void
rhs( const std::vector< std::size_t >& dedge,
     const std::vector< tk::real >& deint,
     const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::size_t >& bedge,
     const std::vector< tk::real >& beint,
     const std::vector< std::uint8_t >& bpsym,
     const std::vector< std::uint8_t >& besym,
     const tk::UnsMesh::Coords& coord,
     const tk::Fields& G,
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  // zero right hand side for all components
  R.fill( 0.0 );

  // advection: domain-edge integral
  advdom( coord, dedge, deint, G, U, R );

  // advection: boundary integrals
  advbnd( coord, bpoin, bpint, bedge, beint, bpsym, besym, G, U, R );

  // source
  src( coord, v, t, tp, R );
}

tk::real
dt( const std::vector< tk::real >& vol, const tk::Fields& U )
// *****************************************************************************
//  Compute the minimum time step size (for unsteady time stepping)
//! \param[in] vol Nodal volume (with contributions from other chares)
//! \param[in] U Solution vector at recent time step
//! \return Minimum time step size
// *****************************************************************************
{
  using inciter::g_inputdeck;

  tk::real mindt = std::numeric_limits< tk::real >::max();
  for (std::size_t p=0; p<U.nunk(); ++p) {
    auto r  = U(p,0,0);
    auto u  = U(p,1,0) / r;
    auto v  = U(p,2,0) / r;
    auto w  = U(p,3,0) / r;
    auto re = U(p,4,0);
    auto vel = tk::length( u, v, w );
    auto pr = eos::pressure( r, u, v, w, re );
    auto c = eos::soundspeed( r, pr );
    auto L = std::cbrt( vol[p] );
    auto euler_dt = L / std::max( vel+c, 1.0e-8 );
    mindt = std::min( mindt, euler_dt );
  }

  mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  return mindt;
}

void
dt( const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp )
// *****************************************************************************
//  Compute a time step size for each mesh node (for steady time stepping)
//! \param[in] vol Nodal volume (with contributions from other chares)
//! \param[in] U Solution vector at recent time step
//! \param[in,out] dtp Time step size for each mesh node
// *****************************************************************************
{
  using inciter::g_inputdeck;

  for (std::size_t i=0; i<U.nunk(); ++i) {
    // compute cubic root of element volume as the characteristic length
    const auto L = std::cbrt( vol[i] );
    // access solution at node p at recent time step
    const auto u = U[i];
    // compute pressure
    auto ie = u[4]/u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0]/u[0];
    auto p = eos::pressure( u[0], ie );
    //auto p = eos::pressure( u[0], u[1]/u[0], u[2]/u[0], u[3]/u[0], u[4] );
    if (p < 0) p = 0.0;
    auto c = eos::soundspeed( u[0], p );
    // characteristic velocity
    auto v = std::sqrt((u[1]*u[1] + u[2]*u[2] + u[3]*u[3])/u[0]/u[0]) + c;
    // compute dt for node
    dtp[i] = L / v * g_inputdeck.get< tag::discr, tag::cfl >();
  }
}

} // physics::
