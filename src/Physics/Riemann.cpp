// *****************************************************************************
/*!
  \file      src/Physics/Riemann.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Riemann, MUSCL, limiting for edge-based continuous Galerkin
*/
// *****************************************************************************

#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Riemann.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace riemann {

static const tk::real muscl_eps = 1.0e-9;
static const tk::real muscl_const = 1.0/3.0;

using inciter::g_cfg;

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
  tk::real vw[3] = { x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[5], delta2[5], delta3[5];
  tk::real ls[5] = { rL, uL, vL, wL, eL },
           rs[5] = { rR, uR, vR, wR, eR },
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

  // force first order if the reconstructions for density or internal energy
  // would have allowed negative values
  if (ls[0] < delta1[0] || ls[4] < delta1[4]) memcpy( url, ls, sizeof ls );
  if (rs[0] < -delta3[0] || rs[4] < -delta3[4]) memcpy( urr, rs, sizeof rs );

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

  // MUSCL reconstruction of edge-end-point of primitive variables
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

static void
primitive( std::size_t ncomp, std::size_t i, const tk::Fields& U, tk::real u[] )
// *****************************************************************************
//! Compute primitive flow variables from conserved ones
//! \param[in] ncomp Number of scalar components: 5 + scalars
//! \param[in] i Index to read conserved variables from
//! \param[in] U Solution vector to read conserved variables from
//! \param[in,out] u Computed primitive variables
// *****************************************************************************
{
  u[0] = U(i,0);
  u[1] = U(i,1) / u[0];
  u[2] = U(i,2) / u[0];
  u[3] = U(i,3) / u[0];
  u[4] = U(i,4) / u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
  for (std::size_t c=5; c<ncomp; ++c) u[c] = U(i,c);
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
//! \param[in] U Solution vector at recent time step
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
    tk::real u[4][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    primitive( ncomp, N[3], U, u[3] );
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[6];
        const auto d = dsupint[0].data();
        f[0] = d[(e*6+0)*3+j] * (u[1][c] + u[0][c]);
        f[1] = d[(e*6+1)*3+j] * (u[2][c] + u[1][c]);
        f[2] = d[(e*6+2)*3+j] * (u[0][c] + u[2][c]);
        f[3] = d[(e*6+3)*3+j] * (u[3][c] + u[0][c]);
        f[4] = d[(e*6+4)*3+j] * (u[3][c] + u[1][c]);
        f[5] = d[(e*6+5)*3+j] * (u[3][c] + u[2][c]);
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
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[3];
        const auto d = dsupint[1].data();
        f[0] = d[(e*3+0)*3+j] * (u[1][c] + u[0][c]);
        f[1] = d[(e*3+1)*3+j] * (u[2][c] + u[1][c]);
        f[2] = d[(e*3+2)*3+j] * (u[0][c] + u[2][c]);
        G(N[0],c*3+j) = G(N[0],c*3+j) - f[0] + f[2];
        G(N[1],c*3+j) = G(N[1],c*3+j) + f[0] - f[1];
        G(N[2],c*3+j) = G(N[2],c*3+j) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real u[2][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    const auto d = dsupint[2].data() + e*3;
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = d[j] * (u[1][c] + u[0][c]);
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
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto n = tk::cross( ba, ca );
    n[0] /= 12.0;
    n[1] /= 12.0;
    n[2] /= 12.0;
    for (std::size_t c=0; c<ncomp; ++c) {
      auto uab = (u[0][c] + u[1][c])/4.0;
      auto ubc = (u[1][c] + u[2][c])/4.0;
      auto uca = (u[2][c] + u[0][c])/4.0;
      tk::real g[] = { uab + uca + u[0][c],
                       uab + ubc + u[1][c],
                       ubc + uca + u[2][c] };
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

static void
rusanov( const tk::UnsMesh::Coords& coord,
         const tk::Fields& G,
         const tk::real dsupint[],
         std::size_t p,
         std::size_t q,
         const tk::real L[],
         const tk::real R[],
         tk::real f[],
         std::size_t symL,
         std::size_t symR )
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
//! \param[in] symL Non-zero if left edge end-point is on a symmetry boundary
//! \param[in] symR Non-zero if right edge end-point is on a symmetry boundary
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

  auto ncomp = G.nprop() / 3;

  // will work on copies of physics variables
  tk::real l[ncomp], r[ncomp];
  memcpy( l, L, sizeof l );
  memcpy( r, R, sizeof r );

  // MUSCL reconstruction in edge-end points for flow variables
  muscl( p, q, coord, G, l[0], l[1], l[2], l[3], l[4],
                         r[0], r[1], r[2], r[3], r[4] );

  // pressure
  auto pL = eos::pressure( l[0]*l[4] );
  auto pR = eos::pressure( r[0]*r[4] );

  // dualface-normal velocities
  auto nx = dsupint[0];
  auto ny = dsupint[1];
  auto nz = dsupint[2];
  auto vnL = symL ? 0.0 : (l[1]*nx + l[2]*ny + l[3]*nz);
  auto vnR = symR ? 0.0 : (r[1]*nx + r[2]*ny + r[3]*nz);

  // back to conserved variables
  l[4] = (l[4] + 0.5*(l[1]*l[1] + l[2]*l[2] + l[3]*l[3])) * l[0];
  l[1] *= l[0];
  l[2] *= l[0];
  l[3] *= l[0];
  r[4] = (r[4] + 0.5*(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])) * r[0];
  r[1] *= r[0];
  r[2] *= r[0];
  r[3] *= r[0];

  // dissipation
  auto len = tk::length( nx, ny, nz );
  auto sl = std::abs(vnL) + eos::soundspeed(l[0],pL)*len;
  auto sr = std::abs(vnR) + eos::soundspeed(r[0],pR)*len;
  auto fw = std::max( sl, sr );

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
      tk::real f[],
      std::size_t symL,
      std::size_t symR )
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
//! \param[in] symL Non-zero if left edge end-point is on a symmetry boundary
//! \param[in] symR Non-zero if right edge end-point is on a symmetry boundary
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

  auto ncomp = G.nprop() / 3;

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
  auto qL = symL ? 0.0 : (l[1]*nx + l[2]*ny + l[3]*nz);
  auto qR = symR ? 0.0 : (r[1]*nx + r[2]*ny + r[3]*nz);

  // pressure
  auto pL = eos::pressure( l[0]*l[4] );
  auto pR = eos::pressure( r[0]*r[4] );

  // back to conserved variables
  l[4] = (l[4] + 0.5*(l[1]*l[1] + l[2]*l[2] + l[3]*l[3])) * l[0];
  l[1] *= l[0];
  l[2] *= l[0];
  l[3] *= l[0];
  r[4] = (r[4] + 0.5*(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])) * r[0];
  r[1] *= r[0];
  r[2] *= r[0];
  r[3] *= r[0];

  // sound speed
  auto cL = eos::soundspeed(l[0],pL);
  auto cR = eos::soundspeed(r[0],pR);

  // left and right wave speeds
  auto sL = fmin( qL - cL, qR - cR );
  auto sR = fmax( qL + cL, qR + cR );

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
  const auto stab2 = g_cfg.get< tag::stab2 >();
  if (stab2) {
    auto stab2coef = g_cfg.get< tag::stab2coef >();
    auto sl = std::abs(qL) + cL;
    auto sr = std::abs(qR) + cR;
    auto fws = stab2coef * std::max(sl,sr) * len;
    f[0] -= fws * (l[0] - r[0]);
    f[1] -= fws * (l[1] - r[1]);
    f[2] -= fws * (l[2] - r[2]);
    f[3] -= fws * (l[3] - r[3]);
    f[4] -= fws * (l[4] - r[4]);
  }

  if (ncomp == 5) return;

  // MUSCL reconstruction in edge-end points for scalars
  muscl( p, q, coord, G, l, r );

  // scalar
  auto sw = std::max( std::abs(qL), std::abs(qR) ) * len;
  for (std::size_t c=5; c<ncomp; ++c) {
    f[c] = (l[c]*qL + r[c]*qR)*len + sw*(r[c] - l[c]);
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
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
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // number of transported scalars
  auto ncomp = U.nprop();

  // configure flux function
  auto cfgflux = [](){
    const auto& cfg = g_cfg.get< tag::flux >();
    if (cfg == "rusanov")
      return rusanov;
    else if (cfg == "hllc")
      return hllc;
    else
      Throw( "Flux not configured" );
  };

  auto flux = cfgflux();

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
    // primitive variables
    tk::real u[4][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    primitive( ncomp, N[3], U, u[3] );
    // edge fluxes
    tk::real f[6][ncomp];
    const auto d = dsupint[0].data();
    flux( coord, G, d+(e*6+0)*3, N[0], N[1], u[0], u[1], f[0], 0, 0 );
    flux( coord, G, d+(e*6+1)*3, N[1], N[2], u[1], u[2], f[1], 0, 0 );
    flux( coord, G, d+(e*6+2)*3, N[2], N[0], u[2], u[0], f[2], 0, 0 );
    flux( coord, G, d+(e*6+3)*3, N[0], N[3], u[0], u[3], f[3], 0, 0 );
    flux( coord, G, d+(e*6+4)*3, N[1], N[3], u[1], u[3], f[4], 0, 0 );
    flux( coord, G, d+(e*6+5)*3, N[2], N[3], u[2], u[3], f[5], 0, 0 );
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
    // primitive variables
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    // edge fluxes
    tk::real f[3][ncomp];
    const auto d = dsupint[1].data();
    flux( coord, G, d+(e*3+0)*3, N[0], N[1], u[0], u[1], f[0], 0, 0 );
    flux( coord, G, d+(e*3+1)*3, N[1], N[2], u[1], u[2], f[1], 0, 0 );
    flux( coord, G, d+(e*3+2)*3, N[2], N[0], u[2], u[0], f[2], 0, 0 );
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
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    // edge fluxes
    tk::real f[ncomp];
    const auto d = dsupint[2].data();
    flux( coord, G, d+e*3, N[0], N[1], u[0], u[1], f, 0, 0 );
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

    auto rA  = U(N[0],0);
    auto ruA = U(N[0],1);
    auto rvA = U(N[0],2);
    auto rwA = U(N[0],3);
    auto reA = U(N[0],4);

    auto rB  = U(N[1],0);
    auto ruB = U(N[1],1);
    auto rvB = U(N[1],2);
    auto rwB = U(N[1],3);
    auto reB = U(N[1],4);

    auto rC  = U(N[2],0);
    auto ruC = U(N[2],1);
    auto rvC = U(N[2],2);
    auto rwC = U(N[2],3);
    auto reC = U(N[2],4);

    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto [nx,ny,nz] = tk::cross( ba, ca );
    nx /= 12.0;
    ny /= 12.0;
    nz /= 12.0;

    tk::real p, vn, f[ncomp][3];
    const auto sym = besym.data() + e*3;

    p = eos::pressure( reA - 0.5*(ruA*ruA + rvA*rvA + rwA*rwA)/rA );
    vn = sym[0] ? 0.0 : (nx*ruA + ny*rvA + nz*rwA)/rA;
    // flow
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + p*nx;
    f[2][0] = rvA*vn + p*ny;
    f[3][0] = rwA*vn + p*nz;
    f[4][0] = (reA + p)*vn;
    // scalar
    for (std::size_t c=5; c<ncomp; ++c) f[c][0] = U(N[0],c)*vn;

    p = eos::pressure( reB - 0.5*(ruB*ruB + rvB*rvB + rwB*rwB)/rB );
    vn = sym[1] ? 0.0 : (nx*ruB + ny*rvB + nz*rwB)/rB;
    // flow
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + p*nx;
    f[2][1] = rvB*vn + p*ny;
    f[3][1] = rwB*vn + p*nz;
    f[4][1] = (reB + p)*vn;
    // scalar
    for (std::size_t c=5; c<ncomp; ++c) f[c][1] = U(N[1],c)*vn;

    p = eos::pressure( reC - 0.5*(ruC*ruC + rvC*rvC + rwC*rwC)/rC );
    vn = sym[2] ? 0.0 : (nx*ruC + ny*rvC + nz*rwC)/rC;
    // flow
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + p*nx;
    f[2][2] = rvC*vn + p*ny;
    f[3][2] = rwC*vn + p*nz;
    f[4][2] = (reC + p)*vn;
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
//! \param[in] U Unknowns/solution vector in mesh nodes
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

} // riemann::
