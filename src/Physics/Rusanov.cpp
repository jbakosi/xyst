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
  tk::real vw[3] = { x[q]-x[p], y[q]-y[p], z[q]-z[p] };

  tk::real delta1[5], delta2[5], delta3[5];
  tk::real ls[5] = { rL, uL, vL, wL, eL },
           rs[5] = { rR, uR, vR, wR, eR },
           url[5], urr[5];
  memcpy( url, ls, sizeof ls );
  memcpy( urr, rs, sizeof rs );

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<5; ++c) {

    auto g1 = G(p,c*3+0,0)*vw[0] + G(p,c*3+1,0)*vw[1] + G(p,c*3+2,0)*vw[2];
    auto g2 = G(q,c*3+0,0)*vw[0] + G(q,c*3+1,0)*vw[1] + G(q,c*3+2,0)*vw[2];

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

  // MUSCL reconstruction of edge-end-point primitive variables
  for (std::size_t c=0; c<ns; ++c) {
    auto g = (5+c)*3;
    auto g1 = G(p,g+0,0)*vw[0] + G(p,g+1,0)*vw[1] + G(p,g+2,0)*vw[2];
    auto g2 = G(q,g+0,0)*vw[0] + G(q,g+1,0)*vw[1] + G(q,g+2,0)*vw[2];

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
  u[0] = U(i,0,0);
  u[1] = U(i,1,0) / u[0];
  u[2] = U(i,2,0) / u[0];
  u[3] = U(i,3,0) / u[0];
  u[4] = U(i,4,0) / u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
  for (std::size_t c=5; c<ncomp; ++c) u[c] = U(i,c,0);
}

static void
advedge( const tk::UnsMesh::Coords& coord,
         const tk::Fields& G,
         const tk::real dsupint[],
         std::size_t p,
         std::size_t q,
         const tk::real L[],
         const tk::real R[],
         tk::real f[],
         std::size_t symL = 0,
         std::size_t symR = 0 )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] coord Mesh node coordinates
//! \param[in] G Nodal gradients
//! \param[in] dsupint Domain superedge integral for this edge
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] L Left physics state variables
//! \param[in,out] R Rigth physics state variables
//! \param[in,out] f Flux computed
//! \param[in] symL Non-zero if left edge end-points is  on a symmetry boundary
//! \param[in] symR Non-zero if right edge end-points is  on a symmetry boundary
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
  auto pL = eos::pressure( l[0], l[4] );
  auto pR = eos::pressure( r[0], r[4] );

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

void
grad( const std::vector< std::size_t >& bpoin,
      const std::vector< tk::real >& bpint,
      const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< std::size_t >, 2 >& bsupedge,
      const std::array< std::vector< tk::real >, 2 >& bsupint,
      const tk::Fields& U,
      tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of primitive variables in all points
//! \param[in] bpoin Streamable boundary point local ids
//! \param[in] bpint Streamable boundary point integrals
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] bsupedge Boundary superedges
//! \param[in] bsupint Boundary superedge integrals
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
        f[0] = d[(e*6+0)*3+j] * ( u[1][c] + u[0][c] );
        f[1] = d[(e*6+1)*3+j] * ( u[2][c] + u[1][c] );
        f[2] = d[(e*6+2)*3+j] * ( u[0][c] + u[2][c] );
        f[3] = d[(e*6+3)*3+j] * ( u[3][c] + u[0][c] );
        f[4] = d[(e*6+4)*3+j] * ( u[3][c] + u[1][c] );
        f[5] = d[(e*6+5)*3+j] * ( u[3][c] + u[2][c] );
        G(N[0],c*3+j,0) = G(N[0],c*3+j,0) - f[0] + f[2] - f[3];
        G(N[1],c*3+j,0) = G(N[1],c*3+j,0) + f[0] - f[1] - f[4];
        G(N[2],c*3+j,0) = G(N[2],c*3+j,0) + f[1] - f[2] - f[5];
        G(N[3],c*3+j,0) = G(N[3],c*3+j,0) + f[3] + f[4] + f[5];
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
        f[0] = d[(e*3+0)*3+j] * ( u[1][c] + u[0][c] );
        f[1] = d[(e*3+1)*3+j] * ( u[2][c] + u[1][c] );
        f[2] = d[(e*3+2)*3+j] * ( u[0][c] + u[2][c] );
        G(N[0],c*3+j,0) = G(N[0],c*3+j,0) - f[0] + f[2];
        G(N[1],c*3+j,0) = G(N[1],c*3+j,0) + f[0] - f[1];
        G(N[2],c*3+j,0) = G(N[2],c*3+j,0) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real u[2][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = dsupint[2][e*3+j] * ( u[1][c] + u[0][c] );
        G(N[0],c*3+j,0) -= f;
        G(N[1],c*3+j,0) += f;
      }
    }
  }

  // boundary integrals

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto p = bpoin[b];
    tk::real u[ncomp];
    primitive( ncomp, p, U, u );
    for (std::size_t c=0; c<ncomp; ++c) {
      G(p,c*3+0,0) += bpint[b*3+0] * u[c];
      G(p,c*3+1,0) += bpint[b*3+1] * u[c];
      G(p,c*3+2,0) += bpint[b*3+2] * u[c];
    }
  }

  // boundary edge contributions: triangle superedges
  for (std::size_t e=0; e<bsupedge[0].size()/6; ++e) {
    const auto N = bsupedge[0].data() + e*6;
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[3];
        const auto b = bsupint[0].data();
        f[0] = b[(e*3+0)*3+j] * ( u[1][c] + u[0][c] );
        f[1] = b[(e*3+1)*3+j] * ( u[2][c] + u[1][c] );
        f[2] = b[(e*3+2)*3+j] * ( u[0][c] + u[2][c] );
        G(N[0],c*3+j,0) = G(N[0],c*3+j,0) - f[0] + f[2];
        G(N[1],c*3+j,0) = G(N[1],c*3+j,0) + f[0] - f[1];
        G(N[2],c*3+j,0) = G(N[2],c*3+j,0) + f[1] - f[2];
      }
    }
  }

  // boundary edge contributions: edges
  for (std::size_t e=0; e<bsupedge[1].size()/4; ++e) {
    const auto N = bsupedge[1].data() + e*4;
    tk::real u[2][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    for (std::size_t c=0; c<ncomp; ++c) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = bsupint[1][e*3+j] * ( u[1][c] + u[0][c] );
        G(N[0],c*3+j,0) -= f;
        G(N[1],c*3+j,0) += f;
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
adv( const tk::UnsMesh::Coords& coord,
     const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::uint8_t >& bpsym,
     const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< std::size_t >, 2 >& bsupedge,
     const std::array< std::vector< tk::real >, 2 >& bsupint,
     const tk::Fields& G,
     const tk::Fields& U,
     // cppcheck-suppress constParameter
     tk::Fields& R )
// *****************************************************************************
//! Compute integrals for advection
//! \param[in] coord Mesh node coordinates
//! \param[in] bpoin Boundary point local ids
//! \param[in] bpint Boundary point integrals
//! \param[in] bpsym Boundary point symmetry BC flags
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] bsupedge Boundary superedges
//! \param[in] bsupint Boundary superedge integrals
//! \param[in] G Nodal gradients
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  // number of transported scalars
  auto ncomp = U.nprop();

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // domain integral

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
    advedge( coord, G, d+(e*6+0)*3, N[0], N[1], u[0], u[1], f[0] );
    advedge( coord, G, d+(e*6+1)*3, N[1], N[2], u[1], u[2], f[1] );
    advedge( coord, G, d+(e*6+2)*3, N[2], N[0], u[2], u[0], f[2] );
    advedge( coord, G, d+(e*6+3)*3, N[0], N[3], u[0], u[3], f[3] );
    advedge( coord, G, d+(e*6+4)*3, N[1], N[3], u[1], u[3], f[4] );
    advedge( coord, G, d+(e*6+5)*3, N[2], N[3], u[2], u[3], f[5] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c] - f[3][c];
      R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c] - f[4][c];
      R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c] - f[5][c];
      R(N[3],c,0) = R(N[3],c,0) + f[3][c] + f[4][c] + f[5][c];
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
    advedge( coord, G, d+(e*3+0)*3, N[0], N[1], u[0], u[1], f[0] );
    advedge( coord, G, d+(e*3+1)*3, N[1], N[2], u[1], u[2], f[1] );
    advedge( coord, G, d+(e*3+2)*3, N[2], N[0], u[2], u[0], f[2] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c];
      R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c];
      R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c];
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
    advedge( coord, G, dsupint[2].data()+e*3, N[0], N[1], u[0], u[1], f );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) -= f[c];
      R(N[1],c,0) += f[c];
    }
  }

  // boundary integrals

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto p = bpoin[b];
    tk::real u[ncomp];
    primitive( ncomp, p, U, u );
    auto pr = eos::pressure( u[0], u[4] );
    // boundary-normal velocity
    auto nx = bpint[b*3+0];
    auto ny = bpint[b*3+1];
    auto nz = bpint[b*3+2];
    auto vn = bpsym[b] ? 0.0 : (nx*u[1] + ny*u[2] + nz*u[3]);
    // flow fluxes
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

  // boundary edge contributions: triangle superedges
  for (std::size_t e=0; e<bsupedge[0].size()/6; ++e) {
    const auto N = bsupedge[0].data() + e*6;
    // primitive variables
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    // edge fluxes
    tk::real f[3][ncomp];
    const auto b = bsupint[0].data();
    advedge( coord, G, b+(e*3+0)*3, N[0], N[1], u[0], u[1], f[0], N[3], N[4] );
    advedge( coord, G, b+(e*3+1)*3, N[1], N[2], u[1], u[2], f[1], N[4], N[5] );
    advedge( coord, G, b+(e*3+2)*3, N[2], N[0], u[2], u[0], f[2], N[5], N[3] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c];
      R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c];
      R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c];
    }
  }

  // boundary edge contributions: edges
  for (std::size_t e=0; e<bsupedge[1].size()/4; ++e) {
    const auto N = bsupedge[1].data() + e*4;
    tk::real u[2][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    // edge fluxes
    tk::real f[ncomp];
    const auto b = bsupint[1].data();
    advedge( coord, G, b+e*3, N[0], N[1], u[0], u[1], f, N[2], N[3] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) -= f[c];
      R(N[1],c,0) += f[c];
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
  using inciter::g_inputdeck;

  auto src = problems::SRC();
  if (!src) return;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t p=0; p<R.nunk(); ++p) {
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) t = tp[p];
    auto s = src( x[p], y[p], z[p], t );
    for (std::size_t c=0; c<s.size(); ++c) R(p,c,0) -= s[c] * v[p];
  }
}

void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< std::size_t >, 2 >& bsupedge,
     const std::array< std::vector< tk::real >, 2 >& bsupint,
     const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::uint8_t >& bpsym,
     const std::array< std::vector< tk::real >, 3 >& coord,
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
//! \param[in] bsupedge Boundary superedges
//! \param[in] bsupint Boundary superedge integrals
//! \param[in] bpoin Boundary point local ids
//! \param[in] bpint Boundary point integrals
//! \param[in] bpsym Boundary point symmetry BC flags
//! \param[in] coord Mesh node coordinates
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

  // zero right hand side for all components
  R.fill( 0.0 );

  // advection
  adv( coord, bpoin, bpint, bpsym, dsupedge, dsupint, bsupedge, bsupint,
       G, U, R );

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
    auto re = U(p,4,0) / r - 0.5*(u*u + v*v + w*w);
    auto vel = tk::length( u, v, w );
    auto pr = eos::pressure( r, re );
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
  auto cfl = g_inputdeck.get< tag::discr, tag::cfl >();

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
