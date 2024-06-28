// *****************************************************************************
/*!
  \file      src/Physics/Chorin.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     ChoCG: Projection-based solver for incompressible flow
*/
// *****************************************************************************

#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Chorin.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace chorin {

using inciter::g_cfg;

void
div( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const tk::Fields& U,
     std::vector< tk::real >& D )
// *****************************************************************************
//  Compute divergence of a vector in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Vector whose divergence to compute
//! \param[in,out] D Nodal divergence of vector in all points
//! \return Divergence of a vector in all mesh points
// *****************************************************************************
{

  // flux function to compute velocity divergence for edge p,q
  auto flux = [&]( const tk::real d[], std::size_t p, std::size_t q ) {
    return d[0] * ( U(p,0) + U(q,0) ) +
           d[1] * ( U(p,1) + U(q,1) ) +
           d[2] * ( U(p,2) + U(q,2) );
  };

  std::fill( begin(D), end(D), 0.0 );

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
    const auto d = dsupint[0].data();
    // edge fluxes
    tk::real f[6] = {
      flux( d+(e*6+0)*3, N[0], N[1] ),
      flux( d+(e*6+1)*3, N[1], N[2] ),
      flux( d+(e*6+2)*3, N[2], N[0] ),
      flux( d+(e*6+3)*3, N[0], N[3] ),
      flux( d+(e*6+4)*3, N[1], N[3] ),
      flux( d+(e*6+5)*3, N[2], N[3] ) };
    // edge flux contributions
    D[ N[0] ] = D[ N[0] ] - f[0] + f[2] - f[3];
    D[ N[1] ] = D[ N[1] ] + f[0] - f[1] - f[4];
    D[ N[2] ] = D[ N[2] ] + f[1] - f[2] - f[5];
    D[ N[3] ] = D[ N[3] ] + f[3] + f[4] + f[5];
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    // edge fluxes
    tk::real f[3] = {
      flux( d+(e*3+0)*3, N[0], N[1] ),
      flux( d+(e*3+1)*3, N[1], N[2] ),
      flux( d+(e*3+2)*3, N[2], N[0] ) };
    // edge flux contributions
    D[ N[0] ] = D[ N[0] ] - f[0] + f[2];
    D[ N[1] ] = D[ N[1] ] + f[0] - f[1];
    D[ N[2] ] = D[ N[2] ] + f[1] - f[2];
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data();
    // edge flux
    tk::real f = flux( d+e*3, N[0], N[1] );
    // edge flux contributions
    D[ N[0] ] -= f;
    D[ N[1] ] += f;
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
    auto [nx,ny,nz] = tk::cross( ba, ca );
    nx /= 12.0;
    ny /= 12.0;
    nz /= 12.0;
    tk::real f[3] = { nx*U(N[0],0) + ny*U(N[0],1) + nz*U(N[0],2),
                      nx*U(N[1],0) + ny*U(N[1],1) + nz*U(N[1],2),
                      nx*U(N[2],0) + ny*U(N[2],1) + nz*U(N[2],2) };
    auto fab = (f[0] + f[1])/4.0;
    auto fbc = (f[1] + f[2])/4.0;
    auto fca = (f[2] + f[0])/4.0;
    D[ N[0] ] += fab + fca + f[0];
    D[ N[1] ] += fab + fbc + f[1];
    D[ N[2] ] += fbc + fca + f[2];
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
      const std::vector< tk::real >& U,
      tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of a scalar in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Scalar whose gradient to compute
//! \param[in,out] G Nodal gradient of scalar in all points
//! \return Gradients of a scalar in all mesh points
// *****************************************************************************
{
  Assert( G.nunk() == U.size(), "Size mismatch" );
  Assert( G.nprop() == 3, "Size mismatch" );
  G.fill( 0.0 );

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
    const auto d = dsupint[0].data();
    tk::real u[] = { U[N[0]], U[N[1]], U[N[2]], U[N[3]] };
    for (std::size_t j=0; j<3; ++j) {
      tk::real f[] = {
        d[(e*6+0)*3+j] * (u[1] + u[0]),
        d[(e*6+1)*3+j] * (u[2] + u[1]),
        d[(e*6+2)*3+j] * (u[0] + u[2]),
        d[(e*6+3)*3+j] * (u[3] + u[0]),
        d[(e*6+4)*3+j] * (u[3] + u[1]),
        d[(e*6+5)*3+j] * (u[3] + u[2]) };
      G(N[0],j) = G(N[0],j) - f[0] + f[2] - f[3];
      G(N[1],j) = G(N[1],j) + f[0] - f[1] - f[4];
      G(N[2],j) = G(N[2],j) + f[1] - f[2] - f[5];
      G(N[3],j) = G(N[3],j) + f[3] + f[4] + f[5];
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    tk::real u[] = { U[N[0]], U[N[1]], U[N[2]] };
    for (std::size_t j=0; j<3; ++j) {
      tk::real f[] = {
        d[(e*3+0)*3+j] * (u[1] + u[0]),
        d[(e*3+1)*3+j] * (u[2] + u[1]),
        d[(e*3+2)*3+j] * (u[0] + u[2]) };
      G(N[0],j) = G(N[0],j) - f[0] + f[2];
      G(N[1],j) = G(N[1],j) + f[0] - f[1];
      G(N[2],j) = G(N[2],j) + f[1] - f[2];
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*3;
    tk::real u[] = { U[N[0]], U[N[1]] };
    for (std::size_t j=0; j<3; ++j) {
      tk::real f = d[j] * (u[1] + u[0]);
      G(N[0],j) -= f;
      G(N[1],j) += f;
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
    tk::real u[] = { U[N[0]], U[N[1]], U[N[2]] };
    auto uab = (u[0] + u[1])/4.0;
    auto ubc = (u[1] + u[2])/4.0;
    auto uca = (u[2] + u[0])/4.0;
    tk::real g[] = { uab + uca + u[0],
                     uab + ubc + u[1],
                     ubc + uca + u[2] };
    for (std::size_t j=0; j<3; ++j) {
      G(N[0],j) += g[j] * n[j];
      G(N[1],j) += g[j] * n[j];
      G(N[2],j) += g[j] * n[j];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
advedge( const tk::real supint[],
         const tk::Fields& U,
         const std::vector< tk::real >& P,
         const std::array< std::vector< tk::real >, 3 >& coord,
         tk::real t,
         tk::real dt,
         const std::vector< tk::real >& tp,
         const std::vector< tk::real >& dtp,
         std::size_t p,
         std::size_t q,
         tk::real f[],
         const std::function< std::vector< tk::real >
                 ( tk::real, tk::real, tk::real, tk::real ) >& src )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] supint Edge integral
//! \param[in] U Velocity and transported scalars at recent time step
//! \param[in] P Pressure
//! \param[in] coord Mesh node coordinates
//! \param[in] t Physical time
//! \param[in] dt Physical time step size
//! \param[in] tp Phisical time step size for each mesh node (if steady state)
//! \param[in] dtp Time step size for each mesh node (if steady state)
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] f Flux computed
//! \param[in] src Function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
  const auto steady = g_cfg.get< tag::steady >();
  const auto ncomp = U.nprop();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // edge vector
  auto dx = x[p] - x[q];
  auto dy = y[p] - y[q];
  auto dz = z[p] - z[q];
  auto dl = dx*dx + dy*dy + dz*dz;
  dx /= dl;
  dy /= dl;
  dz /= dl;

  // left state
  auto uL = U(p,0);
  auto vL = U(p,1);
  auto wL = U(p,2);
  auto pL = P[p];
  auto dnL = uL*dx + vL*dy + wL*dz;

  // right state
  auto uR = U(q,0);
  auto vR = U(q,1);
  auto wR = U(q,2);
  auto pR = P[q];
  auto dnR = uR*dx + vR*dy + wR*dz;

  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // Taylor-Galerkin first half step

  if (steady) dt = (dtp[p] + dtp[q])/2.0;

  tk::real ue[ncomp];

  // flow
  auto dp = pL - pR;
  ue[0] = 0.5*(uL + uR - dt*(uL*dnL - uR*dnR + dp*dx));
  ue[1] = 0.5*(vL + vR - dt*(vL*dnL - vR*dnR + dp*dy));
  ue[2] = 0.5*(wL + wR - dt*(wL*dnL - wR*dnR + dp*dz));
  // scalar
  for (std::size_t c=3; c<ncomp; ++c) {
    ue[c] = 0.5*(U(p,c) + U(q,c) - dt*(U(p,c)*dnL - U(q,c)*dnR));
  }

  // source
  if (src) {
    if (steady) t = (tp[p] + tp[q])/2.0;
    auto coef = dt/4.0;
    auto sL = src( x[p], y[p], z[p], t );
    auto sR = src( x[q], y[q], z[q], t );
    // flow + scalar
    for (std::size_t c=0; c<ncomp; ++c) {
      ue[c] += coef*(sL[c] + sR[c]);
    }
  }

  // Taylor-Galerkin second half step

  auto ruh = ue[0];
  auto rvh = ue[1];
  auto rwh = ue[2];
  auto ph = (P[p] + P[q])/2.0;
  auto vn = ruh*nx + rvh*ny + rwh*nz;

  // flow
  f[0] = 2.0*(ruh*vn + ph*nx);
  f[1] = 2.0*(rvh*vn + ph*ny);
  f[2] = 2.0*(rwh*vn + ph*nz);
  // scalar
  for (std::size_t c=3; c<ncomp; ++c) {
    f[c] = 2.0*ue[c]*vn;
  }

  // source
  if (src) {
    auto coef = -5.0/3.0*supint[3];
    auto xe = (x[p] + x[q])/2.0;
    auto ye = (y[p] + y[q])/2.0;
    auto ze = (z[p] + z[q])/2.0;
    auto se = src( xe, ye, ze, t+dt/2.0 );
    // flow + scalar
    for (std::size_t c=0; c<ncomp; ++c) {
      f[ncomp+c] = coef*se[c];
    }
  }

  // artificial viscosity

  const auto stab2 = g_cfg.get< tag::stab2 >();
  if (!stab2) return;

  auto stab2coef = g_cfg.get< tag::stab2coef >();
  auto vnL = uL*nx + vL*ny + wL*nz;
  auto vnR = uR*nx + vR*ny + wR*nz;
  auto sl = std::abs(vnL);
  auto sr = std::abs(vnR);
  auto fw = stab2coef * std::max( sl, sr );

  // flow
  f[0] -= fw*(uL - uR);
  f[1] -= fw*(vL - vR);
  f[2] -= fw*(wL - wR);
  // scalar
  for (std::size_t c=3; c<ncomp; ++c) {
    f[c] -= fw*(U(p,c) - U(q,c));
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
advdom( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
        const std::array< std::vector< tk::real >, 3 >& dsupint,
        const std::array< std::vector< tk::real >, 3 >& coord,
        tk::real t,
        tk::real dt,
        const std::vector< tk::real >& tp,
        const std::vector< tk::real >& dtp,
        const tk::Fields& U,
        const std::vector< tk::real >& P,
        // cppcheck-suppress constParameter
        tk::Fields& R )
// *****************************************************************************
//! Compute domain integral for advection
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] t Physical time
//! \param[in] dt Physical time step size
//! \param[in] tp Phisical time step size for each mesh node (if steady state)
//! \param[in] dtp Time step size for each mesh node (if steady state)
//! \param[in] U Velocity and transported scalars at recent time step
//! \param[in] P Pressure
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  auto ncomp = U.nprop();
  auto src = problems::SRC();

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
    tk::real f[6][ncomp*2];
    const auto d = dsupint[0].data();
    advedge( d+(e*6+0)*4, U, P, coord, t, dt, tp, dtp, N[0], N[1], f[0], src );
    advedge( d+(e*6+1)*4, U, P, coord, t, dt, tp, dtp, N[1], N[2], f[1], src );
    advedge( d+(e*6+2)*4, U, P, coord, t, dt, tp, dtp, N[2], N[0], f[2], src );
    advedge( d+(e*6+3)*4, U, P, coord, t, dt, tp, dtp, N[0], N[3], f[3], src );
    advedge( d+(e*6+4)*4, U, P, coord, t, dt, tp, dtp, N[1], N[3], f[4], src );
    advedge( d+(e*6+5)*4, U, P, coord, t, dt, tp, dtp, N[2], N[3], f[5], src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c] - f[3][c];
      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c] - f[4][c];
      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c] - f[5][c];
      R(N[3],c) = R(N[3],c) + f[3][c] + f[4][c] + f[5][c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c) += f[0][nc] + f[2][nc] + f[3][nc];
        R(N[1],c) += f[0][nc] + f[1][nc] + f[4][nc];
        R(N[2],c) += f[1][nc] + f[2][nc] + f[5][nc];
        R(N[3],c) += f[3][nc] + f[4][nc] + f[5][nc];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    tk::real f[3][ncomp*2];
    const auto d = dsupint[1].data();
    advedge( d+(e*3+0)*4, U, P, coord, t, dt, tp, dtp, N[0], N[1], f[0], src );
    advedge( d+(e*3+1)*4, U, P, coord, t, dt, tp, dtp, N[1], N[2], f[1], src );
    advedge( d+(e*3+2)*4, U, P, coord, t, dt, tp, dtp, N[2], N[0], f[2], src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c];
      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c];
      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c) += f[0][nc] + f[2][nc];
        R(N[1],c) += f[0][nc] + f[1][nc];
        R(N[2],c) += f[1][nc] + f[2][nc];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real f[ncomp*2];
    const auto d = dsupint[2].data();
    advedge( d+e*4, U, P, coord, t, dt, tp, dtp, N[0], N[1], f, src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) -= f[c];
      R(N[1],c) += f[c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c) += f[nc];
        R(N[1],c) += f[nc];
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
advbnd( const std::vector< std::size_t >& triinpoel,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::vector< std::uint8_t >& besym,
        const tk::Fields& U,
        const std::vector< tk::real >& P,
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integral for advection
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] coord Mesh node coordinates
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] U Solution vector at recent time step
//! \param[in] P Pressure
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

    auto uA = U(N[0],0);
    auto vA = U(N[0],1);
    auto wA = U(N[0],2);

    auto uB = U(N[1],0);
    auto vB = U(N[1],1);
    auto wB = U(N[1],2);

    auto uC = U(N[2],0);
    auto vC = U(N[2],1);
    auto wC = U(N[2],2);

    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto [nx,ny,nz] = tk::cross( ba, ca );
    nx /= 12.0;
    ny /= 12.0;
    nz /= 12.0;

    tk::real p, vn, f[ncomp][3];
    const auto sym = besym.data() + e*3;

    p = P[N[0]];
    vn = sym[0] ? 0.0 : nx*uA + ny*vA + nz*wA;
    // flow
    f[0][0] = uA*vn + p*nx;
    f[1][0] = vA*vn + p*ny;
    f[2][0] = wA*vn + p*nz;
    // scalar
    for (std::size_t c=3; c<ncomp; ++c) f[c][0] = U(N[0],c)*vn;

    p = P[N[1]];
    vn = sym[1] ? 0.0 : nx*uB + ny*vB + nz*wB;
    // flow
    f[0][1] = uB*vn + p*nx;
    f[1][1] = vB*vn + p*ny;
    f[2][1] = wB*vn + p*nz;
    // scalar
    for (std::size_t c=3; c<ncomp; ++c) f[c][1] = U(N[1],c)*vn;

    p = P[N[2]];
    vn = sym[2] ? 0.0 : nx*uC + ny*vC + nz*wC;
    // flow
    f[0][2] = uC*vn + p*nx;
    f[1][2] = vC*vn + p*ny;
    f[2][2] = wC*vn + p*nz;
    // scalar
    for (std::size_t c=3; c<ncomp; ++c) f[c][2] = U(N[2],c)*vn;

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

void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::uint8_t >& besym,
     const tk::Fields& U,
     const std::vector< tk::real >& P,
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     const std::vector< tk::real >& dtp,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] U Solution vector of primitive variables at recent time step
//! \param[in] P Pressure
//! \param[in] t Physical time
//! \param[in] dt Physical time size
//! \param[in] tp Physical time for each mesh node
//! \param[in] dtp Time step size for each mesh node (if steady state)
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  R.fill( 0.0 );
  advdom( dsupedge, dsupint, coord, t, dt, tp, dtp, U, P, R );
  advbnd( triinpoel, coord, besym, U, P, R );
}

} // chorin::
