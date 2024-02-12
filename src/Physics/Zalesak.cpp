// *****************************************************************************
/*!
  \file      src/Physics/Zalesak.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Zalesak, FCT limiting for edge-based continuous Galerkin
*/
// *****************************************************************************

#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Zalesak.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace zalesak {

using inciter::g_cfg;

static void
advedge( const tk::real supint[],
         const tk::Fields& U,
         const tk::Fields& S,
         const tk::Fields& G,
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
//! \param[in] U Solution vector to read conserved variables from
//! \param[in] S Stabilization coefficients at nodes
//! \param[in] G Nodal gradients
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
  auto rL  = U(p,0,0);
  auto ruL = U(p,1,0);
  auto rvL = U(p,2,0);
  auto rwL = U(p,3,0);
  auto reL = U(p,4,0);
  auto pL = eos::pressure( reL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)/rL );
  auto dnL = (ruL*dx + rvL*dy + rwL*dz)/rL;

  // right state
  auto rR  = U(q,0,0);
  auto ruR = U(q,1,0);
  auto rvR = U(q,2,0);
  auto rwR = U(q,3,0);
  auto reR = U(q,4,0);
  auto pR = eos::pressure( reR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)/rR );
  auto dnR = (ruR*dx + rvR*dy + rwR*dz)/rR;

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
  ue[0] = 0.5*(rL + rR - dt*(rL*dnL - rR*dnR));
  ue[1] = 0.5*(ruL + ruR - dt*(ruL*dnL - ruR*dnR + dp*dx));
  ue[2] = 0.5*(rvL + rvR - dt*(rvL*dnL - rvR*dnR + dp*dy));
  ue[3] = 0.5*(rwL + rwR - dt*(rwL*dnL - rwR*dnR + dp*dz));
  ue[4] = 0.5*(reL + reR - dt*((reL+pL)*dnL - (reR+pR)*dnR));

  // scalar
  for (std::size_t c=5; c<ncomp; ++c) {
    ue[c] = 0.5*(U(p,c,0) + U(q,c,0) - dt*(U(p,c,0)*dnL - U(q,c,0)*dnR));
  }

  // source
  if (src) {
    if (steady) t = (tp[p] + tp[q])/2.0;
    auto coef = dt/4.0;
    auto sL = src( x[p], y[p], z[p], t );
    auto sR = src( x[q], y[q], z[q], t );
    for (std::size_t c=0; c<ncomp; ++c) {
      ue[c] += coef*(sL[c] + sR[c]);
    }
  }

  // Taylor-Galerkin second half step

  auto rh  = ue[0];
  auto ruh = ue[1];
  auto rvh = ue[2];
  auto rwh = ue[3];
  auto reh = ue[4];
  auto ph = eos::pressure( reh - 0.5*(ruh*ruh + rvh*rvh + rwh*rwh)/rh );
  auto vn = (ruh*nx + rvh*ny + rwh*nz)/rh;

  // flow
  f[0] = 2.0*rh*vn;
  f[1] = 2.0*(ruh*vn + ph*nx);
  f[2] = 2.0*(rvh*vn + ph*ny);
  f[3] = 2.0*(rwh*vn + ph*nz);
  f[4] = 2.0*(reh + ph)*vn;

  // scalar
  for (std::size_t c=5; c<ncomp; ++c) {
    f[c] = 2.0*ue[c]*vn;
  }

  // source
  if (src) {
    auto coef = -5.0/3.0*supint[3];
    auto xe = (x[p] + x[q])/2.0;
    auto ye = (y[p] + y[q])/2.0;
    auto ze = (z[p] + z[q])/2.0;
    auto se = src( xe, ye, ze, t+dt/2.0 );
    for (std::size_t c=0; c<ncomp; ++c) {
      f[ncomp+c] = coef*se[c];
    }
  }

  // artificial viscosity

  const auto stab2 = g_cfg.get< tag::stab2 >();
  const auto stab4 = g_cfg.get< tag::stab4 >();
  if (!stab2 && !stab4) return;

  auto stab2coef = g_cfg.get< tag::stab2coef >();
  auto b = 0.0;
  tk::real d2rL = 0.0, d2ruL = 0.0, d2rvL = 0.0, d2rwL = 0.0, d2reL = 0.0;
  tk::real d2rR = 0.0, d2ruR = 0.0, d2rvR = 0.0, d2rwR = 0.0, d2reR = 0.0;
  if (stab4) {
    stab2coef = 1.0;
    b = std::max( S(p,0,0), S(q,0,0) );
    d2rL  = G(p,3+0,0);
    d2ruL = G(p,3+1,0);
    d2rvL = G(p,3+2,0);
    d2rwL = G(p,3+3,0);
    d2reL = G(p,3+4,0);
    d2rR  = G(q,3+0,0);
    d2ruR = G(q,3+1,0);
    d2rvR = G(q,3+2,0);
    d2rwR = G(q,3+3,0);
    d2reR = G(q,3+4,0);
  }

  auto vnL = (ruL*nx + rvL*ny + rwL*nz)/rL;
  auto vnR = (ruR*nx + rvR*ny + rwR*nz)/rR;
  auto len = tk::length( nx, ny, nz );
  auto cL = eos::soundspeed( std::max(rL,1.0e-8), std::max(pL,0.0) );
  auto cR = eos::soundspeed( std::max(rR,1.0e-8), std::max(pR,0.0) );
  auto sl = std::abs(vnL) + cL*len;
  auto sr = std::abs(vnR) + cR*len;
  auto fw = stab2coef * std::max( sl, sr );

  auto s2c = 1.0 - b;
  auto s4c = 0.25 * b * dl;
  f[0] -= fw*(s2c*(rL - rR) + s4c*(d2rR - d2rL));
  f[1] -= fw*(s2c*(ruL - ruR) + s4c*(d2ruR - d2ruL));
  f[2] -= fw*(s2c*(rvL - rvR) + s4c*(d2rvR - d2rvL));
  f[3] -= fw*(s2c*(rwL - rwR) + s4c*(d2rwR - d2rwL));
  f[4] -= fw*(s2c*(reL - reR) + s4c*(d2reR - d2reL));
  for (std::size_t c=5; c<ncomp; ++c) {
    tk::real d2cL = 0.0, d2cR = 0.0;
    if (stab4) {
      d2cL = G(p,3+c,0);
      d2cR = G(q,3+c,0);
    }
    f[c] -= fw*(s2c*(U(p,c,0) - U(q,c,0)) + s4c*(d2cR - d2cL));
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
        const tk::Fields& S,
        const tk::Fields& G,
        // cppcheck-suppress constParameter
        tk::Fields& R )
// *****************************************************************************
//! Compute domain-edge integrals for advection
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] t Physical time
//! \param[in] dt Physical time step size
//! \param[in] tp Phisical time step size for each mesh node (if steady state)
//! \param[in] dtp Time step size for each mesh node (if steady state)
//! \param[in] U Solution vector at recent time step
//! \param[in] S Stabilization coefficients at nodes
//! \param[in] G Nodal gradients
//! \param[in,out] R Right-hand side vector
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
    advedge( d+(e*6+0)*4, U,S,G, coord, t, dt, tp, dtp, N[0], N[1], f[0], src );
    advedge( d+(e*6+1)*4, U,S,G, coord, t, dt, tp, dtp, N[1], N[2], f[1], src );
    advedge( d+(e*6+2)*4, U,S,G, coord, t, dt, tp, dtp, N[2], N[0], f[2], src );
    advedge( d+(e*6+3)*4, U,S,G, coord, t, dt, tp, dtp, N[0], N[3], f[3], src );
    advedge( d+(e*6+4)*4, U,S,G, coord, t, dt, tp, dtp, N[1], N[3], f[4], src );
    advedge( d+(e*6+5)*4, U,S,G, coord, t, dt, tp, dtp, N[2], N[3], f[5], src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c] - f[3][c];
      R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c] - f[4][c];
      R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c] - f[5][c];
      R(N[3],c,0) = R(N[3],c,0) + f[3][c] + f[4][c] + f[5][c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c,0) += f[0][nc] + f[2][nc] + f[3][nc];
        R(N[1],c,0) += f[0][nc] + f[1][nc] + f[4][nc];
        R(N[2],c,0) += f[1][nc] + f[2][nc] + f[5][nc];
        R(N[3],c,0) += f[3][nc] + f[4][nc] + f[5][nc];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    tk::real f[3][ncomp*2];
    const auto d = dsupint[1].data();
    advedge( d+(e*3+0)*4, U,S,G, coord, t, dt, tp, dtp, N[0], N[1], f[0], src );
    advedge( d+(e*3+1)*4, U,S,G, coord, t, dt, tp, dtp, N[1], N[2], f[1], src );
    advedge( d+(e*3+2)*4, U,S,G, coord, t, dt, tp, dtp, N[2], N[0], f[2], src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c];
      R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c];
      R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c,0) += f[0][nc] + f[2][nc];
        R(N[1],c,0) += f[0][nc] + f[1][nc];
        R(N[2],c,0) += f[1][nc] + f[2][nc];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real f[ncomp*2];
    const auto d = dsupint[2].data();
    advedge( d+e*4, U,S,G, coord, t, dt, tp, dtp, N[0], N[1], f, src );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) -= f[c];
      R(N[1],c,0) += f[c];
      if (src) {
        auto nc = ncomp + c;
        R(N[0],c,0) += f[nc];
        R(N[1],c,0) += f[nc];
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
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integrals for advection
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

    auto rA  = U(N[0],0,0);
    auto ruA = U(N[0],1,0);
    auto rvA = U(N[0],2,0);
    auto rwA = U(N[0],3,0);
    auto reA = U(N[0],4,0);

    auto rB  = U(N[1],0,0);
    auto ruB = U(N[1],1,0);
    auto rvB = U(N[1],2,0);
    auto rwB = U(N[1],3,0);
    auto reB = U(N[1],4,0);

    auto rC  = U(N[2],0,0);
    auto ruC = U(N[2],1,0);
    auto rvC = U(N[2],2,0);
    auto rwC = U(N[2],3,0);
    auto reC = U(N[2],4,0);

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
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + p*nx;
    f[2][0] = rvA*vn + p*ny;
    f[3][0] = rwA*vn + p*nz;
    f[4][0] = (reA + p)*vn;
    for (std::size_t c=5; c<ncomp; ++c) f[c][0] = U(N[0],0,0)*vn;

    p = eos::pressure( reB - 0.5*(ruB*ruB + rvB*rvB + rwB*rwB)/rB );
    vn = sym[1] ? 0.0 : (nx*ruB + ny*rvB + nz*rwB)/rB;
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + p*nx;
    f[2][1] = rvB*vn + p*ny;
    f[3][1] = rwB*vn + p*nz;
    f[4][1] = (reB + p)*vn;
    for (std::size_t c=5; c<ncomp; ++c) f[c][1] = U(N[1],0,0)*vn;

    p = eos::pressure( reC - 0.5*(ruC*ruC + rvC*rvC + rwC*rwC)/rC );
    vn = sym[2] ? 0.0 : (nx*ruC + ny*rvC + nz*rwC)/rC;
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + p*nx;
    f[2][2] = rvC*vn + p*ny;
    f[3][2] = rwC*vn + p*nz;
    f[4][2] = (reC + p)*vn;
    for (std::size_t c=5; c<ncomp; ++c) f[c][2] = U(N[2],0,0)*vn;

    for (std::size_t c=0; c<ncomp; ++c) {
      auto fab = (f[c][0] + f[c][1])/4.0;
      auto fbc = (f[c][1] + f[c][2])/4.0;
      auto fca = (f[c][2] + f[c][0])/4.0;
      R(N[0],c,0) += fab + fca + f[c][0];
      R(N[1],c,0) += fab + fbc + f[c][1];
      R(N[2],c,0) += fbc + fca + f[c][2];
    }
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
//  Compute nodal gradients in all points
//! \param[in] bpoin Streamable boundary point local ids
//! \param[in] bpint Streamable boundary point integrals
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] bsupedge Boundary superedges
//! \param[in] bsupint Boundary superedge integrals
//! \param[in] U Solution vector at recent time step
//! \param[in,out] G Nodal gradients computed
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
  Assert( G.nprop() == 3 + ncomp, "Size mismatch" );
  G.fill( 0.0 );

  // domain integral

  // domain edge contributions: tetrahedron superedges
  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
    const auto N = dsupedge[0].data() + e*4;
    tk::real p[4];
    for (std::size_t a=0; a<4; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }
    tk::real f[6];
    const auto d = dsupint[0].data();
    // pressure gradient
    for (std::size_t j=0; j<3; ++j) {
      f[0] = d[(e*6+0)*4+j] * (p[1] + p[0]);
      f[1] = d[(e*6+1)*4+j] * (p[2] + p[1]);
      f[2] = d[(e*6+2)*4+j] * (p[0] + p[2]);
      f[3] = d[(e*6+3)*4+j] * (p[3] + p[0]);
      f[4] = d[(e*6+4)*4+j] * (p[3] + p[1]);
      f[5] = d[(e*6+5)*4+j] * (p[3] + p[2]);
      G(N[0],j,0) = G(N[0],j,0) - f[0] + f[2] - f[3];
      G(N[1],j,0) = G(N[1],j,0) + f[0] - f[1] - f[4];
      G(N[2],j,0) = G(N[2],j,0) + f[1] - f[2] - f[5];
      G(N[3],j,0) = G(N[3],j,0) + f[3] + f[4] + f[5];
    }
    // Laplacian of conserved quantities
    for (std::size_t c=0; c<ncomp; ++c) {
      f[0] = d[(e*6+0)*4+3] * (U(N[0],c,0) - U(N[1],c,0));
      f[1] = d[(e*6+1)*4+3] * (U(N[1],c,0) - U(N[2],c,0));
      f[2] = d[(e*6+2)*4+3] * (U(N[2],c,0) - U(N[0],c,0));
      f[3] = d[(e*6+3)*4+3] * (U(N[0],c,0) - U(N[3],c,0));
      f[4] = d[(e*6+4)*4+3] * (U(N[1],c,0) - U(N[3],c,0));
      f[5] = d[(e*6+5)*4+3] * (U(N[2],c,0) - U(N[3],c,0));
      G(N[0],3+c,0) = G(N[0],3+c,0) - f[0] + f[2] - f[3];
      G(N[1],3+c,0) = G(N[1],3+c,0) + f[0] - f[1] - f[4];
      G(N[2],3+c,0) = G(N[2],3+c,0) + f[1] - f[2] - f[5];
      G(N[3],3+c,0) = G(N[3],3+c,0) + f[3] + f[4] + f[5];
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    tk::real p[3];
    for (std::size_t a=0; a<3; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }
    tk::real f[3];
    const auto d = dsupint[1].data();
    // pressure gradient
    for (std::size_t j=0; j<3; ++j) {
      f[0] = d[(e*3+0)*4+j] * (p[1] + p[0]);
      f[1] = d[(e*3+1)*4+j] * (p[2] + p[1]);
      f[2] = d[(e*3+2)*4+j] * (p[0] + p[2]);
      G(N[0],j,0) = G(N[0],j,0) - f[0] + f[2];
      G(N[1],j,0) = G(N[1],j,0) + f[0] - f[1];
      G(N[2],j,0) = G(N[2],j,0) + f[1] - f[2];
    }
    // Laplacian of conserved quantities
    for (std::size_t c=0; c<ncomp; ++c) {
      f[0] = d[(e*3+0)*4+3] * (U(N[0],c,0) - U(N[1],c,0));
      f[1] = d[(e*3+1)*4+3] * (U(N[1],c,0) - U(N[2],c,0));
      f[2] = d[(e*3+2)*4+3] * (U(N[2],c,0) - U(N[0],c,0));
      G(N[0],3+c,0) = G(N[0],3+c,0) - f[0] + f[2];
      G(N[1],3+c,0) = G(N[1],3+c,0) + f[0] - f[1];
      G(N[2],3+c,0) = G(N[2],3+c,0) + f[1] - f[2];
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real p[2];
    for (std::size_t a=0; a<2; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }
    const auto d = dsupint[2].data() + e*4;
    // pressure gradient
    for (std::size_t j=0; j<3; ++j) {
      tk::real f = d[j] * (p[1] + p[0]);
      G(N[0],j,0) -= f;
      G(N[1],j,0) += f;
    }
    // Laplacian of conserved quantities
    for (std::size_t c=0; c<ncomp; ++c) {
      tk::real f = d[3] * (U(N[0],c,0) - U(N[1],c,0));
      G(N[0],3+c,0) -= f;
      G(N[1],3+c,0) += f;
    }
  }

  // boundary integrals

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto i = bpoin[b];
    auto r  = U(i,0,0);
    auto ru = U(i,1,0);
    auto rv = U(i,2,0);
    auto rw = U(i,3,0);
    auto re = U(i,4,0);
    auto p = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    // pressure gradient
    G(i,0,0) += bpint[b*3+0] * p;
    G(i,1,0) += bpint[b*3+1] * p;
    G(i,2,0) += bpint[b*3+2] * p;
  }

  // boundary edge contributions: triangle superedges
  for (std::size_t e=0; e<bsupedge[0].size()/6; ++e) {
    const auto N = bsupedge[0].data() + e*6;
    tk::real p[3];
    for (std::size_t a=0; a<3; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }
    // pressure gradient
    for (std::size_t j=0; j<3; ++j) {
      tk::real f[3];
      const auto b = bsupint[0].data();
      f[0] = b[(e*3+0)*3+j] * (p[1] + p[0]);
      f[1] = b[(e*3+1)*3+j] * (p[2] + p[1]);
      f[2] = b[(e*3+2)*3+j] * (p[0] + p[2]);
      G(N[0],j,0) = G(N[0],j,0) - f[0] + f[2];
      G(N[1],j,0) = G(N[1],j,0) + f[0] - f[1];
      G(N[2],j,0) = G(N[2],j,0) + f[1] - f[2];
    }
  }

  // boundary edge contributions: edges
  for (std::size_t e=0; e<bsupedge[1].size()/4; ++e) {
    const auto N = bsupedge[1].data() + e*4;
    tk::real p[2];
    for (std::size_t a=0; a<2; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }
    // pressure gradient
    const auto b = bsupint[1].data() + e*3;
    for (std::size_t j=0; j<3; ++j) {
      tk::real f = b[j] * (p[1] + p[0]);
      G(N[0],j,0) -= f;
      G(N[1],j,0) += f;
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
stab( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const tk::Fields& U,
      const tk::Fields& G,
      tk::Fields& S )
// *****************************************************************************
//  Compute stabilization coefficient in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] coord Mesh node coordinates
//! \param[in] U Unknowns/solution vector in mesh nodes
//! \param[in] G Nodal gradients
//! \param[in,out] S Stabilization coefficients at nodes
//! \return Gradients in all mesh points
// *****************************************************************************
{
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

  using std::abs;
  using std::max;

  Assert( S.nunk() == U.nunk(), "Size mismatch" );
  Assert( S.nprop() == 1, "Size mismatch" );
  S.fill( 0.0 );

  // Lambda to compute the pressure sensor function on edge {p,q}
  const auto eps = std::numeric_limits< tk::real >::epsilon()*1.0e+3;
  auto beta = [&]( std::size_t p, std::size_t q,
                   tk::real pL, tk::real pR,
                   tk::real dpL[3], tk::real dpR[3] )
  {
    auto dx = x[p] - x[q];
    auto dy = y[p] - y[q];
    auto dz = z[p] - z[q];
    auto dpx = dpL[0] + dpR[0];
    auto dpy = dpL[1] + dpR[1];
    auto dpz = dpL[2] + dpR[2];
    auto dp = 0.5*(dx*dpx + dy*dpy + dz*dpz);
    auto ps = pL - pR;
    auto e = abs(ps - dp);
    auto d = abs(ps) + abs(dp) + eps;
    return 1.0 - e/d;
  };

  // domain edge contributions: tetrahedron superedges
  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
    const auto N = dsupedge[0].data() + e*4;
    tk::real p[4], dp[4][3];
    for (std::size_t a=0; a<4; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
      dp[a][0] = G(N[a],0,0);
      dp[a][1] = G(N[a],1,0);
      dp[a][2] = G(N[a],2,0);
    }
    tk::real b[6];  // pressure sensor function in edges
    b[0] = beta( 0, 1, p[0], p[1], dp[0], dp[1] );
    b[1] = beta( 1, 2, p[1], p[2], dp[1], dp[2] );
    b[2] = beta( 2, 0, p[2], p[0], dp[2], dp[0] );
    b[3] = beta( 0, 3, p[0], p[3], dp[0], dp[3] );
    b[4] = beta( 1, 3, p[1], p[3], dp[1], dp[3] );
    b[5] = beta( 2, 3, p[2], p[3], dp[2], dp[3] );
    // store in points max of connecting edges
    S(N[0],0,0) = max( S(N[0],0,0), max(b[0],max(b[2],b[3])) );
    S(N[1],0,0) = max( S(N[1],0,0), max(b[0],max(b[1],b[4])) );
    S(N[2],0,0) = max( S(N[2],0,0), max(b[1],max(b[2],b[5])) );
    S(N[3],0,0) = max( S(N[3],0,0), max(b[3],max(b[4],b[5])) );
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    tk::real p[3], dp[3][3];
    for (std::size_t a=0; a<3; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
      dp[a][0] = G(N[a],0,0);
      dp[a][1] = G(N[a],1,0);
      dp[a][2] = G(N[a],2,0);
    }
    tk::real b[3];  // pressure sensor function in edges
    b[0] = beta( 0, 1, p[0], p[1], dp[0], dp[1] );
    b[1] = beta( 1, 2, p[1], p[2], dp[1], dp[2] );
    b[2] = beta( 2, 0, p[2], p[0], dp[2], dp[0] );
    // store in points max of connecting edges
    S(N[0],0,0) = max( S(N[0],0,0), max(b[0],b[2]) );
    S(N[1],0,0) = max( S(N[1],0,0), max(b[0],b[1]) );
    S(N[2],0,0) = max( S(N[2],0,0), max(b[1],b[2]) );
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real p[2], dp[2][3];
    for (std::size_t a=0; a<2; ++a) {
      auto r  = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      auto re = U(N[a],4,0);
      p[a] = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
      dp[a][0] = G(N[a],0,0);
      dp[a][1] = G(N[a],1,0);
      dp[a][2] = G(N[a],2,0);
    }
    // pressure sensor function in edge
    auto b = beta( 0, 1, p[0], p[1], dp[0], dp[1] );
    // store in points max of connecting edges
    S(N[0],0,0) = max( S(N[0],0,0), b );
    S(N[1],0,0) = max( S(N[1],0,0), b );
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
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     const std::vector< tk::real >& dtp,
     const tk::Fields& U,
     const tk::Fields& S,
     const tk::Fields& G,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] t Physical time
//! \param[in] dt Physical time size
//! \param[in] tp Phisical time step size for each mesh node (if steady state)
//! \param[in] dtp Time step size for each mesh node (if steady state)
//! \param[in] U Unknowns/solution vector in mesh nodes
//! \param[in] S Stabilization coefficients
//! \param[in] G Nodal gradients
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  R.fill( 0.0 );
  advdom( dsupedge, dsupint, coord, t, dt, tp, dtp, U, S, G, R );
  advbnd( triinpoel, coord, besym, U, R );
}

} // zalesak::
