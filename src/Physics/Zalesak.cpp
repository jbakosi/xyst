// *****************************************************************************
/*!
  \file      src/Physics/Zalesak.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
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
         const std::array< std::vector< tk::real >, 3 >& coord,
         tk::real dt,
         tk::real t,
         std::size_t p,
         std::size_t q,
         tk::real f[],
         const std::function< std::vector< tk::real >
                 ( tk::real, tk::real, tk::real, tk::real ) >& src )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] supint Edge integral
//! \param[in] U Solution vector to read conserved variables from
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time step size
//! \param[in] t Physical time
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] f Flux computed
//! \param[in] src Function to call to evaluate a problem-sepcific source term
// *****************************************************************************
{
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

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  // Taylor-Galerkin first half step

  tk::real ue[ncomp];
  auto dp = pL - pR;
  ue[0] = 0.5*(rL + rR - dt*(rL*dnL - rR*dnR));
  ue[1] = 0.5*(ruL + ruR - dt*(ruL*dnL - ruR*dnR + dp*dx));
  ue[2] = 0.5*(rvL + rvR - dt*(rvL*dnL - rvR*dnR + dp*dy));
  ue[3] = 0.5*(rwL + rwR - dt*(rwL*dnL - rwR*dnR + dp*dz));
  ue[4] = 0.5*(reL + reR - dt*((reL+pL)*dnL - (reR+pR)*dnR));

  for (std::size_t c=5; c<ncomp; ++c) {
    ue[c] = 0.5*(U(p,c,0) + U(q,c,0) - dt*(U(p,c,0)*dnL - U(q,c,0)*dnR));
  }

  if (src) {
    auto coef = dt/4.0;
    auto sL = src( x[p], y[p], z[p], t );
    auto sR = src( x[q], y[q], z[q], t );
    for (std::size_t c=0; c<ncomp; ++c) {
      ue[c] += coef*(sL[c] + sR[c]);
    }
  }

  // Taylor-Galerkin second half step

  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  auto rh  = ue[0];
  auto ruh = ue[1];
  auto rvh = ue[2];
  auto rwh = ue[3];
  auto reh = ue[4];
  auto ph = eos::pressure( reh - 0.5*(ruh*ruh + rvh*rvh + rwh*rwh)/rh );
  auto vn = (ruh*nx + rvh*ny + rwh*nz)/rh;
  f[0] = dt*rh*vn;
  f[1] = dt*(ruh*vn + ph*nx);
  f[2] = dt*(rvh*vn + ph*ny);
  f[3] = dt*(rwh*vn + ph*nz);
  f[4] = dt*(reh + ph)*vn;

  for (std::size_t c=5; c<ncomp; ++c) {
    f[c] = dt*ue[c]*vn;
  }

  if (src) {
    auto coef = -5.0/3.0*dt*supint[3];
    auto xe = (x[p] + x[q])/2.0;
    auto ye = (y[p] + y[q])/2.0;
    auto ze = (z[p] + z[q])/2.0;
    auto se = src( xe, ye, ze, t+dt/2.0 );
    for (std::size_t c=0; c<ncomp; ++c) {
      f[ncomp+c] = coef*se[c];
    }
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
        tk::real dt,
        tk::real t,
        const tk::Fields& U,
        // cppcheck-suppress constParameter
        tk::Fields& R )
// *****************************************************************************
//! Compute domain-edge integrals for advection
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time step size
//! \param[in] t Physical time
//! \param[in] U Solution vector at recent time step
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
    advedge( d+(e*6+0)*4, U, coord, dt, t, N[0], N[1], f[0], src );
    advedge( d+(e*6+1)*4, U, coord, dt, t, N[1], N[2], f[1], src );
    advedge( d+(e*6+2)*4, U, coord, dt, t, N[2], N[0], f[2], src );
    advedge( d+(e*6+3)*4, U, coord, dt, t, N[0], N[3], f[3], src );
    advedge( d+(e*6+4)*4, U, coord, dt, t, N[1], N[3], f[4], src );
    advedge( d+(e*6+5)*4, U, coord, dt, t, N[2], N[3], f[5], src );
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
    advedge( d+(e*3+0)*4, U, coord, dt, t, N[0], N[1], f[0], src );
    advedge( d+(e*3+1)*4, U, coord, dt, t, N[1], N[2], f[1], src );
    advedge( d+(e*3+2)*4, U, coord, dt, t, N[2], N[0], f[2], src );
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
    advedge( d+e*4, U, coord, dt, t, N[0], N[1], f, src );
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
        tk::real dt,
        const tk::Fields& U,
        tk::Fields& R )
// *****************************************************************************
//! Compute boundary integrals for advection
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] coord Mesh node coordinates
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] dt Physical time size
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector
// *****************************************************************************
{
  auto ncomp = U.nprop();

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

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

    tk::real p, vn, f[5][3];
    const auto sym = besym.data() + e*3;

    p = eos::pressure( reA - 0.5*(ruA*ruA + rvA*rvA + rwA*rwA)/rA );
    vn = sym[0] ? 0.0 : (nx*ruA + ny*rvA + nz*rwA)/rA;
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + p*nx;
    f[2][0] = rvA*vn + p*ny;
    f[3][0] = rwA*vn + p*nz;
    f[4][0] = (reA + p)*vn;

    p = eos::pressure( reB - 0.5*(ruB*ruB + rvB*rvB + rwB*rwB)/rB );
    vn = sym[1] ? 0.0 : (nx*ruB + ny*rvB + nz*rwB)/rB;
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + p*nx;
    f[2][1] = rvB*vn + p*ny;
    f[3][1] = rwB*vn + p*nz;
    f[4][1] = (reB + p)*vn;

    p = eos::pressure( reC - 0.5*(ruC*ruC + rvC*rvC + rwC*rwC)/rC );
    vn = sym[2] ? 0.0 : (nx*ruC + ny*rvC + nz*rwC)/rC;
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + p*nx;
    f[2][2] = rvC*vn + p*ny;
    f[3][2] = rwC*vn + p*nz;
    f[4][2] = (reC + p)*vn;

    for (std::size_t c=0; c<ncomp; ++c) {
      auto fab = (f[c][0] + f[c][1])/4.0;
      auto fbc = (f[c][1] + f[c][2])/4.0;
      auto fca = (f[c][2] + f[c][0])/4.0;
      R(N[0],c,0) += dt*(fab + fca + f[c][0]);
      R(N[1],c,0) += dt*(fab + fbc + f[c][1]);
      R(N[2],c,0) += dt*(fbc + fca + f[c][2]);
    }
  }
}

void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::uint8_t >& besym,
     const tk::Fields& U,
     tk::real t,
     tk::real dt,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Unknowns/solution vector in mesh nodes
//! \param[in] besym Boundary element symmetry BC flags
//! \param[in] t Physical time
//! \param[in] dt Physical time size
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

  // advection + src
  advdom( dsupedge, dsupint, coord, dt, t, U, R );
  advbnd( triinpoel, coord, besym, dt, U, R );
}

} // zalesak::