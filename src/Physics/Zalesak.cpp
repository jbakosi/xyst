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

#include <iostream>     // NOT NEEDED
#include <iomanip>     // NOT NEEDED

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
         const std::vector< tk::real >& vol,
         const std::size_t N[],
         tk::real f[] )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] supint Edge integral
//! \param[in] U Solution vector to read conserved variables from
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time size
//! \param[in,out] f Flux computed
// *****************************************************************************
{
  const auto ncomp = U.nprop();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto p = N[0];
  auto q = N[1];
  auto volL = 1.0;//vol[p]/4.0;
  auto volR = 1.0;//vol[q]/4.0;

  // integral coefficient
  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  // edge vector
  auto dx = supint[3];//x[p] - x[q];
  auto dy = supint[4];//y[p] - y[q];
  auto dz = supint[5];//z[p] - z[q];
  //auto dl = dx*dx + dy*dy + dz*dz;
  //dx /= supint[3];//dl;
  //dy /= supint[3];//dl;
  //dz /= supint[3];//dl;

  // left state
  auto rL  = U(p,0,0);
  auto ruL = U(p,1,0);
  auto rvL = U(p,2,0);
  auto rwL = U(p,3,0);
  auto reL = U(p,4,0);
  auto pL = eos::pressure( reL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)/rL );
  auto dnL = (ruL*dx + rvL*dy + rwL*dz)/rL/volL;

  // right state
  auto rR  = U(q,0,0);
  auto ruR = U(q,1,0);
  auto rvR = U(q,2,0);
  auto rwR = U(q,3,0);
  auto reR = U(q,4,0);
  auto pR = eos::pressure( reR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)/rR );
  auto dnR = (ruR*dx + rvR*dy + rwR*dz)/rR/volR;

  // flow fluxes
  auto dp = pL/volL - pR/volR;
  auto rh  = 0.5*(rL + rR - dt*(rL*dnL - rR*dnR));
  auto ruh = 0.5*(ruL + ruR - dt*(ruL*dnL - ruR*dnR + dp*dx));
  auto rvh = 0.5*(rvL + rvR - dt*(rvL*dnL - rvR*dnR + dp*dy));
  auto rwh = 0.5*(rwL + rwR - dt*(rwL*dnL - rwR*dnR + dp*dz));
  auto reh = 0.5*(reL + reR - dt*((reL+pL)*dnL - (reR+pR)*dnR));
  auto ph = eos::pressure( reh - 0.5*(ruh*ruh + rvh*rvh + rwh*rwh)/rh );
  auto vn = (ruh*nx + rvh*ny + rwh*nz)/rh;
  f[0] = rh*vn;
  f[1] = ruh*vn + ph*nx;
  f[2] = rvh*vn + ph*ny;
  f[3] = rwh*vn + ph*nz;
  f[4] = (reh + ph)*vn;

//if (p==13994 || q==13994) std::cout << "ed: " << dx << ", " << dy << ", " << dz << ": " << dnL << ", " << dnR << " > " << f[0] << '\n';

  // scalar fluxes
  for (std::size_t c=5; c<ncomp; ++c) {
    auto hc = 0.5*(U(p,c,0) + U(q,c,0) - dt*(U(p,c,0)*dnL - U(q,c,0)*dnR));
    f[c] = hc*vn;
  }
}

static void
adv( const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::uint8_t >& bpsym,
     const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< std::size_t >, 2 >& bsupedge,
     const std::array< std::vector< tk::real >, 2 >& bsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< tk::real >& vol,
     tk::real dt,
     const tk::Fields& U,
     // cppcheck-suppress constParameter
     tk::Fields& R )
// *****************************************************************************
//! Compute integrals for advection
//! \param[in] bpoin Boundary point local ids
//! \param[in] bpint Boundary point integrals
//! \param[in] bpsym Boundary point symmetry BC flags
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] bsupedge Boundary superedges
//! \param[in] bsupint Boundary superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time size
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
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
  //for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
  //  const auto N = dsupedge[0].data() + e*4;
  //  // edge fluxes
  //  tk::real f[6][ncomp];
  //  const auto d = dsupint[0].data();
  //  advedge( d+(e*6+0)*4, U, coord, dt, N[0], N[1], f[0] );
  //  advedge( d+(e*6+1)*4, U, coord, dt, N[1], N[2], f[1] );
  //  advedge( d+(e*6+2)*4, U, coord, dt, N[2], N[0], f[2] );
  //  advedge( d+(e*6+3)*4, U, coord, dt, N[0], N[3], f[3] );
  //  advedge( d+(e*6+4)*4, U, coord, dt, N[1], N[3], f[4] );
  //  advedge( d+(e*6+5)*4, U, coord, dt, N[2], N[3], f[5] );
  //  // edge flux contributions
  //  for (std::size_t c=0; c<ncomp; ++c) {
  //    R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c] - f[3][c];
  //    R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c] - f[4][c];
  //    R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c] - f[5][c];
  //    R(N[3],c,0) = R(N[3],c,0) + f[3][c] + f[4][c] + f[5][c];
  //  }
  //}

  //// domain edge contributions: triangle superedges
  //for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
  //  const auto N = dsupedge[1].data() + e*3;
  //  // edge fluxes
  //  tk::real f[3][ncomp];
  //  const auto d = dsupint[1].data();
  //  advedge( d+(e*3+0)*4, U, coord, dt, N[0], N[1], f[0] );
  //  advedge( d+(e*3+1)*4, U, coord, dt, N[1], N[2], f[1] );
  //  advedge( d+(e*3+2)*4, U, coord, dt, N[2], N[0], f[2] );
  //  // edge flux contributions
  //  for (std::size_t c=0; c<ncomp; ++c) {
  //    R(N[0],c,0) = R(N[0],c,0) - f[0][c] + f[2][c];
  //    R(N[1],c,0) = R(N[1],c,0) + f[0][c] - f[1][c];
  //    R(N[2],c,0) = R(N[2],c,0) + f[1][c] - f[2][c];
  //  }
  //}

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    // edge fluxes
    tk::real f[ncomp];
    const auto d = dsupint[2].data();
    advedge( d+e*4, U, coord, dt, vol, N, f );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
//if (c==0 && N[0]==13994) std::cout << "d: v: " << R(N[0],c,0) << ", -: " << f[c] << '\n';
//if (c==0 && N[1]==13994) std::cout << "d: v: " << R(N[1],c,0) << ", +: " << f[c] << '\n';
      R(N[0],c,0) -= f[c];
      R(N[1],c,0) += f[c];
    }
  }

  // boundary integrals

  // boundary point contributions
  for (std::size_t b=0; b<bpoin.size(); ++b) {
    auto i = bpoin[b];
    auto r = U(i,0,0);
    auto ru = U(i,1,0);
    auto rv = U(i,2,0);
    auto rw = U(i,3,0);
    auto re = U(i,4,0);
    auto p = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    const auto n = bpint.data() + b*3;
    auto nx = n[0];
    auto ny = n[1];
    auto nz = n[2];
    auto vn = bpsym[b] ? 0.0 : (ru*nx + rv*ny + rw*nz)/r;
//if (i==13994) std::cout << "p: + " << U(i,0,0)*vn << '\n';
    R(i,0,0) += r*vn;
    R(i,1,0) += ru*vn + p*nx;
    R(i,2,0) += rv*vn + p*ny;
    R(i,3,0) += rw*vn + p*nz;
    R(i,4,0) += (re + p)*vn;
    for (std::size_t c=5; c<ncomp; ++c) {
      R(i,c,0) += U(i,c,0)*vn;
    }
  }
 
  //// boundary edge contributions: triangle superedges
  //for (std::size_t e=0; e<bsupedge[0].size()/6; ++e) {
  //  const auto N = bsupedge[0].data() + e*6;
  //  // edge fluxes
  //  tk::real f[3][ncomp];
  //  const auto b = bsupint[0].data();
  //  advedge( b+(e*3+0)*3, U, coord, dt, N[0], N[1], f[0], N[3], N[4] );
  //  advedge( b+(e*3+1)*3, U, coord, dt, N[1], N[2], f[1], N[4], N[5] );
  //  advedge( b+(e*3+2)*3, U, coord, dt, N[2], N[0], f[2], N[5], N[3] );
  //  // edge flux contributions
  //  for (std::size_t c=0; c<ncomp; ++c) {
  //    R(N[0],c,0) += f[0][c] + f[2][c];
  //    R(N[1],c,0) += f[0][c] + f[1][c];
  //    R(N[2],c,0) += f[1][c] + f[2][c];
  //  }
  //}
 
  // boundary edge contributions: edges
  for (std::size_t e=0; e<bsupedge[1].size()/4; ++e) {
    const auto N = bsupedge[1].data() + e*4;
    const auto n = bsupint[1].data() + e*3;
    auto nx = n[0];
    auto ny = n[1];
    auto nz = n[2];
 
    auto p = N[0];
    auto rL  = U(p,0,0);
    auto ruL = U(p,1,0);
    auto rvL = U(p,2,0);
    auto rwL = U(p,3,0);
    auto reL = U(p,4,0);
    auto pL = eos::pressure( reL - 0.5*(ruL*ruL + rvL*rvL + rwL*rwL)/rL );
    auto vnL = N[2] ? 0.0 : (nx*ruL + ny*rvL + nz*rwL)/rL;
 
    auto q = N[1];
    auto rR  = U(q,0,0);
    auto ruR = U(q,1,0);
    auto rvR = U(q,2,0);
    auto rwR = U(q,3,0);
    auto reR = U(q,4,0);
    auto pR = eos::pressure( reR - 0.5*(ruR*ruR + rvR*rvR + rwR*rwR)/rR );
    auto vnR = N[3] ? 0.0 : (nx*ruR + ny*rvR + nz*rwR)/rR;
 
    tk::real f[ncomp];
    auto p2 = pL + pR;
    f[0] = rL*vnL + rR*vnR;
    f[1] = ruL*vnL + ruR*vnR + p2*nx;
    f[2] = rvL*vnL + rvR*vnR + p2*ny;
    f[3] = rwL*vnL + rwR*vnR + p2*nz;
    f[4] = (reL + pL)*vnL + (reR + pR)*vnR;
 
    for (std::size_t c=5; c<ncomp; ++c) {
      f[c] = U(p,c,0)*vnL + U(q,c,0)*vnR;
    }

    for (std::size_t c=0; c<ncomp; ++c) {
//if (c==0 && N[0]==13994) std::cout << "b: - " << f[c] << '\n';
//if (c==0 && N[1]==13994) std::cout << "b: + " << f[c] << '\n';
      R(p,c,0) -= f[c];
      R(q,c,0) += f[c];
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
        tk::real dt,
        const tk::Fields& U,
        tk::Fields& R )
// *****************************************************************************
// *****************************************************************************
{
  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t t=0; t<triinpoel.size()/3; ++t) {
    std::size_t N[3] =
      { triinpoel[t*3+0], triinpoel[t*3+1], triinpoel[t*3+2] };

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
    nx /= 6.0;
    ny /= 6.0;
    nz /= 6.0;

    // boundary flux
    tk::real u, v, w, p, vn;
    u = ruA / rA;
    v = rvA / rA;
    w = rwA / rA;
    p = eos::pressure( rA, reA / rA - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    auto flu11 = rA*vn;
    auto flu21 = ruA*vn + p*nx;
    auto flu31 = rvA*vn + p*ny;
    auto flu41 = rwA*vn + p*nz;
    auto flu51 = (reA + p)*vn;
    u = ruB / rB;
    v = rvB / rB;
    w = rwB / rB;
    p = eos::pressure( rB, reB / rB - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    auto flu12 = rB*vn;
    auto flu22 = ruB*vn + p*nx;
    auto flu32 = rvB*vn + p*ny;
    auto flu42 = rwB*vn + p*nz;
    auto flu52 = (reB + p)*vn;
    u = ruC / rC;
    v = rvC / rC;
    w = rwC / rC;
    p = eos::pressure( rC, reC / rC - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    auto flu13 = rC*vn;
    auto flu23 = ruC*vn + p*nx;
    auto flu33 = rvC*vn + p*ny;
    auto flu43 = rwC*vn + p*nz;
    auto flu53 = (reC + p)*vn;

    auto flua1 = flu11 + flu12 + flu13;
    auto flua2 = flu21 + flu22 + flu23;
    auto flua3 = flu31 + flu32 + flu33;
    auto flua4 = flu41 + flu42 + flu43;
    auto flua5 = flu51 + flu52 + flu53;

    R(N[0],0,0) += (5.0*flu11 + flua1)/8.0;
    R(N[0],1,0) += (5.0*flu21 + flua2)/8.0;
    R(N[0],2,0) += (5.0*flu31 + flua3)/8.0;
    R(N[0],3,0) += (5.0*flu41 + flua4)/8.0;
    R(N[0],4,0) += (5.0*flu51 + flua5)/8.0;

    R(N[1],0,0) += (5.0*flu12 + flua1)/8.0;
    R(N[1],1,0) += (5.0*flu22 + flua2)/8.0;
    R(N[1],2,0) += (5.0*flu32 + flua3)/8.0;
    R(N[1],3,0) += (5.0*flu42 + flua4)/8.0;
    R(N[1],4,0) += (5.0*flu52 + flua5)/8.0;

    R(N[2],0,0) += (5.0*flu13 + flua1)/8.0;
    R(N[2],1,0) += (5.0*flu23 + flua2)/8.0;
    R(N[2],2,0) += (5.0*flu33 + flua3)/8.0;
    R(N[2],3,0) += (5.0*flu43 + flua4)/8.0;
    R(N[2],4,0) += (5.0*flu53 + flua5)/8.0;
  }
}

static void
advbnd2( const std::vector< std::size_t >& triinpoel,
         const std::array< std::vector< tk::real >, 3 >& coord,
         tk::real dt,
         const tk::Fields& U,
         tk::Fields& R )
// *****************************************************************************
// *****************************************************************************
{
  // access node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // boundary integrals: compute fluxes in edges
  std::vector< tk::real > bflux( triinpoel.size() * 5 * 2 );

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    std::size_t N[3] = { triinpoel[e*3+0], triinpoel[e*3+1], triinpoel[e*3+2] };
    // access solution at element nodes
    tk::real rA  = U(N[0],0,0);
    tk::real rB  = U(N[1],0,0);
    tk::real rC  = U(N[2],0,0);
    tk::real ruA = U(N[0],1,0);
    tk::real ruB = U(N[1],1,0);
    tk::real ruC = U(N[2],1,0);
    tk::real rvA = U(N[0],2,0);
    tk::real rvB = U(N[1],2,0);
    tk::real rvC = U(N[2],2,0);
    tk::real rwA = U(N[0],3,0);
    tk::real rwB = U(N[1],3,0);
    tk::real rwC = U(N[2],3,0);
    tk::real reA = U(N[0],4,0);
    tk::real reB = U(N[1],4,0);
    tk::real reC = U(N[2],4,0);

    const std::array< tk::real, 3 >
      ba{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] },
      ca{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] };
    auto [nx,ny,nz] = tk::cross( ba, ca );
    auto A = sqrt( nx*nx + ny*ny + nz*nz );
    nx /= A;
    ny /= A;
    nz /= A;
    A /= 2.0;

    tk::real u, v, w, p, vn, f[5][3];

    u = ruA / rA;
    v = rvA / rA;
    w = rwA / rA;
    p = eos::pressure( rA, reA / rA - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    f[0][0] = rA*vn;
    f[1][0] = ruA*vn + p*nx;
    f[2][0] = rvA*vn + p*ny;
    f[3][0] = rwA*vn + p*nz;
    f[4][0] = (reA + p)*vn;

    u = ruB / rB;
    v = rvB / rB;
    w = rwB / rB;
    p = eos::pressure( rB, reB / rB - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    f[0][1] = rB*vn;
    f[1][1] = ruB*vn + p*nx;
    f[2][1] = rvB*vn + p*ny;
    f[3][1] = rwB*vn + p*nz;
    f[4][1] = (reB + p)*vn;

    u = ruC / rC;
    v = rvC / rC;
    w = rwC / rC;
    p = eos::pressure( rC, reC / rC - 0.5*(u*u + v*v + w*w) );
    vn = nx*u + ny*v + nz*w;
    f[0][2] = rC*vn;
    f[1][2] = ruC*vn + p*nx;
    f[2][2] = rvC*vn + p*ny;
    f[3][2] = rwC*vn + p*nz;
    f[4][2] = (reC + p)*vn;

    // store flux in boundary elements
    for (std::size_t c=0; c<5; ++c) {
      auto eb = (e*5+c)*6;
      auto Bab = A/24 * (f[c][0] + f[c][1]);
      bflux[eb+0] = Bab + A/6 * f[c][0];
      bflux[eb+1] = Bab;
      Bab = A/24 * (f[c][1] + f[c][2]);
      bflux[eb+2] = Bab + A/6 * f[c][1];
      bflux[eb+3] = Bab;
      Bab = A/24 * (f[c][2] + f[c][0]);
      bflux[eb+4] = Bab + A/6 * f[c][2];
      bflux[eb+5] = Bab;
    }
  }

  // boundary integrals: sum flux contributions to points
  for (std::size_t e=0; e<triinpoel.size()/3; ++e)
    for (std::size_t c=0; c<5; ++c) {
      auto eb = (e*5+c)*6;
      R(triinpoel[e*3+0],c,0) += bflux[eb+0] + bflux[eb+5];
      R(triinpoel[e*3+1],c,0) += bflux[eb+1] + bflux[eb+2];
      R(triinpoel[e*3+2],c,0) += bflux[eb+3] + bflux[eb+4];
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
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     const std::vector< tk::real >& vol,
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     tk::Fields& R,
     const std::vector< std::size_t >& triinpoel )
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
//! \param[in] U Unknowns/solution vector in mesh nodes
//! \param[in] v Nodal mesh volumes without contributions from other chares
//! \param[in] t Physical time
//! \param[in] dt Physical time size
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
  //adv( bpoin, bpint, bpsym, dsupedge, dsupint, bsupedge, bsupint, coord, vol, dt,
  //     U, R );

  //advbnd( triinpoel, coord, dt, U, R );
  //advbnd2( triinpoel, coord, dt, U, R );

  // source
  src( coord, v, t, tp, R );
}

} // zalesak::
