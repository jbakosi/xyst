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
         std::size_t p,
         std::size_t q,
         tk::real f[],
         std::size_t symL = 0,
         std::size_t symR = 0 )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] supint Edge integral
//! \param[in] U Solution vector to read conserved variables from
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time size
//! \param[in] p Left node id of edge-end
//! \param[in] q Right node id of edge-end
//! \param[in,out] f Flux computed
//! \param[in] symL Non-zero if left edge end-point is on a symmetry boundary
//! \param[in] symR Non-zero if right edge end-point is on a symmetry boundary
// *****************************************************************************
{
  auto ncomp = U.nprop();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // internal energy
  auto el = ( U(p,4,0) - 0.5 * ( U(p,1,0)*U(p,1,0)
                               + U(p,2,0)*U(p,2,0)
                               + U(p,3,0)*U(p,3,0) ) ) / U(p,0,0);
  auto er = ( U(q,4,0) - 0.5 * ( U(q,1,0)*U(q,1,0)
                               + U(q,2,0)*U(q,2,0)
                               + U(q,3,0)*U(q,3,0) ) ) / U(q,0,0);

  // pressure
  auto pL = eos::pressure( U(p,0,0), el );
  auto pR = eos::pressure( U(q,0,0), er );

  // edge vector
  auto dx = x[p] - x[q];
  auto dy = y[p] - y[q];
  auto dz = z[p] - z[q] ;
  auto dl = dx*dx + dy*dy + dz*dz;
  auto dnL = symL ? 0.0 :
             (dx*U(p,1,0) + dy*U(p,2,0) + dz*U(p,3,0)) / U(p,0,0) / dl;
  auto dnR = symR ? 0.0 :
             (dx*U(q,1,0) + dy*U(q,2,0) + dz*U(q,3,0)) / U(q,0,0) / dl;

  // integral coefficient
  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  // flow fluxes
  auto dp = pR - pL;
  auto rh  = 0.5*(U(p,0,0) + U(q,0,0)) + 0.5*dt*(U(q,0,0)*dnR - U(p,0,0)*dnL);
  auto ruh = 0.5*(U(p,1,0) + U(q,1,0))
           + 0.5*dt*(U(q,1,0)*dnR - U(p,1,0)*dnL + dp*dx/dl);
  auto rvh = 0.5*(U(p,2,0) + U(q,2,0))
           + 0.5*dt*(U(q,2,0)*dnR - U(p,2,0)*dnL + dp*dy/dl);
  auto rwh = 0.5*(U(p,3,0) + U(q,3,0))
           + 0.5*dt*(U(q,3,0)*dnR - U(p,3,0)*dnL + dp*dz/dl);
  auto reh = 0.5*(U(p,4,0) + U(q,4,0))
           + 0.5*dt*((U(q,4,0) + pR)*dnR - (U(p,4,0)+pL)*dnL);
  auto ph = eos::pressure( rh, (reh-0.5*(ruh*ruh + rvh*rvh + rwh*rwh)/rh)/rh );
  auto vn = (ruh*nx + rvh*ny + rwh*nz)/rh;
  f[0] = 2.0*rh*vn;
  f[1] = 2.0*(ruh*vn + ph*nx);
  f[2] = 2.0*(rvh*vn + ph*ny);
  f[3] = 2.0*(rwh*vn + ph*nz);
  f[4] = 2.0*(reh + ph)*vn;

  // scalar fluxes
  for (std::size_t c=5; c<ncomp; ++c) {
    auto hc = 0.5*(U(p,c,0) + U(q,c,0)) + 0.5*dt*(U(q,c,0)*dnR - U(p,c,0)*dnL);
    f[c] = 2.0*hc*vn;
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
  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
    const auto N = dsupedge[0].data() + e*4;
    // edge fluxes
    tk::real f[6][ncomp];
    const auto d = dsupint[0].data();
    advedge( d+(e*6+0)*4, U, coord, dt, N[0], N[1], f[0] );
    advedge( d+(e*6+1)*4, U, coord, dt, N[1], N[2], f[1] );
    advedge( d+(e*6+2)*4, U, coord, dt, N[2], N[0], f[2] );
    advedge( d+(e*6+3)*4, U, coord, dt, N[0], N[3], f[3] );
    advedge( d+(e*6+4)*4, U, coord, dt, N[1], N[3], f[4] );
    advedge( d+(e*6+5)*4, U, coord, dt, N[2], N[3], f[5] );
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
    // edge fluxes
    tk::real f[3][ncomp];
    const auto d = dsupint[1].data();
    advedge( d+(e*3+0)*4, U, coord, dt, N[0], N[1], f[0] );
    advedge( d+(e*3+1)*4, U, coord, dt, N[1], N[2], f[1] );
    advedge( d+(e*3+2)*4, U, coord, dt, N[2], N[0], f[2] );
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
    // edge fluxes
    tk::real f[ncomp];
    const auto d = dsupint[2].data();
    advedge( d+e*4, U, coord, dt, N[0], N[1], f );
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
    auto e = ( U(p,4,0) - 0.5 * ( U(p,1,0)*U(p,1,0)
                                + U(p,2,0)*U(p,2,0)
                                + U(p,3,0)*U(p,3,0) ) ) / U(p,0,0);
    auto pr = eos::pressure( U(p,0,0), e );
    // boundary-normal velocity
    auto nx = bpint[b*3+0];
    auto ny = bpint[b*3+1];
    auto nz = bpint[b*3+2];
    auto vn = bpsym[b] ? 0.0 :
              (nx*U(p,1,0) + ny*U(p,2,0) + nz*U(p,3,0)) / U(p,0,0 );
    // flow fluxes
    R(p,0,0) += U(p,0,0)*vn;
    R(p,1,0) += U(p,1,0)*vn + pr*nx;
    R(p,2,0) += U(p,2,0)*vn + pr*ny;
    R(p,3,0) += U(p,3,0)*vn + pr*nz;
    R(p,4,0) += (U(p,4,0) + pr)*vn;
    // scalar fluxes
    for (std::size_t c=5; c<ncomp; ++c) {
      R(p,c,0) += U(p,c,0) * vn;
    }
  }

  // boundary edge contributions: triangle superedges
  for (std::size_t e=0; e<bsupedge[0].size()/6; ++e) {
    const auto N = bsupedge[0].data() + e*6;
    // edge fluxes
    tk::real f[3][ncomp];
    const auto b = bsupint[0].data();
    advedge( b+(e*3+0)*3, U, coord, dt, N[0], N[1], f[0], N[3], N[4] );
    advedge( b+(e*3+1)*3, U, coord, dt, N[1], N[2], f[1], N[4], N[5] );
    advedge( b+(e*3+2)*3, U, coord, dt, N[2], N[0], f[2], N[5], N[3] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) += f[0][c] + f[2][c];
      R(N[1],c,0) += f[0][c] + f[1][c];
      R(N[2],c,0) += f[1][c] + f[2][c];
    }
  }

  // boundary edge contributions: edges
  for (std::size_t e=0; e<bsupedge[1].size()/4; ++e) {
    const auto N = bsupedge[1].data() + e*4;
    // edge fluxes
    tk::real f[ncomp];
    const auto b = bsupint[1].data();
    advedge( b+e*3, U, coord, dt, N[0], N[1], f, N[2], N[3] );
    // edge flux contributions
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c,0) += f[c];
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
     tk::real t,
     tk::real dt,
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
  adv( bpoin, bpint, bpsym, dsupedge, dsupint, bsupedge, bsupint, coord, dt,
       U, R );

  // source
  src( coord, v, t, tp, R );
}

} // zalesak::
