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
advedge( std::size_t ncomp,
         const tk::real dsupint[],
         tk::real dt,
         const tk::real L[],
         const tk::real R[],
         tk::real f[],
         std::size_t symL = 0,
         std::size_t symR = 0 )
// *****************************************************************************
//! Compute advection fluxes on a single edge
//! \param[in] ncomp Number of scalar components to solve for
//! \param[in] dsupint Domain superedge integral for this edge
//! \param[in] dt Physical time size
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

  // will work with copies of physics variables
  tk::real l[ncomp], r[ncomp];
  memcpy( l, L, sizeof l );
  memcpy( r, R, sizeof r );

  // pressure
  auto pL = eos::pressure( l[0], l[4] );
  auto pR = eos::pressure( r[0], r[4] );

  // dualface-normal velocities
  auto nx = dsupint[0];
  auto ny = dsupint[1];
  auto nz = dsupint[2];
  auto vnL = symL ? 0.0 : (l[1]*nx + l[2]*ny + l[3]*nz);
  auto vnR = symR ? 0.0 : (r[1]*nx + r[2]*ny + r[3]*nz);
  auto len = tk::length( nx, ny, nz );

  // back to conserved variables
  l[4] = (l[4] + 0.5*(l[1]*l[1] + l[2]*l[2] + l[3]*l[3])) * l[0];
  l[1] *= l[0];
  l[2] *= l[0];
  l[3] *= l[0];
  r[4] = (r[4] + 0.5*(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])) * r[0];
  r[1] *= r[0];
  r[2] *= r[0];
  r[3] *= r[0];

  // flow fluxes
  tk::real h[ncomp];
  h[0] = 0.5*(l[0] + r[0] + dt/len*(r[0]*vnR - l[0]*vnL));
  h[1] = 0.5*(l[1] + r[1] + dt/len*(r[1]*vnR - l[1]*vnL + (pR-pL)*nx));
  h[2] = 0.5*(l[2] + r[2] + dt/len*(r[2]*vnR - l[2]*vnL + (pR-pL)*ny));
  h[3] = 0.5*(l[3] + r[3] + dt/len*(r[3]*vnR - l[3]*vnL + (pR-pL)*nz));
  h[4] = 0.5*(l[4] + r[4] + dt/len*((r[4] + pR)*vnR - (l[4] + pL)*vnL));
  auto rh = h[0];
  auto uh = h[1] / rh;
  auto vh = h[2] / rh;
  auto wh = h[3] / rh;
  auto eh = h[4] / rh - 0.5*(uh*uh + vh*vh + wh*wh);
  auto ph = eos::pressure( rh, eh );
  auto vnh = uh*nx + vh*ny + wh*nz;
  f[0] = 2.0*h[0]*vnh;
  f[1] = 2.0*(h[1]*vnh + ph*nx);
  f[2] = 2.0*(h[2]*vnh + ph*ny);
  f[3] = 2.0*(h[3]*vnh + ph*nz);
  f[4] = 2.0*(h[4] + ph)*vnh;

  if (ncomp == 5) return;

  // scalar fluxes
  for (std::size_t c=5; c<ncomp; ++c) {
    tk::real hc = 0.5*(l[c] + r[c] + dt/len*(r[c]*vnR - l[c]*vnL));
    f[c] = 2.0*hc*vnh;
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
adv( const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::uint8_t >& bpsym,
     const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< std::size_t >, 2 >& bsupedge,
     const std::array< std::vector< tk::real >, 2 >& bsupint,
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
//! \param[in] dt Physical time size
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
    tk::real u[4][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    primitive( ncomp, N[3], U, u[3] );
    // edge fluxes
    tk::real f[6][ncomp];
    const auto d = dsupint[0].data();
    advedge( ncomp, d+(e*6+0)*3, dt, u[0], u[1], f[0] );
    advedge( ncomp, d+(e*6+1)*3, dt, u[1], u[2], f[1] );
    advedge( ncomp, d+(e*6+2)*3, dt, u[2], u[0], f[2] );
    advedge( ncomp, d+(e*6+3)*3, dt, u[0], u[3], f[3] );
    advedge( ncomp, d+(e*6+4)*3, dt, u[1], u[3], f[4] );
    advedge( ncomp, d+(e*6+5)*3, dt, u[2], u[3], f[5] );
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
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    // edge fluxes
    tk::real f[3][ncomp];
    const auto d = dsupint[1].data();
    advedge( ncomp, d+(e*3+0)*3, dt, u[0], u[1], f[0] );
    advedge( ncomp, d+(e*3+1)*3, dt, u[1], u[2], f[1] );
    advedge( ncomp, d+(e*3+2)*3, dt, u[2], u[0], f[2] );
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
    const auto d = dsupint[2].data();
    advedge( ncomp, d+e*3, dt, u[0], u[1], f );
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
    tk::real u[3][ncomp];
    primitive( ncomp, N[0], U, u[0] );
    primitive( ncomp, N[1], U, u[1] );
    primitive( ncomp, N[2], U, u[2] );
    // edge fluxes
    tk::real f[3][ncomp];
    const auto b = bsupint[0].data();
    advedge( ncomp, b+(e*3+0)*3, dt, u[0], u[1], f[0], N[3], N[4] );
    advedge( ncomp, b+(e*3+1)*3, dt, u[1], u[2], f[1], N[4], N[5] );
    advedge( ncomp, b+(e*3+2)*3, dt, u[2], u[0], f[2], N[5], N[3] );
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
    advedge( ncomp, b+e*3, dt, u[0], u[1], f, N[2], N[3] );
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
  adv( bpoin, bpint, bpsym, dsupedge, dsupint, bsupedge, bsupint, dt, U, R );

  // source
  src( coord, v, t, tp, R );
}

} // zalesak::
