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

static void
advdom( const tk::UnsMesh::Coords& /*coord*/,
        const std::array< std::vector< std::size_t >, 3 >& /*dsupedge*/,
        const std::array< std::vector< tk::real >, 3 >& /*dsupint*/,
        const tk::Fields& /*G*/,
        const tk::Fields& /*U*/,
        // cppcheck-suppress constParameter
        tk::Fields& /*R*/ )
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
//  // number of transported scalars
//  auto ncomp = U.nprop();
//
//  #if defined(__clang__)
//    #pragma clang diagnostic push
//    #pragma clang diagnostic ignored "-Wvla"
//    #pragma clang diagnostic ignored "-Wvla-extension"
//  #elif defined(STRICT_GNUC)
//    #pragma GCC diagnostic push
//    #pragma GCC diagnostic ignored "-Wvla"
//  #endif
//
//  // domain edge contributions: tetrahedron superedges
//  for (std::size_t e=0; e<dsupedge[0].size()/4; ++e) {
//    const auto N = dsupedge[0].data() + e*4;
//    tk::real u[4][ncomp];
//    for (std::size_t c=0; c<ncomp; ++c) {
//      u[0][c] = U(N[0],c);
//      u[1][c] = U(N[1],c);
//      u[2][c] = U(N[2],c);
//      u[3][c] = U(N[3],c);
//    }
//    // edge fluxes
//    tk::real f[6][ncomp];
//    const auto d = dsupint[0].data();
//    flux( coord, G, d+(e*6+0)*3, N[0], N[1], u[0], u[1], f[0] );
//    flux( coord, G, d+(e*6+1)*3, N[1], N[2], u[1], u[2], f[1] );
//    flux( coord, G, d+(e*6+2)*3, N[2], N[0], u[2], u[0], f[2] );
//    flux( coord, G, d+(e*6+3)*3, N[0], N[3], u[0], u[3], f[3] );
//    flux( coord, G, d+(e*6+4)*3, N[1], N[3], u[1], u[3], f[4] );
//    flux( coord, G, d+(e*6+5)*3, N[2], N[3], u[2], u[3], f[5] );
//    // edge flux contributions
//    for (std::size_t c=0; c<ncomp; ++c) {
//      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c] - f[3][c];
//      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c] - f[4][c];
//      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c] - f[5][c];
//      R(N[3],c) = R(N[3],c) + f[3][c] + f[4][c] + f[5][c];
//    }
//  }
//
//  // domain edge contributions: triangle superedges
//  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
//    const auto N = dsupedge[1].data() + e*3;
//    tk::real u[3][ncomp];
//    for (std::size_t c=0; c<ncomp; ++c) {
//      u[0][c] = U(N[0],c);
//      u[1][c] = U(N[1],c);
//      u[2][c] = U(N[2],c);
//    }
//    // edge fluxes
//    tk::real f[3][ncomp];
//    const auto d = dsupint[1].data();
//    flux( coord, G, d+(e*3+0)*3, N[0], N[1], u[0], u[1], f[0] );
//    flux( coord, G, d+(e*3+1)*3, N[1], N[2], u[1], u[2], f[1] );
//    flux( coord, G, d+(e*3+2)*3, N[2], N[0], u[2], u[0], f[2] );
//    // edge flux contributions
//    for (std::size_t c=0; c<ncomp; ++c) {
//      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c];
//      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c];
//      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c];
//    }
//  }
//
//  // domain edge contributions: edges
//  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
//    const auto N = dsupedge[2].data() + e*2;
//    tk::real u[2][ncomp];
//    for (std::size_t c=0; c<ncomp; ++c) {
//      u[0][c] = U(N[0],c);
//      u[1][c] = U(N[1],c);
//    }
//    // edge fluxes
//    tk::real f[ncomp];
//    const auto d = dsupint[2].data();
//    flux( coord, G, d+e*3, N[0], N[1], u[0], u[1], f );
//    // edge flux contributions
//    for (std::size_t c=0; c<ncomp; ++c) {
//      R(N[0],c) -= f[c];
//      R(N[1],c) += f[c];
//    }
//  }
//
//  #if defined(__clang__)
//    #pragma clang diagnostic pop
//  #elif defined(STRICT_GNUC)
//    #pragma GCC diagnostic pop
//  #endif
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

} // chorin::
