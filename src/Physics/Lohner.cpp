// *****************************************************************************
/*!
  \file      src/Physics/Lohner.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LohCG: Artificial compressibility solver for incompressible flow
*/
// *****************************************************************************

#include "Vector.hpp"
#include "Around.hpp"
#include "DerivedData.hpp"
#include "EOS.hpp"
#include "Lohner.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace lohner {

static const tk::real muscl_eps = 1.0e-9;
static const tk::real muscl_const = 1.0/3.0;

using inciter::g_cfg;

void
div( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const tk::Fields& U,
     std::vector< tk::real >& D,
     std::size_t pos )
// *****************************************************************************
//  Compute divergence of a vector in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Vector whose divergence to compute
//! \param[in,out] D Vector added to
//! \param[in] pos Position at which the three vector components are in U
// *****************************************************************************
{
  Assert( U.nunk() == D.size(), "Size mismatch" );

  // Lambda to compute the velocity divergence contribution for edge p-q
  auto div = [&]( const tk::real d[], std::size_t p, std::size_t q ){
    return d[0] * (U(p,pos+0) + U(q,pos+0)) +
           d[1] * (U(p,pos+1) + U(q,pos+1)) +
           d[2] * (U(p,pos+2) + U(q,pos+2));
  };

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
    tk::real f[] = {
      div( d+(e*6+0)*4, N[0], N[1] ),
      div( d+(e*6+1)*4, N[1], N[2] ),
      div( d+(e*6+2)*4, N[2], N[0] ),
      div( d+(e*6+3)*4, N[0], N[3] ),
      div( d+(e*6+4)*4, N[1], N[3] ),
      div( d+(e*6+5)*4, N[2], N[3] ) };
    // edge flux contributions
    D[N[0]] = D[N[0]] - f[0] + f[2] - f[3];
    D[N[1]] = D[N[1]] + f[0] - f[1] - f[4];
    D[N[2]] = D[N[2]] + f[1] - f[2] - f[5];
    D[N[3]] = D[N[3]] + f[3] + f[4] + f[5];
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    // edge fluxes
    tk::real f[] = {
      div( d+(e*3+0)*4, N[0], N[1] ),
      div( d+(e*3+1)*4, N[1], N[2] ),
      div( d+(e*3+2)*4, N[2], N[0] ) };
    // edge flux contributions
    D[N[0]] = D[N[0]] - f[0] + f[2];
    D[N[1]] = D[N[1]] + f[0] - f[1];
    D[N[2]] = D[N[2]] + f[1] - f[2];
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data();
    // edge flux
    tk::real f = div( d+e*4, N[0], N[1] );
    // edge flux contributions
    D[N[0]] -= f;
    D[N[1]] += f;
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    auto uxA = U(N[0],pos+0);
    auto uyA = U(N[0],pos+1);
    auto uzA = U(N[0],pos+2);
    auto uxB = U(N[1],pos+0);
    auto uyB = U(N[1],pos+1);
    auto uzB = U(N[1],pos+2);
    auto uxC = U(N[2],pos+0);
    auto uyC = U(N[2],pos+1);
    auto uzC = U(N[2],pos+2);
    auto ux = (6.0*uxA + uxB + uxC)/8.0;
    auto uy = (6.0*uyA + uyB + uyC)/8.0;
    auto uz = (6.0*uzA + uzB + uzC)/8.0;
    D[N[0]] += ux*n[0] + uy*n[1] + uz*n[2];
    ux = (uxA + 6.0*uxB + uxC)/8.0;
    uy = (uyA + 6.0*uyB + uyC)/8.0;
    uz = (uzA + 6.0*uzB + uzC)/8.0;
    D[N[1]] += ux*n[0] + uy*n[1] + uz*n[2];
    ux = (uxA + uxB + 6.0*uxC)/8.0;
    uy = (uyA + uyB + 6.0*uyC)/8.0;
    uz = (uzA + uzB + 6.0*uzC)/8.0;
    D[N[2]] += ux*n[0] + uy*n[1] + uz*n[2];
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
//  Compute gradients of scalar in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Scalar whose gradient to compute
//! \param[in,out] G Nodal gradient of scalar in all points
// *****************************************************************************
{
  if (G.nprop() == 0) return;

  Assert( G.nunk() == U.size(), "Size mismatch" );
  Assert( G.nprop() > 2, "Size mismatch" );
  Assert( G.nprop() % 3 == 0, "Size mismatch" );

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
        d[(e*6+0)*4+j] * (u[1] + u[0]),
        d[(e*6+1)*4+j] * (u[2] + u[1]),
        d[(e*6+2)*4+j] * (u[0] + u[2]),
        d[(e*6+3)*4+j] * (u[3] + u[0]),
        d[(e*6+4)*4+j] * (u[3] + u[1]),
        d[(e*6+5)*4+j] * (u[3] + u[2]) };
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
        d[(e*3+0)*4+j] * (u[1] + u[0]),
        d[(e*3+1)*4+j] * (u[2] + u[1]),
        d[(e*3+2)*4+j] * (u[0] + u[2]) };
      G(N[0],j) = G(N[0],j) - f[0] + f[2];
      G(N[1],j) = G(N[1],j) + f[0] - f[1];
      G(N[2],j) = G(N[2],j) + f[1] - f[2];
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*4;
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
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    auto uA = U[N[0]];
    auto uB = U[N[1]];
    auto uC = U[N[2]];
    auto f = (6.0*uA + uB + uC)/8.0;
    G(N[0],0) += f * n[0];
    G(N[0],1) += f * n[1];
    G(N[0],2) += f * n[2];
    f = (uA + 6.0*uB + uC)/8.0;
    G(N[1],0) += f * n[0];
    G(N[1],1) += f * n[1];
    G(N[1],2) += f * n[2];
    f = (uA + uB + 6.0*uC)/8.0;
    G(N[2],0) += f * n[0];
    G(N[2],1) += f * n[1];
    G(N[2],2) += f * n[2];
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

void
vgrad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
       const std::array< std::vector< tk::real >, 3 >& dsupint,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::vector< std::size_t >& triinpoel,
       const tk::Fields& U,
       tk::Fields& G )
// *****************************************************************************
//  Compute velocity gradients in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Velocity whose gradient to compute
//! \param[in,out] G Nodal velocity gradients (9 components) in all points
// *****************************************************************************
{
  Assert( G.nunk() == U.nunk(), "Size mismatch" );
  Assert( G.nprop() == 9, "Size mismatch" );

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
    for (std::size_t i=0; i<3; ++i) {
      tk::real u[] = { U(N[0],i+1), U(N[1],i+1), U(N[2],i+1), U(N[3],i+1) };
      auto i3 = i*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = { d[(e*6+0)*4+j] * (u[1] + u[0]),
                         d[(e*6+1)*4+j] * (u[2] + u[1]),
                         d[(e*6+2)*4+j] * (u[0] + u[2]),
                         d[(e*6+3)*4+j] * (u[3] + u[0]),
                         d[(e*6+4)*4+j] * (u[3] + u[1]),
                         d[(e*6+5)*4+j] * (u[3] + u[2]) };
        G(N[0],i3+j) = G(N[0],i3+j) - f[0] + f[2] - f[3];
        G(N[1],i3+j) = G(N[1],i3+j) + f[0] - f[1] - f[4];
        G(N[2],i3+j) = G(N[2],i3+j) + f[1] - f[2] - f[5];
        G(N[3],i3+j) = G(N[3],i3+j) + f[3] + f[4] + f[5];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    for (std::size_t i=0; i<3; ++i) {
      tk::real u[] = { U(N[0],i+1), U(N[1],i+1), U(N[2],i+1) };
      auto i3 = i*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = { d[(e*3+0)*4+j] * (u[1] + u[0]),
                         d[(e*3+1)*4+j] * (u[2] + u[1]),
                         d[(e*3+2)*4+j] * (u[0] + u[2]) };
        G(N[0],i3+j) = G(N[0],i3+j) - f[0] + f[2];
        G(N[1],i3+j) = G(N[1],i3+j) + f[0] - f[1];
        G(N[2],i3+j) = G(N[2],i3+j) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*4;
    for (std::size_t i=0; i<3; ++i) {
      tk::real u[] = { U(N[0],i+1), U(N[1],i+1) };
      auto i3 = i*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = d[j] * (u[1] + u[0]);
        G(N[0],i3+j) -= f;
        G(N[1],i3+j) += f;
      }
    }
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    for (std::size_t i=0; i<3; ++i) {
      tk::real u[] = { U(N[0],i+1), U(N[1],i+1), U(N[2],i+1) };
      auto i3 = i*3;
      auto f = (6.0*u[0] + u[1] + u[2])/8.0;
      G(N[0],i3+0) += f * n[0];
      G(N[0],i3+1) += f * n[1];
      G(N[0],i3+2) += f * n[2];
      f = (u[0] + 6.0*u[1] + u[2])/8.0;
      G(N[1],i3+0) += f * n[0];
      G(N[1],i3+1) += f * n[1];
      G(N[1],i3+2) += f * n[2];
      f = (u[0] + u[1] + 6.0*u[2])/8.0;
      G(N[2],i3+0) += f * n[0];
      G(N[2],i3+1) += f * n[1];
      G(N[2],i3+2) += f * n[2];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static tk::real
flux( const tk::Fields& U,
      const tk::Fields& G,
      std::size_t i,
      std::size_t j,
      std::size_t p,
      std::size_t q )
// *****************************************************************************
//! Compute momentum flux over edge of points p-q
//! \param[in] U Velocity vector
//! \param[in] G Velocity gradients
//! \param[in] i Tensor component, 1st index
//! \param[in] j Tensor component, 2nd index
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \return Momentum flux contribution for edge p-q
// *****************************************************************************
{
  auto inv = U(p,i)*U(p,j) + U(q,i)*U(q,j);

  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto mu = g_cfg.get< tag::mat_dyn_viscosity >();
  if (mu < eps) return -inv;

  auto vis = G(p,i*3+j) + G(p,j*3+i) + G(q,i*3+j) + G(q,j*3+i);
  if (i == j) {
    vis -= 2.0/3.0 * ( G(p,0) + G(p,4) + G(p,8) + G(q,0) + G(q,4) + G(q,8) );
  }
  return mu*vis - inv;
}

static tk::real
flux( const tk::Fields& U,
      const tk::Fields& G,
      std::size_t i,
      std::size_t j,
      std::size_t p )
// *****************************************************************************
//! Compute momentum flux in point p
//! \param[in] U Velocity vector
//! \param[in] G Velocity gradients
//! \param[in] i Tensor component, 1st index
//! \param[in] j Tensor component, 2nd index
//! \param[in] p Node index of point
//! \return Momentum flux contribution for point p
// *****************************************************************************
{
  auto inv = U(p,i)*U(p,j);

  auto eps = std::numeric_limits< tk::real >::epsilon();
  auto mu = g_cfg.get< tag::mat_dyn_viscosity >();
  if (mu < eps) return -inv;

  auto vis = G(p,i*3+j) + G(p,j*3+i);
  if (i == j) {
    vis -= 2.0/3.0 * ( G(p,0) + G(p,4) + G(p,8) );
  }
  return mu*vis - inv;
}

void
flux( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const tk::Fields& U,
      const tk::Fields& G,
      tk::Fields& F )
// *****************************************************************************
//  Compute momentum flux in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Velocity field
//! \param[in] G Velocity gradients, dui/dxj, 9 components
//! \param[in,out] F Momentum flux, Fi = d[ sij - uiuj ] / dxj, where
//!    s_ij = mu[ dui/dxj + duj/dxi - 2/3 duk/dxk delta_ij ]
// *****************************************************************************
{
  Assert( F.nunk() == U.nunk(), "Size mismatch" );
  Assert( F.nprop() == 3, "Size mismatch" );
  Assert( G.nunk() == U.nunk(), "Size mismatch" );
  Assert( G.nprop() == 9, "Size mismatch" );

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
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = { d[(e*6+0)*4+j] * flux(U,G,i,j,N[1],N[0]),
                         d[(e*6+1)*4+j] * flux(U,G,i,j,N[2],N[1]),
                         d[(e*6+2)*4+j] * flux(U,G,i,j,N[0],N[2]),
                         d[(e*6+3)*4+j] * flux(U,G,i,j,N[3],N[0]),
                         d[(e*6+4)*4+j] * flux(U,G,i,j,N[3],N[1]),
                         d[(e*6+5)*4+j] * flux(U,G,i,j,N[3],N[2]) };
        F(N[0],i) = F(N[0],i) - f[0] + f[2] - f[3];
        F(N[1],i) = F(N[1],i) + f[0] - f[1] - f[4];
        F(N[2],i) = F(N[2],i) + f[1] - f[2] - f[5];
        F(N[3],i) = F(N[3],i) + f[3] + f[4] + f[5];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = { d[(e*3+0)*4+j] * flux(U,G,i,j,N[1],N[0]),
                         d[(e*3+1)*4+j] * flux(U,G,i,j,N[2],N[1]),
                         d[(e*3+2)*4+j] * flux(U,G,i,j,N[0],N[2]), };
        F(N[0],i) = F(N[0],i) - f[0] + f[2];
        F(N[1],i) = F(N[1],i) + f[0] - f[1];
        F(N[2],i) = F(N[2],i) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*4;
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = d[j] * flux(U,G,i,j,N[1],N[0]);
        F(N[0],i) -= f;
        F(N[1],i) += f;
      }
    }
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    for (std::size_t i=0; i<3; ++i) {
      auto fxA = flux(U,G,i,0,N[0]);
      auto fyA = flux(U,G,i,1,N[0]);
      auto fzA = flux(U,G,i,2,N[0]);
      auto fxB = flux(U,G,i,0,N[1]);
      auto fyB = flux(U,G,i,1,N[1]);
      auto fzB = flux(U,G,i,2,N[1]);
      auto fxC = flux(U,G,i,0,N[2]);
      auto fyC = flux(U,G,i,1,N[2]);
      auto fzC = flux(U,G,i,2,N[2]);
      auto fx = (6.0*fxA + fxB + fxC)/8.0;
      auto fy = (6.0*fyA + fyB + fyC)/8.0;
      auto fz = (6.0*fzA + fzB + fzC)/8.0;
      F(N[0],i) += fx*n[0] + fy*n[1] + fz*n[2];
      fx = (fxA + 6.0*fxB + fxC)/8.0;
      fy = (fyA + 6.0*fyB + fyC)/8.0;
      fz = (fzA + 6.0*fzB + fzC)/8.0;
      F(N[1],i) += fx*n[0] + fy*n[1] + fz*n[2];
      fx = (fxA + fxB + 6.0*fxC)/8.0;
      fy = (fyA + fyB + 6.0*fyC)/8.0;
      fz = (fzA + fzB + 6.0*fzC)/8.0;
      F(N[2],i) += fx*n[0] + fy*n[1] + fz*n[2];
    }
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
      const tk::Fields& U,
      tk::Fields& G )
// *****************************************************************************
//  Compute gradients of scalar in all points
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Unknown vector whose gradient to compute
//! \param[in,out] G Nodal gradients in all points
// *****************************************************************************
{
  if (G.nprop() == 0) return;

  Assert( G.nunk() == U.nunk(), "Size mismatch" );
  Assert( G.nprop() > 2, "Size mismatch" );
  Assert( G.nprop() % 3 == 0, "Size mismatch" );
  Assert( G.nprop() == (U.nprop()-1)*3, "Size mismatch" );

  const auto ncomp = U.nprop();

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
    for (std::size_t c=1; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c), U(N[2],c), U(N[3],c) };
      auto g = (c-1)*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = {
          d[(e*6+0)*4+j] * (u[1] + u[0]),
          d[(e*6+1)*4+j] * (u[2] + u[1]),
          d[(e*6+2)*4+j] * (u[0] + u[2]),
          d[(e*6+3)*4+j] * (u[3] + u[0]),
          d[(e*6+4)*4+j] * (u[3] + u[1]),
          d[(e*6+5)*4+j] * (u[3] + u[2]) };
        G(N[0],g+j) = G(N[0],g+j) - f[0] + f[2] - f[3];
        G(N[1],g+j) = G(N[1],g+j) + f[0] - f[1] - f[4];
        G(N[2],g+j) = G(N[2],g+j) + f[1] - f[2] - f[5];
        G(N[3],g+j) = G(N[3],g+j) + f[3] + f[4] + f[5];
      }
    }
  }

  // domain edge contributions: triangle superedges
  for (std::size_t e=0; e<dsupedge[1].size()/3; ++e) {
    const auto N = dsupedge[1].data() + e*3;
    const auto d = dsupint[1].data();
    for (std::size_t c=1; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c), U(N[2],c) };
      auto g = (c-1)*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f[] = {
          d[(e*3+0)*4+j] * (u[1] + u[0]),
          d[(e*3+1)*4+j] * (u[2] + u[1]),
          d[(e*3+2)*4+j] * (u[0] + u[2]) };
        G(N[0],g+j) = G(N[0],g+j) - f[0] + f[2];
        G(N[1],g+j) = G(N[1],g+j) + f[0] - f[1];
        G(N[2],g+j) = G(N[2],g+j) + f[1] - f[2];
      }
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    const auto d = dsupint[2].data() + e*4;
    for (std::size_t c=1; c<ncomp; ++c) {
      tk::real u[] = { U(N[0],c), U(N[1],c) };
      auto g = (c-1)*3;
      for (std::size_t j=0; j<3; ++j) {
        tk::real f = d[j] * (u[1] + u[0]);
        G(N[0],g+j) -= f;
        G(N[1],g+j) += f;
      }
    }
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    for (std::size_t c=1; c<ncomp; ++c) {
      auto g = (c-1)*3;
      auto uA = U(N[0],c);
      auto uB = U(N[1],c);
      auto uC = U(N[2],c);
      auto f = (6.0*uA + uB + uC)/8.0;
      G(N[0],g+0) += f * n[0];
      G(N[0],g+1) += f * n[1];
      G(N[0],g+2) += f * n[2];
      f = (uA + 6.0*uB + uC)/8.0;
      G(N[1],g+0) += f * n[0];
      G(N[1],g+1) += f * n[1];
      G(N[1],g+2) += f * n[2];
      f = (uA + uB + 6.0*uC)/8.0;
      G(N[2],g+0) += f * n[0];
      G(N[2],g+1) += f * n[1];
      G(N[2],g+2) += f * n[2];
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
adv_damp2( const tk::real supint[],
           const tk::Fields& U,
           const tk::Fields&,
           const std::array< std::vector< tk::real >, 3 >&,
           std::size_t p,
           std::size_t q,
           tk::real f[] )
// *****************************************************************************
//! Compute advection fluxes on a single edge using 2nd-order damping
//! \param[in] supint Edge integral
//! \param[in] U Velocity and transported scalars at recent time step
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] f Flux computed
// *****************************************************************************
{
  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  // left state
  auto pL = U(p,0);
  auto uL = U(p,1);
  auto vL = U(p,2);
  auto wL = U(p,3);
  auto vnL = uL*nx + vL*ny + wL*nz;

  // right state
  auto pR = U(q,0);
  auto uR = U(q,1);
  auto vR = U(q,2);
  auto wR = U(q,3);
  auto vnR = uR*nx + vR*ny + wR*nz;

  auto s = g_cfg.get< tag::soundspeed >();
  auto s2 = s*s;
  auto d = supint[3] * g_cfg.get< tag::mat_dyn_viscosity >();

  // flow
  auto pf = pL + pR;
  f[0] = (vnL + vnR)*s2;
  f[1] = uL*vnL + uR*vnR + pf*nx - d*(uR - uL);
  f[2] = vL*vnL + vR*vnR + pf*ny - d*(vR - vL);
  f[3] = wL*vnL + wR*vnR + pf*nz - d*(wR - wL);

  // artificial viscosity
  const auto stab2 = g_cfg.get< tag::stab2 >();
  if (!stab2) return;

  auto len = tk::length( nx, ny, nz );
  auto sl = std::abs(vnL) + s*len;
  auto sr = std::abs(vnR) + s*len;
  auto aw = g_cfg.get< tag::stab2coef >() * std::max(sl,sr);
  f[0] += aw * (pR - pL)*s2;
  f[1] += aw * (uR - uL);
  f[2] += aw * (vR - vL);
  f[3] += aw * (wR - wL);
}

static void
adv_damp4( const tk::real supint[],
           const tk::Fields& U,
           const tk::Fields& G,
           const std::array< std::vector< tk::real >, 3 >& coord,
           std::size_t p,
           std::size_t q,
           tk::real f[] )
// *****************************************************************************
//! Compute advection fluxes on a single edge using 4th-order damping
//! \param[in] supint Edge integral
//! \param[in] U Velocity and transported scalars at recent time step
//! \param[in] G Gradients of velocity and transported scalars
//! \param[in] coord Mesh node coordinates
//! \param[in] p Left node index of edge
//! \param[in] q Right node index of edge
//! \param[in,out] f Flux computed
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

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  tk::real uL[] = { U(p,1), U(p,2), U(p,3) };
  tk::real uR[] = { U(q,1), U(q,2), U(q,3) };

  // MUSCL reconstruction in edge-end points
  for (std::size_t c=0; c<ncomp-1; ++c) {
    auto g = c*3;
    auto g1 = G(p,g+0)*dx + G(p,g+1)*dy + G(p,g+2)*dz;
    auto g2 = G(q,g+0)*dx + G(q,g+1)*dy + G(q,g+2)*dz;
    auto delta2 = uR[c] - uL[c];
    auto delta1 = 2.0 * g1 - delta2;
    auto delta3 = 2.0 * g2 - delta2;

    // van Leer limiter
    auto rL = (delta2 + muscl_eps) / (delta1 + muscl_eps);
    auto rR = (delta2 + muscl_eps) / (delta3 + muscl_eps);
    auto rLinv = (delta1 + muscl_eps) / (delta2 + muscl_eps);
    auto rRinv = (delta3 + muscl_eps) / (delta2 + muscl_eps);
    auto phiL = (std::abs(rL) + rL) / (std::abs(rL) + 1.0);
    auto phiR = (std::abs(rR) + rR) / (std::abs(rR) + 1.0);
    auto phi_L_inv = (std::abs(rLinv) + rLinv) / (std::abs(rLinv) + 1.0);
    auto phi_R_inv = (std::abs(rRinv) + rRinv) / (std::abs(rRinv) + 1.0);
    // update unknowns with reconstructed unknowns
    uL[c] += 0.25*(delta1*(1.0-muscl_const)*phiL +
                   delta2*(1.0+muscl_const)*phi_L_inv);
    uR[c] -= 0.25*(delta3*(1.0-muscl_const)*phiR +
                   delta2*(1.0+muscl_const)*phi_R_inv);
  }

  auto nx = supint[0];
  auto ny = supint[1];
  auto nz = supint[2];

  // normal velocities
  auto vnL = uL[0]*nx + uL[1]*ny + uL[2]*nz;
  auto vnR = uR[0]*nx + uR[1]*ny + uR[2]*nz;

  auto s = g_cfg.get< tag::soundspeed >();
  auto s2 = s*s;
  auto d = supint[3] * g_cfg.get< tag::mat_dyn_viscosity >();

  // flow
  auto pf = U(p,0) + U(q,0);
  f[0] = (vnL + vnR)*s2;
  f[1] = uL[0]*vnL + uR[0]*vnR + pf*nx - d*(uR[0] - uL[0]);
  f[2] = uL[1]*vnL + uR[1]*vnR + pf*ny - d*(uR[1] - uL[1]);
  f[3] = uL[2]*vnL + uR[2]*vnR + pf*nz - d*(uR[2] - uL[2]);

  // artificial viscosity
  const auto stab2 = g_cfg.get< tag::stab2 >();
  if (!stab2) return;

  auto len = tk::length( nx, ny, nz );
  auto sl = std::abs(vnL) + s*len;
  auto sr = std::abs(vnR) + s*len;
  auto aw = g_cfg.get< tag::stab2coef >() * std::max(sl,sr);
  f[0] += aw * (U(q,0) - U(p,0))*s2;
  f[1] += aw * (uR[0] - uL[0]);
  f[2] += aw * (uR[1] - uL[1]);
  f[3] += aw * (uR[2] - uL[2]);

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

static void
adv( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const tk::Fields& U,
     const tk::Fields& G,
     // cppcheck-suppress constParameter
     tk::Fields& R )
// *****************************************************************************
//! Add advection to rhs
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Velocity and transported scalars at recent time step
//! \param[in] G Gradients of velocity and transported scalars
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  const auto ncomp = U.nprop();

  // configure advection
  auto adv = [](){
    const auto& flux = g_cfg.get< tag::flux >();
         if (flux == "damp2") return adv_damp2;
    else if (flux == "damp4") return adv_damp4;
    else Throw( "Flux not correctly configured" );
  }();

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
    tk::real f[6][ncomp];
    const auto d = dsupint[0].data();
    adv( d+(e*6+0)*4, U, G, coord, N[0], N[1], f[0] );
    adv( d+(e*6+1)*4, U, G, coord, N[1], N[2], f[1] );
    adv( d+(e*6+2)*4, U, G, coord, N[2], N[0], f[2] );
    adv( d+(e*6+3)*4, U, G, coord, N[0], N[3], f[3] );
    adv( d+(e*6+4)*4, U, G, coord, N[1], N[3], f[4] );
    adv( d+(e*6+5)*4, U, G, coord, N[2], N[3], f[5] );
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
    tk::real f[3][ncomp];
    const auto d = dsupint[1].data();
    adv( d+(e*3+0)*4, U, G, coord, N[0], N[1], f[0] );
    adv( d+(e*3+1)*4, U, G, coord, N[1], N[2], f[1] );
    adv( d+(e*3+2)*4, U, G, coord, N[2], N[0], f[2] );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) = R(N[0],c) - f[0][c] + f[2][c];
      R(N[1],c) = R(N[1],c) + f[0][c] - f[1][c];
      R(N[2],c) = R(N[2],c) + f[1][c] - f[2][c];
    }
  }

  // domain edge contributions: edges
  for (std::size_t e=0; e<dsupedge[2].size()/2; ++e) {
    const auto N = dsupedge[2].data() + e*2;
    tk::real f[ncomp];
    const auto d = dsupint[2].data();
    adv( d+e*4, U, G, coord, N[0], N[1], f );
    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) -= f[c];
      R(N[1],c) += f[c];
    }
  }

  // boundary integral

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto s = g_cfg.get< tag::soundspeed >();
  auto s2 = s * s;

  for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
    const auto N = triinpoel.data() + e*3;
    tk::real n[3];
    tk::crossdiv( x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]],
                  x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]], 6.0,
                  n[0], n[1], n[2] );
    tk::real f[ncomp][3];

    auto p = U(N[0],0);
    auto u = U(N[0],1);
    auto v = U(N[0],2);
    auto w = U(N[0],3);
    auto vn = n[0]*u + n[1]*v + n[2]*w;
    f[0][0] = vn * s2;
    f[1][0] = u*vn + p*n[0];
    f[2][0] = v*vn + p*n[1];
    f[3][0] = w*vn + p*n[2];

    p = U(N[1],0);
    u = U(N[1],1);
    v = U(N[1],2);
    w = U(N[1],3);
    vn = n[0]*u + n[1]*v + n[2]*w;
    f[0][1] = vn * s2;
    f[1][1] = u*vn + p*n[0];
    f[2][1] = v*vn + p*n[1];
    f[3][1] = w*vn + p*n[2];

    p = U(N[2],0);
    u = U(N[2],1);
    v = U(N[2],2);
    w = U(N[2],3);
    vn = n[0]*u + n[1]*v + n[2]*w;
    f[0][2] = vn * s2;
    f[1][2] = u*vn + p*n[0];
    f[2][2] = v*vn + p*n[1];
    f[3][2] = w*vn + p*n[2];

    for (std::size_t c=0; c<ncomp; ++c) {
      R(N[0],c) += (6.0*f[c][0] + f[c][1] + f[c][2])/8.0;
      R(N[1],c) += (f[c][0] + 6.0*f[c][1] + f[c][2])/8.0;
      R(N[2],c) += (f[c][0] + f[c][1] + 6.0*f[c][2])/8.0;
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
     const tk::Fields& U,
     const tk::Fields& G,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] dsupedge Domain superedges
//! \param[in] dsupint Domain superedge integrals
//! \param[in] coord Mesh node coordinates
//! \param[in] triinpoel Boundary face connectivity
//! \param[in] U Solution vector of primitive variables at recent time step
//! \param[in] G Gradients of velocity and transported scalars
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );
  Assert( R.nunk() == coord[0].size(),
          "Number of unknowns and/or number of components in right-hand "
          "side vector incorrect" );

  R.fill( 0.0 );
  adv( dsupedge, dsupint, coord, triinpoel, U, G, R );
}

} // lohner::
