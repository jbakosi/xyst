// *****************************************************************************
/*!
  \file      src/Physics/Kozak.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     KozCG: Taylor-Galerkin, FCT, element-based continuous Galerkin
*/
// *****************************************************************************

#include "Vector.hpp"
#include "EOS.hpp"
#include "Kozak.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace kozak {

using inciter::g_cfg;

void
rhs( const std::vector< std::size_t >& inpoel,
     const std::array< std::vector< tk::real >, 3 >& coord,
     tk::real dt,
     tk::real t,
     const std::vector< tk::real >& tp,
     const tk::Fields& U,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
//! \param[in] inpoel Tetrahedron connectivity
//! \param[in] coord Mesh node coordinates
//! \param[in] U Unknowns/solution vector in mesh nodes
//! \param[in] t Physical time
//! \param[in] dt Physical time size
//! \param[in] tp Physical time for each mesh node
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Size mismatch" );
  Assert( R.nunk() == coord[0].size(), "Size mismatch" );

  // zero right hand side for all components
  R.fill( 0.0 );

  const auto ncomp = U.nprop();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto src = problems::SRC();

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wvla"
    #pragma clang diagnostic ignored "-Wvla-extension"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wvla"
  #endif

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {

    // Element gradients and Jacobian
    const auto N = inpoel.data() + e*4;
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::cross( ca, da );
    grad[2] = tk::cross( da, ba );
    grad[3] = tk::cross( ba, ca );
    for (std::size_t i=0; i<3; ++i) {
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    }

    // Taylor-Galerkin first half step

    tk::real p[4];
    for (std::size_t a=0; a<4; ++a) {
      auto  r = U(N[a],0,0);
      auto ru = U(N[a],1,0);
      auto rv = U(N[a],2,0);
      auto rw = U(N[a],3,0);
      p[a] = eos::pressure( U(N[a],4,0) - 0.5*(ru*ru + rv*rv + rw*rw)/r );
    }

    tk::real ue[ncomp];
    for (std::size_t c=0; c<ncomp; ++c) {
      ue[c] = (U(N[0],c,0) + U(N[1],c,0) + U(N[2],c,0) + U(N[3],c,0))/4.0;
    }

    auto coef = dt/J/2.0;
    for (std::size_t j=0; j<3; ++j) {
      for (std::size_t a=0; a<4; ++a) {
        auto cg = coef * grad[a][j];
        auto uj = U(N[a],j+1,0) / U(N[a],0,0);
        ue[0] -= cg * U(N[a],j+1,0);
        ue[1] -= cg * U(N[a],1,0) * uj;
        ue[2] -= cg * U(N[a],2,0) * uj;
        ue[3] -= cg * U(N[a],3,0) * uj;
        ue[j+1] -= cg * p[a];
        ue[4] -= cg * (U(N[a],4,0) + p[a]) * uj;
        for (std::size_t c=5; c<ncomp; ++c) {
          ue[c] -= cg * U(N[a],c,0) * uj;
        }
      }
    }

    if (src) {
      coef = dt/8.0;
      for (std::size_t a=0; a<4; ++a) {
        if (g_cfg.get< tag::steady >()) t = tp[N[a]];
        auto s = src( x[N[a]], y[N[a]], z[N[a]], t );
        for (std::size_t c=0; c<ncomp; ++c) {
          ue[c] += coef * s[c];
        }
      }
    }

    // Taylor-Galerkin: second half step

    auto  r = ue[0];
    auto ru = ue[1];
    auto rv = ue[2];
    auto rw = ue[3];
    auto pr = eos::pressure( ue[4] - 0.5*(ru*ru + rv*rv + rw*rw)/r );

    coef = dt/6.0;
    for (std::size_t j=0; j<3; ++j) {
      auto uj = ue[j+1] / ue[0];
      for (std::size_t a=0; a<4; ++a) {
        auto cg = coef * grad[a][j];
        R(N[a],0,0) += cg * ue[j+1];
        R(N[a],1,0) += cg * ue[1] * uj;
        R(N[a],2,0) += cg * ue[2] * uj;
        R(N[a],3,0) += cg * ue[3] * uj;
        R(N[a],j+1,0) += cg * pr;
        R(N[a],4,0) += cg * (ue[4] + pr) * uj;
        for (std::size_t c=5; c<ncomp; ++c) {
          R(N[a],c,0) += cg * ue[c] * uj;
        }
      }
    }

    if (src) {
      auto xe = (x[N[0]] + x[N[1]] + x[N[2]] + x[N[3]])/4.0;
      auto ye = (y[N[0]] + y[N[1]] + y[N[2]] + y[N[3]])/4.0;
      auto ze = (z[N[0]] + z[N[1]] + z[N[2]] + z[N[3]])/4.0;
      if (g_cfg.get< tag::steady >()) {
        t = (tp[N[0]] + tp[N[1]] + tp[N[2]] + tp[N[3]])/4.0;
      }
      auto se = src( xe, ye, ze, t+dt/2.0 );
      coef = dt*J/24.0;
      for (std::size_t a=0; a<4; ++a) {
        for (std::size_t c=0; c<ncomp; ++c) {
          R(N[a],c,0) += coef * se[c];
        }
      }
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif
}

} // kozak::
