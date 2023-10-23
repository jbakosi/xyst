// *****************************************************************************
/*!
  \file      src/Physics/Kozak.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Kozak, FCT limiting for element-based continuous Galerkin
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

static void
adv( const std::vector< std::size_t >& inpoel,
     const std::array< std::vector< tk::real >, 3 >& coord,
     tk::real dt,
     const tk::Fields& U,
     // cppcheck-suppress constParameter
     tk::Fields& R )
// *****************************************************************************
//! Compute integrals for advection
//! \param[in] coord Mesh node coordinates
//! \param[in] dt Physical time size
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto ncomp = U.nprop();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
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
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

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
        ue[0] -= coef * grad[a][j] * U(N[a],j+1,0);
        auto uj = U(N[a],j+1,0) / U(N[a],0,0);
        for (std::size_t i=0; i<3; ++i) {
          ue[i+1] -= coef * grad[a][j] * U(N[a],i+1,0) * uj;
        }
        ue[j+1] -= coef * grad[a][j] * p[a];
        ue[4] -= coef * grad[a][j] * (U(N[a],4,0) + p[a]) * uj;
        for (std::size_t c=5; c<ncomp; ++c) {
          ue[c] -= coef * grad[a][j] * U(N[a],c,0) * uj;
        }
      }
    }

    auto  r = ue[0];
    auto ru = ue[1]/r;
    auto rv = ue[2]/r;
    auto rw = ue[3]/r;
    auto pr = eos::pressure( ue[4] - 0.5*(ru*ru + rv*rv + rw*rw)/r );

    coef = dt/6.0;
    for (std::size_t j=0; j<3; ++j) {
      auto uj = ue[j+1] / ue[0];
      for (std::size_t a=0; a<4; ++a) {
        R(N[a],0,0) += coef * grad[a][j] * ue[j+1];
        for (std::size_t i=0; i<3; ++i) {
          R(N[a],i+1,0) += coef * grad[a][j] * ue[i+1] * uj;
        }
        R(N[a],j+1,0) += coef * grad[a][j] * pr;
        R(N[a],4,0) += coef * grad[a][j] * (ue[4] + pr) * uj;
        for (std::size_t c=5; c<ncomp; ++c) {
          R(N[a],c,0) += coef * grad[a][j] * ue[c] * uj;
        }
      }
    }
  }
}

static void
src( const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< tk::real >& v,
     tk::real dt,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R )
// *****************************************************************************
//  Compute source integral
//! \param[in] coord Mesh node coordinates
//! \param[in] v Nodal mesh volumes without contributions from other chares
//! \param[in] dt Physical time size
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
    for (std::size_t c=0; c<s.size(); ++c) R(p,c,0) += dt * s[c] * v[p];
  }
}

void
rhs( const std::vector< std::size_t >& inpoel,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< tk::real >& v,
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     const tk::Fields& U,
     tk::Fields& R )
// *****************************************************************************
//  Compute right hand side
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
  adv( inpoel, coord, dt, U, R );

  // source
  src( coord, v, dt, t, tp, R );
}

} // kozak::
