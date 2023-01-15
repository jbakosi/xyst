// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/PointSrc.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for transport equations
  \details   This file defines a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.h
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.h for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************

#include "PointSrc.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::TransportProblemPointSrc;

std::vector< tk::real >
TransportProblemPointSrc::initialize( ncomp_t, ncomp_t ncomp,
  tk::real, tk::real, tk::real, tk::real )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,z,t) for all components
//! \return Values of all components evaluated at (x,y,t)
// *****************************************************************************
{
  std::vector< tk::real > u( ncomp, 0.0 );
  return u;
}

void
TransportProblemPointSrc::src(
  ncomp_t system,
  ncomp_t offset, 
  tk::real t,
  const std::array< std::vector< tk::real >, 3 >& coord,
  tk::Fields& U )
// *****************************************************************************
//  Add source
//! \param[in] system Equation system index
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] t Physical time
//! \param[in,out] U Solution to update with source
//! \return Values of all components evaluated at (x,y,t)
// *****************************************************************************
{
  const auto& src = g_inputdeck.get< tag::param, eq, tag::source >()[ system ];

  auto sx = src[0];
  auto sy = src[1];
  auto sz = src[2];
  auto sr = src[3];
  auto st = src[4];

  if (t < st) return;

  // access node cooordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t p=0; p<U.nunk(); ++p) {
    if (((sx-x[p])*(sx-x[p]) +
         (sy-y[p])*(sy-y[p]) +
         (sz-z[p])*(sz-z[p])) < sr*sr)
    {
      U(p,0,offset) = 1.0;
    }
  }
}

void
TransportProblemPointSrc::errchk( ncomp_t system, ncomp_t ) const
// *****************************************************************************
//  Do error checking on PDE parameters
//! \param[in] system Equation system index, i.e., which transport equation
//!   system we operate on among the systems of PDEs
// *****************************************************************************
{
  using tag::param;

  const auto& s = g_inputdeck.get< param, eq, tag::source >()[ system ];
  ErrChk( s.size() % 5 == 0,
    "Wrong number of advection-diffusion PDE parameters 'source'" );
}
