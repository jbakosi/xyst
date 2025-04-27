// *****************************************************************************
/*!
  \file      src/Physics/Problems.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem-specific functions. Initial conditions, source terms.
*/
// *****************************************************************************
#pragma once

#include "EOS.hpp"
#include "Fields.hpp"

namespace problems {

//! Query user config and assign function to set initial conditions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real, std::size_t ) >
IC();

//! Query user config and assign function to query analytic solutions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real, std::size_t ) >
SOL();

//! Set inital conditions
void
initialize(
  const std::array< std::vector< tk::real >, 3 >& coord,
  tk::Fields& U,
  tk::real t,
  std::size_t meshid,
  const std::vector< std::unordered_set<std::size_t> >& boxnodes = {} );

//! Query user config and assign function to set pressure initial conditions
std::function< tk::real( tk::real, tk::real, tk::real, std::size_t ) >
PRESSURE_IC();

//! Query user config and assign function to query analytic pressure solutions
std::function< tk::real( tk::real, tk::real, tk::real, std::size_t ) >
PRESSURE_SOL();

//! Assign function to query pressure gradient at a point
std::function< std::array< tk::real, 3 >( tk::real, tk::real, tk::real ) >
PRESSURE_GRAD();

//! Assign function to set pressure solve right hand side
std::function< tk::real( tk::real, tk::real, tk::real ) >
PRESSURE_RHS();

//! Set pressure initial condition
tk::real
initialize( tk::real x, tk::real y, tk::real z, std::size_t meshid );

//! Query user config and assign function to add a source term
std::function< std::vector< tk::real >
  ( tk::real, tk::real, tk::real, tk::real, std::size_t ) >
SRC();

//  Query user config and assign function to apply source to numerical solution
std::function< void( const std::array< std::vector< tk::real >, 3 >&,
                     tk::real,
                     tk::Fields& ) >
PHYS_SRC();

} // problems::
