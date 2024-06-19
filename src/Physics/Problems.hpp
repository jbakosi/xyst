// *****************************************************************************
/*!
  \file      src/Physics/Problems.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
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
             ( tk::real, tk::real, tk::real, tk::real ) >
IC();

//! Query user config and assign function to query analytic solutions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
SOL();

//! Set inital conditions
void
initialize(
  const std::array< std::vector< tk::real >, 3 >& coord,
  tk::Fields& U,
  tk::real t,
  const std::vector< std::unordered_set<std::size_t> >& boxnodes = {} );

//! Query user config and assign function to set pressure initial conditions
std::function< tk::real( tk::real, tk::real, tk::real ) >
PRESSURE_IC();

//! Query user config and assign function to query analytic pressure solutions
std::function< tk::real( tk::real, tk::real, tk::real ) >
PRESSURE_SOL();

//! Set pressure right hand side
void
pressure_rhs( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< tk::real >& vol,
              std::vector< tk::real >& r );

//! Set pressure initial condition
tk::real
initialize( tk::real x, tk::real y, tk::real z );

//! Query user config and assign function to add a source term
std::function< std::vector< tk::real >
  ( tk::real, tk::real, tk::real, tk::real ) >
SRC();

//  Query user config and assign function to apply source to numerical solution
std::function< void( const std::array< std::vector< tk::real >, 3 >&,
                     tk::real,
                     tk::Fields& ) >
PHYS_SRC();

} // problems::
