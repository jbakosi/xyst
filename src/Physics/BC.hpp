// *****************************************************************************
/*!
  \file      src/Physics/BC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace physics {

//! Set Dirichlet boundary conditions
void
dirbc(
  tk::Fields& U,
  tk::real t,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const std::vector< std::unordered_set< std::size_t > >& boxnodes,
  const std::vector< std::size_t >& dirbcmasks );

//! Set symmetry boundary conditions
void
symbc( tk::Fields& U,
       const std::vector< std::size_t >& symbcnodes,
       const std::vector< tk::real >& symbcnorms,
       std::size_t pos );

//! Set farfield boundary conditions
void
farbc( tk::Fields& U,
       const std::vector< std::size_t >& farbcnodes,
       const std::vector< tk::real >& farbcnorms );

//! Set pressure boundary conditions at nodes
void
prebc( tk::Fields& U,
       const std::vector< std::size_t >& prebcnodes,
       const std::vector< tk::real >& prebcvals );

} // physics::
