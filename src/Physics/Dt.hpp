// *****************************************************************************
/*!
  \file      src/Physics/Dt.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Compute dt for next time step for edge-based continuous Galerkin
*/
// *****************************************************************************
#pragma once

#include <vector>

#include "Fields.hpp"

namespace physics {

//! Compute minimum time step size
tk::real
dt( const std::vector< tk::real >& vol, const tk::Fields& U );

//! Compute time step size for each mesh node (for steady time stepping)
void
dt( const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp );

} // physics::
