// *****************************************************************************
/*!
  \file      src/Physics/Box.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Initial condition box related functionality
*/
// *****************************************************************************
#pragma once

#include <vector>
#include <unordered_set>

#include "Types.hpp"

namespace problems {

//! Determine nodes that lie inside user-defined IC box(es)
std::vector< std::unordered_set< std::size_t > >
boxnodes( const std::array< std::vector< tk::real >, 3 >& coord );

//! Evaluate solution in user-defined IC box
void
box( std::size_t p,
     std::vector< tk::real >& u,
     const std::vector< std::unordered_set< std::size_t > >& boxnodes );

} // problems::
