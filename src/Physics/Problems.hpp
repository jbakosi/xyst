// *****************************************************************************
/*!
  \file      src/Physics/Problems.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem-specific functions. Initial conditions, source terms.
*/
// *****************************************************************************
#pragma once

#include "EOS.hpp"
#include "Fields.hpp"

namespace problems {

//! Determine nodes that lie inside user-defined IC box(es)
std::vector< std::unordered_set< std::size_t > >
boxnodes( const std::array< std::vector< tk::real >, 3 >& coord );

//! Query user config and assign function to set initial conditions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
IC();

//! Query user config and assign function to query analytic solutions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
SOL();

//! Query user config and assign function to add a source term
std::function< std::vector< tk::real >
  ( tk::real, tk::real, tk::real, tk::real ) >
SRC();

//  Query user config and assign function to apply source to numerical solution
std::function< void( const std::array< std::vector< tk::real >, 3 >&,
                     tk::real,
                     tk::Fields& ) >
PHYS_SRC();

//! Set inital conditions
void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t,
            const std::vector< std::unordered_set<std::size_t> >& inbox = {} );

} // problems::
