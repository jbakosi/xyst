// *****************************************************************************
/*!
  \file      src/Physics/Problems.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \bried     Problem-specific functions. Initial conditions, source terms.
*/
// *****************************************************************************
#pragma once

#include "Inciter/InputDeck/InputDeck.hpp"
#include "EOS.hpp"
#include "Fields.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

namespace problems {

//! Query user config and assign function to set initial conditions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
IC();

//! Query user config and assign function to query analytic solutions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
SOL();

//! Query user config and assign function to add a source term
std::function< void( tk::real, tk::real, tk::real, tk::real,
  tk::real&, tk::real&, tk::real&, tk::real&, tk::real&, tk::real& ) >
SRC();

//! Initalize the compressible flow equations, prepare for time integration
void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t );

} // problems::
