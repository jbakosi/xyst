// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Types for Incitier's parsers
  \details   Types for Incitier's parsers. This file defines the components of
    the agged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during parsing (both command-line argument parsing and control file
    parsing).
*/
// *****************************************************************************
#pragma once

#include "Tags.hpp"
#include "TaggedTuple.hpp"
#include "Base/Types.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

using namespace tao;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::control,   kw::control::info::expect::type //!< Control filename
  , tag::input,     kw::input::info::expect::type   //!< Input filename
  , tag::output,    kw::output::info::expect::type  //!< Output filename
    //! Refined output (output field data on a refined mesh)
  , tag::refined,   kw::refined::info::expect::type
    //! Diagnostics filename
  , tag::diag,      kw::diagnostics_cmd::info::expect::type
  , tag::particles, std::string                     //!< Particles filename
  , tag::restart,   kw::restart::info::expect::type //!< Restart dirname
> >;

} // ctr::
} // inciter::
