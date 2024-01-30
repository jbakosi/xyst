// *****************************************************************************
/*!
  \file      src/Inciter/History.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Types for collecting history output
  \details   Types for collecting history output.
*/
// *****************************************************************************
#pragma once

#include "TaggedTuple.hpp"

namespace tag {
struct id;
struct elem;
struct fn;
} // tag::

namespace inciter {

//! History point data
using HistData = tk::TaggedTuple< brigand::list<
    tag::id,    std::string               //!< Point identifier
  , tag::elem,  std::size_t               //!< Host elem id
  , tag::fn,    std::array< tk::real, 4 > //!< Shapefunctions evaluated at point
> >;

} // inciter::
