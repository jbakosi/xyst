// *****************************************************************************
/*!
  \file      src/Inciter/IntegralReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging integrals across PEs
  \details   Custom Charm++ reducer for merging integrals across PEs.
*/
// *****************************************************************************
#pragma once

#include <vector>
#include <map>
#include <memory>
#include <utility>

#include "NoWarning/charm++.hpp"

#include "Types.hpp"

namespace inciter {
namespace integrals {

//! Serialize to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid,
           const std::vector< std::map< int, tk::real > >& d );

//! Charm++ custom reducer for merging duuring reduction across PEs
CkReductionMsg*
mergeIntegrals( int nmsg, CkReductionMsg **msgs );

} // integrals::
} // inciter::
