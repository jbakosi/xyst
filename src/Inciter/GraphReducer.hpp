// *****************************************************************************
/*!
  \file      src/Inciter/GraphReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging mesh graphs across PEs
*/
// *****************************************************************************
#pragma once

#include <unordered_map>
#include <vector>

#include "NoWarning/charm++.hpp"

namespace tk {

//! Serialize to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_map< std::size_t,
                                     std::vector< std::size_t > >& d );

//! Charm++ custom reducer for merging during reduction across PEs
CkReductionMsg*
mergeGraph( int nmsg, CkReductionMsg **msgs );

} // tk::
