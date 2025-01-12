// *****************************************************************************
/*!
  \file      src/Inciter/PartsReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging mesh part assignments across PEs
*/
// *****************************************************************************
#pragma once

#include <unordered_map>

#include "NoWarning/charm++.hpp"

namespace tk {

//! Serialize to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_map< std::size_t, std::size_t >& d );

//! Charm++ custom reducer for merging during reduction across PEs
CkReductionMsg*
mergeParts( int nmsg, CkReductionMsg **msgs );

} // tk::
