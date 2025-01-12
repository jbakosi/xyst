// *****************************************************************************
/*!
  \file      src/Partition/ZoltanGraph.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library's graph partitioners
*/
// *****************************************************************************
#pragma once

#include <vector>
#include <unordered_map>

#include "Types.hpp"
#include "Exception.hpp"

namespace inciter {

//! Partition mesh using Zoltan with a graph partitioner
std::unordered_map< std::size_t, std::size_t >
graphPartMesh( const std::vector< std::size_t >& ginpoel,
               const std::unordered_map< std::size_t,
                                         std::vector< std::size_t > >& graph,
               const std::vector< std::string >& zoltan_params,
               int npart );

} // inciter::
