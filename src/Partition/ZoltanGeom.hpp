// *****************************************************************************
/*!
  \file      src/Partition/ZoltanGeom.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library's geometric partitioners
*/
// *****************************************************************************
#pragma once

#include <vector>
#include <string>
#include <array>
#include <unordered_map>

#include "Types.hpp"

namespace inciter {

//! Partition mesh using Zoltan with a geometric partitioner
std::vector< std::size_t >
geomPartMesh( const char* alg,
              const std::vector< std::string >& zoltan_params,
              const std::vector< std::size_t >& inpoel,
              const std::array< std::vector< tk::real >, 3 >& coord,
              int npart );

} // inciter::
