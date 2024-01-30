// *****************************************************************************
/*!
  \file      src/Partition/ZoltanInterOp.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
// *****************************************************************************
#pragma once

#include <vector>
#include <string>
#include <array>
#include <numeric>

#include "Types.hpp"
#include "Exception.hpp"

namespace inciter {

//! Partition mesh using Zoltan with a geometric partitioner
std::vector< std::size_t >
geomPartMesh( const char* alg,
              const std::vector< std::string >& zoltan_params,
              const std::vector< std::size_t >& inpoel,
              const std::array< std::vector< tk::real >, 3 >& coord,
              int npart );

//! Partition mesh using Zoltan with a geometric partitioner
std::vector< std::size_t >
graphPartMesh( const std::vector< std::size_t >& ginpoel,
               const std::vector< std::string >& zoltan_params,
               int npart );

//! Partition mesh using Zoltan2 with a geometric partitioner
std::vector< std::size_t >
geomPartMeshZ2( const char* alg,
                const std::vector< std::string >& zoltan_params,
                const std::vector< std::size_t >& inpoel,
                const std::array< std::vector< tk::real >, 3 >& coord,
                int npart );

//! Partition mesh using Zoltan with a geometric or graph partitioner
std::vector< std::size_t >
partMesh( const std::string& alg,
          const std::vector< std::string >& zoltan_params,
          const std::vector< std::size_t >& inpoel,
          const std::vector< std::size_t >& ginpoel,
          const std::array< std::vector< tk::real >, 3 >& coord,
          int npart );

} // inciter::
