// *****************************************************************************
/*!
  \file      src/ZoltanInterOp/ZoltanInterOp.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
// *****************************************************************************
#pragma once

#include "Options/PartitioningAlgorithm.hpp"

//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan with a geometric partitioner
std::vector< std::size_t >
geomPartMesh( tk::ctr::PartitioningAlgorithmType alg,
              const std::vector< std::size_t >& inpoel,
              const std::array< std::vector< tk::real >, 3 >& coord,
              int npart );

//! Partition mesh using Zoltan with a geometric partitioner
std::vector< std::size_t >
graphPartMesh( const std::vector< std::size_t >& ginpoel, int npart );

//! Partition mesh using Zoltan with a geometric or graph partitioner
std::vector< std::size_t >
partMesh( tk::ctr::PartitioningAlgorithmType alg,
          const std::vector< std::size_t >& inpoel,
          const std::vector< std::size_t >& ginpoel,
          const std::array< std::vector< tk::real >, 3 >& coord,
          int npart );

} // zoltan::