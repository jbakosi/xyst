// *****************************************************************************
/*!
  \file      src/LoadBalance/ZoltanInterOp.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library's partitioners.
*/
// *****************************************************************************

#include "ZoltanInterOp.hpp"

namespace zoltan {

std::vector< std::size_t >
partMesh( tk::ctr::PartitioningAlgorithmType alg,
          const std::vector< std::size_t >& inpoel,
          const std::vector< std::size_t >& /*ginpoel*/,
          const std::array< std::vector< tk::real >, 3 >& coord,
          int npart )
// *****************************************************************************
//  Partition mesh using Zoltan with a geometric or graph partitioner
//! \param[in] alg Partitioning algorithm type
//! \param[in] inpoel Mesh connectivity with local ids
//! \param[in] ginpoel Mesh connectivity with global ids
//! \param[in] coord Node coordinates
//! \param[in] npart Number of desired partitions
//! \return Array of chare ownership IDs mapping elements to chares
// *****************************************************************************
{
  std::vector< std::size_t > chare;

  if (alg == tk::ctr::PartitioningAlgorithmType::PHG)
    //chare = graphPartMesh( ginpoel, npart );
    Throw( "Unimplemented" );
  else
    chare = geomPartMesh( alg, inpoel, coord, npart );

  return chare;
}

} // zoltan::
