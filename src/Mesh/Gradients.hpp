// *****************************************************************************
/*!
  \file      src/Mesh/Gradients.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Functions computing gradients on unstructured meshes for tetrahedra
  \details   Functions computing gradients using linear finite element shape
             functions on unstructured meshes for tetrahedra.
*/
// *****************************************************************************
#pragma once

#include <array>
#include <stddef.h>
#include <vector>
#include <utility>

#include "Fields.hpp"
#include "Types.hpp"

namespace tk {

//! Compute gradient at a mesh node
std::array< tk::real, 3 >
nodegrad( std::size_t node,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup,
          const tk::Fields& U,
          uint64_t c );

//! Compute gradient at a mesh edge
std::array< tk::real, 3 >
edgegrad( const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::vector< std::size_t >& esued,
          const tk::Fields& U,
          uint64_t c );

} // tk::
