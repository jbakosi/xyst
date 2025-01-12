// *****************************************************************************
/*!
  \file      src/Physics/Lohner.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LohCG: Artificial compressibility solver for incompressible flow
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace lohner {

//! Compute divergence of a vector in all points
void
div( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const tk::Fields& U,
     std::vector< tk::real >& D,
     std::size_t pos = 0 );

//! Compute gradient of scalar in all points
void
grad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const std::vector< tk::real >& U,
      tk::Fields& G );

//! Compute velocity gradients in all points
void
vgrad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
       const std::array< std::vector< tk::real >, 3 >& dsupint,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::vector< std::size_t >& triinpoel,
       const tk::Fields& U,
       tk::Fields& G );

//! Compute momentum flux in all points
void
flux( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const tk::Fields& U,
      const tk::Fields& G,
      tk::Fields& F );

//! Compute gradient of scalar in all points
void
grad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const tk::Fields& U,
      tk::Fields& G );

//! Compute right hand side
void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const tk::Fields& U,
     const tk::Fields& G,
     tk::Fields& R );

} // lohner::
