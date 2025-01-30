// *****************************************************************************
/*!
  \file      src/Physics/Chorin.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     ChoCG: Projection-based solver for incompressible flow
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace chorin {

//! Compute divergence of a vector in all points
void
div( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     tk::real dt,
     const std::vector< tk::real >& P,
     const tk::Fields& G,
     const tk::Fields& U,
     std::vector< tk::real >& D,
     bool stab );

//! Compute gradient of scalar in all points
void
grad( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& triinpoel,
      const std::vector< tk::real >& U,
      tk::Fields& G );

//  Compute velocity+scalar gradients in all points
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

//! Compute right hand side
void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& P,
     const tk::Fields& U,
     const tk::Fields& G,
     const tk::Fields& W,
     tk::Fields& R );

} // chorin::
