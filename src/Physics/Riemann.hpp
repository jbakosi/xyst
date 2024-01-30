// *****************************************************************************
/*!
  \file      src/Physics/Riemann.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Riemann, MUSCL, limiting for edge-based continuous Galerkin
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace riemann {

//! Compute nodal gradients of primitive variables in all points
void
grad( const std::vector< std::size_t >& bpoin,
      const std::vector< tk::real >& bpint,
      const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< std::size_t >, 2 >& bsupedge,
      const std::array< std::vector< tk::real >, 2 >& bsupint,
      const tk::Fields& U,
      tk::Fields& G );

//! Compute right hand side
void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< std::size_t >, 2 >& bsupedge,
     const std::array< std::vector< tk::real >, 2 >& bsupint,
     const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::uint8_t >& bpsym,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const tk::Fields& G,
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R );

} // riemann::
