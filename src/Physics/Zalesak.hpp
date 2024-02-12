// *****************************************************************************
/*!
  \file      src/Physics/Zalesak.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Zalesak, FCT limiting for edge-based continuous Galerkin
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace zalesak {

//! Compute nodal gradients in all points
void
grad( const std::vector< std::size_t >& bpoin,
      const std::vector< tk::real >& bpint,
      const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& dsupint,
      const std::array< std::vector< std::size_t >, 2 >& bsupedge,
      const std::array< std::vector< tk::real >, 2 >& bsupint,
      const tk::Fields& U,
      tk::Fields& G );

//! Compute stabilization coefficient in all points
void
stab( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const tk::Fields& U,
      const tk::Fields& G,
      tk::Fields& S );

//! Compute right hand side
void
rhs( const std::array< std::vector< std::size_t >, 3 >& dsupedge,
     const std::array< std::vector< tk::real >, 3 >& dsupint,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::uint8_t >& besym,
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     const std::vector< tk::real >& dtp,
     const tk::Fields& U,
     const tk::Fields& S,
     const tk::Fields& G,
     tk::Fields& R );

} // zalesak::
