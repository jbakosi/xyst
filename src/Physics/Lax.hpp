// *****************************************************************************
/*!
  \file      src/Physics/Lax.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     LaxCG: Time-derivative preconditioning for all Ma
  \see       Luo, Baum, Lohner, "Extension of Harten-Lax-van Leer Scheme for
             Flows at All Speeds", AIAA Journal, Vol. 43, No. 6, 2005
  \see       Weiss & Smith, "Preconditioning Applied to Variable and Constant
             Density Time-Accurate Flows on Unstructured Meshes", AIAA Journal,
             Vol. 33, No. 11, 1995, pp. 2050-2057.
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace lax {

//! Compute reference velocitity of the preconditioned system
tk::real
refvel( tk::real r, tk::real p, tk::real v );

//! Compute nodal gradients of primitive variables in all points
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
     const std::vector< std::uint8_t >& besym,
     const tk::Fields& G,
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R );

} // lax::
