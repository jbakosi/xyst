// *****************************************************************************
/*!
  \file      src/Physics/Kozak.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     KozCG: Taylor-Galerkin, FCT, element-based continuous Galerkin
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"

namespace kozak {

//! Compute right hand side
void
rhs( const std::vector< std::size_t >& inpoel,
     const std::array< std::vector< tk::real >, 3 >& coord,
     tk::real t,
     tk::real dt,
     const std::vector< tk::real >& tp,
     const std::vector< tk::real >& dtp,
     const tk::Fields& U,
     tk::Fields& R );

} // kozak::
