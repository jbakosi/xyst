// *****************************************************************************
/*!
  \file      src/Physics/Operators.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible single-material flow using continuous Galerkin
  \details   Physics operators governing compressible single-material flow using
             continuous Galerkin finite elements.
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "Table.hpp"

namespace physics {

//! Compute nodal gradients of primitive variables in all points
void
grad( const std::vector< std::size_t >& dedge,
      const std::vector< tk::real >& deint,
      const std::vector< std::size_t >& bpoin,
      const std::vector< tk::real >& bpint,
      const std::vector< std::size_t >& bedge,
      const std::vector< tk::real >& beint,
      const tk::Fields& U,
      tk::Fields& G );

//! Compute right hand side
void
rhs( const std::vector< std::size_t >& dedge,
     const std::vector< tk::real >& deint,
     const std::vector< std::size_t >& bpoin,
     const std::vector< tk::real >& bpint,
     const std::vector< std::size_t >& bedge,
     const std::vector< tk::real >& beint,
     const std::vector< std::uint8_t >& bpsym,
     const std::vector< std::uint8_t >& besym,
     const tk::UnsMesh::Coords& coord,
     const tk::Fields& G,
     const tk::Fields& U,
     const std::vector< tk::real >& v,
     tk::real t,
     const std::vector< tk::real >& tp,
     tk::Fields& R );

//! Compute minimum time step size
tk::real
dt( const std::vector< tk::real >& vol, const tk::Fields& U );

//! Compute time step size for each mesh node (for steady time stepping)
void
dt( const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp );

} // physics::
