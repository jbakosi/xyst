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

namespace inciter {

//! Determine nodes that lie inside the user-defined IC box
void
ICBoxNodes( const tk::UnsMesh::Coords& coord,
            std::vector< std::unordered_set< std::size_t > >& inbox );

//! Query user config and assign function to set initial conditions
std::function< std::vector< tk::real >
             ( tk::real, tk::real, tk::real, tk::real ) >
IC();

//! Initalize the compressible flow equations, prepare for time integration
void
initialize( const std::array< std::vector< tk::real >, 3 >& coord,
            tk::Fields& U,
            tk::real t );

//! Compute nodal gradients of primitive variables along chare-boundary
void
bndgrad( const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& bndel,
         const std::vector< std::size_t >& gid,
         const std::unordered_map< std::size_t, std::size_t >& bid,
         const tk::Fields& U,
         tk::Fields& G );

//! Compute right hand side
void
rhs( tk::real t,
     const std::array< std::vector< tk::real >, 3 >& coord,
     const std::vector< std::size_t >& inpoel,
     const std::vector< std::size_t >& triinpoel,
     const std::vector< std::size_t >& gid,
     const std::unordered_map< std::size_t, std::size_t >& bid,
     const std::unordered_map< std::size_t, std::size_t >& lid,
     const std::vector< tk::real >& dfn,
     const std::pair< std::vector< std::size_t >,
                      std::vector< std::size_t > >& esup,
     const std::pair< std::vector< std::size_t >,
                      std::vector< std::size_t > >& psup,
     const std::vector< int >& symbctri,
     const std::unordered_set< std::size_t >& spongenodes,
     const std::vector< tk::real >& vol,
     const std::vector< std::size_t >& edgenode,
     const std::vector< std::size_t >& edgeid,
     const std::vector< tk::real >& tp,
     const tk::Fields& bG,
     const tk::Fields& U,
     tk::Fields& R );

//! Compute the minimum time step size (for unsteady time stepping)
tk::real
dt( const std::array< std::vector< tk::real >, 3 >& coord,
    const std::vector< std::size_t >& inpoel,
    tk::real t,
    const tk::Fields& U );

//! Compute a time step size for each mesh node (for steady time stepping)
void
dt( uint64_t,
    const std::vector< tk::real >& vol,
    const tk::Fields& U,
    std::vector< tk::real >& dtp );

//! Query Dirichlet boundary condition value on a given side set
std::map< std::size_t, std::vector< std::pair< bool, tk::real > > >
dirbc( tk::real t,
       tk::real deltat,
       const std::vector< tk::real >& tp,
       const std::vector< tk::real >& dtp,
       const std::pair< const int, std::vector< std::size_t > >& ss,
       const std::array< std::vector< tk::real >, 3 >& coord );

//! Set symmetry boundary conditions at nodes
void
symbc( tk::Fields& U,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::unordered_map< int,
         std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
       const std::unordered_set< std::size_t >& nodes );

//! Set farfield boundary conditions at nodes
void
farfieldbc(
  tk::Fields& U,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const std::unordered_map< int,
    std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
  const std::unordered_set< std::size_t >& nodes );

//! Apply sponge conditions at sponge nodes
void
sponge( tk::Fields& U,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::unordered_set< std::size_t >& nodes );

//! Apply user defined time dependent BCs
void
timedepbc( tk::real t,
  tk::Fields& U,
  const std::vector< std::unordered_set< std::size_t > >& nodes,
  const std::vector< tk::Table<5> >& timedepfn );

} // physics::
