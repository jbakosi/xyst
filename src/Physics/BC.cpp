// *****************************************************************************
/*!
  \file      src/Physics/BC.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions
*/
// *****************************************************************************

#include "BC.hpp"
#include "EOS.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace physics {

using inciter::g_cfg;

void
dirbc( tk::Fields& U,
       tk::real t,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::vector< std::size_t >& dirbcmasks )
// *****************************************************************************
//  Set symmetry boundary conditions at nodes
//! \param[in] t Physical time at which to evaluate BCs
//! \param[in] U Solution vector at recent time step
//! \param[in] coord Mesh node coordinates
//! \param[in] dirbcmasks Nodes and component masks for Dirichlet BCs
// *****************************************************************************
{
  if (g_cfg.get< tag::bc_dir >().empty()) return;

  auto ncomp = U.nprop();
  auto nmask = ncomp + 1;

  Assert( dirbcmasks.size() % nmask == 0, "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  auto ic = problems::IC();

  for (std::size_t i=0; i<dirbcmasks.size()/nmask; ++i) {
    auto p = dirbcmasks[i*nmask+0];
    auto u = ic( x[p], y[p], z[p], t );
    for (std::size_t c=0; c<ncomp; ++c)
      if (dirbcmasks[i*nmask+1+c])
        U(p,c,0) = u[c];
  }
}

void
symbc( tk::Fields& U,
       const std::vector< std::size_t >& symbcnodes,
       const std::vector< tk::real >& symbcnorms )
// *****************************************************************************
//  Set symmetry boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] symbcnodes Node ids at which to set symmetry BCs
//! \param[in] symbcnorms Normals at nodes at which to set symmetry BCs
// *****************************************************************************
{
  if (g_cfg.get< tag::bc_sym >().empty()) return;

  Assert( symbcnodes.size()*3 == symbcnorms.size(), "Size mismatch" );

  for (std::size_t i=0; i<symbcnodes.size(); ++i) {
    auto p  = symbcnodes[i];
    auto nx = symbcnorms[i*3+0];
    auto ny = symbcnorms[i*3+1];
    auto nz = symbcnorms[i*3+2];
    auto rvn = U(p,1,0)*nx + U(p,2,0)*ny + U(p,3,0)*nz;
    U(p,1,0) -= rvn * nx;
    U(p,2,0) -= rvn * ny;
    U(p,3,0) -= rvn * nz;
  }
}

void
farbc( tk::Fields& U,
       const std::vector< std::size_t >& farbcnodes,
       const std::vector< tk::real >& farbcnorms )
// *****************************************************************************
//  Set farfield boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] farbcnodes Nodes ids at which to set farfield BCs
//! \param[in] farbcnorms Normals at nodes at which to set farfield BCs
// *****************************************************************************
{
  if (g_cfg.get< tag::bc_far >().empty()) return;

  Assert( farbcnodes.size()*3 == farbcnorms.size(), "Size mismatch" );

  // cppcheck-suppress unreadVariable
  tk::real fr = g_cfg.get< tag::bc_far_density >();

  const auto& fue = g_cfg.get< tag::bc_far_velocity >();
  ErrChk( !fue.empty(), "No farfield velocity specified" );
  // cppcheck-suppress unreadVariable
  tk::real fu = fue[0];
  // cppcheck-suppress unreadVariable
  tk::real fv = fue[1];
  // cppcheck-suppress unreadVariable
  tk::real fw = fue[2];

  tk::real fp = g_cfg.get< tag::bc_far_pressure >();

  for (std::size_t i=0; i<farbcnodes.size(); ++i) {
    auto p  = farbcnodes[i];
    auto nx = farbcnorms[i*3+0];
    auto ny = farbcnorms[i*3+1];
    auto nz = farbcnorms[i*3+2];
    auto& r  = U(p,0,0);
    auto& ru = U(p,1,0);
    auto& rv = U(p,2,0);
    auto& rw = U(p,3,0);
    auto& re = U(p,4,0);
    //auto vn = (ru*nx + rv*ny + rw*nz)/r;
    auto vn = fu*nx + fv*ny + fw*nz;
    //auto a = eos::soundspeed( r,
    //           eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r ) );
    auto a = eos::soundspeed( fr, fp );
    auto M = vn / a;
    if (M <= -1.0) {
      // supersonic inflow, all characteristics from outside
      r  = fr;
      ru = fr * fu;
      rv = fr * fv;
      rw = fr * fw;
      re = eos::totalenergy( fr, fu, fv, fw, fp );
    } else if (M > -1.0 && M < 0.0) {
      // subsonic inflow: 1 outgoing and 4 incoming characteristics,
      // pressure from inside, rest from outside
      auto pr = eos::pressure( re - 0.5*(ru*ru + rv*rv + rw*rw)/r );
      r  = fr;
      ru = fr * fu;
      rv = fr * fv;
      rw = fr * fw;
      re = eos::totalenergy( fr, fu, fv, fw, pr );
    } else if (M >= 0.0 && M < 1.0) {
      // subsonic outflow: 1 incoming and 4 outgoing characteristics,
      // pressure from outside, rest from inside
      re = eos::totalenergy( r, ru/r, rv/r, rw/r, fp );
    }
  }
}

void
prebc( tk::Fields& U,
       const std::vector< std::size_t >& prebcnodes,
       const std::vector< tk::real >& prebcvals )
// *****************************************************************************
//  Set pressure boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] prebcnodes Node ids at which to set pressure BCs
//! \param[in] prebcvals Density and pressure values at pressure BC nodes
// *****************************************************************************
{
  Assert( prebcnodes.size()*2 == prebcvals.size(), "Size mismatch" );

  for (std::size_t i=0; i<prebcnodes.size(); ++i) {
    auto p = prebcnodes[i];
    U(p,0,0) = prebcvals[i*2+0];
    U(p,4,0) = eos::totalenergy( U(p,0,0), U(p,1,0)/U(p,0,0),
                 U(p,2,0)/U(p,0,0), U(p,3,0)/U(p,0,0), prebcvals[i*2+1] );
  }
}

} // physics::
