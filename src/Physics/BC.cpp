// *****************************************************************************
/*!
  \file      src/Physics/BC.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions
*/
// *****************************************************************************

#include "BC.hpp"
#include "Tags.hpp"
#include "EOS.hpp"
#include "Problems.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

namespace physics {

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
  auto ncomp = U.nprop();
  auto nmask = ncomp + 1;

  Assert( dirbcmasks.size() % nmask == 0, "Dirichlet BC masks size mismatch" );

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
  using inciter::g_inputdeck;

  const auto& sbc =
    g_inputdeck.get< tag::param, tag::compflow, tag::bc, tag::symmetry >();
  if (sbc.empty()) return;

  Assert( symbcnodes.size()*3 == symbcnorms.size(), "Size mismaatch" );

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
  using inciter::g_inputdeck;

  const auto& compflow = g_inputdeck.get< tag::param, tag::compflow >();

  const auto& fbc = compflow.get< tag::bc, tag::farfield >();
  if (fbc.empty()) return;

  Assert( farbcnodes.size() == farbcnorms.size()*3, "Size mismaatch" );

  const auto& fre = compflow.get< tag::farfield_density >();
  ErrChk( !fre.empty(), "No farfield density specified" );
  tk::real fr = fre[0];

  const auto& fue = compflow.get< tag::farfield_velocity >();
  ErrChk( !fue.empty(), "No farfield velocity specified" );
  tk::real fu = fue[0][0];
  tk::real fv = fue[0][1];
  tk::real fw = fue[0][2];

  const auto& fpe = compflow.get< tag::farfield_pressure >();
  ErrChk( !fpe.empty(), "No farfield pressure specified" );
  tk::real fp = fpe[0];

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
    //auto a = eos::soundspeed( r, eos::pressure( r, ru/r, rv/r, rw/r, re ) );
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
      auto pr = eos::pressure( r, ru/r, rv/r, rw/r, re );
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

} // physics::
