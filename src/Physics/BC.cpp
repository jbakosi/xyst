// *****************************************************************************
/*!
  \file      src/Physics/BC.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions
*/
// *****************************************************************************

#include "BC.hpp"
#include "EOS.hpp"
#include "Box.hpp"
#include "Problems.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace physics {

using inciter::g_cfg;

void
dirbc( std::size_t meshid,
       tk::Fields& U,
       tk::real t,
       const std::array< std::vector< tk::real >, 3 >& coord,
       const std::vector< std::unordered_set< std::size_t > >& boxnodes,
       const std::vector< std::size_t >& dirbcmask,
       const std::vector< double >& dirbcval )
// *****************************************************************************
//  Set Dirichlet boundary conditions at nodes
//! \param[in] meshid Mesh id to use
//! \param[in] U Solution vector at recent time step
//! \param[in] t Physical time at which to evaluate BCs
//! \param[in] coord Mesh node coordinates
//! \param[in] boxnodes List of nodes at which box user ICs are set
//! \param[in] dirbcmask Nodes and component masks for Dirichlet BCs
//! \param[in] dirbcval Nodes and component values for Dirichlet BCs
// *****************************************************************************
{
  auto ncomp = U.nprop();
  auto nmask = ncomp + 1;

  Assert( dirbcmask.size() % nmask == 0, "Size mismatch" );
  Assert( dirbcval.size() % nmask == 0, "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  auto ic = problems::IC();

  for (std::size_t i=0; i<dirbcmask.size()/nmask; ++i) {
    auto p = dirbcmask[i*nmask+0];      // local node id
    auto u = ic( x[p], y[p], z[p], t, meshid ); // evaluate solution/ic
    problems::box( p, u, boxnodes );    // overwrite with box value
    for (std::size_t c=0; c<ncomp; ++c) {
      auto mask = dirbcmask[i*nmask+1+c];
      if (mask == 1) {                              // mask == 1: IC+box value
        U(p,c) = u[c];
      } else if (mask == 2 && !dirbcval.empty()) {  // mask == 2: BC value
        U(p,c) = dirbcval[i*nmask+1+c];
      }
    }
  }
}

void
dirbcp( std::size_t meshid,
        tk::Fields& U,
        const std::array< std::vector< tk::real >, 3 >& coord,
        const std::vector< std::size_t >& dirbcmaskp,
        const std::vector< double >& dirbcvalp )
// *****************************************************************************
//  Set pressure Dirichlet boundary conditions at nodes
//! \param[in] meshid Mesh id to use
//! \param[in] U Solution vector at recent time step
//! \param[in] coord Mesh node coordinates
//! \param[in] dirbcmaskp Nodes and component masks for Dirichlet BCs
//! \param[in] dirbcvalp Nodes and component values for Dirichlet BCs
// *****************************************************************************
{
  std::size_t nmask = 1 + 1;

  Assert( dirbcmaskp.size() % nmask == 0, "Size mismatch" );
  Assert( dirbcvalp.size() % nmask == 0, "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  auto ic = problems::PRESSURE_IC();

  for (std::size_t i=0; i<dirbcmaskp.size()/nmask; ++i) {
    auto p = dirbcmaskp[i*nmask+0];      // local node id
    auto mask = dirbcmaskp[i*nmask+1];
    if (mask == 1) {                               // mask == 1: IC value
      U(p,0) = ic( x[p], y[p], z[p], meshid );
    } else if (mask == 2 && !dirbcvalp.empty()) {  // mask == 2: BC value
      U(p,0) = dirbcvalp[i*nmask+1];
    }
  }
}

void
symbc( tk::Fields& U,
       const std::vector< std::size_t >& symbcnodes,
       const std::vector< tk::real >& symbcnorms,
       std::size_t pos )
// *****************************************************************************
//  Set symmetry boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] symbcnodes Node ids at which to set symmetry BCs
//! \param[in] symbcnorms Normals at nodes at which to set symmetry BCs
//! \param[in] pos Position at which the three velocity components are in U
// *****************************************************************************
{
  Assert( symbcnodes.size()*3 == symbcnorms.size(), "Size mismatch" );

  for (std::size_t i=0; i<symbcnodes.size(); ++i) {
    auto p = symbcnodes[i];
    auto n = symbcnorms.data() + i*3;
    auto& u = U(p,pos+0);
    auto& v = U(p,pos+1);
    auto& w = U(p,pos+2);
    auto vn = u*n[0] + v*n[1] + w*n[2];
    u -= vn * n[0];
    v -= vn * n[1];
    w -= vn * n[2];
  }
}

void
noslipbc( tk::Fields& U,
          const std::vector< std::size_t >& noslipbcnodes,
          std::size_t pos )
// *****************************************************************************
//  Set noslip boundary conditions at nodes
//! \param[in] U Solution vector at recent time step
//! \param[in] noslipbcnodes Node ids at which to set noslip BCs
//! \param[in] pos Position at which the three velocity components are in U
// *****************************************************************************
{
  for (auto p : noslipbcnodes) U(p,pos+0) = U(p,pos+1) = U(p,pos+2) = 0.0;
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
  const auto& tf = g_cfg.get< tag::bc_far >();
  if (tf.get< tag::sidesets >().empty()) return;

  Assert( farbcnodes.size()*3 == farbcnorms.size(), "Size mismatch" );

  // cppcheck-suppress unreadVariable
  tk::real fr = tf.get< tag::density >();

  const auto& fue = tf.get< tag::velocity >();
  ErrChk( !fue.empty(), "No farfield velocity specified" );
  // cppcheck-suppress unreadVariable
  tk::real fu = fue[0];
  // cppcheck-suppress unreadVariable
  tk::real fv = fue[1];
  // cppcheck-suppress unreadVariable
  tk::real fw = fue[2];

  tk::real fp = tf.get< tag::pressure >();

  for (std::size_t i=0; i<farbcnodes.size(); ++i) {
    auto p  = farbcnodes[i];
    auto nx = farbcnorms[i*3+0];
    auto ny = farbcnorms[i*3+1];
    auto nz = farbcnorms[i*3+2];
    auto& r  = U(p,0);
    auto& ru = U(p,1);
    auto& rv = U(p,2);
    auto& rw = U(p,3);
    auto& re = U(p,4);
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
    U(p,0) = prebcvals[i*2+0];
    U(p,4) = eos::totalenergy( U(p,0), U(p,1)/U(p,0), U(p,2)/U(p,0),
                               U(p,3)/U(p,0), prebcvals[i*2+1] );
  }
}

} // physics::
