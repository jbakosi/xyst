// *****************************************************************************
/*!
  \file      src/Physics/Box.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Initial conditions box related functionality
*/
// *****************************************************************************

#include "Box.hpp"
#include "EOS.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // ::inciter

namespace problems {

using inciter::g_cfg;

std::vector< std::unordered_set< std::size_t > >
boxnodes( const std::array< std::vector< tk::real >, 3 >& coord )
// *****************************************************************************
//  Determine nodes that lie inside user-defined IC box(es)
//! \param[in] coord Mesh node coordinates
//! \return inbox List of nodes at which box user ICs are set for each IC box
// *****************************************************************************
{
  const auto& icbox = g_cfg.get< tag::ic, tag::boxes >();

  if (icbox.empty()) return {};

  std::vector< std::unordered_set< std::size_t > > inbox;

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  std::size_t bcnt = 0;
  for (const auto& b : icbox) {
    inbox.emplace_back();
    const auto& bx = b.get< tag::box_x >();
    const auto& by = b.get< tag::box_y >();
    const auto& bz = b.get< tag::box_z >();
    std::vector< tk::real > box{ bx[0], bx[1], by[0], by[1], bz[0], bz[1] };

    const auto eps = std::numeric_limits< tk::real >::epsilon();
    // Determine which nodes lie in the IC box
    if ( std::any_of( begin(box), end(box), [=](auto p)
                      { return abs(p) > eps; } ) )
    {
      std::array< tk::real, 3 > b_min{{box[0], box[2], box[4]}};
      std::array< tk::real, 3 > b_max{{box[1], box[3], box[5]}};
      for (std::size_t i=0; i<x.size(); ++i) {
        std::array< tk::real, 3 > node{{ x[i], y[i], z[i] }};
        if ( node[0]>b_min[0] && node[0]<b_max[0] &&
             node[1]>b_min[1] && node[1]<b_max[1] &&
             node[2]>b_min[2] && node[2]<b_max[2] )
        {
          inbox[bcnt].insert( i );
        }
      }
    }
    ++bcnt;
  }

  return inbox;
}

void
box( std::size_t p,
     std::vector< tk::real >& u,
     const std::vector< std::unordered_set< std::size_t > >& boxnodes )
// *****************************************************************************
//  Evaluate solution in user-defined IC box
//! \param[in] p Local mesh node id at which to evaluate
//! \param[in,out] u Solution to overwrite with box value
//! \param[in] boxnodes Nodes at which box user ICs are set (for each box IC)
// *****************************************************************************
{
  const auto& icbox = g_cfg.get< tag::ic, tag::boxes >();
  if (icbox.empty()) return;

  auto large = std::numeric_limits< double >::max();

  for (std::size_t j=0; j<icbox.size(); ++j) {
    const auto& b = icbox[j];
    if (boxnodes.size() > j && boxnodes[j].find(p) != boxnodes[j].end()) {
      auto boxr = b.get< tag::box_density >();
      const auto& boxv = b.get< tag::box_velocity >();
      auto boxp = b.get< tag::box_pressure >();
      auto boxe = b.get< tag::box_energy >();
      auto boxt = b.get< tag::box_temperature >();

      tk::real r = 0.0, ru = 0.0, rv = 0.0, rw = 0.0, re = 0.0;
      if (std::abs(boxr - large) > 1.0e-12 && boxr > 0.0) {
        r = boxr;
      }
      if (std::abs(boxv[0] - large) > 1.0e-12 &&
          std::abs(boxv[1] - large) > 1.0e-12 &&
          std::abs(boxv[2] - large) > 1.0e-12)
      {
        ru = r * boxv[0];
        rv = r * boxv[1];
        rw = r * boxv[2];
      }
      if (std::abs(boxp - large) > 1.0e-12 && boxp> 0.0) {
        re = eos::totalenergy( r, ru/r, rv/r, rw/r, boxp);
      }
      if (std::abs(boxe - large) > 1.0e-12 && boxe > 0.0) {
        auto ux = ru/r, uy = rv/r, uz = rw/r;
        auto ke = 0.5*(ux*ux + uy*uy + uz*uz);
        re = r * (boxe + ke);
      }
      if (std::abs(boxt - large) > 1.0e-12 && boxt > 0.0) {
        auto cv = g_cfg.get< tag::mat_spec_heat_const_vol >();
        if (std::abs(cv - large) > 1.0e-12 && cv > 0.0) {
          re = r * boxt * cv;
        }
      }
      u[0] = r;
      u[1] = ru;
      u[2] = rv;
      u[3] = rw;
      u[4] = re;
    }
  }
}

} // problems::
