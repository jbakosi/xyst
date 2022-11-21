// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************
#ifndef FieldOutput_h
#define FieldOutput_h

#include "Types.hpp"
#include "Fields.hpp"
#include "ContainerUtil.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Collect field output names from numerical solution based on user input
std::vector< std::string >
numericFieldNames();

//! Collect field output from numerical solution based on user input
std::vector< std::vector< tk::real > >
numericFieldOutput( const tk::Fields& U );

//! Collect field output names from analytic solutions based on user input
//! \tparam PDE Partial differential equation type
//! \param[in] eq PDE whose analytic solution field names to query
//! \param[in,out] f Output field names augmented
template< class PDE >
void
analyticFieldNames( const PDE& eq,
                    std::vector< std::string >& f )
{
  for (const auto& v : g_inputdeck.get< tag::cmd, tag::io, tag::outvar >())
    if (v.analytic()) tk::concat( eq.analyticFieldNames(), f );
}

//! Collect field output from analytic solutions based on user input
//! \tparam PDE Partial differential equation type
//! \param[in] eq PDE whose analytic solution to output
//! \param[in] x x coordinates at which to evaluate the analytic solution
//! \param[in] y y coordinates at which to evaluate the analytic solution
//! \param[in] z z coordinates at which to evaluate the analytic solution
//! \param[in] t Physical time at which to evaluate the analytic solution
//! \param[in,out] f Output fields augmented by analytic solutions requested
template< class PDE >
void
analyticFieldOutput( const PDE& eq,
                     const std::vector< tk::real >& x,
                     const std::vector< tk::real >& y,
                     const std::vector< tk::real >& z,
                     tk::real t,
                     std::vector< std::vector< tk::real > >& f )
{
  for (const auto& v : g_inputdeck.get< tag::cmd, tag::io, tag::outvar >()) {
    if (v.analytic()) {
      auto ncomp = eq.analyticSolution( x[0], y[0], z[0], t ).size();
      f.resize( f.size() + ncomp, std::vector< tk::real >( x.size() ) );
      for (std::size_t i=0; i<x.size(); ++i) {
        auto s = eq.analyticSolution( x[i], y[i], z[i], t );
        for (std::size_t j=0; j<ncomp; ++j) f[f.size()-ncomp+j][i] = s[j];
      }
    }
  }
}

} // inciter::

#endif // FieldOutput_h
