// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************

#include "FieldOutput.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

std::vector< std::string >
numericFieldNames()
// *****************************************************************************
// Collect field output names from numerical solution based on user input
//! \return Output field names requested by user
// *****************************************************************************
{
  std::vector< std::string > f;
  for (const auto& v : g_inputdeck.get< tag::cmd, tag::io, tag::outvar >()) {
    if (not v.analytic()) {
      std::stringstream s;
      if (v.alias.empty()) s << v; else s << v.alias;
      f.push_back( s.str() );
    }
  }

  return f;
}

std::vector< std::vector< tk::real > >
numericFieldOutput( const tk::Fields& U )
// *****************************************************************************
// Collect field output from numerical solution based on user input
//! \param[in] U Solution data to extract from
//! \return Output fields requested by user
// *****************************************************************************
{
  // Get offset map
  const auto& offset =
    g_inputdeck.get< tag::component >().offsetmap( g_inputdeck );

  std::vector< std::vector< tk::real > > f;
  for (const auto& v : g_inputdeck.get< tag::cmd, tag::io, tag::outvar >()) {
    auto o = tk::cref_find( offset, v.var );
    if (v.name.empty()) {        // depvar-based direct access
      f.push_back( U.extract( v.field, o ) );
    } else if (!v.analytic()) {  // human-readable non-analytic via custom fn
      Assert( v.getvar, "getvar() not configured for " + v.name );
      f.push_back( v.getvar( U, o, 1 ) );
    }
  }

  return f;
}

} // inciter::
