// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************

#include <regex>

#include <brigand/algorithms/for_each.hpp>

#include "InciterPrint.hpp"

using inciter::InciterPrint;

void
InciterPrint::inthead( const std::string& t,
                       const std::string& name,
                       const std::string& legend,
                       const std::string& head ) const
// *****************************************************************************
//  Print time integration header
//! \param[in] t Section title
//! \param[in] name Section name
//! \param[in] legend Legend to print
//! \param[in] head Head to append
// *****************************************************************************
{
  section( t, name );
  std::string l( legend );
  l = std::regex_replace( l, std::regex("\n"), "\n" + m_item_indent );
  raw( m_item_indent + l + head );
}

void InciterPrint::refvar( const std::vector< std::string >& rvar,
                           const std::vector< std::size_t >& refidx ) const
// *****************************************************************************
// Print mesh refinement variables and their indices in the unknown vector
//! \param[in] rvar Refinement variable name list
//! \param[in] refidx Refinement variable index (location in data array) list
// *****************************************************************************
{
  Assert( rvar.size() == refidx.size(), "Size mismatch" );

  if (rvar.empty()) return;

  std::string c;
  for (std::size_t i=0; i<rvar.size(); ++i)
    c += rvar[i] + '[' + std::to_string(refidx[i]) + "] ";
  auto name = kw::amr_refvar::name() + " & id(s)";
  name[0] = static_cast< char >( std::toupper( name[0] ) );
  item( name, c );
}

void InciterPrint::edgeref( const std::vector< std::size_t >& edgenodes ) const
// *****************************************************************************
// Print initial mesh refinement edge-node pairs
// *****************************************************************************
{
   if (edgenodes.empty()) return;

   std::string c;
   for (auto i : edgenodes) c += std::to_string(i) + ' ';
   auto name = kw::amr_edgelist::name();
   name[0] = static_cast< char >( std::toupper( name[0] ) );
   item( name, c );
}
