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
