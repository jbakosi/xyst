// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's input deck file parser
  \details   This file defines the input deck, i.e., control file, parser for
    the computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************

#include <ostream>
#include <type_traits>

#include "XystConfig.hpp"

#include "NoWarning/pegtl.hpp"

#include "Print.hpp"
#include "Tags.hpp"
#include "Inciter/Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/InputDeck/Parser.hpp"
#include "Inciter/InputDeck/Grammar.hpp"

using inciter::InputDeckParser;

InputDeckParser::InputDeckParser( const ctr::CmdLine& cmdline,
                                  ctr::InputDeck& inputdeck ) :
  FileParser( cmdline.get< tag::io, tag::control >() )
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line stack
//! \param[in,out] inputdeck Input deck stack where data is stored during
//!    parsing
// *****************************************************************************
{
  // Create InputDeck (a tagged tuple) to store parsed input
  ctr::InputDeck id( cmdline );

  // Parse input file and populate the underlying tagged tuple
  tao::pegtl::file_input<> in( m_filename );
  tao::pegtl::parse< deck::read_file, tk::grm::action >( in, id );

  // Echo errors and warnings accumulated during parsing
  diagnostics( tk::Print(), id.get< tag::error >() );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  inputdeck = std::move( id );
}
