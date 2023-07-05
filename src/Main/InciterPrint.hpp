// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef InciterPrint_h
#define InciterPrint_h

#include <iostream>
#include <string>

#include "NoWarning/format.hpp"

#include "Print.hpp"
#include "ContainerUtil.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! InciterPrint : tk::Print
class InciterPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[in] screen Screen output filename
    //! \param[in,out] str Verbose stream
    //! \param[in] mode Open mode for screen output file, see
    //!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    explicit InciterPrint( const std::string& screen,
                           std::ostream& str = std::clog,
                           std::ios_base::openmode mode = std::ios_base::out,
                           std::ostream& qstr = std::cout ) :
      Print( screen, str, mode, qstr ) {}

    //! Print control option: 'group : option'
    template< typename Option, typename... tags >
    void Item() const {
      Option opt;
      m_stream << m_item_name_value_fmt
                  % m_item_indent % opt.group()
                  % opt.name( g_inputdeck.get< tags... >() );
    }

    //! Print list of codes of vector-valued option
    //! \tparam Option Option type
    //! \tparam T Enum type for option type
    //! \param[in] v Vector of option types (enums) whose code vector to print
    template< typename Option, typename T >
    void ItemVec( const std::vector< T >& v ) const {
      Option opt;
      std::string codes;
      for (auto e : v) codes += opt.code(e);
      item( opt.group(), codes );
    }

    //! Print time integration header
    void inthead( const std::string& t, const std::string& name,
                  const std::string& legend, const std::string& head ) const;

    //! Print initial mesh refinement edge-node pairs
    void edgeref( const std::vector< std::size_t >& edgenodes ) const;
};

} // inciter::

#endif // InciterPrint_h
