// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/BC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary condition options for inciter
  \details   Boundary condition options for inciter
*/
// *****************************************************************************
#ifndef InciterBCOptions_h
#define InciterBCOptions_h

#include <brigand/sequences/list.hpp>

#include "TaggedTuple.hpp"
#include "Toggle.hpp"
#include "Keywords.hpp"

namespace inciter {
namespace ctr {

//! Boundary condition types
enum class BCType : uint8_t { SYM
                            , FARFIELD
                            };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, BCType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class BC : public tk::Toggle< BCType > {

  public:
    // List valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::bc_sym
                                  , kw::bc_inlet
                                  , kw::bc_outlet
                                  , kw::bc_extrapolate
                                  , kw::bc_farfield
                                  >;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit BC() :
      tk::Toggle< BCType >( "Boundary condition",
        //! Enums -> names
        {   { BCType::SYM, kw::bc_sym::name() }
          , { BCType::FARFIELD, kw::bc_farfield::name() }
        },
        //! keywords -> Enums
        {   { kw::bc_sym::string(), BCType::SYM }
          , { kw::bc_farfield::string(), BCType::FARFIELD }
        } ) {}
};

} // ctr::
} // inciter::

#endif // InciterBCOptions_h
