// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Scheme.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Discretization scheme options for inciter
  \details   Discretization scheme options for inciter
*/
// *****************************************************************************
#ifndef SchemeOptions_h
#define SchemeOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Scheme types
enum class SchemeType : uint8_t { DiagCG
                                , ALECG
                                };

//! Pack/Unpack SchemeType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, SchemeType& e ) { PUP::pup( p, e ); }

//! \brief Scheme options: outsource to base templated on enum type
class Scheme : public tk::Toggle< SchemeType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::diagcg
                                  , kw::alecg
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Scheme() :
      tk::Toggle< SchemeType >(
        //! Group, i.e., options, name
        kw::scheme::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { SchemeType::DiagCG, kw::diagcg::name() },
          { SchemeType::ALECG, kw::alecg::name() }
        },
        //! keywords -> Enums
        { { kw::diagcg::string(), SchemeType::DiagCG }
        , { kw::alecg::string(), SchemeType::ALECG }
        } ) {}
};

} // ctr::
} // inciter::

#endif // SchemeOptions_h
