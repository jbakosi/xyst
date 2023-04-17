// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem options for inciter
  \details   Problem options for inciter
*/
// *****************************************************************************
#ifndef ProblemOptions_h
#define ProblemOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Problem types
enum class ProblemType : uint8_t { USER_DEFINED
                                 , VORTICAL_FLOW
                                 , TAYLOR_GREEN
                                 , NONLIN_ENER_GROWTH
                                 , RAYLEIGH_TAYLOR
                                 , SOD
                                 , ROTATED_SOD
                                 , SEDOV
                                 , SLOT_CYL
                                 , POINT_SRC
                                 , SHEAR_DIFF
                                 };

//! Pack/Unpack ProblemType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, ProblemType& e ) { PUP::pup( p, e ); }

//! \brief Problem options: outsource to base templated on enum type
class Problem : public tk::Toggle< ProblemType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = tk::unique_codes< kw::user_defined
                                     , kw::vortical_flow
                                     , kw::taylor_green
                                     , kw::nonlin_ener_growth
                                     , kw::rayleigh_taylor
                                     , kw::sod
                                     , kw::rotated_sod
                                     , kw::sedov
                                     , kw::slot_cyl
                                     , kw::point_src
                                     , kw::shear_diff
                                     >::list;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Problem() :
      tk::Toggle< ProblemType >(
        //! Group, i.e., options, name
        kw::problem::name(),
        //! Enums -> names
        { { ProblemType::USER_DEFINED, kw::user_defined::name() }
        , { ProblemType::VORTICAL_FLOW, kw::vortical_flow::name() }
        , { ProblemType::TAYLOR_GREEN, kw::taylor_green::name() }
        , { ProblemType::NONLIN_ENER_GROWTH, kw::nonlin_ener_growth::name()}
        , { ProblemType::RAYLEIGH_TAYLOR, kw::rayleigh_taylor::name() }
        , { ProblemType::SOD, kw::sod::name() }
        , { ProblemType::ROTATED_SOD, kw::rotated_sod::name() }
        , { ProblemType::SEDOV, kw::sedov::name() }
        , { ProblemType::SLOT_CYL, kw::slot_cyl::name() }
        , { ProblemType::POINT_SRC, kw::point_src::name() }
        , { ProblemType::SHEAR_DIFF, kw::shear_diff::name() }
        },
        //! keywords -> Enums
        { { kw::user_defined::string(), ProblemType::USER_DEFINED }
        , { kw::vortical_flow::string(), ProblemType::VORTICAL_FLOW }
        , { kw::taylor_green::string(), ProblemType::TAYLOR_GREEN }
        , { kw::nonlin_ener_growth::string(), ProblemType::NONLIN_ENER_GROWTH }
        , { kw::rayleigh_taylor::string(), ProblemType::RAYLEIGH_TAYLOR }
        , { kw::sod::string(), ProblemType::SOD }
        , { kw::rotated_sod::string(), ProblemType::ROTATED_SOD }
        , { kw::sedov::string(), ProblemType::SEDOV }
        , { kw::slot_cyl::string(), ProblemType::SLOT_CYL }
        , { kw::point_src::string(), ProblemType::POINT_SRC }
        , { kw::shear_diff::string(), ProblemType::SHEAR_DIFF }
        } )
    {
      brigand::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the problem option requested
    //! \return Policy code of the option
    const std::string& code( ProblemType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for problem \"") << p <<
                "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( brigand::type_<U> ) {
        static_assert( tk::HasTypedef_code_v< typename U::info >,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< ProblemType, std::string > policy {
        { ProblemType::USER_DEFINED, *kw::user_defined::code() }
      , { ProblemType::VORTICAL_FLOW, *kw::vortical_flow::code() }
      , { ProblemType::TAYLOR_GREEN, *kw::taylor_green::code() }
      , { ProblemType::NONLIN_ENER_GROWTH, *kw::nonlin_ener_growth::code() }
      , { ProblemType::RAYLEIGH_TAYLOR, *kw::rayleigh_taylor::code() }
      , { ProblemType::SOD, *kw::sod::code() }
      , { ProblemType::ROTATED_SOD, *kw::rotated_sod::code() }
      , { ProblemType::SEDOV, *kw::sedov::code() }
      , { ProblemType::SLOT_CYL, *kw::slot_cyl::code() }
      , { ProblemType::POINT_SRC, *kw::point_src::code() }
      , { ProblemType::SHEAR_DIFF, *kw::shear_diff::code() }
    };
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
