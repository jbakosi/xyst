// *****************************************************************************
/*!
  \file      src/Control/InciterInputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's input deck definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the control file parsing of the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#pragma once

#include "InciterCmdLine.hpp"

namespace tag {
DEFTAG( cmd );
DEFTAG( nstep );
DEFTAG( ttyi );
DEFTAG( term );
DEFTAG( cfl );
DEFTAG( t0 );
DEFTAG( dt );
DEFTAG( reorder );
DEFTAG( steady );
DEFTAG( residual );
DEFTAG( rescomp );
DEFTAG( problem );
DEFTAG( problem_ncomp );
DEFTAG( problem_alpha );
DEFTAG( problem_kappa );
DEFTAG( problem_beta );
DEFTAG( problem_r0 );
DEFTAG( problem_p0 );
DEFTAG( problem_ce );
DEFTAG( problem_src );
DEFTAG( location );
DEFTAG( radius );
DEFTAG( release_time );
DEFTAG( part );
DEFTAG( fieldout );
DEFTAG( fieldout_iter );
DEFTAG( fieldout_time );
DEFTAG( fieldout_range );
DEFTAG( histout );
DEFTAG( histout_iter );
DEFTAG( histout_time );
DEFTAG( histout_range );
DEFTAG( histout_precision );
DEFTAG( histout_format );
DEFTAG( integout );
DEFTAG( integout_iter );
DEFTAG( integout_time );
DEFTAG( integout_range );
DEFTAG( integout_precision );
DEFTAG( integout_format );
DEFTAG( diag_iter );
DEFTAG( diag_precision );
DEFTAG( diag_format );
DEFTAG( ic );
DEFTAG( x );
DEFTAG( y );
DEFTAG( z );
DEFTAG( ic_density );
DEFTAG( ic_pressure );
DEFTAG( ic_energy );
DEFTAG( ic_temperature );
DEFTAG( ic_velocity );
DEFTAG( bc_dir );
DEFTAG( bc_sym );
DEFTAG( bc_far );
DEFTAG( bc_far_density );
DEFTAG( bc_far_pressure );
DEFTAG( bc_far_velocity );
DEFTAG( bc_pre );
DEFTAG( bc_pre_density );
DEFTAG( bc_pre_pressure );
DEFTAG( mat_spec_heat_ratio );
DEFTAG( mat_spec_heat_const_vol );
DEFTAG( mat_heat_conductivity );
DEFTAG( href_t0 );
DEFTAG( href_dt );
DEFTAG( href_dtfreq );
DEFTAG( href_maxlevels );
DEFTAG( href_error );
DEFTAG( href_init );
DEFTAG( href_refvar );
} // tag::

namespace inciter {
namespace ctr {

//! Member data for tagged tuple
using InputDeckMembers = brigand::list<
    tag::cmd, CmdLine
  , tag::nstep, uint64_t
  , tag::ttyi, uint64_t
  , tag::term, double
  , tag::cfl, double
  , tag::t0, double
  , tag::dt, double
  , tag::reorder, bool
  , tag::steady, bool
  , tag::residual, double
  , tag::rescomp, uint64_t
  , tag::problem, std::string
  , tag::problem_ncomp, uint64_t
  , tag::problem_alpha, double
  , tag::problem_kappa, double
  , tag::problem_beta, std::vector< double >
  , tag::problem_r0, double
  , tag::problem_p0, double
  , tag::problem_ce, double
  , tag::problem_src, tk::TaggedTuple< brigand::list<
                        tag::location, std::vector< double >
                          , tag::radius, double
                          , tag::release_time, double
                      > >
  , tag::part, std::string
  , tag::fieldout, std::vector< int >
  , tag::fieldout_iter, uint64_t
  , tag::fieldout_time, double
  , tag::fieldout_range, std::vector< std::vector< double > >
  , tag::histout, std::vector< std::vector< double > >
  , tag::histout_iter, uint64_t
  , tag::histout_time, double
  , tag::histout_range, std::vector< std::vector< double > >
  , tag::histout_precision, std::streamsize
  , tag::histout_format, std::string
  , tag::integout, std::vector< int >
  , tag::integout_iter, uint64_t
  , tag::integout_time, double
  , tag::integout_range, std::vector< std::vector< double > >
  , tag::integout_precision, std::streamsize
  , tag::integout_format, std::string
  , tag::diag_iter, uint64_t
  , tag::diag_precision, std::streamsize
  , tag::diag_format, std::string
  , tag::ic, std::vector<
               tk::TaggedTuple< brigand::list<
                   tag::x,              std::vector< double >
                 , tag::y,              std::vector< double >
                 , tag::z,              std::vector< double >
                 , tag::ic_density,     double
                 , tag::ic_pressure,    double
                 , tag::ic_energy,      double
                 , tag::ic_temperature, double
                 , tag::ic_velocity,    std::vector< double >
               > >
             >
  , tag::ic_density, double
  , tag::ic_pressure, double
  , tag::ic_energy, double
  , tag::ic_temperature, double
  , tag::ic_velocity,  std::vector< double >
  , tag::bc_dir, std::vector< std::vector< int > >
  , tag::bc_sym, std::vector< int >
  , tag::bc_far, std::vector< int >
  , tag::bc_far_density, double
  , tag::bc_far_pressure, double
  , tag::bc_far_velocity, std::vector< double >
  , tag::bc_pre, std::vector< std::vector< int > >
  , tag::bc_pre_density, std::vector< double >
  , tag::bc_pre_pressure, std::vector< double >
  , tag::mat_spec_heat_ratio, double
  , tag::mat_spec_heat_const_vol, double
  , tag::mat_heat_conductivity, double
  , tag::href_t0, bool
  , tag::href_dt, bool
  , tag::href_dtfreq, uint64_t
  , tag::href_maxlevels, uint64_t
  , tag::href_refvar, std::vector< uint64_t >
  , tag::href_error, std::string
  , tag::href_init, std::vector< std::string >
>;

//! \brief InputDeck : Control< specialized to Inciter >, see Types.h,
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/Inciter/Types.h
class InputDeck : public tk::TaggedTuple< InputDeckMembers > {

  public:
    //! Parse control file
    void parse( const CmdLine& cmdline );

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< InputDeckMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i InputDeck object reference
    friend void operator|( PUP::er& p, InputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // inciter::
