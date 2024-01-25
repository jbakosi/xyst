// *****************************************************************************
/*!
  \file      src/Control/InciterConfig.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter confguration (parsed from cmdline and control file)
*/
// *****************************************************************************
#pragma once

#include <getopt.h>

#include "TaggedTuple.hpp"

namespace tag {
DEFTAG( commit );
DEFTAG( input );
DEFTAG( control );
DEFTAG( output );
DEFTAG( diag );
DEFTAG( diag_iter );
DEFTAG( diag_precision );
DEFTAG( diag_format );
DEFTAG( checkpoint );
DEFTAG( quiescence );
DEFTAG( virt );
DEFTAG( nonblocking );
DEFTAG( benchmark );
DEFTAG( feedback );
DEFTAG( lbfreq );
DEFTAG( lbtime );
DEFTAG( rsfreq );
DEFTAG( nstep );
DEFTAG( ttyi );
DEFTAG( term );
DEFTAG( cfl );
DEFTAG( t0 );
DEFTAG( dt );
DEFTAG( reorder );
DEFTAG( part );
DEFTAG( zoltan_params );
DEFTAG( solver );
DEFTAG( stab2 );
DEFTAG( stab2coef );
DEFTAG( stab4 );
DEFTAG( fct );
DEFTAG( fctdif );
DEFTAG( fctclip );
DEFTAG( fctsys );
DEFTAG( deactivate );
DEFTAG( deatol );
DEFTAG( deadif );
DEFTAG( deafreq );
DEFTAG( deasys );
DEFTAG( deatime );
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
DEFTAG( freezeflow );
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
//! Inciter control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using ConfigMembers = brigand::list<
    tag::commit, std::string
  , tag::input, std::string
  , tag::control, std::string
  , tag::output, std::string
  , tag::diag, std::string
  , tag::diag_iter, uint64_t
  , tag::diag_precision, std::streamsize
  , tag::diag_format, std::string
  , tag::checkpoint, std::string
  , tag::quiescence, bool
  , tag::virt, double
  , tag::nonblocking, bool
  , tag::benchmark, bool
  , tag::feedback, bool
  , tag::lbfreq, uint64_t
  , tag::lbtime, double
  , tag::rsfreq, uint64_t
  , tag::nstep, uint64_t
  , tag::ttyi, uint64_t
  , tag::term, double
  , tag::cfl, double
  , tag::t0, double
  , tag::dt, double
  , tag::reorder, bool
  , tag::part, std::string
  , tag::zoltan_params, std::vector< std::string >
  , tag::solver, std::string
  , tag::stab2, bool
  , tag::stab2coef, double
  , tag::stab4, bool
  , tag::fct, bool
  , tag::fctdif, double
  , tag::fctclip, bool
  , tag::fctsys, std::vector< uint64_t >
  , tag::deactivate, bool
  , tag::deatol, double
  , tag::deadif, double
  , tag::deafreq, uint64_t
  , tag::deasys, std::vector< uint64_t >
  , tag::deatime, double
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
                      , tag::freezeflow, double
                      > >
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

//! Config is a TaggedTuple specialized to Inciter
class Config : public tk::TaggedTuple< ConfigMembers > {

  public:
    //! Contructor: parse inciter command line
    void cmdline( int argc, char** argv );

    //! Parse control file
    void control();

    /** @name Pack/Unpack: Serialize Config object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< ConfigMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c Config object reference
    friend void operator|( PUP::er& p, Config& c ) { c.pup(p); }
    //@}

  private:
    //! Echo help on command line arguments
    void help( char** argv );
};

} // ctr::
} // inciter::
