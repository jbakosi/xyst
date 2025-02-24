// *****************************************************************************
/*!
  \file      src/Control/InciterConfig.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
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
DEFTAG( rk );
DEFTAG( theta );
DEFTAG( t0 );
DEFTAG( dt );
DEFTAG( turkel );
DEFTAG( soundspeed );
DEFTAG( velinf );
DEFTAG( presure );
DEFTAG( tol );
DEFTAG( verbose );
DEFTAG( hydrostat );
DEFTAG( pc );
DEFTAG( mom_iter );
DEFTAG( mom_tol );
DEFTAG( mom_verbose );
DEFTAG( mom_pc );
DEFTAG( reorder );
DEFTAG( part );
DEFTAG( zoltan_params );
DEFTAG( solver );
DEFTAG( stab );
DEFTAG( stab2 );
DEFTAG( stab2coef );
DEFTAG( fct );
DEFTAG( fctdif );
DEFTAG( fctclip );
DEFTAG( fctsys );
DEFTAG( fctfreeze );
DEFTAG( deactivate );
DEFTAG( deatol );
DEFTAG( deadif );
DEFTAG( deafreq );
DEFTAG( deasys );
DEFTAG( deatime );
DEFTAG( flux );
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
DEFTAG( freezetime );
DEFTAG( fieldout );
DEFTAG( fieldout_ );
DEFTAG( histout );
DEFTAG( histout_ );
DEFTAG( integout );
DEFTAG( integout_ );
DEFTAG( iter );
DEFTAG( time );
DEFTAG( range );
DEFTAG( sidesets );
DEFTAG( points );
DEFTAG( integrals );
DEFTAG( precision );
DEFTAG( format );
DEFTAG( ic );
DEFTAG( ic_ );
DEFTAG( density );
DEFTAG( pressure );
DEFTAG( pressure_ );
DEFTAG( energy );
DEFTAG( temperature );
DEFTAG( velocity );
DEFTAG( boxes );
DEFTAG( box_x );
DEFTAG( box_y );
DEFTAG( box_z );
DEFTAG( box_density );
DEFTAG( box_pressure );
DEFTAG( box_energy );
DEFTAG( box_temperature );
DEFTAG( box_velocity );
DEFTAG( bc_dir );
DEFTAG( bc_dirval );
DEFTAG( bc_dir_ );
DEFTAG( bc_sym );
DEFTAG( bc_sym_ );
DEFTAG( bc_noslip );
DEFTAG( bc_noslip_ );
DEFTAG( bc_far );
DEFTAG( bc_far_ );
DEFTAG( bc_pre );
DEFTAG( bc_pre_ );
DEFTAG( mat_spec_heat_ratio );
DEFTAG( mat_spec_heat_const_vol );
DEFTAG( mat_spec_gas_const );
DEFTAG( mat_heat_conductivity );
DEFTAG( mat_dyn_viscosity );
DEFTAG( mat_dyn_diffusivity );
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
  , tag::input, std::vector< std::string >
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
  , tag::rk, uint64_t
  , tag::theta, double
  , tag::t0, double
  , tag::dt, double
  , tag::turkel, double
  , tag::soundspeed, double
  , tag::velinf, std::vector< double >
  , tag::pressure, tk::TaggedTuple< brigand::list<
                     tag::iter,      uint64_t
                   , tag::tol,       double
                   , tag::verbose,   uint64_t
                   , tag::hydrostat, uint64_t
                   , tag::pc,        std::string
                   , tag::bc_dir,    std::vector< std::vector< int > >
                   , tag::bc_dirval, std::vector< std::vector< double > >
                   , tag::bc_sym,    std::vector< int >
                   > >
  , tag::pressure_, std::vector<
                      tk::TaggedTuple< brigand::list<
                        tag::iter,      uint64_t
                      , tag::tol,       double
                      , tag::verbose,   uint64_t
                      , tag::hydrostat, uint64_t
                      , tag::pc,        std::string
                      , tag::bc_dir,    std::vector< std::vector< int > >
                      , tag::bc_dirval, std::vector< std::vector< double > >
                      , tag::bc_sym,    std::vector< int >
                      > >
                    >
  , tag::mom_iter, uint64_t
  , tag::mom_tol, double
  , tag::mom_verbose, uint64_t
  , tag::mom_pc, std::string
  , tag::reorder, bool
  , tag::part, std::string
  , tag::zoltan_params, std::vector< std::string >
  , tag::solver, std::string
  , tag::stab, bool
  , tag::stab2, bool
  , tag::stab2coef, double
  , tag::fct, bool
  , tag::fctdif, double
  , tag::fctclip, bool
  , tag::fctsys, std::vector< uint64_t >
  , tag::fctfreeze, double
  , tag::deactivate, bool
  , tag::deatol, double
  , tag::deadif, double
  , tag::deafreq, uint64_t
  , tag::deasys, std::vector< uint64_t >
  , tag::deatime, double
  , tag::flux, std::string
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
  , tag::freezeflow, double
  , tag::freezetime, double
  , tag::fieldout, tk::TaggedTuple< brigand::list<
                     tag::sidesets, std::vector< int >
                   , tag::iter,     uint64_t
                   , tag::time,     double
                   , tag::range,    std::vector< std::vector< double > >
                   > >
  , tag::fieldout_, std::vector<
                      tk::TaggedTuple< brigand::list<
                        tag::sidesets, std::vector< int >
                      , tag::iter,     uint64_t
                      , tag::time,     double
                      , tag::range,    std::vector< std::vector< double > >
                      > >
                    >
  , tag::histout, tk::TaggedTuple< brigand::list<
                    tag::points,    std::vector< std::vector< double > >
                  , tag::iter,      uint64_t
                  , tag::time,      double
                  , tag::range,     std::vector< std::vector< double > >
                  , tag::precision, std::streamsize
                  , tag::format,    std::string
                  > >
  , tag::histout_, std::vector<
                     tk::TaggedTuple< brigand::list<
                       tag::points,    std::vector< std::vector< double > >
                     , tag::iter,      uint64_t
                     , tag::time,      double
                     , tag::range,     std::vector< std::vector< double > >
                     , tag::precision, std::streamsize
                     , tag::format,    std::string
                     > >
                   >
  , tag::integout, tk::TaggedTuple< brigand::list<
                    tag::sidesets, std::vector< int >
                  , tag::integrals, std::vector< std::string >
                  , tag::iter, uint64_t
                  , tag::time, double
                  , tag::range, std::vector< std::vector< double > >
                  , tag::precision, std::streamsize
                  , tag::format, std::string
                  > >
  , tag::integout_, std::vector<
                      tk::TaggedTuple< brigand::list<
                        tag::sidesets, std::vector< int >
                      , tag::integrals, std::vector< std::string >
                      , tag::iter, uint64_t
                      , tag::time, double
                      , tag::range, std::vector< std::vector< double > >
                      , tag::precision, std::streamsize
                      , tag::format, std::string
                      > >
                   >
  , tag::ic, tk::TaggedTuple< brigand::list<
               tag::density,     double
             , tag::pressure,    double
             , tag::energy,      double
             , tag::temperature, double
             , tag::velocity,    std::vector< double >
             , tag::boxes, std::vector<
                 tk::TaggedTuple< brigand::list<
                     tag::box_x,           std::vector< double >
                   , tag::box_y,           std::vector< double >
                   , tag::box_z,           std::vector< double >
                   , tag::box_density,     double
                   , tag::box_pressure,    double
                   , tag::box_energy,      double
                   , tag::box_temperature, double
                   , tag::box_velocity,    std::vector< double >
                 > >
               >
             > >
  , tag::ic_, std::vector<
                tk::TaggedTuple< brigand::list<
                  tag::density,     double
                , tag::pressure,    double
                , tag::energy,      double
                , tag::temperature, double
                , tag::velocity,    std::vector< double >
                , tag::boxes, std::vector<
                    tk::TaggedTuple< brigand::list<
                        tag::box_x,           std::vector< double >
                      , tag::box_y,           std::vector< double >
                      , tag::box_z,           std::vector< double >
                      , tag::box_density,     double
                      , tag::box_pressure,    double
                      , tag::box_energy,      double
                      , tag::box_temperature, double
                      , tag::box_velocity,    std::vector< double >
                    > >
                  >
                > >
              >
  , tag::bc_dir, std::vector< std::vector< int > >
  , tag::bc_dirval, std::vector< std::vector< double > >
  , tag::bc_dir_, std::vector< std::vector< std::vector< int > > >
  , tag::bc_sym, std::vector< int >
  , tag::bc_sym_, std::vector< std::vector< int > >
  , tag::bc_noslip, std::vector< int >
  , tag::bc_noslip_, std::vector< std::vector< int > >
  , tag::bc_far, tk::TaggedTuple< brigand::list<
                     tag::sidesets, std::vector< int >
                   , tag::density,  double
                   , tag::pressure, double
                   , tag::velocity, std::vector< double >
                 > >
  , tag::bc_far_, std::vector<
                    tk::TaggedTuple< brigand::list<
                      tag::sidesets, std::vector< int >
                    , tag::density,  double
                    , tag::pressure, double
                    , tag::velocity, std::vector< double >
                    > >
                  >
  , tag::bc_pre, tk::TaggedTuple< brigand::list<
                   tag::sidesets, std::vector< std::vector< int > >
                 , tag::density,  std::vector< double >
                 , tag::pressure, std::vector< double >
                 > >
  , tag::bc_pre_, std::vector<
                    tk::TaggedTuple< brigand::list<
                      tag::sidesets, std::vector< std::vector< int > >
                    , tag::density,  std::vector< double >
                    , tag::pressure, std::vector< double >
                    > >
                  >
  , tag::mat_spec_heat_ratio, double
  , tag::mat_spec_heat_const_vol, double
  , tag::mat_spec_gas_const, double
  , tag::mat_heat_conductivity, double
  , tag::mat_dyn_viscosity, double
  , tag::mat_dyn_diffusivity, double
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
