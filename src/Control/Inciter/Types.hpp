// *****************************************************************************
/*!
  \file      src/Control/Inciter/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Types for Incitier's parsers
  \details   Types for Incitier's parsers. This file defines the components of
    the agged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during parsing (both command-line argument parsing and control file
    parsing).
*/
// *****************************************************************************

#pragma once

#include "Tags.hpp"
#include "TaggedTuple.hpp"
#include "Base/Types.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/Options/AMRInitial.hpp"
#include "Inciter/Options/AMRError.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

using namespace tao;

//! Adaptive-mesh refinement options
using amr = tk::TaggedTuple< brigand::list<
    tag::amr,     bool                            //!< AMR on/off
  , tag::t0ref,   bool                            //!< AMR before t<0 on/off
  , tag::dtref,   bool                            //!< AMR during t>0 on/off
  , tag::dtref_uniform, bool                      //!< Force dtref uniform-only
  , tag::dtfreq,  kw::amr_dtfreq::info::expect::type //!< Refinement frequency
  , tag::maxlevels, kw::amr_maxlevels::info::expect::type //!< Max refine levels
  , tag::init,    std::vector< AMRInitialType >   //!< List of initial AMR types
  , tag::refvar,  std::vector< std::string >      //!< List of refinement vars
  , tag::id,      std::vector< std::size_t >      //!< List of refvar indices
  , tag::error,   AMRErrorType                    //!< Error estimator for AMR
  , tag::tolref,  tk::real                        //!< Refine tolerance
  , tag::tolderef, tk::real                       //!< De-refine tolerance
  //! List of edges-node pairs
  , tag::edge,    std::vector< kw::amr_edgelist::info::expect::type >
  //! Refinement tagging edges with end-point coordinates lower than x coord
  , tag::xminus,  kw::amr_xminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than x coord
  , tag::xplus,  kw::amr_xplus::info::expect::type
  //! Refinement tagging edges with end-point coordinates lower than y coord
  , tag::yminus,  kw::amr_yminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than y coord
  , tag::yplus,  kw::amr_yplus::info::expect::type
  //! Refinement tagging edges with end-point coordinates lower than z coord
  , tag::zminus,  kw::amr_zminus::info::expect::type
  //! Refinement tagging edges with end-point coordinates higher than z coord
  , tag::zplus,  kw::amr_zplus::info::expect::type
> >;

//! Discretization parameters storage
using discretization = tk::TaggedTuple< brigand::list<
    tag::nstep,  kw::nstep::info::expect::type  //!< Number of time steps
  , tag::term,   kw::term::info::expect::type   //!< Time to terminate
  , tag::t0,     kw::t0::info::expect::type     //!< Starting time
  , tag::dt,     kw::dt::info::expect::type     //!< Size of time step
  , tag::cfl,    kw::cfl::info::expect::type    //!< CFL coefficient
  , tag::pelocal_reorder, bool                  //!< PE-locality reordering
  , tag::steady_state, bool                     //!< March to steady state
  , tag::residual, kw::residual::info::expect::type //!< Convergence residual
  , tag::rescomp, kw::rescomp::info::expect::type //!< Convergence residual comp
  , tag::partitioner, tk::ctr::PartitioningAlgorithmType //!< Mesh partitioner
> >;

//! ASCII output floating-point precision in digits
using precision = tk::TaggedTuple< brigand::list<
    //! Diagnostics output precision
    tag::diag, kw::precision::info::expect::type
    //! History output precision
  , tag::history, kw::precision::info::expect::type
    //! Integral output precision
  , tag::integral, kw::precision::info::expect::type
> >;

//! ASCII output floating-point format
using floatformat = tk::TaggedTuple< brigand::list<
    tag::diag,    tk::ctr::TxtFloatFormatType  //!< Diagnostics output format
  , tag::history, tk::ctr::TxtFloatFormatType  //!< History output format
  , tag::integral,tk::ctr::TxtFloatFormatType  //!< Integral output format
> >;

//! Output intervals in units of iteration count
using interval_iter = tk::TaggedTuple< brigand::list<
    //! TTY output interval
    tag::tty,     kw::ttyi::info::expect::type
    //! Field output interval
  , tag::field,   kw::interval_iter::info::expect::type
    //! History output interval
  , tag::history, kw::interval_iter::info::expect::type
    //! Integral output interval
  , tag::integral,kw::interval_iter::info::expect::type
    //! Diags output interval
  , tag::diag,    kw::interval_iter::info::expect::type
> >;

//! Output intervals in units of physics time
using interval_time = tk::TaggedTuple< brigand::list<
    //! Field output interval
    tag::field,   kw::interval_time::info::expect::type
    //! History output interval
  , tag::history, kw::interval_time::info::expect::type
    //! Integral output interval
  , tag::integral, kw::interval_time::info::expect::type
> >;

//! Output time ranges in units of physics time
using time_range = tk::TaggedTuple< brigand::list<
    //! \brief Field output configuration: outer vector: multiple ranges, inner
    //!        vector: mintime, maxtime, dt
    tag::field,   std::vector<
                    std::vector< kw::time_range::info::expect::type > >
    //! \brief History output configuration: outer vector: multiple ranges,
    //!        inner vector: mintime, maxtime, dt
  , tag::history, std::vector<
                    std::vector< kw::time_range::info::expect::type > >
    //! \brief Integral output configuration: outer vector: multiple ranges,
    //!        inner vector: mintime, maxtime, dt
  , tag::integral, std::vector<
                     std::vector< kw::time_range::info::expect::type > >
> >;

//! Output configuration parameters
using output_parameters = tk::TaggedTuple< brigand::list<
    //! Output intervals in units of iteration count
    tag::iter,  interval_iter
    //! Output intervals in units of physics time
  , tag::time,  interval_time
    //! Output time ranges in units of physics time
  , tag::range, time_range
> >;

//! History output parameters storage
using history = tk::TaggedTuple< brigand::list<
    tag::point,   std::vector< std::vector< kw::point::info::expect::type > >
  , tag::id,      std::vector< std::string >     //!< Point identifiers
> >;

//! Source input parameters storage
using Sources = tk::TaggedTuple< brigand::list<
    tag::point,   std::vector< std::vector< kw::point::info::expect::type > >
  , tag::id,      std::vector< std::string >     //!< Point identifiers
> >;

//! Surface I/O parameters storage
using surface = tk::TaggedTuple< brigand::list<
    //! List of side sets to save as field output
    tag::field,    std::vector< kw::sideset::info::expect::type >
    //! List of side sets at which to save integral output
  , tag::integral, std::vector< kw::sideset::info::expect::type >
> >;

//! IO parameters storage
using ios = tk::TaggedTuple< brigand::list<
    tag::nrestart,  int                             //!< Number of restarts
  , tag::control,   kw::control::info::expect::type //!< Control filename
  , tag::input,     kw::input::info::expect::type   //!< Input filename
  , tag::output,    kw::output::info::expect::type  //!< Output filename
    //! Refined output (output field data on a refined mesh)
  , tag::refined,   kw::refined::info::expect::type
  , tag::screen,    kw::screen::info::expect::type  //!< Screen output filename
  , tag::surface,   surface
    //! Diagnostics filename
  , tag::diag,      kw::diagnostics_cmd::info::expect::type
  , tag::particles, std::string                     //!< Particles filename
  , tag::restart,   kw::restart::info::expect::type //!< Restart dirname
> >;

//! Box, given by coordinates, specifying physics variables
using box = tk::TaggedTuple< brigand::list<
    tag::xmin,          kw::xmin::info::expect::type
  , tag::xmax,          kw::xmax::info::expect::type
  , tag::ymin,          kw::ymin::info::expect::type
  , tag::ymax,          kw::ymax::info::expect::type
  , tag::zmin,          kw::zmin::info::expect::type
  , tag::zmax,          kw::zmax::info::expect::type
  , tag::mass,          kw::mass::info::expect::type
  , tag::density,       kw::density::info::expect::type
  , tag::velocity,      std::vector< kw::velocity::info::expect::type >
  , tag::pressure,      kw::pressure::info::expect::type
  , tag::energy,        kw::energy::info::expect::type
  , tag::energy_content,kw::energy_content::info::expect::type
  , tag::temperature,   kw::temperature::info::expect::type
> >;

//! Initial condition configuration
using ic = tk::TaggedTuple< brigand::list<
    tag::density,       std::vector<
                          std::vector< kw::density::info::expect::type > >
  , tag::velocity,      std::vector<
                          std::vector< kw::velocity::info::expect::type > >
  , tag::pressure,      std::vector<
                          std::vector< kw::pressure::info::expect::type > >
  , tag::energy,        std::vector<
                          std::vector< kw::energy::info::expect::type > >
  , tag::temperature,   std::vector<
                          std::vector< kw::temperature::info::expect::type > >
  , tag::box,           std::vector< std::vector< box > >
> >;

//! Boundary conditions configuration (list of side sets for each eq system)
using bc = tk::TaggedTuple< brigand::list<
    tag::dirichlet,     std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::symmetry,      std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::farfield,      std::vector< std::vector<
                          kw::sideset::info::expect::type > >
  , tag::pressure,      std::vector< std::vector<
                          kw::sideset::info::expect::type > >
> >;

//! Mesh assignment and configuration
using mesh = tk::TaggedTuple< brigand::list<
    tag::id,          std::vector< std::size_t >
  , tag::filename,    std::vector< std::string >
  , tag::location,    std::vector<
                        std::vector< kw::location::info::expect::type > >
  , tag::orientation, std::vector<
                        std::vector< kw::orientation::info::expect::type > >
  , tag::reference,   std::vector< char >
> >;

//! Compressible flow equation parameters storage
using CompFlowPDEParameters = tk::TaggedTuple< brigand::list<
    tag::depvar,        std::vector< char >
  , tag::mesh,          mesh
  , tag::farfield_pressure, std::vector< kw::pressure::info::expect::type >
  , tag::farfield_density,  std::vector< kw::density::info::expect::type >
  , tag::farfield_velocity, std::vector< std::vector<
                              kw::velocity::info::expect::type > >
  , tag::pressure_pressure, std::vector< kw::pressure::info::expect::type >
  , tag::pressure_density,  std::vector< kw::density::info::expect::type >
  , tag::bc,            bc
  , tag::ic,            ic
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::alpha,         std::vector< kw::pde_alpha::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::beta,          std::vector< kw::pde_beta::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betax,         std::vector< kw::pde_betax::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betay,         std::vector< kw::pde_betay::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::betaz,         std::vector< kw::pde_betaz::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::r0,            std::vector< kw::pde_r0::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::ce,            std::vector< kw::pde_ce::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::kappa,         std::vector< kw::pde_kappa::info::expect::type >
    //! Parameter vector (for specific, e.g., verification problems)
  , tag::p0,            std::vector< kw::pde_p0::info::expect::type >
    //! Ratio of spec heats
  , tag::gamma,         std::vector<
                          std::vector< kw::mat_gamma::info::expect::type > >
    //! Dynamic viscosity
  , tag::mu,            std::vector<
                          std::vector< kw::mat_mu::info::expect::type > >
    //! Spec. heat at const vol.
  , tag::cv,            std::vector<
                          std::vector< kw::mat_cv::info::expect::type > >
    //! Heat conductivity
  , tag::k,             std::vector<
                          std::vector< kw::mat_k::info::expect::type > >
  , tag::diffusivity,   std::vector< std::vector<
                          kw::pde_diffusivity::info::expect::type > >
  , tag::lambda,        std::vector< std::vector<
                          kw::pde_lambda::info::expect::type > >
  , tag::u0,            std::vector< std::vector<
                          kw::pde_u0::info::expect::type > >
  , tag::source,        std::vector< std::vector<
                          kw::pde_source::info::expect::type > >
> >;

//! Parameters storage
using parameters = tk::TaggedTuple< brigand::list<
    tag::compflow,  CompFlowPDEParameters
> >;

//! PEGTL location/position type to use throughout all of Inciter's parsers
using Location = pegtl::position;

} // ctr::
} // inciter::
