// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/InputDeck.hpp
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
#ifndef InciterInputDeck_h
#define InciterInputDeck_h

#include <limits>
#include <iomanip>
#include <iostream>

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/set.hpp"

#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/Components.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {
namespace ctr {

//! Member data for tagged tuple
using InputDeckMembers = brigand::list<
    tag::cmd,           CmdLine
  , tag::title,         kw::title::info::expect::type
  , tag::amr,           amr
  , tag::discr,         discretization
  , tag::problem,       ProblemType
  , tag::prec,          precision
  , tag::flformat,      floatformat
  , tag::component,     ncomps
  , tag::sys,           std::map< tk::ctr::ncomp_t, tk::ctr::ncomp_t >
  , tag::output,        output_parameters
  , tag::param,         parameters
  , tag::diag,          diagnostics
  , tag::error,         std::vector< std::string >
  , tag::history,       history
>;

//! \brief InputDeck : Control< specialized to Inciter >, see Types.h,
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/Inciter/Types.h
class InputDeck : public tk::TaggedTuple< InputDeckMembers > {

  public:
    //! \brief Inciter input deck keywords
    //! \see tk::grm::use and its documentation
    using keywords = brigand::set< kw::title
                                 , kw::nstep
                                 , kw::term
                                 , kw::t0
                                 , kw::dt
                                 , kw::ttyi
                                 , kw::transport
                                 , kw::end
                                 , kw::shear_diff
                                 , kw::point_src
                                 , kw::slot_cyl
                                 , kw::problem
                                 , kw::field_output
                                 , kw::refined
                                 , kw::interval_iter
                                 , kw::interval_time
                                 , kw::time_range
                                 , kw::partitioning
                                 , kw::algorithm
                                 , kw::rcb
                                 , kw::rib
                                 , kw::hsfc
                                 , kw::phg
                                 , kw::inciter
                                 , kw::ncomp
                                 , kw::pde_diffusivity
                                 , kw::pde_lambda
                                 , kw::pde_u0
                                 , kw::pde_source
                                 , kw::bc_dirichlet
                                 , kw::sideset
                                 , kw::compflow
                                 , kw::ic
                                 , kw::box
                                 , kw::materialid
                                 , kw::mass
                                 , kw::density
                                 , kw::velocity
                                 , kw::position
                                 , kw::acceleration
                                 , kw::fntype
                                 , kw::fn
                                 , kw::move
                                 , kw::pressure
                                 , kw::energy
                                 , kw::energy_content
                                 , kw::temperature
                                 , kw::xmin
                                 , kw::xmax
                                 , kw::ymin
                                 , kw::ymax
                                 , kw::zmin
                                 , kw::zmax
                                 , kw::txt_float_format
                                 , kw::txt_float_default
                                 , kw::txt_float_fixed
                                 , kw::txt_float_scientific
                                 , kw::precision
                                 , kw::diagnostics
                                 , kw::history_output
                                 , kw::mesh
                                 , kw::filename
                                 , kw::location
                                 , kw::orientation
                                 , kw::reference
                                 , kw::material
                                 , kw::id
                                 , kw::eos
                                 , kw::stiffenedgas
                                 , kw::jwl
                                 , kw::mat_gamma
                                 , kw::mat_pstiff
                                 , kw::mat_mu
                                 , kw::mat_cv
                                 , kw::mat_k
                                 , kw::physics
                                 , kw::advection
                                 , kw::advdiff
                                 , kw::navierstokes
                                 , kw::euler
                                 , kw::user_defined
                                 , kw::vortical_flow
                                 , kw::pde_alpha
                                 , kw::pde_beta
                                 , kw::pde_p0
                                 , kw::cfl
                                 , kw::vortmult
                                 , kw::mj
                                 , kw::elem
                                 , kw::node
                                 , kw::depvar
                                 , kw::nonlinear_energy_growth
                                 , kw::pde_betax
                                 , kw::pde_betay
                                 , kw::pde_betaz
                                 , kw::pde_ce
                                 , kw::pde_kappa
                                 , kw::pde_r0
                                 , kw::rayleigh_taylor
                                 , kw::taylor_green
                                 , kw::filetype
                                 , kw::exodusii
                                 , kw::root
                                 , kw::error
                                 , kw::l2
                                 , kw::linf
                                 , kw::pelocal_reorder
                                 , kw::steady_state
                                 , kw::residual
                                 , kw::rescomp
                                 , kw::amr
                                 , kw::amr_t0ref
                                 , kw::amr_dtref
                                 , kw::amr_dtref_uniform
                                 , kw::amr_dtfreq
                                 , kw::amr_maxlevels
                                 , kw::amr_initial
                                 , kw::amr_uniform
                                 , kw::amr_uniform_derefine
                                 , kw::amr_initial_conditions
                                 , kw::amr_coords
                                 , kw::amr_error
                                 , kw::amr_jump
                                 , kw::amr_hessian
                                 , kw::amr_refvar
                                 , kw::amr_tolref
                                 , kw::amr_tolderef
                                 , kw::amr_edgelist
                                 , kw::amr_xminus
                                 , kw::amr_xplus
                                 , kw::amr_yminus
                                 , kw::amr_yplus
                                 , kw::amr_zminus
                                 , kw::amr_zplus
                                 , kw::bc_sym
                                 , kw::bc_farfield
                                 , kw::point
                                 , kw::radius
                                 , kw::rotated_sod
                                 , kw::sod
                                 , kw::sedov
                                 >;

    //! Set of tags to ignore when printing this InputDeck
    using ignore = CmdLine::ignore;

    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< tag::cmd >() = cl;
      // Default discretization parameters
      get< tag::discr, tag::nstep >() =
         std::numeric_limits< kw::nstep::info::expect::type >::max();
      get< tag::discr, tag::term >() =
         std::numeric_limits< kw::term::info::expect::type >::max();
      get< tag::discr, tag::t0 >() = 0.0;
      get< tag::discr, tag::dt >() = 0.0;
      get< tag::discr, tag::cfl >() = 0.0;
      get< tag::discr, tag::pelocal_reorder >() = false;
      get< tag::discr, tag::steady_state >() = false;
      get< tag::discr, tag::residual >() = 1.0e-8;
      get< tag::discr, tag::rescomp >() = 1;
      // Default AMR settings
      get< tag::amr, tag::amr >() = false;
      get< tag::amr, tag::t0ref >() = false;
      get< tag::amr, tag::dtref >() = false;
      get< tag::amr, tag::dtref_uniform >() = false;
      get< tag::amr, tag::dtfreq >() = 3;
      get< tag::amr, tag::maxlevels >() = 2;
      get< tag::amr, tag::error >() = AMRErrorType::JUMP;
      get< tag::amr, tag::tolref >() = 0.2;
      get< tag::amr, tag::tolderef >() = 0.05;
      auto rmax =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max() / 100;
      get< tag::amr, tag::xminus >() = rmax;
      get< tag::amr, tag::xplus >() = -rmax;
      get< tag::amr, tag::yminus >() = rmax;
      get< tag::amr, tag::yplus >() = -rmax;
      get< tag::amr, tag::zminus >() = rmax;
      get< tag::amr, tag::zplus >() = -rmax;
      // Default txt floating-point output precision in digits
      get< tag::prec, tag::diag >() = std::cout.precision();
      get< tag::prec, tag::history >() = std::cout.precision();
      // Default intervals
      get< tag::output, tag::iter, tag::tty >() = 1;
      get< tag::output, tag::iter, tag::diag >() = 1;
      get< tag::output, tag::iter, tag::field >() =
        std::numeric_limits< kw::interval_iter::info::expect::type >::max();
      get< tag::output, tag::iter, tag::history >() =
        std::numeric_limits< kw::interval_iter::info::expect::type >::max();
      // Initialize help: fill own keywords
      const auto& ctrinfoFill = tk::ctr::Info( get< tag::cmd, tag::ctrinfo >() );
      brigand::for_each< keywords >( ctrinfoFill );
    }

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

    //! Extract surface side set ids along which user wants to save solution
    //! \return Unique set of surface side set ids along which user wants to
    //!   save solution field variables
    //! \note This returns an ordered set so the order of the set ids are
    //!   always the same.
    std::set< int > outsets() const {
      std::set< int > ids;
      for (const auto& s : get< tag::cmd, tag::io, tag::surface >())
        ids.insert( s );
      return ids;
    }

    //! Extract list of mesh filenames (each assigned to a solver)
    std::vector< std::string > mesh() const {
      //using PDETypes = parameters::Keys;
      std::vector< std::string > meshes;
      //brigand::for_each< PDETypes >( Meshes( *this, meshes ) );
      return meshes;
    }

    //! Extract list of dependent variables (each configuring a solver)
    std::vector< char > depvar() const {
      //using PDETypes = parameters::Keys;
      std::vector< char > depvar;
      //brigand::for_each< PDETypes >( Depvar( *this, depvar ) );
      return depvar;
    }

  private:
    //! Function object to extract the mesh filenames assigned to solvers
    //! \details This is instantiated for all PDE types at compile time. It goes
    //!   through all configured solvers (equation system configuration blocks)
    //!   and builds a list of all mesh filenames associated to all solvers in
    //!   the input file.
    struct Meshes {
      const InputDeck& inputdeck;
      std::vector< std::string >& filenames;
      explicit Meshes( const InputDeck& i, std::vector< std::string >& f )
        : inputdeck(i), filenames(f) {}
      template< typename eq > void operator()( brigand::type_<eq> ) {
        const auto& eq_mesh_filename =
           inputdeck.get< tag::param, eq, tag::mesh, tag::filename >();
        for (const auto& f : eq_mesh_filename) filenames.push_back( f );
      }
    };

    //! Function object to extract the dependent variables assigned to solvers
    //! \details This is instantiated for all PDE types at compile time. It goes
    //!   through all configured solvers (equation system configuration blocks)
    //!   and builds a list of all dependent variables associated to all solvers
    //!   in the input file.
    struct Depvar {
      const InputDeck& inputdeck;
      std::vector< char >& depvar;
      explicit Depvar( const InputDeck& i, std::vector< char >& d ) :
        inputdeck(i), depvar(d) {}
      template< typename eq > void operator()( brigand::type_<eq> ) {
        const auto& eq_depvar = inputdeck.get< tag::param, eq, tag::depvar >();
        for (const auto& d : eq_depvar) depvar.push_back( d );
      }
    };
};

} // ctr::
} // inciter::

#endif // InciterInputDeck_h
