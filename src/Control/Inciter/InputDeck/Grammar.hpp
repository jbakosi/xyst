// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's input deck grammar definition
  \details   Inciter's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterInputDeckGrammar_h
#define InciterInputDeckGrammar_h

#include <limits>
#include <cmath>

#include "CommonGrammar.hpp"
#include "Keywords.hpp"
#include "ContainerUtil.hpp"
#include "Centering.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  //! \brief Specialization of tk::grm::use for Inciter's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::InputDeck::keywords >;

  // Inciter's InputDeck state

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  static tk::TaggedTuple< brigand::list<
             tag::transport,   std::size_t
           , tag::compflow,    std::size_t
         > > neq;

  //! \brief Parser-lifetime storage for point names
  //! \details Used to track the point names registered so that parsing new ones
  //!    can be required to be unique.
  static std::set< std::string > pointnames;

} // ::deck
} // ::inciter

namespace tk {
namespace grm {

  using namespace tao;

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // Inciter's InputDeck actions

  //! Rule used to trigger action
  template< class eq > struct register_inciter_eq : pegtl::success {};
  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  template< class eq >
  struct action< register_inciter_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {
      using inciter::deck::neq;
      ++neq.get< eq >();
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_mesh : pegtl::success {};
  //! \brief Check mesh ... end block for correctness
  template< class eq >
  struct action< check_mesh< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      auto& mesh = stack.template get< tag::param, eq, tag::mesh >();
      auto& mesh_ref = mesh.template get< tag::reference >();
      // if no mesh reference given by user
      if (mesh_ref.empty() || mesh_ref.size() != neq.get< eq >()) {
        // put in '-', meaning no reference
        mesh_ref.push_back('-');
        auto& location = mesh.template get< tag::location >();
        // if no location, put in the origin
        if (location.size() != neq.get< eq >())
          location.push_back( { 0.0, 0.0, 0.0 } );
        else    // reference was not given, but location was, error out
          Message< Stack, ERROR, MsgKey::LOC_NOMESHREF >( stack, in );
        auto& orientation = mesh.template get< tag::orientation >();
        if (orientation.size() != neq.get< eq >())
          orientation.push_back( { 0.0, 0.0, 0.0 } );
        else    // reference was not given, but orientation was, error out
          Message< Stack, ERROR, MsgKey::ORI_NOMESHREF >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_compflow : pegtl::success {};
  //! \brief Set defaults and do error checking on the compressible flow
  //!   equation block
  //! \details This is error checking that only the compressible flow equation
  //!   block must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a compflow...end block has
  //!   just finished.
  template< class eq >
  struct action< check_compflow< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using tag::param;
      using ProblemType = inciter::ctr::ProblemType;

      // Error out if no dependent variable has been selected
      auto& depvar = stack.template get< param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );

      // Setup number of scalar components based on problem specified
      auto problem = stack.template get< tag::problem >();
      std::size_t nc = 5;
      if (problem == ProblemType::SLOT_CYL ||
          problem == ProblemType::POINT_SRC) {
         ++nc;
      }
      stack.template get< tag::component, eq >().push_back( nc );

      // Error check Dirichlet boundary condition block for all transport eq
      // configurations
      const auto& bc =
        stack.template get< param, eq, tag::bc, tag::dirichlet >();
      for (const auto& s : bc)
        if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_transport : pegtl::success {};
  //! \brief Set defaults and do error checking on the transport equation block
  //! \details This is error checking that only the transport equation block
  //!   must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a transport...end block has
  //!   just finished.
  template< class eq >
  struct action< check_transport< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using tag::param;

      // Error out if no dependent variable has been selected
      auto& depvar = stack.template get< param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );

      // If no number of components has been selected, default to 1
      auto& ncomp = stack.template get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        ncomp.push_back( 1 );

      // Error check Dirichlet boundary condition block for all transport eq
      // configurations
      const auto& bc =
        stack.template get< param, eq, tag::bc, tag::dirichlet >();
      for (const auto& s : bc)
        if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class Option, typename...tags >
  struct store_inciter_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for inciter.
  template< class Option, typename... tags >
  struct action< store_inciter_option< Option, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      store_option< Stack, inciter::deck::use, Option, inciter::ctr::InputDeck,
                    Input, tags... >
                  ( stack, in, inciter::g_inputdeck_defaults );
    }
  };

  //! Function object to count the number of meshes assigned to solvers
  //! \details This is instantiated for all PDE types at compile time. It goes
  //!   through all configured solvers (equation system configuration blocks)
  //!   and counts the number of mesh filenames configured.
  template< typename Stack >
  struct count_meshes {
    const Stack& stack;
    std::size_t& count;
    explicit
    count_meshes( const Stack& s, std::size_t& c ) : stack(s), count(c) {}
    template< typename eq > void operator()( brigand::type_<eq> ) {
      count +=
        stack.template get< tag::param, eq, tag::mesh, tag::filename >().size();
    }
  };

  //! Function object to assign mesh ids to solvers
  //! \details This is instantiated for all PDE types at compile time. It goes
  //!   through all configured solvers (equation system configuration blocks)
  //!   and assigns a new mesh id to all solvers configured in the input file.
  template< typename Stack >
  struct assign_meshid {
    Stack& stack;
    std::size_t& meshid;
    explicit assign_meshid( Stack& s, std::size_t& m ) : stack(s), meshid(m) {}
    template< typename eq > void operator()( brigand::type_<eq> ) {
      const auto& eq_mesh_filename =
        stack.template get< tag::param, eq, tag::mesh, tag::filename >();
      auto& id = stack.template get< tag::param, eq, tag::mesh, tag::id >();
      for (std::size_t i=0; i<eq_mesh_filename.size(); ++i)
        id.push_back( meshid++ );
    }
  };

  //! Function object to do error checking on output time ranges
  template< typename Stack, typename Input >
  struct range_errchk {
    Stack& stack;
    const Input& input;
    explicit range_errchk( Stack& s, const Input& in ) : stack(s), input(in) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      for (const auto& r : stack.template get< tag::output, tag::range, U >())
        if ( r.size() != 3 or r[0] > r[1] or r[2] < 0.0 or r[2] > r[1]-r[0] )
          Message< Stack, ERROR, MsgKey::BADRANGE >( stack, input );
    }
  };

  //! Rule used to trigger action
  struct check_inciter : pegtl::success {};
  //! \brief Do error checking on the inciter block
  template<> struct action< check_inciter > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using inciter::g_inputdeck_defaults;

      // Error out if no dt policy has been selected
      const auto& dt = stack.template get< tag::discr, tag::dt >();
      const auto& cfl = stack.template get< tag::discr, tag::cfl >();
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) <
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) <
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, ERROR, MsgKey::NODT >( stack, in );

      // If both dt and cfl are given, warn that dt wins over cfl
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, WARNING, MsgKey::MULDT >( stack, in );

      // Do error checking on time history points
      const auto& hist = stack.template get< tag::history, tag::point >();
      if (std::any_of( begin(hist), end(hist),
           [](const auto& p){ return p.size() != 3; } ) )
      {
        Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
      }

      // Do error checking on residual eq system component index
      const auto rc = stack.template get< tag::discr, tag::rescomp >();
      const auto& ncomps = stack.template get< tag::component >();
      if (rc < 1 || rc > ncomps.nprop())
        Message< Stack, ERROR, MsgKey::LARGECOMP >( stack, in );

      // Do error checking on output time range configuration parameters: they
      // all must be a 3 reals: mintime, maxtime, and dt with maxtime >
      // mintime, and dt<maxtime-mintime.
      brigand::for_each< inciter::ctr::time_range::Keys >
                       ( range_errchk< Stack, Input >( stack, in ) );

      // Do error checking on time history point names (this is a programmer
      // error if triggers, hence assert)
      Assert(
        (stack.template get< tag::history, tag::id >().size() == hist.size()),
        "Number of history points and ids must equal" );

      using PDETypes = inciter::ctr::parameters::Keys;

      // If at least a mesh filename is assigned to a solver, all solvers must
      // have a mesh filename assigned
      std::size_t nmesh = 0;
      brigand::for_each< PDETypes >( count_meshes< Stack >( stack, nmesh ) );
      if (nmesh > 0 && nmesh != depvars.size())
        Message< Stack, ERROR, MsgKey::MULTIMESH >( stack, in );

      // Now that the inciter ... end block is finished, assign mesh ids to
      // solvers configured
      std::size_t meshid = 0;
      brigand::for_each< PDETypes >( assign_meshid< Stack >( stack, meshid ) );
      Assert( meshid == nmesh, "Not all meshes configured have mesh ids" );
    }
  };

  //! Rule used to trigger action
  template< typename Feature >
  struct enable : pegtl::success {};
  //! Enable adaptive mesh refinement (AMR)
  template< typename Feature >
  struct action< enable< Feature > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template get< Feature, Feature >() = true;
    }
  };

  //! Rule used to trigger action
  struct compute_refvar_idx : pegtl::success {};
  //! Compute indices of refinement variables
  //! \details This functor computes the indices in the unknown vector for all
  //!   refinement variables in the system of systems of dependent variables
  //!   after the refvar...end block has been parsed in the amr...end block.
  //!   After basic error checking, the vector at stack.get<tag::amr,tag::id>()
  //!   is filled.
  template<>
  struct action< compute_refvar_idx > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // reference variables just parsed by refvar...end block
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      // get ncomponents object from this input deck
      const auto& ncomps = stack.template get< tag::component >();
      // compute offset map associating offsets to dependent variables
      auto offsetmap = ncomps.offsetmap( stack );
      // compute number of components associated to dependent variabels
      auto ncompmap = ncomps.ncompmap( stack );
      // reference variable index vector to fill
      auto& refidx = stack.template get< tag::amr, tag::id >();
      // Compute indices for all refvars
      for (const auto& v : refvar) {    // for all reference variables parsed
        // depvar is the first char of a refvar
        auto depvar = v[0];
        // the field ID is optional and is the rest of the depvar string
        std::size_t f = (v.size()>1 ? std::stoul(v.substr(1)) : 1) - 1;
        // field ID must be less than or equal to the number of scalar
        // components configured for the eq system for this dependent variable
        if (f >= tk::cref_find( ncompmap, depvar ))
          Message< Stack, ERROR, MsgKey::NOSUCHCOMPONENT >( stack, in );
        // get offset for depvar
        auto eqsys_offset = tk::cref_find( offsetmap, depvar );
        // the index is the eq offset + field ID
        auto idx = eqsys_offset + f;
        // save refvar index in system of all systems
        refidx.push_back( idx );
      }
    }
  };

  //! Rule used to trigger action
  struct check_amr_errors : pegtl::success {};
  //! Do error checking for the amr...end block
  //! \details This is error checking that only the amr...end block
  //!   must satisfy. Besides error checking this can also set defaults
  //!   as this block is called when parsing of a amr...end block has
  //!   just finished.
  template<>
  struct action< check_amr_errors > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Error out if refvar size does not equal refidx size (programmer error)
      Assert( (stack.template get< tag::amr, tag::refvar >().size() ==
               stack.template get< tag::amr, tag::id >().size()),
              "The size of refvar and refidx vectors must equal" );
      const auto& initref = stack.template get< tag::amr, tag::init >();
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      const auto& edgelist = stack.template get< tag::amr, tag::edge >();
      // Error out if initref edge list is not divisible by 2 (user error)
      if (edgelist.size() % 2 == 1)
        Message< Stack, ERROR, MsgKey::T0REFODD >( stack, in );
      // Warn if initial AMR will be a no-op
      if ( stack.template get< tag::amr, tag::t0ref >() && initref.empty() )
        Message< Stack, WARNING, MsgKey::T0REFNOOP >( stack, in );
      // Error out if timestepping AMR will be a no-op (user error)
      if ( stack.template get< tag::amr, tag::dtref >() && refvar.empty() )
        Message< Stack, ERROR, MsgKey::DTREFNOOP >( stack, in );
      // Error out if mesh refinement frequency is zero (programmer error)
      Assert( (stack.template get< tag::amr, tag::dtfreq >() > 0),
              "Mesh refinement frequency must be positive" );
    }
  };

  //! Rule used to trigger action
  struct check_pref_errors : pegtl::success {};
  //! Do error checking for the pref...end block
  template<>
  struct action< check_pref_errors > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto& tolref = stack.template get< tag::pref, tag::tolref >();
      if (tolref < 0.0 || tolref > 1.0)
        Message< Stack, ERROR, MsgKey::PREFTOL >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct match_pointname : pegtl::success {};
  //! \brief Match PDF name to the registered ones
  //! \details This is used to check the set of PDF names dependent previously
  //!    registered to make sure all are unique.
  template<>
  struct action< match_pointname > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::pointnames;
      // find matched name in set of registered ones
      if (pointnames.find( in.string() ) == pointnames.end()) {
        pointnames.insert( in.string() );
        stack.template get< tag::history, tag::id >().push_back( in.string() );
      }
      else  // error out if name matched var is already registered
        Message< Stack, ERROR, MsgKey::POINTEXISTS >( stack, in );
    }
  };

} // ::grm
} // ::tk

namespace inciter {

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  using namespace tao;

  // Inciter's InputDeck grammar

  //! Error checks after an equation...end block has been parsed
  template< class eq, template< class > class eqchecker >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           tk::grm::register_inciter_eq< eq >,
           // check mesh ... end block
           tk::grm::check_mesh< eq >,
           // do error checking on this block
           eqchecker< eq > > {};

  //! Discretization parameters
  struct discretization :
         pegtl::sor<
           tk::grm::discrparam< use, kw::nstep, tag::nstep >,
           tk::grm::discrparam< use, kw::term, tag::term >,
           tk::grm::discrparam< use, kw::t0, tag::t0 >,
           tk::grm::discrparam< use, kw::dt, tag::dt >,
           tk::grm::discrparam< use, kw::cfl, tag::cfl >,
           tk::grm::discrparam< use, kw::residual, tag::residual >,
           tk::grm::discrparam< use, kw::rescomp, tag::rescomp >,
           tk::grm::process< use< kw::pelocal_reorder >,
                             tk::grm::Store< tag::discr, tag::pelocal_reorder >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::steady_state >,
                             tk::grm::Store< tag::discr, tag::steady_state >,
                             pegtl::alpha >,
           tk::grm::interval_iter< use< kw::ttyi >,
                                   tag::output, tag::iter, tag::tty >
         > {};

  //! PDE parameter vector
  template< class keyword, class eq, class param, class... xparams >
  struct pde_parameter_vector :
         tk::grm::parameter_vector< use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq, param, xparams... > {};

   //! Match box parameter
  template< class eq, typename keyword, typename target >
  struct box_parameter :
           pegtl::if_must<
             tk::grm::readkw< typename use<keyword>::pegtl_string >,
             tk::grm::scan<
               pegtl::sor< tk::grm::number,
                           tk::grm::msg< tk::grm::ERROR,
                                         tk::grm::MsgKey::MISSING > >,
               tk::grm::Back_back_store< target,
                 tag::param, eq, tag::ic, tag::box > > > {};

   //! Match box parameter and store deep
  template< class eq, typename keyword, typename target, typename subtarget >
  struct box_deep_parameter :
           pegtl::if_must<
             tk::grm::readkw< typename use<keyword>::pegtl_string >,
             tk::grm::scan<
               pegtl::sor< tk::grm::number,
                           tk::grm::msg< tk::grm::ERROR,
                                         tk::grm::MsgKey::MISSING > >,
               tk::grm::Back_back_deep_store< target, subtarget,
                 tag::param, eq, tag::ic, tag::box > > > {};

   //! Match box parameter vector
  template< class eq, typename keyword, typename target >
  struct box_vector :
         tk::grm::vector< use< keyword >,
                          tk::grm::Back_back_store_back< target,
                            tag::param, eq, tag::ic, tag::box >,
                          use< kw::end > > {};

   //! Match box parameter vector and store deep
  template< class eq, typename keyword, typename target, typename subtarget >
  struct box_deep_vector :
         tk::grm::vector< use< keyword >,
                          tk::grm::Back_back_deep_store_back< target, subtarget,
                            tag::param, eq, tag::ic, tag::box >,
                          use< kw::end > > {};


   //! Match box option
  template< class eq, typename Option, typename keyword, typename target,
            typename subtarget >
  struct box_option :
         tk::grm::process<
           use< keyword >,
           tk::grm::back_back_deep_store_option< target, subtarget, use,
             Option, tag::param, eq, tag::ic, tag::box >,
           pegtl::alpha > {};

  //! put in PDE parameter for equation matching keyword
  template< typename eq, typename keyword, typename param,
            class kw_type = tk::grm::number >
  struct parameter :
         tk::grm::process< use< keyword >,
                           tk::grm::Store_back< tag::param, eq, param >,
                           kw_type > {};

  //! put in PDE bool parameter for equation matching keyword into vector< int >
  template< typename eq, typename keyword, typename p >
  struct parameter_bool :
         tk::grm::process< use< keyword >,
                           tk::grm::Store_back_bool< tag::param, eq, p >,
                           pegtl::alpha > {};

  //! Dirichlet boundary conditions block
  struct bc_dirichlet :
         pegtl::if_must<
           tk::grm::readkw< use< kw::bc_dirichlet >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector<
               use,
               use< kw::sideset >,
               tk::grm::Store_back_back,
               tk::grm::start_vector,
               tk::grm::check_vector,
               tag::compflow, tag::bc, tag::dirichlet > > > {};

  //! Symmetry boundary conditions block
  struct bc_sym :
         pegtl::if_must<
           tk::grm::readkw< typename use< kw::bc_sym >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector<
               use,
               use< kw::sideset >,
               tk::grm::Store_back_back,
               tk::grm::start_vector,
               tk::grm::check_vector,
               tag::compflow, tag::bc, tag::symmetry > > > {};

  //! Farfield boundary conditions block
  struct bc_farfield :
         pegtl::if_must<
           tk::grm::readkw< typename use< kw::bc_farfield >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             parameter< tag::compflow, kw::pressure, tag::farfield_pressure >,
             parameter< tag::compflow, kw::density, tag::farfield_density >,
             pde_parameter_vector< kw::velocity, tag::compflow,
                                   tag::farfield_velocity >,
             tk::grm::parameter_vector<
               use,
               use< kw::sideset >,
               tk::grm::Store_back_back,
               tk::grm::start_vector,
               tk::grm::check_vector,
               tag::compflow, tag::bc, tag::farfield > > > {};

  //! Pressure boundary conditions block
  struct bc_pressure :
         pegtl::if_must<
           tk::grm::readkw< typename use< kw::bc_pressure >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             parameter< tag::compflow, kw::pressure, tag::pressure_pressure >,
             parameter< tag::compflow, kw::density, tag::pressure_density >,
             tk::grm::parameter_vector<
               use,
               use< kw::sideset >,
               tk::grm::Store_back_back,
               tk::grm::start_vector,
               tk::grm::check_vector,
               tag::compflow, tag::bc, tag::pressure > > > {};

  //! edgelist ... end block
  struct edgelist :
         tk::grm::vector< use< kw::amr_edgelist >,
                          tk::grm::Store_back< tag::amr, tag::edge >,
                          use< kw::end >,
                          tk::grm::check_vector< tag::amr, tag::edge > > {};

  //! xminus configuring coordinate-based edge tagging for mesh refinement
  template< typename keyword, typename Tag >
  struct half_world :
         tk::grm::control< use< keyword >, pegtl::digit, tk::grm::Store,
                           tag::amr, Tag > {};

  //! coords ... end block
  struct coords :
           pegtl::if_must<
             tk::grm::readkw< use< kw::amr_coords >::pegtl_string >,
             tk::grm::block< use< kw::end >,
                             half_world< kw::amr_xminus, tag::xminus >,
                             half_world< kw::amr_xplus, tag::xplus >,
                             half_world< kw::amr_yminus, tag::yminus >,
                             half_world< kw::amr_yplus, tag::yplus >,
                             half_world< kw::amr_zminus, tag::zminus >,
                             half_world< kw::amr_zplus, tag::zplus > > > {};

  //! initial conditions box block
  template< class eq >
  struct box :
         pegtl::if_must<
           tk::grm::readkw< use< kw::box >::pegtl_string >,
           tk::grm::start_vector_back< tag::param, eq, tag::ic, tag::box >,
           tk::grm::block< use< kw::end >
             , box_parameter< eq, kw::xmin, tag::xmin >
             , box_parameter< eq, kw::xmax, tag::xmax >
             , box_parameter< eq, kw::ymin, tag::ymin >
             , box_parameter< eq, kw::ymax, tag::ymax >
             , box_parameter< eq, kw::zmin, tag::zmin >
             , box_parameter< eq, kw::zmax, tag::zmax >
             , box_parameter< eq, kw::density, tag::density >
             , box_parameter< eq, kw::pressure, tag::pressure >
             , box_parameter< eq, kw::temperature, tag::temperature >
             , box_parameter< eq, kw::energy_content, tag::energy_content >
             , box_parameter< eq, kw::energy, tag::energy >
             , box_parameter< eq, kw::mass, tag::mass >
             , box_vector< eq, kw::velocity, tag::velocity >
             > > {};

  //! initial conditions block for compressible flow
  template< class eq >
  struct ic :
         pegtl::if_must<
           tk::grm::readkw< use< kw::ic >::pegtl_string >,
           tk::grm::block< use< kw::end >,
             pegtl::sor<
               pde_parameter_vector< kw::density, eq,
                                     tag::ic, tag::density >,
               pde_parameter_vector< kw::velocity, eq,
                                     tag::ic, tag::velocity >,
               pde_parameter_vector< kw::pressure, eq,
                                     tag::ic, tag::pressure >,
               pde_parameter_vector< kw::temperature, eq,
                                     tag::ic, tag::temperature >,
               pde_parameter_vector< kw::energy, eq,
                                     tag::ic, tag::energy > >,
               pegtl::seq< box< eq > > > > {};

  //! put in material property for equation matching keyword
  template< typename eq, typename keyword, typename property >
  struct material_property :
         pde_parameter_vector< keyword, eq, property > {};

  //! Material properties block for compressible flow
  template< class eq >
  struct material_properties :
         pegtl::seq<
          pegtl::if_must<
            tk::grm::readkw< use< kw::material >::pegtl_string >,
            tk::grm::block< use< kw::end >,
             material_property< eq, kw::mat_gamma, tag::gamma >,
             material_property< eq, kw::mat_mu, tag::mu >,
             material_property< eq, kw::mat_cv, tag::cv >,
             material_property< eq, kw::mat_k, tag::k >
            > > > {};

  //! Mesh ... end block
  template< class eq >
  struct mesh :
         pegtl::if_must<
           tk::grm::readkw< use< kw::mesh >::pegtl_string >,
           tk::grm::block< use< kw::end >,
             tk::grm::filename< use, tag::param, eq, tag::mesh, tag::filename >
           , pde_parameter_vector< kw::location, eq, tag::mesh, tag::location >
           , pde_parameter_vector< kw::orientation, eq,
                                   tag::mesh, tag::orientation >
           , tk::grm::process<
               use< kw::reference >,
               tk::grm::Store_back< tag::param, eq, tag::mesh, tag::reference >,
               pegtl::alpha >
           > > {};

  //! flow
  struct compflow :
         pegtl::if_must<
           tk::grm::readkw< use< kw::compflow >::pegtl_string >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic, tag::box >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::compflow,
                                            tag::depvar >,
                           mesh< tag::compflow >,
                           ic< tag::compflow >,
                           material_properties< tag::compflow >,
                           parameter< tag::compflow, kw::pde_alpha,
                                      tag::alpha >,
                           parameter< tag::compflow, kw::pde_p0,
                                      tag::p0 >,
                           parameter< tag::compflow, kw::pde_betax,
                                      tag::betax >,
                           parameter< tag::compflow, kw::pde_betay,
                                      tag::betay >,
                           parameter< tag::compflow, kw::pde_betaz,
                                      tag::betaz >,
                           parameter< tag::compflow, kw::pde_beta,
                                      tag::beta >,
                           parameter< tag::compflow, kw::pde_r0,
                                      tag::r0 >,
                           parameter< tag::compflow, kw::pde_ce,
                                      tag::ce >,
                           parameter< tag::compflow, kw::pde_kappa,
                                      tag::kappa >,
                           pde_parameter_vector< kw::pde_diffusivity,
                                                 tag::compflow,
                                                 tag::diffusivity >,
                           pde_parameter_vector< kw::pde_lambda,
                                                 tag::compflow,
                                                 tag::lambda >,
                           pde_parameter_vector< kw::pde_u0,
                                                 tag::compflow,
                                                 tag::u0 >,
                           pde_parameter_vector< kw::pde_source,
                                                 tag::compflow,
                                                 tag::source >,
                           bc_dirichlet,
                           bc_sym,
                           bc_farfield,
                           bc_pressure
                         >,
           check_errors< tag::compflow, tk::grm::check_compflow > > {};

  //! scalar transport
  struct transport :
         pegtl::if_must<
           tk::grm::readkw< use< kw::transport >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::transport,
                                            tag::depvar >,
                           mesh< tag::transport >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::transport >,
                           pde_parameter_vector< kw::pde_diffusivity,
                                                 tag::transport,
                                                 tag::diffusivity >,
                           pde_parameter_vector< kw::pde_lambda,
                                                 tag::transport,
                                                 tag::lambda >,
                           pde_parameter_vector< kw::pde_u0,
                                                 tag::transport,
                                                 tag::u0 >,
                           pde_parameter_vector< kw::pde_source,
                                                 tag::transport,
                                                 tag::source >
                         >,
           check_errors< tag::transport, tk::grm::check_transport > > {};

  //! partitioning ... end block
  struct partitioning :
         pegtl::if_must<
           tk::grm::readkw< use< kw::partitioning >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::process<
                             use< kw::algorithm >,
                             tk::grm::store_inciter_option<
                               tk::ctr::PartitioningAlgorithm,
                               tag::discr, tag::partitioner >,
                             pegtl::alpha > > > {};

  //! refinement variable(s) (refvar) ... end block
  struct refvars :
         pegtl::if_must<
           tk::grm::vector< use< kw::amr_refvar >,
                            tk::grm::match_depvar<
                              tk::grm::Store_back< tag::amr, tag::refvar > >,
                            use< kw::end >,
                            tk::grm::check_vector< tag::amr, tag::refvar >,
                            tk::grm::fieldvar< pegtl::alpha > >,
           tk::grm::compute_refvar_idx > {};

  //! adaptive mesh refinement (AMR) amr...end block
  struct amr :
         pegtl::if_must<
           tk::grm::readkw< use< kw::amr >::pegtl_string >,
           // enable AMR if amr...end block encountered
           tk::grm::enable< tag::amr >,
           tk::grm::block< use< kw::end >,
                           refvars,
                           edgelist,
                           coords,
                           tk::grm::process<
                             use< kw::amr_initial >,
                             tk::grm::store_back_option< use,
                                                         ctr::AMRInitial,
                                                         tag::amr,
                                                         tag::init >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::amr_error >,
                             tk::grm::store_inciter_option<
                               ctr::AMRError,
                               tag::amr, tag::error >,
                             pegtl::alpha >,
                           tk::grm::control< use< kw::amr_tolref >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::amr,
                                             tag::tolref >,
                           tk::grm::control< use< kw::amr_tolderef >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::amr,
                                             tag::tolderef >,
                           tk::grm::process< use< kw::amr_t0ref >,
                             tk::grm::Store< tag::amr, tag::t0ref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtref_uniform >,
                             tk::grm::Store< tag::amr, tag::dtref_uniform >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtref >,
                             tk::grm::Store< tag::amr, tag::dtref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtfreq >,
                             tk::grm::Store< tag::amr, tag::dtfreq >,
                             pegtl::digit >,
                           tk::grm::process< use< kw::amr_maxlevels >,
                             tk::grm::Store< tag::amr, tag::maxlevels >,
                             pegtl::digit >
                         >,
           tk::grm::check_amr_errors > {};

  //! field_output ... end block
  struct field_output :
         pegtl::if_must<
           tk::grm::readkw< use< kw::field_output >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::interval_iter< use< kw::interval_iter >,
                                     tag::output, tag::iter, tag::field >,
             tk::grm::interval_time< use< kw::interval_time >,
                                     tag::output, tag::time, tag::field >,
             tk::grm::time_range< use, kw::time_range,
                                  tag::output, tag::range, tag::field >,
             tk::grm::process<
               use< kw::refined >,
               tk::grm::Store< tag::cmd, tag::io, tag::refined >,
               pegtl::alpha >,
             pegtl::if_must<
               tk::grm::vector<
                 use< kw::sideset >,
                 tk::grm::Store_back< tag::cmd, tag::io, tag::surface,
                                      tag::field >,
                 use< kw::end > > >
           > > {};

  //! history_output ... end block
  struct history_output :
         pegtl::if_must<
           tk::grm::readkw< use< kw::history_output >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::interval_iter< use< kw::interval_iter >,
               tag::output, tag::iter, tag::history >,
             tk::grm::interval_time< use< kw::interval_time >,
               tag::output, tag::time, tag::history >,
             tk::grm::time_range< use, kw::time_range,
                                  tag::output, tag::range, tag::history >,
             tk::grm::precision< use, tag::history >,
             tk::grm::process<
               use< kw::txt_float_format >,
               tk::grm::store_inciter_option< tk::ctr::TxtFloatFormat,
                                              tag::flformat,
                                              tag::history >,
               pegtl::alpha >,
             pegtl::if_must<
               tk::grm::readkw< use< kw::point >::pegtl_string >,
               tk::grm::act< pegtl::identifier, tk::grm::match_pointname >,
               pegtl::seq<
                 tk::grm::start_vector< tag::history, tag::point >,
                 tk::grm::block<
                   use< kw::end >,
                   tk::grm::scan< tk::grm::number,
                     tk::grm::Store_back_back< tag::history, tag::point > > >
               > >
           > > {};

  //! integral_output ... end block
  struct integral_output :
         pegtl::if_must<
           tk::grm::readkw< use< kw::integral_output >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::interval_iter< use< kw::interval_iter >,
                                     tag::output, tag::iter, tag::integral >,
             tk::grm::interval_time< use< kw::interval_time >,
                                     tag::output, tag::time, tag::integral >,
             tk::grm::time_range< use, kw::time_range,
                                  tag::output, tag::range, tag::integral >,
             tk::grm::precision< use, tag::integral >,
             tk::grm::process<
               use< kw::txt_float_format >,
               tk::grm::store_inciter_option< tk::ctr::TxtFloatFormat,
                                              tag::flformat,
                                              tag::integral >,
               pegtl::alpha >,
             pegtl::if_must<
               tk::grm::vector<
                 use< kw::sideset >,
                 tk::grm::Store_back< tag::cmd, tag::io, tag::surface,
                                      tag::integral >,
                 use< kw::end > > >
           > > {};

  //! 'inciter' block
  struct inciter :
         pegtl::if_must<
           tk::grm::readkw< use< kw::inciter >::pegtl_string >,
           pegtl::sor<
             pegtl::seq< tk::grm::block<
                           use< kw::end >,
                           discretization,
                           compflow,
                           transport,
                           tk::grm::process<
                             use< kw::problem >,
                             tk::grm::store_inciter_option<
                               ctr::Problem, tag::problem >,
                             pegtl::alpha >,
                           amr,
                           partitioning,
                           field_output,
                           history_output,
                           integral_output,
                           tk::grm::diagnostics<
                             use,
                             tk::grm::store_inciter_option > >,
                         tk::grm::check_inciter >,
            tk::grm::msg< tk::grm::MsgType::ERROR,
                          tk::grm::MsgKey::UNFINISHED > > > {};

  //! \brief All keywords
  struct keywords :
         pegtl::sor< tk::grm::title< use >, inciter > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< keywords, tk::grm::ignore > {};

} // deck::
} // inciter::

#endif // InciterInputDeckGrammar_h
