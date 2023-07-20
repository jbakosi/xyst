// *****************************************************************************
/*!
  \file      src/Control/Keywords.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Definition of all keywords
  \details   This file contains the definition of all keywords, including those
    of command-line argument parsers as well as input, i.e., control, file
    parsers. All keywords are shared among all parsers of all executables.

    All keywords are case-sensitive.

    The information contained in this file is used to build data structures for
    on-screen help on the command-line arguments and control file keywords,
    available via the --help command-line arguments.

    A note on design: Defining structs that have static member functions
    returning a std::string is a way of storing C++-style strings at
    compile-time (which is not possible in a straightforward manner at this
    time). This could also be done with C-style const char* as well. The
    '*_info' structs store these strings, which then is used to specialize the
    _kw::keyword_ template, defined in Control/Keyword.hpp. Specializing the
    _keyword_ template also requires a specification of the precise string of
    characters that make up a keyword eventually matched by the parsers. Since
    these are all template arguments, the construction of keywords, their help,
    as well as all grammars, are entirely assembled at compile-time. Since the
    '*_info' struct member functions are static, they can be called without
    instantiating an object and thus available at compile-time.

    The definition of an '*_info' struct _requires_ at least the name, short,
    and long descriptions, defined by member functions _name()_,
    _shortDescription()_, and _longDescription()_, respectively. Everything else
    is optional. However, when adding a new keyword it is highly recommended to
    define all of the _optional_ members if they make sense for the given
    keyword. If an expect value type is also given, that can be hooked up into
    where it is used.

    The most general definition of a keyword is as follows:

    \code{.cpp}
      // Keyword info definition
      struct keyword_info {

        // Required very short name, usually a single word or (for e.g.
        // policies) a single character. This can be the keyword itself, but
        // does not have to be. This field is used as an id of the option or
        // setting.
        static std::string name() { return "Name"; }

        // Required short keyword description
        static std::string shortDescription() { return "Short description"; }

        // Required detailed keyword description. This returns a string literal,
        // since this is usually multi-line and it is less work to maintain this
        // way.
        static std::string longDescription() { return
          R"(Longer, possibly multi-line description of the keyword. Example
          usage of the keyword is welcome here. Don't worry about formatting:
          when this field is printed, extra spaces will be removed and line
          breaks will be inserted.)";
        }

        // Optional keyword alias. See also kw::Alias and
        // <tpl_install_dir>/include/pegtl/pegtl/constants.hh for examples
        // of what can be passed as template arguments (basically single
        // characters). The one below defines the character 'c' as the alias.
        // Aliases are single character long. The idea of an alias is to have a
        // long as well as a short keyword for the same functionality.
        // Currently, this is only hooked up for command-line arguments and not
        // for control-file keywords, which is intentional. Command-line
        // arguments are a lot less than control file keywords and are more
        // frequently typed by the user. Thus command-line argument aliases are
        // user-friendly. There are many control file keywords and aliases would
        // only cause confusion. Defining an alias for a command-line argument
        // enables the command-line parser to match on '--longer_keyword' as
        // well as on '-c'. Depending on whether the alias typedef is defined
        // for a keyword or not, the correct grammar is automatically generated
        // at compile-time, matching on both the longer keyword as well as on
        // the alias. Defining an alias for a control file keyword can be done
        // but has no effect in a control file parser.
        using alias = Alias< c >;

        // Optional single-character (policy) code. See also kw::Code and
        // <tpl_install_dir/include/pegtl/pegtl/constants.hh for examples
        // of what can be passed as template arguments (basically single
        // characters). The one below defines the character 'C' as the (policy)
        // code. This code is used for abbreviating policies used to configure
        // various orthogonal behaviors of classes using policy-based design.
        using code = Code< C >;

        // Optional expected data for the keyword - bundled to struct expect.
        // This struct is entirely optional within a keyword definition.
        // However, if it is defined, it must at least define the static member
        // function description() which returns the description of the type the
        // keyword expects. It may also optionally define the following fields:
        //
        //    - type - defining the expected type
        //    - lower - lower bound of the expected value
        //    - upper - upper bound of the expected value
        //    - choices - valid choices for the expected value
        //
        struct expect {

          // If this struct is defined, required expected type description, max
          // 10 characters long
          static std::string description() { return "int"; }

          // Optional expected type
          using type = std::streamsize;

          // Optional expected value lower bound
          static const type lower = 1;

          // Optional expected value upper bound
          static const type upper = 10;

          // Optional expected valid choices description, here giving
          // information on the expected type and the valid bounds. Note that
          // this can be any string, but if they exist, it is a good idea give
          // at least some information on the bounds, as below, since the bounds
          // are NOT displayed in the help for a keyword. This decision keeps
          // the bounds specifications generic since they can be any type. As a
          // result, the help structures, defined in HelpFactory.hpp, are simpler
          // as they do not have to be parameterized by the type of the bounds,
          // which greatly simplifies that code.
          static std::string choices() {
            return "integer [" + std::to_string(lower) + "..." +
                   std::to_string(upper) + ']';
          }

        };

      };

      // Keyword definition, passing the above info struct as the first
      // template argument, and the rest of the template arguments are the
      // characters of the keyword to be matched by the parser. Remember: all
      // keywords are case-sensitive. This one contrived example below defines
      // keyword 'kw', matching the keyword 'KeYwOrD'.
      using kw = keyword< keyword_info, K,e,Y,w,O,r,D >;
    \endcode
  \see Control/Keyword.hpp
  \see Control/HelpFactory.hpp
*/
// *****************************************************************************
#ifndef Keywords_h
#define Keywords_h

#include <limits>

#undef I
#include <pegtl/contrib/alphabet.hpp>

#include "Types.hpp"
#include "Keyword.hpp"
#include "XystBuildConfig.hpp"

//! Keywords used by all input deck and command line parsers
namespace kw {

using namespace tao::pegtl::alphabet;

struct title_info {
  static std::string name() { return "title"; }
  static std::string shortDescription() { return "Set analysis title"; }
  static std::string longDescription() { return
    R"(The analysis title may be specified in the input file using the 'title'
    keyword. The 'title' keyword must be followed by a double-quoted string
    specifying the analysis title. Example: title "Example problem".
    Specifying a title is optional.)";
  }
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using title = keyword< title_info, TAOCPP_PEGTL_STRING("title") >;

struct end_info {
  static std::string name() { return "end"; }
  static std::string shortDescription() { return "End of an input block"; }
  static std::string longDescription() { return
    R"(The end of a block is given by the 'end' keyword in the input file.)";
  }
};
using end = keyword< end_info, TAOCPP_PEGTL_STRING("end") >;

struct help_info {
  static std::string name() { return "help"; }
  static std::string shortDescription() { return
    R"(Display one-liner help on all command-line arguments)"; }
  static std::string longDescription() { return
    R"(Get a short one-liner help on all command-line arguments from an
    executable. It also triggers the help from the Charm++ runtime system and in
    addition to that of the executable, it also lists command-line arguments
    from Converse Machine, Tracing, Load Balancer, Record/Replay, and Charm++
    command-line parameters.)";
  }
  using alias = Alias< h >;
};
using help = keyword< help_info, TAOCPP_PEGTL_STRING("help") >;

struct txt_float_default_info {
  static std::string name() { return "default"; }
  static std::string shortDescription() { return
   "Select the default ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'default' floating-point output format for ASCII floating-point real
    number output. Example: "format default", which selects the default
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_default =
  keyword< txt_float_default_info, TAOCPP_PEGTL_STRING("default") >;

struct txt_float_scientific_info {
  static std::string name() { return "scientific"; }
  static std::string shortDescription() { return
   "Select the scientific ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'scientific' floating-point output format for ASCII floating-point real
    number output. Example: "format scientific", which selects the
    scientific floating-point output. Valid options are 'default', 'fixed',
    and 'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_scientific =
  keyword< txt_float_scientific_info, TAOCPP_PEGTL_STRING("scientific") >;

struct txt_float_fixed_info {
  static std::string name() { return "fixed"; }
  static std::string shortDescription() { return
   "Select the fixed ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'fixed' floating-point output format for ASCII floating-point real
    number output. Example: \"format fixed\", which selects the fixed
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_fixed = keyword< txt_float_fixed_info, TAOCPP_PEGTL_STRING("fixed") >;

struct txt_float_format_info {
  static std::string name() { return "float format"; }
  static std::string shortDescription() { return
    "Specify the ASCII floating-point output format"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    floating-point output format for ASCII floating-point number output.
    Example: "format scientific", which selects the scientific
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + txt_float_default::string() + "\' | \'"
                  + txt_float_scientific::string() + "\' | \'"
                  + txt_float_fixed::string() + '\'';
    }
  };
};
using txt_float_format =
  keyword< txt_float_format_info, TAOCPP_PEGTL_STRING("format") >;

struct precision_info {
  static std::string name() { return "precision"; }
  static std::string shortDescription() { return
    R"(Precision in digits for ASCII floating-point output)"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    the precision in digits for ASCII floating-point real number output.
    Example: "precision 10", which selects ten digits for floating-point
    output, e.g., 3.141592654. The number of digits must be larger than zero
    and lower than the maximum representable digits for the given floating-point
    type. For more info on setting the precision in C++, see
    http://en.cppreference.com/w/cpp/io/manip/setprecision, and
    http://en.cppreference.com/w/cpp/types/numeric_limits/digits10)";
  }
  struct expect {
    using type = std::streamsize;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< tk::real >::digits10 + 1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using precision = keyword< precision_info, TAOCPP_PEGTL_STRING("precision") >;

struct elem_info {
  static std::string name() { return "elem"; }
  static std::string shortDescription() { return
    "Specify elem-centering for output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select elem-centering for variable output. In
    walker for example, this is used to configure probability values on the
    sample space grid for file output of probability density functions (PDFs).
    Example: "centering elem", which selects element-centered values. Valid
    options are 'elem' and 'node', denoting cell-centered and point-centered
    output, respectively. In inciter this keyword is used in output variable
    specification blocks, prefixing variable names by either 'node' or 'elem',
    to specify their centering for output to file.)"; }
};
using elem = keyword< elem_info, TAOCPP_PEGTL_STRING("elem") >;

struct node_info {
  static std::string name() { return "node"; }
  static std::string shortDescription() { return
    "Specify node-centering for output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select node-centering for variable output. In
    walker for example, this is used to configure probability values on the
    sample space grid for file output of probability density functions (PDFs).
    Example: "centering elem", which selects element-centered values. Valid
    options are 'elem' and 'node', denoting cell-centered and point-centered
    output, respectively. In inciter this keyword is used in output variable
    specification blocks, prefixing variable names by either 'node' or 'elem',
    to specify their centering for output to file.)"; }
};
using node = keyword< node_info, TAOCPP_PEGTL_STRING("node") >;

struct centering_info {
  static std::string name() { return "centering"; }
  static std::string shortDescription() { return
    "Specify data-centering for PDF output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the data
    centering of the probability value output on the sample space grid for
    file output of probability density functions (PDFs). Example: "centering
    elem", which selects element-centered values. Valid options are 'elem'
    and 'node', denoting cell-centered and point-centered output,
    respectively.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + elem::string() + "\' | \'"
                  + node::string() + '\'';
    }
  };
};
using pdf_centering = keyword< centering_info, TAOCPP_PEGTL_STRING("centering") >;

struct nstep_info {
  static std::string name() { return "nstep"; }
  static std::string shortDescription() { return
    "Set number of time steps to take"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of time steps to take in a
    simulation. The number of time steps are used in conjunction with the
    maximmum time specified by keyword 'term': the simulation stops whichever is
    reached first. Both 'nstep' and 'term' can be left unspecified, in which
    case their default values are used. See also 'term'.)";
  }
  struct expect {
    using type = uint64_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using nstep = keyword< nstep_info, TAOCPP_PEGTL_STRING("nstep") >;

struct term_info {
  static std::string name() { return "term"; }
  static std::string shortDescription() { return
    "Set maximum non-dimensional time to simulate"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the termination time in a simulation. The
    termination time and number of time steps, specified by 'nstep', are used in
    conjunction to determine when to stop a simulation: whichever is
    reached first. Both 'nstep' and 'term' can be left unspecified, in which
    case their default values are used. See also 'nstep'.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using term = keyword< term_info, TAOCPP_PEGTL_STRING("term") >;

struct t0_info {
  static std::string name() { return "t0"; }
  static std::string shortDescription() { return
    "Set starting non-dimensional time"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the starting time in a simulation.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using t0 = keyword< t0_info, TAOCPP_PEGTL_STRING("t0") >;

struct dt_info {
  static std::string name() { return "dt"; }
  static std::string shortDescription() { return
    "Select constant time step size"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the time step size that used as a
    constant during simulation. Setting 'cfl' and 'dt' are mutually
    exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using dt = keyword< dt_info, TAOCPP_PEGTL_STRING("dt") >;

struct cfl_info {
  static std::string name() { return "CFL"; }
  static std::string shortDescription() { return
    "Set the Courant-Friedrichs-Lewy (CFL) coefficient"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the CFL coefficient for
    variable-time-step-size simulations. Setting 'cfl' and 'dt' are mutually
    exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using cfl = keyword< cfl_info, TAOCPP_PEGTL_STRING("cfl") >;

struct ncomp_info {
  static std::string name() { return "ncomp"; }
  static std::string shortDescription() { return
    "Set number of scalar components for a system of differential equations"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of scalar
    components of a vector. 'ncomp' means "number of components". It is also
    used for specifying the number of scalar components of a transporter scalar
    (see also the keywords 'transport').)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using ncomp = keyword< ncomp_info,  TAOCPP_PEGTL_STRING("ncomp") >;

struct ttyi_info {
  static std::string name() { return "ttyi"; }
  static std::string shortDescription() { return
    "Set screen output interval"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the interval in time steps for screen
    output during a simulation.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using ttyi = keyword< ttyi_info, TAOCPP_PEGTL_STRING("ttyi") >;

struct pari_info {
  static std::string name() { return "pari"; }
  static std::string shortDescription() { return
    "Set particles output  interval"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the interval in time steps for particles
    output during a simulation.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using pari = keyword< pari_info, TAOCPP_PEGTL_STRING("pari") >;

struct interval_iter_info {
  static std::string name() { return "interval"; }
  static std::string shortDescription() { return
    "Set interval (in units of iteration count)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify an interval in units of iteration count
    (i.e., number of time steps). This must be used within a relevant block.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using interval_iter =
  keyword< interval_iter_info, TAOCPP_PEGTL_STRING("interval") >;

struct interval_time_info {
  static std::string name() { return "time_interval"; }
  static std::string shortDescription() { return
    "Set interval (in units of physics time)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify an interval in units of physics time.
    This must be used within a relevant block.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = std::numeric_limits< tk::real >::epsilon();
    static std::string description() { return "real"; }
  };
};
using interval_time =
  keyword< interval_time_info, TAOCPP_PEGTL_STRING("time_interval") >;

struct time_range_info {
  static std::string name() { return "time_range"; }
  static std::string shortDescription() { return
    "Configure physics time range for output (in units of physics time)"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure field-, or history-output, specifying
    a start time, a stop time, and an output frequency in physics time units.
    Example: 'time_range 0.2 0.3 0.001 end', which specifies that from t=0.2 to
    t=0.3 output should happen at physics time units of dt=0.001. This must be
    used within a relevant block.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "3 reals"; }
  };
};
using time_range =
  keyword< time_range_info, TAOCPP_PEGTL_STRING("time_range") >;

struct statistics_info {
  static std::string name() { return "statistics"; }
  static std::string shortDescription() { return
    "Start of statistics input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested output for statistical moments.
    Example: "statistics <Y> <yy> end", which requests the first two moments of
    the flutcutating variable 'Y'. For more info on the structure of the
    statistics ... end block, see doc/pages/statistics_output.dox.)";
  }
};
using statistics = keyword< statistics_info, TAOCPP_PEGTL_STRING("statistics") >;

struct history_output_info {
  static std::string name() { return "history_output"; }
  static std::string shortDescription() { return
    "Start of history_output input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested history output.)";
  }
};
using history_output =
  keyword< history_output_info, TAOCPP_PEGTL_STRING("history_output") >;

struct integral_output_info {
  static std::string name() { return "integral_output"; }
  static std::string shortDescription() { return
    "Start of integral_output input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested integral output.)";
  }
};
using integral_output =
  keyword< integral_output_info, TAOCPP_PEGTL_STRING("integral_output") >;

struct field_output_info {
  static std::string name() { return "field_output"; }
  static std::string shortDescription() { return
    "Start of field_output input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    list and settings of requested field output.)";
  }
};
using field_output =
  keyword< field_output_info, TAOCPP_PEGTL_STRING("field_output") >;

struct velocity_info {
  static std::string name() { return "velocity"; }
  static std::string shortDescription() { return "Specify velocity"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure a velocity vector, used for, e.g.,
    boundary or initial conditions or as a keyword that selects velocity in some
    other context-specific way, e.g., 'velocity' as opposed to 'position'.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using velocity = keyword< velocity_info, TAOCPP_PEGTL_STRING("velocity") >;

struct acceleration_info {
  static std::string name() { return "acceleration"; }
  static std::string shortDescription() { return "Specify acceleration"; }
  static std::string longDescription() { return
    R"(This keyword is used as a keyword that selects acceleration in some
    other context-specific way, e.g., as opposed to 'velocity' or 'position'.)";
  }
};
using acceleration =
  keyword< acceleration_info, TAOCPP_PEGTL_STRING("acceleration") >;

struct materialid_info {
  static std::string name() { return "materialid"; }
  static std::string shortDescription() { return "Specify material id"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure the material id within a box as a part
    of the initialization.)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using materialid = keyword< materialid_info,
  TAOCPP_PEGTL_STRING("materialid") >;

struct density_info {
  static std::string name() { return "density"; }
  static std::string shortDescription() { return "Specify density"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure a density, used for, e.g., boundary or
    initial conditions.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using density = keyword< density_info, TAOCPP_PEGTL_STRING("density") >;

struct pressure_info {
  static std::string name() { return "pressure"; }
  static std::string shortDescription() { return "Specify pressure"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure a pressure, used for, e.g., boundary or
    initial conditions or as a keyword that selects pressure in some other
    context-specific way.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pressure = keyword< pressure_info, TAOCPP_PEGTL_STRING("pressure") >;

struct energy_info {
  static std::string name() { return "energy"; }
  static std::string shortDescription() { return
    "Specify energy per unit mass"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure energy per unit mass, used for, e.g.,
    boundary or initial conditions.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using energy = keyword< energy_info, TAOCPP_PEGTL_STRING("energy") >;

struct temperature_info {
  static std::string name() { return "temperature"; }
  static std::string shortDescription() { return "Specify temperature"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure temperature, used for, e.g.,
    boundary or initial conditions.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using temperature =
  keyword< temperature_info, TAOCPP_PEGTL_STRING("temperature") >;

struct xmin_info {
  static std::string name() { return "xmin"; }
  static std::string shortDescription() { return "Minimum x coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a minimum x coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using xmin = keyword< xmin_info, TAOCPP_PEGTL_STRING("xmin") >;

struct xmax_info {
  static std::string name() { return "xmax"; }
  static std::string shortDescription() { return "Maximum x coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a maximum x coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using xmax = keyword< xmax_info, TAOCPP_PEGTL_STRING("xmax") >;

struct ymin_info {
  static std::string name() { return "ymin"; }
  static std::string shortDescription() { return "Minimum y coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a minimum y coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using ymin = keyword< ymin_info, TAOCPP_PEGTL_STRING("ymin") >;

struct ymax_info {
  static std::string name() { return "ymax"; }
  static std::string shortDescription() { return "Maximum y coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a maximum y coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using ymax = keyword< ymax_info, TAOCPP_PEGTL_STRING("ymax") >;

struct zmin_info {
  static std::string name() { return "zmin"; }
  static std::string shortDescription() { return "Minimum z coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a minimum z coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using zmin = keyword< zmin_info, TAOCPP_PEGTL_STRING("zmin") >;

struct zmax_info {
  static std::string name() { return "zmax"; }
  static std::string shortDescription() { return "Maximum z coordinate"; }
  static std::string longDescription() { return
    R"(This keyword used to configure a maximum z coordinate, e.g., to specify
    a box.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using zmax = keyword< zmax_info, TAOCPP_PEGTL_STRING("zmax") >;

struct box_info {
  static std::string name() { return "box"; }
  static std::string shortDescription() { return
    R"(Introduce a box ... end block used to assign initial conditions)"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a box ... end block used to assign
    initial conditions within a box given by spatial coordinates. Example:
    box x- 0.5 x+ 1.5 y- -0.5 y+ 0.5 z- -0.5 z+ 0.5 density 1.2 end pressure
    1.4 end end", which specifies a box with extends within which the density
    will be set to 1.2 and the pressure to be 1.4. Besides the box dimensions,
    the following physics keywords are allowed in a box ... end block:)"
    + std::string("\'")
    + density::string()+ "\', \'"
    + velocity::string() + "\', \'"
    + energy::string() + "\', \'"
    + temperature::string() + "\', \'"
    + pressure::string() + "\'."; }
};
using box = keyword< box_info, TAOCPP_PEGTL_STRING("box") >;

struct ic_info {
  static std::string name() { return "ic"; }
  static std::string shortDescription() { return
    R"(Introduce an ic...end block used to configure initial conditions)"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an ic...end block used to set initial
    conditions. Keywords allowed in a ic ... end block: )" + std::string("\'")
    + materialid::string()+ "\', \'"
    + density::string()+ "\', \'"
    + velocity::string() + "\', \'"
    + pressure::string() + "\', \'"
    + energy::string() + "\', \'"
    + temperature::string() + "\', \'"
    + box::string() + "\'.";
  }
};
using ic = keyword< ic_info, TAOCPP_PEGTL_STRING("ic") >;

struct depvar_info {
  static std::string name() { return "depvar"; }
  static std::string shortDescription() { return
    "Select dependent variable (in a relevant block)"; }
  static std::string longDescription() { return
    R"(Dependent variable, e.g, in differential equations.)"; }
  struct expect {
    using type = char;
    static std::string description() { return "character"; }
  };
};
using depvar = keyword< depvar_info, TAOCPP_PEGTL_STRING("depvar") >;

struct control_info {
  static std::string name() { return "control"; }
  static std::string shortDescription()
  { return "Specify the control file name [REQUIRED]"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the control file from which
    detailed user input is parsed.)";
  }
  using alias = Alias< c >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using control = keyword< control_info, TAOCPP_PEGTL_STRING("control") >;

struct benchmark_info {
  static std::string name() { return "benchmark"; }
  static std::string shortDescription() { return "Select benchmark mode"; }
  static std::string longDescription() { return
    R"(This keyword is used to select benchmark mode. In benchmark mode no large
       file output is performed, overriding the configuration in the control
       file.)";
  }
  using alias = Alias< b >;
};

using benchmark = keyword< benchmark_info, TAOCPP_PEGTL_STRING("benchmark") >;

struct nonblocking_info {
  static std::string name() { return "nonblocking"; }
  static std::string shortDescription()
  { return "Select non-blocking migration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select non-blocking, instead of the default
       blocking, migration. WARNING: This feature is experimental, not well
       tested, and may not always work as expected.)";
  }
  using alias = Alias< n >;
};

using nonblocking =
  keyword< nonblocking_info, TAOCPP_PEGTL_STRING("nonblocking") >;

struct lbfreq_info {
  static std::string name() { return "Load balancing frequency"; }
  static std::string shortDescription()
  { return "Set load-balancing frequency during time stepping"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the frequency of load-balancing during
       time stepping. The default is 1, which means that load balancing is
       initiated every time step. Note, however, that this does not necessarily
       mean that load balancing will be performed by the runtime system every
       time step, only that the Charm++ load-balancer is initiated. For more
       information, see the Charm++ manual.)";
  }
  using alias = Alias< l >;
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max()-1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using lbfreq = keyword< lbfreq_info, TAOCPP_PEGTL_STRING("lbfreq") >;

struct rsfreq_info {
  static std::string name() { return "Checkpoint/restart frequency"; }
  static std::string shortDescription()
  { return "Set checkpoint/restart frequency during time stepping"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the frequency of dumping checkpoint/restart
       files during time stepping. The default is 1000, which means that
       checkpoint/restart files are dumped at every 1000th time step.)";
  }
  using alias = Alias< r >;
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max()-1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using rsfreq = keyword< rsfreq_info, TAOCPP_PEGTL_STRING("rsfreq") >;

struct feedback_info {
  static std::string name() { return "feedback"; }
  static std::string shortDescription() { return "Enable on-screen feedback"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable more detailed on-screen feedback on
       particular tasks and sub-tasks as they happen. This is useful for large
       problems and debugging.)";
  }
  using alias = Alias< f >;
};
using feedback = keyword< feedback_info, TAOCPP_PEGTL_STRING("feedback") >;

struct version_info {
  static std::string name() { return "Show version"; }
  static std::string shortDescription() { return "Show version information"; }
  static std::string longDescription() { return
    R"(This keyword is used to display version information for the
       executable/tool on the standard output and exit successfully.)";
  }
  using alias = Alias< V >;
};
using version = keyword< version_info, TAOCPP_PEGTL_STRING("version") >;

struct raise_signal_info {
  static std::string name() { return "Raise signal (for testing)"; }
  static std::string shortDescription() { return "Raise signal"; }
  static std::string longDescription() { return
    R"(This keyword is used to test signal handling.)";
  }
  using alias = Alias< I >;
};
using raise_signal = keyword< raise_signal_info,
                              TAOCPP_PEGTL_STRING("raise_signal") >;

struct trace_info {
  static std::string name() { return "trace"; }
  static std::string shortDescription()
  { return "Disable call and stack trace"; }
  static std::string longDescription() { return R"(This keyword can be used to
    disable the on-screen call trace and stack trace after an exception is
    thrown. Trace output is on by default and in some cases, the call and
    stack trace can be huge and not very helpful, hence this command line
    option.)"; }
  using alias = Alias< t >;
};
using trace = keyword< trace_info, TAOCPP_PEGTL_STRING("trace") >;

struct quiescence_info {
  static std::string name() { return "quiescence"; }
  static std::string shortDescription()
  { return "Enable quiescence detection"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable the quiescence detection feature of
       Charm++, used to catch logic errors in the asynchronous control flow,
       resulting in deadlocks. This is useful for automated testing and
       debugging and does have some overhead, so it is off by default.)";
  }
  using alias = Alias< q >;
};
using quiescence =
  keyword< quiescence_info, TAOCPP_PEGTL_STRING("quiescence") >;

struct virtualization_info {
  static std::string name() { return "virtualization"; }
  static std::string shortDescription() { return
    R"(Set degree of virtualization)"; }
  static std::string longDescription() { return
    R"(This option is used to set the degree of virtualization
    (over-decomposition). The virtualization parameter is a real number
    between 0.0 and 1.0, inclusive, which controls the degree of
    virtualization or over-decomposition. Independent of the value of
    virtualization the work is approximately evenly distributed among the
    available processing elements. For zero virtualization (no
    over-decomposition), the work is simply decomposed into
    total_work/numPEs, which yields the smallest number of Charm++ chares and
    the largest chunks of work units. The other extreme is unity
    virtualization, which decomposes the total work into the smallest size
    work units possible, yielding the largest number of Charm++ chares.
    Obviously, the optimum will be between 0.0 and 1.0, depending on the
    problem.)";
  }
  using alias = Alias< u >;
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static constexpr type upper = 1.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "real between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using virtualization =
  keyword< virtualization_info, TAOCPP_PEGTL_STRING("virtualization") >;

struct pdf_info {
  static std::string name() { return "pdf"; }
  static std::string shortDescription() { return
    "Specify the name of the PDF output file"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the output file in which to
    store probability density functions (PDFs) during a simulation.)";
  }
  using alias = Alias< p >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using pdf = keyword< pdf_info, TAOCPP_PEGTL_STRING("pdf") >;

struct stat_info {
  static std::string name() { return "stat"; }
  static std::string shortDescription() { return
    "Specify the name of the statistical moments output file"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the output file in which to
    store statistical moments during a simulation.)";
  }
  using alias = Alias< s >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using stat = keyword< stat_info, TAOCPP_PEGTL_STRING("stat") >;

struct particles_info {
  static std::string name() { return "particles"; }
  static std::string shortDescription() { return
    "Specify the name of the particles position output file"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the output file in which to
    store particles positions during a simulation.)";
  }
  using alias = Alias< x >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using particles = keyword< particles_info, TAOCPP_PEGTL_STRING("particles") >;

struct input_info {
  static std::string name() { return "input"; }
  static std::string shortDescription() { return "Specify the input file"; }
  static std::string longDescription() { return
    R"(This option is used to define the name of input file.)";
  }
  using alias = Alias< i >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using input = keyword< input_info, TAOCPP_PEGTL_STRING("input") >;

struct output_info {
  static std::string name() { return "output"; }
  static std::string shortDescription() { return "Specify the output file"; }
  static std::string longDescription() { return
    R"(This option is used to define the output file name. In MeshConv, this is
    used to specify the output mesh file name. In Inciter this is used to
    specify the output base filename. The base filename is appended by
    ".e-s.<meshid>.<numchares>.<chareid>", where 'e-s' probably stands for
    ExodusII sequence (the output file format), <meshid> counts the number of
    new meshes (this is incremented whenever the mesh is new compared to the
    previous iteration, due to, e.g., mesh refinement), <numchares> is the total
    number of mesh partitions, and <chareid> is the work unit (or mesh
    partition) id.)";
  }
  using alias = Alias< o >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using output = keyword< output_info, TAOCPP_PEGTL_STRING("output") >;

struct refined_info {
  static std::string name() { return "Refined field output"; }
  static std::string shortDescription() { return
    "Turn refined field output on/off"; }
  static std::string longDescription() { return
    R"(This keyword can be used to turn on/off refined field output, which
    refines the mesh and evaluates the solution on the refined mesh for saving
    the solution.)"; }
  struct expect {
    using type = bool;
    static std::string description() { return "string"; }
    static std::string choices() { return "true | false"; }
  };
};
using refined =keyword< refined_info, TAOCPP_PEGTL_STRING("refined") >;

struct restart_info {
  static std::string name() { return "checkpoint/restart directory name"; }
  static std::string shortDescription()
    { return "Specify the directory for restart files"; }
  static std::string longDescription() { return
    R"(This option is used to specify the directory name in which to save
    checkpoint/restart files.)";
  }
  using alias = Alias< R >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using restart = keyword< restart_info, TAOCPP_PEGTL_STRING("restart") >;

struct diagnostics_cmd_info {
  static std::string name() { return "diagnostics"; }
  static std::string shortDescription()
  { return "Specify the diagnostics file name"; }
  static std::string longDescription() { return
    R"(This option is used to define the diagnostics file name.)";
  }
  using alias = Alias< d >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using diagnostics_cmd =
  keyword< diagnostics_cmd_info, TAOCPP_PEGTL_STRING("diagnostics") >;

struct diagnostics_info {
  static std::string name() { return "diagnostics"; }
  static std::string shortDescription()
  { return "Specify the diagnostics file name"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the dagnostics ... end block, used to
    configure diagnostics output. Keywords allowed in this block: )"
    + std::string("\'")
    + interval_iter::string() + "\' | \'"
    + txt_float_format::string() + "\' | \'"
    + precision::string() + "\'.";
  }
};
using diagnostics =
  keyword< diagnostics_info, TAOCPP_PEGTL_STRING("diagnostics") >;

struct reorder_cmd_info {
  static std::string name() { return "reorder"; }
  static std::string shortDescription() { return "Reorder mesh nodes"; }
  static std::string longDescription() { return
    R"(This keyword is used as a command line argument to instruct the mesh
    converter to not only convert but also reorder the mesh nodes using the
    advancing front technique. Reordering is optional in meshconv and
    inciter.)";
  }
  using alias = Alias< r >;
  struct expect {
    using type = bool;
    static std::string description() { return "string"; }
  };
};
using reorder_cmd = keyword< reorder_cmd_info, TAOCPP_PEGTL_STRING("reorder") >;

struct pelocal_reorder_info {
  static std::string name() { return "PE-local reorder"; }
  static std::string shortDescription() { return "PE-local reorder"; }
  static std::string longDescription() { return
    R"(This keyword is used in inciter as a keyword in the inciter...end block
    as "pelocal_reorder true" (or false) to do (or not do) a global distributed
    mesh reordering across all PEs that yields an approximately continous mesh
    node ID order as mesh partitions are assigned to PEs after mesh
    partitioning. This reordering is optional.)";
  }
  struct expect {
    using type = bool;
    static std::string choices() { return "true | false"; }
    static std::string description() { return "string"; }
  };
};
using pelocal_reorder =
  keyword< pelocal_reorder_info, TAOCPP_PEGTL_STRING("pelocal_reorder") >;

struct steady_state_info {
  static std::string name() { return "steady_state"; }
  static std::string shortDescription() { return "March to steady state"; }
  static std::string longDescription() { return
    R"(This keyword is used indicate that local time stepping should be used
       towards a stationary solution.)"; }
  struct expect {
    using type = bool;
    static std::string choices() { return "true | false"; }
    static std::string description() { return "string"; }
  };
};
using steady_state =
  keyword< steady_state_info, TAOCPP_PEGTL_STRING("steady_state") >;

struct residual_info {
  static std::string name() { return "residual"; }
  static std::string shortDescription() { return
    "Set the convergence criterion for the residual to reach"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a convergence criterion for, e.g., local
    time stepping marching to steady state, below which the simulation is
    considered converged.)"; }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 1.0e-14;
    static std::string description() { return "real"; }
  };
};
using residual = keyword< residual_info, TAOCPP_PEGTL_STRING("residual") >;

struct rescomp_info {
  static std::string name() { return "rescomp"; }
  static std::string shortDescription() { return
    "Equation system component index for convergence"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a single integer that is used to denote
    the equation component index in the complete system of equation systems
    configured in an input file to use for the convergence criterion for local
    time stepping marching towards steady state.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using rescomp = keyword< rescomp_info, TAOCPP_PEGTL_STRING("rescomp") >;

struct group_info {
  static std::string name() { return "group"; }
  static std::string shortDescription() { return
    "Select test group(s) to run"; }
  static std::string longDescription() { return
    R"(This option can be used to select one or more test groups to run by
    specifying the full or a partial name of a test group. All tests of a
    selected group will be executed. If this option is not given, all test
    groups are executed by default. Examples: '--group make_list' - run only
    the 'make_list' test group, '--group Parser' - run the test groups that have
    the string 'Parser' in their name, e.g., groups 'Control/FileParser' and
    'Control/StringParser'.)";
  }
  using alias = Alias< g >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using group = keyword< group_info, TAOCPP_PEGTL_STRING("group") >;

struct inciter_info {
  static std::string name() { return "inciter"; }
  static std::string shortDescription() { return
    "Start configuration block for inciter"; }
  static std::string longDescription() { return
    R"(This keyword is used to select inciter. Inciter, is a continuum-realm
    shock hydrodynamics tool, solving a PDE.)";
  }
};
using inciter = keyword< inciter_info, TAOCPP_PEGTL_STRING("inciter") >;

struct user_defined_info {
  using code = Code< U >;
  static std::string name() { return "User-defined"; }
  static std::string shortDescription() { return
    "Select user-defined specification for a problem"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the user-defined specification for an
    option. This could be a 'problem' to be solved by a partial differential
    equation, but can also be a 'user-defined' mesh velocity specification for
    ALE mesh motion.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};

using user_defined =
  keyword< user_defined_info, TAOCPP_PEGTL_STRING("user_defined") >;

struct shear_diff_info {
  using code = Code< S >;
  static std::string name() { return "Shear-diffusion"; }
  static std::string shortDescription() { return
    "Select the shear + diffusion test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the shear diffusion test problem. The
    initial and boundary conditions are specified to set up the test problem
    suitable to exercise and test the advection and diffusion terms of the
    scalar transport equation. Example: "problem shear_diff".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using shear_diff = keyword< shear_diff_info, TAOCPP_PEGTL_STRING("shear_diff") >;

struct point_src_info {
  using code = Code< P >;
  static std::string name() { return "Point source"; }
  static std::string shortDescription() { return
    "Select the point source test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the point source test problem. The
    initial and boundary conditions are specified to set up the test problem
    suitable to exercise and test the advection and diffusion terms of the
    scalar transport equation. Example: "problem point_src".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using point_src = keyword< point_src_info, TAOCPP_PEGTL_STRING("point_src") >;

struct slot_cyl_info {
  using code = Code< Z >;
  static std::string name() { return "Zalesak's slotted cylinder"; }
  static std::string shortDescription() { return
    "Select Zalesak's slotted cylinder test problem"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Zalesak's slotted cylinder test
    problem. The initial and boundary conditions are specified to set up the
    test problem suitable to exercise and test the advection and diffusion
    terms of the scalar transport equation. Example: "problem slot_cyl".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using slot_cyl = keyword< slot_cyl_info, TAOCPP_PEGTL_STRING("slot_cyl") >;

struct vortical_flow_info {
  using code = Code< V >;
  static std::string name() { return "Vortical flow"; }
  static std::string shortDescription() { return
    "Select the vortical flow test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the vortical flow test problem. The
    purpose of this test problem is to test velocity errors generated by spatial
    operators in the presence of 3D vorticity and in particluar the
    superposition of planar and vortical flows, analogous to voritcity
    stretching. Example: "problem vortical_flow. For more details, see Waltz,
    et. al, "Manufactured solutions for the three-dimensional Euler equations
    with relevance to Inertial Confinement Fusion", Journal of Computational
    Physics 267 (2014) 196-209.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using vortical_flow =
  keyword< vortical_flow_info, TAOCPP_PEGTL_STRING("vortical_flow") >;

struct nonlin_ener_growth_info {
  using code = Code< N >;
  static std::string name() { return "Nonlinear energy growth"; }
  static std::string shortDescription() { return
    "Select the nonlinear energy growth test problem ";}
  static std::string longDescription() { return
    R"(This keyword is used to select the nonlinear energy growth test problem.
    The purpose of this test problem is to test nonlinear, time dependent energy
    growth and the subsequent development of pressure gradients due to coupling
    between the internal energy and the equation of state. Example: "problem
    nonlinear_energy_growth". For more details, see Waltz, et. al, "Manufactured
    solutions for the three-dimensional Euler equations with relevance to
    Inertial Confinement Fusion", Journal of Computational Physics 267 (2014)
    196-209.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using nonlin_ener_growth =
  keyword< nonlin_ener_growth_info,
           TAOCPP_PEGTL_STRING("nonlinear_energy_growth") >;

struct rayleigh_taylor_info {
  using code = Code< R >;
  static std::string name() { return "Rayleigh-Taylor"; }
  static std::string shortDescription() { return
    "Select the Rayleigh-Taylor test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Rayleigh-Taylor unstable configuration
    test problem. The purpose of this test problem is to assess time dependent
    fluid motion in the presence of Rayleigh-Taylor unstable conditions, i.e.
    opposing density and pressure gradients. Example: "problem rayleigh_taylor".
    For more details, see Waltz, et. al, "Manufactured solutions for the
    three-dimensional Euler equations with relevance to Inertial Confinement
    Fusion", Journal of Computational Physics 267 (2014) 196-209.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using rayleigh_taylor =
  keyword< rayleigh_taylor_info, TAOCPP_PEGTL_STRING("rayleigh_taylor") >;

struct taylor_green_info {
  using code = Code< T >;
  static std::string name() { return "Taylor-Green"; }
  static std::string shortDescription() { return
    "Select the Taylor-Green test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Taylor-Green vortex test problem. The
    purpose of this problem is to test time accuracy and the correctness of the
    discretization of the viscous term in the Navier-Stokes equation. Example:
    "problem taylor_green". For more details on the flow, see G.I. Taylor, A.E.
    Green, "Mechanism of the Production of Small Eddies from Large Ones", Proc.
    R. Soc. Lond. A 1937 158 499-521; DOI: 10.1098/rspa.1937.0036. Published 3
    February 1937.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using taylor_green =
  keyword< taylor_green_info, TAOCPP_PEGTL_STRING("taylor_green") >;

struct sod_info {
  using code = Code< H >;
  static std::string name() { return "Sod shock-tube"; }
  static std::string shortDescription() { return
    "Select the Sod shock-tube test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Sod shock-tube test problem. The
    purpose of this test problem is to test the correctness of the
    approximate Riemann solver and its shock and interface capturing
    capabilities. Example: "problem sod_shocktube". For more details, see
    G. A. Sod, "A Survey of Several Finite Difference Methods for Systems of
    Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27 (1978)
    1–31.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using sod = keyword< sod_info, TAOCPP_PEGTL_STRING("sod") >;

struct rotated_sod_info {
  using code = Code< O >;
  static std::string name() { return "Rotated Sod shock-tube"; }
  static std::string shortDescription() { return
    "Select the rotated Sod shock-tube test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the rotated Sod shock-tube test problem.
    This the same as Sod shocktube but the geometry is rotated about X, Y, Z
    each by 45 degrees (in that order) so that none of the domain boundary align
    with any of the coordinate directions. The purpose of this test problem is
    to test the correctness of the approximate Riemann solver and its shock and
    interface capturing capabilities in an arbitrarily oriented geometry.
    Example: "problem rotated_sod_shocktube". For more details on the Sod
    problem, see G. A. Sod, "A Survey of Several Finite Difference Methods for
    Systems of Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27
    (1978) 1–31.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using rotated_sod =
  keyword< rotated_sod_info, TAOCPP_PEGTL_STRING("rotated_sod") >;

struct sedov_info {
  using code = Code< B >;
  static std::string name() { return "Sedov blast-wave"; }
  static std::string shortDescription() { return
    "Select the Sedov blast-wave test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Sedov blast-wave test problem. The
    purpose of this test problem is to test the correctness of the
    approximate Riemann solver and its strong shock and interface capturing
    capabilities. Example: "problem sedov".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using sedov = keyword< sedov_info, TAOCPP_PEGTL_STRING("sedov") >;


struct problem_info {
  using code = Code< t >;
  static std::string name() { return "Test problem"; }
  static std::string shortDescription() { return
    "Specify problem configuration for a partial differential equation solver";
  }
  static std::string longDescription() { return
    R"(This keyword is used to specify the problem configuration for a partial
    differential equation solver in the input file.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + user_defined::string() + "\' | \'"
                  + shear_diff::string() + "\' | \'"
                  + point_src::string() + "\' | \'"
                  + slot_cyl::string() + "\' | \'"
                  + vortical_flow::string() + "\' | \'"
                  + nonlin_ener_growth::string() + "\' | \'"
                  + rayleigh_taylor::string() + "\' | \'"
                  + taylor_green::string() + "\' | \'"
                  + sod::string() + "\' | \'"
                  + rotated_sod::string() + '\'';
    }
  };
};
using problem = keyword< problem_info, TAOCPP_PEGTL_STRING("problem") >;

struct navierstokes_info {
  using code = Code< N >;
  static std::string name() { return "Navier-Stokes"; }
  static std::string shortDescription() { return "Specify the Navier-Stokes "
    "(viscous) compressible flow physics configuration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Navier-Stokes (viscous) compressible
    flow physics configuration. Example: "compflow physics navierstokes end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using navierstokes =
  keyword< navierstokes_info, TAOCPP_PEGTL_STRING("navierstokes") >;

struct euler_info {
  using code = Code< E >;
  static std::string name() { return "Euler"; }
  static std::string shortDescription() { return "Specify the Euler (inviscid) "
    "compressible flow physics configuration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Euler (inviscid) compressible
    flow physics configuration. Example: "compflow physics euler end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using euler = keyword< euler_info, TAOCPP_PEGTL_STRING("euler") >;

struct advection_info {
  using code = Code< A >;
  static std::string name() { return "Advection"; }
  static std::string shortDescription() { return
    "Specify the advection physics configuration for a PDE "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the advection physics configuration for a
    PDE. Example: "transport physics advection end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using advection = keyword< advection_info, TAOCPP_PEGTL_STRING("advection") >;

struct advdiff_info {
  using code = Code< D >;
  static std::string name() { return "Advection + diffusion"; }
  static std::string shortDescription() { return
    "Specify the advection + diffusion physics configuration for a PDE "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the advection +diffusion physics
    configuration for a PDE. Example: "transport physics advdiff end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using advdiff = keyword< advdiff_info, TAOCPP_PEGTL_STRING("advdiff") >;

struct physics_info {
  using code = Code< p >;
  static std::string name() { return "Physics configuration"; }
  static std::string shortDescription() { return
    "Specify the physics configuration for a system of PDEs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the physics configuration for a particular
    PDE system. Example: "physics navierstokes", which selects the Navier-Stokes
    equations for solving viscous compressible flow, given within the
    compflow ... end block. Valid options depend on the given block the keyword
    is used.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + advection::string() + "\' | \'"
                  + advdiff::string() + "\' | \'"
                  + navierstokes::string() + "\' | \'"
                  + euler::string() + '\'';
    }
  };
};
using physics = keyword< physics_info, TAOCPP_PEGTL_STRING("physics") >;

struct pde_diffusivity_info {
  static std::string name() { return "diffusivity"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) diffusivity)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "diffusivity 5.0 2.0 3.0 end". The length of the vector depends on the
    particular type of PDE system and is controlled by the preceding keyword
    'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_diffusivity =
  keyword< pde_diffusivity_info, TAOCPP_PEGTL_STRING("diffusivity") >;

struct pde_lambda_info {
  static std::string name() { return "lambda"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) lambda)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "lambda 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of PDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_lambda = keyword< pde_lambda_info, TAOCPP_PEGTL_STRING("lambda") >;

struct pde_u0_info {
  static std::string name() { return "u0"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) u0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "u0 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of PDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_u0 = keyword< pde_u0_info, TAOCPP_PEGTL_STRING("u0") >;

struct pde_source_info {
  static std::string name() { return "source"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) source)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "source 5.0 2.0 3.0 0.1 end". The length of the vector depends on the
    particular
    type of PDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_source = keyword< pde_source_info, TAOCPP_PEGTL_STRING("source") >;

struct pde_alpha_info {
  static std::string name() { return "alpha"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) alpha)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "alpha 5.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_alpha = keyword< pde_alpha_info, TAOCPP_PEGTL_STRING("alpha") >;

struct pde_beta_info {
  static std::string name() { return "beta"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) beta)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "beta 5.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_beta = keyword< pde_beta_info, TAOCPP_PEGTL_STRING("beta") >;

struct pde_p0_info {
  static std::string name() { return "p0"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) p0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "p0 10.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_p0 = keyword< pde_p0_info, TAOCPP_PEGTL_STRING("p0") >;

// nonlinear energy parameters here
struct pde_betax_info {
  static std::string name() { return "betax"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) betax)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "betax 1.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_betax = keyword< pde_betax_info, TAOCPP_PEGTL_STRING("betax") >;

struct pde_betay_info {
  static std::string name() { return "betay"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) betay)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "betay 0.75".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_betay = keyword< pde_betay_info, TAOCPP_PEGTL_STRING("betay") >;

struct pde_betaz_info {
  static std::string name() { return "betaz"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) betaz)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "betaz 0.5".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_betaz = keyword< pde_betaz_info, TAOCPP_PEGTL_STRING("betaz") >;

struct pde_ce_info {
  static std::string name() { return "ce"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) ce)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to parameterize the
    Euler equations solving the manufactured solution test case "non-linear
    energy growth". Example: "ce -1.0". For more information on the test case see
    Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
    equations with relevance to Inertial Confinement Fusion", Journal of
    Computational Physics 267 (2014) 196-209.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_ce = keyword< pde_ce_info, TAOCPP_PEGTL_STRING("ce") >;

struct pde_kappa_info {
  static std::string name() { return "kappa"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) kappa)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "kappa 0.8")"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_kappa = keyword< pde_kappa_info, TAOCPP_PEGTL_STRING("kappa") >;

struct pde_r0_info {
  static std::string name() { return "r0"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) r0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to parameterize the
    Euler equations solving the manufactured solution test case "non-linear
    energy growth". Example: "r0 2.0". For more information on the test case see
    Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
    equations with relevance to Inertial Confinement Fusion", Journal of
    Computational Physics 267 (2014) 196-209.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_r0 = keyword< pde_r0_info, TAOCPP_PEGTL_STRING("r0") >;

struct cweight_info {
  static std::string name() { return "cweight"; }
  static std::string shortDescription() { return
    R"(Set value for central linear weight used by WENO, cweight)"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the central linear weight used for the
    central stencil in the Weighted Essentially Non-Oscillatory (WENO) limiter
    for discontinuous Galerkin (DG) methods. Example:
    "cweight 10.0".)"; }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 1.0;
    static constexpr type upper = 1000.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "real between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "]";
    }
  };
};
using cweight = keyword< cweight_info, TAOCPP_PEGTL_STRING("cweight") >;

struct sideset_info {
  static std::string name() { return "sideset"; }
  static std::string shortDescription() { return
    "Specify configuration for setting BC on a side set";
  }
  static std::string longDescription() { return
    R"(This keyword is used to specify boundary conditions on a side set for a
    solving partial differential equation.)";
  }
  struct expect {
    using type = int;
    static std::string description() { return "strings"; }
  };
};
using sideset = keyword< sideset_info, TAOCPP_PEGTL_STRING("sideset") >;

struct bc_dirichlet_info {
  static std::string name() { return "Dirichlet boundary condition"; }
  static std::string shortDescription() { return
    "Start configuration block describing Dirichlet boundary conditions"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an bc_dirichlet ... end block, used to
    specify the configuration for setting Dirichlet boundary conditions (BC) for
    a partial differential equation. This keyword is used to list multiple side
    sets on which a prescribed Dirichlet BC is then applied. Such prescribed BCs
    at each point in space and time are evaluated using a built-in function,
    e.g., using the method of manufactured solutions.
    Keywords allowed in a bc_dirichlet ... end block: )" + std::string("\'")
    + sideset::string() + "\'. "
    + R"(For an example bc_dirichlet ... end block, see
      doc/html/inicter_example_shear.html.)";
  }
};
using bc_dirichlet =
  keyword< bc_dirichlet_info, TAOCPP_PEGTL_STRING("bc_dirichlet") >;

struct bc_sym_info {
  static std::string name() { return "Symmetry boundary condition"; }
  static std::string shortDescription() { return
    "Start configuration block describing symmetry boundary conditions"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an bc_sym ... end block, used to
    specify the configuration for setting symmetry boundary conditions for a
    partial differential equation. Keywords allowed in a bc_sym ... end
    block: )" + std::string("\'")
    + sideset::string() + "\'. "
    + R"(For an example bc_sym ... end block, see
      doc/html/inicter_example_gausshump.html.)";
  }
};
using bc_sym =
  keyword< bc_sym_info, TAOCPP_PEGTL_STRING("bc_sym") >;

struct point_info {
  static std::string name() { return "point"; }
  static std::string shortDescription() { return "Specify a point"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a point, used, e.g., in specifying a
    point in 3D space for setting a stagnation (velocity vector = 0).  Example
    specification: 'point 0.0 0.1 0.2 end')";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "3 reals"; }
  };
};
using point = keyword< point_info, TAOCPP_PEGTL_STRING("point") >;

struct radius_info {
  static std::string name() { return "radius"; }
  static std::string shortDescription() { return "Specify a radius"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a radius, used, e.g., in specifying a
    point in 3D space for setting a stagnation (velocity vector = 0).  Example
    specification: 'radius 1.0e-5')";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using radius = keyword< radius_info, TAOCPP_PEGTL_STRING("radius") >;

struct bc_farfield_info {
  static std::string name() { return "Farfield boundary condition"; }
  static std::string shortDescription() { return
    "Start configuration block describing farfield boundary conditions"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a bc_farfield ... end block, used
    to specify the configuration for setting farfield boundary conditions
    for the compressible flow equations. Keywords allowed in a bc_farfield
    ... end block: )" + std::string("\'")
    + density::string() + "\', \'"
    + pressure::string() + "\', \'"
    + velocity::string() + "\', \'"
    + sideset::string() + "\'. ";
  }
};
using bc_farfield =
  keyword< bc_farfield_info, TAOCPP_PEGTL_STRING("bc_farfield") >;

struct bc_pressure_info {
  static std::string name() { return "Pressure boundary condition"; }
  static std::string shortDescription() { return
    "Start configuration block describing pressure boundary conditions"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a bc_pressure ... end block, used
    to specify the configuration for setting pressure boundary conditions
    for the compressible flow equations. Keywords allowed in a bc_pressure
    ... end block: )" + std::string("\'")
    + density::string() + "\', \'"
    + pressure::string() + "\', \'"
    + sideset::string() + "\'. ";
  }
};
using bc_pressure =
  keyword< bc_pressure_info, TAOCPP_PEGTL_STRING("bc_pressure") >;

struct id_info {
  static std::string name() { return "id"; }
  static std::string shortDescription() { return "ID"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify an ID, a positive integer.)";
  }
  struct expect {
    using type = uint64_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using id = keyword< id_info, TAOCPP_PEGTL_STRING("id") >;

struct prelax_info {
  static std::string name() { return "Pressure relaxation"; }
  static std::string shortDescription() { return
    "Turn multi-material finite pressure relaxation on/off"; }
  static std::string longDescription() { return
    R"(This keyword is used to turn finite pressure relaxation between multiple
       materials on/off. It is used only for the multi-material solver, and has
       no effect when used for the other PDE types.)";
  }
  struct expect {
    using type = int;
    static std::string description() { return "string"; }
    static std::string choices() { return "1 | 0"; }
  };
};
using prelax = keyword< prelax_info, TAOCPP_PEGTL_STRING("prelax") >;

struct mat_gamma_info {
  static std::string name() { return "gamma"; }
  static std::string shortDescription() { return "ratio of specific heats"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, ratio of specific
       heats.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_gamma = keyword< mat_gamma_info, TAOCPP_PEGTL_STRING("gamma") >;

struct mat_pstiff_info {
  static std::string name() { return "pstiff"; }
  static std::string shortDescription() { return "EoS stiffness parameter"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, stiffness
       parameter in the stiffened gas equation of state.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_pstiff = keyword< mat_pstiff_info, TAOCPP_PEGTL_STRING("pstiff") >;

struct mat_mu_info {
  static std::string name() { return "mu"; }
  static std::string shortDescription() { return "dynamic viscosity"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, dynamic
       viscosity.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_mu = keyword< mat_mu_info, TAOCPP_PEGTL_STRING("mu") >;

struct mat_cv_info {
  static std::string name() { return "cv"; }
  static std::string shortDescription() {
    return "specific heat at constant volume"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, specific heat at
       constant volume.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_cv = keyword< mat_cv_info, TAOCPP_PEGTL_STRING("cv") >;

struct mat_k_info {
  static std::string name() { return "k"; }
  static std::string shortDescription() { return "heat conductivity"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, heat
       conductivity.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_k = keyword< mat_k_info, TAOCPP_PEGTL_STRING("k") >;

struct stiffenedgas_info {
  static std::string name() { return "Stiffened gas"; }
  static std::string shortDescription() { return
    "Select the stiffened gas equation of state"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the stiffened gas equation of state.)"; }
};
using stiffenedgas =
  keyword< stiffenedgas_info, TAOCPP_PEGTL_STRING("stiffenedgas") >;

struct jwl_info {
  static std::string name() { return "JWL"; }
  static std::string shortDescription() { return
    "Select the JWL equation of state"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Jones, Wilkins, Lee equation of
    state.)"; }
};
using jwl = keyword< jwl_info, TAOCPP_PEGTL_STRING("jwl") >;

struct eos_info {
  static std::string name() { return "Equation of state"; }
  static std::string shortDescription() { return
    "Select equation of state (type)"; }
  static std::string longDescription() { return
    R"(This keyword is used to select an equation of state for a material.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + stiffenedgas::string() + "\' | \'"
                  + jwl::string() + '\'';
    }
  };
};
using eos = keyword< eos_info, TAOCPP_PEGTL_STRING("eos") >;

struct material_info {
  static std::string name() { return "Material properties block"; }
  static std::string shortDescription() { return
    "Start configuration block for material properties"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a material ... end block, used to
    specify material properties. Keywords allowed in a material ... end
    block: )" + std::string("\'")
    + id::string()+ "\', \'"
    + eos::string()+ "\', \'"
    + mat_gamma::string()+ "\', \'"
    + mat_pstiff::string()+ "\', \'"
    + mat_mu::string()+ "\', \'"
    + mat_cv::string()+ "\', \'"
    + mat_k::string() + "\'. "
    + R"(For an example material ... end block, see
      doc/html/inicter_example_compflow.html.)";
  }
};
using material = keyword< material_info, TAOCPP_PEGTL_STRING("material") >;

struct transport_info {
  static std::string name() { return "Transport"; }
  static std::string shortDescription() { return
    "Start configuration block for an transport equation"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an transport ... end block, used to
    specify the configuration for a transport equation type. Keywords allowed
    in an transport ... end block: )" + std::string("\'")
    + depvar::string() + "\', \'"
    + ncomp::string() + "\', \'"
    + problem::string() + "\', \'"
    + physics::string() + "\', \'"
    + pde_diffusivity::string() + "\', \'"
    + pde_lambda::string() + "\', \'"
    + bc_dirichlet::string() + "\', \'"
    + bc_sym::string() + "\', \'"
    + pde_u0::string() + "\'. "
    + pde_source::string() + "\'. "
    + R"(For an example transport ... end block, see
      doc/html/inicter_example_transport.html.)";
  }
};
using transport = keyword< transport_info, TAOCPP_PEGTL_STRING("transport") >;

struct compflow_info {
  static std::string name() { return "Compressible single-material flow"; }
  static std::string shortDescription() { return
    "Start configuration block for the compressible flow equations"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the compflow ... end block, used to
    specify the configuration for a system of partial differential equations,
    governing compressible fluid flow. Keywords allowed in an compflow ... end
    block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + physics::string() + "\', \'"
    + problem::string() + "\', \'"
    + material::string() + "\', \'"
    + pde_alpha::string() + "\', \'"
    + pde_p0::string() + "\', \'"
    + pde_betax::string() + "\', \'"
    + pde_betay::string() + "\', \'"
    + pde_betaz::string() + "\', \'"
    + pde_beta::string() + "\', \'"
    + pde_r0::string() + "\', \'"
    + pde_ce::string() + "\', \'"
    + pde_kappa::string() + "\', \'"
    + bc_dirichlet::string() + "\', \'"
    + bc_sym::string() + "\', \'"
    + bc_farfield::string() + "\', \'"
    + R"(For an example compflow ... end block, see
      doc/html/inicter_example_compflow.html.)";
  }
};
using compflow = keyword< compflow_info, TAOCPP_PEGTL_STRING("compflow") >;

struct rcb_info {
  static std::string name() { return "recursive coordinate bisection"; }
  static std::string shortDescription() { return
    "Select recursive coordinate bisection mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the recursive coordinate bisection (RCB)
    mesh partitioner. RCB is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.hpp for other valid options.)"; }
};
using rcb = keyword< rcb_info, TAOCPP_PEGTL_STRING("rcb") >;

struct rib_info {
  static std::string name() { return "recursive inertial bisection"; }
  static std::string shortDescription() { return
    "Select recursive inertial bisection mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the recursive inertial bisection (RIB)
    mesh partitioner. RIB is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.hpp for other valid options.)"; }
};
using rib = keyword< rib_info, TAOCPP_PEGTL_STRING("rib") >;

struct hsfc_info {
  static std::string name() { return "Hilbert space filling curve"; }
  static std::string shortDescription() { return
    "Select Hilbert Space Filling Curve (HSFC) mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Hilbert Space Filling Curve (HSFC)
    mesh partitioner. HSFC is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.hpp for other valid options.)"; }
};
using hsfc = keyword< hsfc_info, TAOCPP_PEGTL_STRING("hsfc") >;

struct phg_info {
  static std::string name() { return "hypergraph"; }
  static std::string shortDescription() { return
    "Select parallel hypergraph mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the parallel hypergraph (PHG)
    mesh partitioner. PHG is a graph-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.hpp for other valid options.)"; }
};
using phg = keyword< phg_info, TAOCPP_PEGTL_STRING("phg") >;

struct algorithm_info {
  static std::string name() { return "algorithm"; }
  static std::string shortDescription() { return
    "Select mesh partitioning algorithm"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a mesh partitioning algorithm. See
    Control/Options/PartitioningAlgorithm.hpp for valid options.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + rcb::string() + "\' | \'"
                  + rib::string() + "\' | \'"
                  + hsfc::string() + '\'';
    }
  };
};
using algorithm = keyword< algorithm_info, TAOCPP_PEGTL_STRING("algorithm") >;

struct partitioning_info {
  static std::string name() { return "partitioning"; }
  static std::string shortDescription() { return
    "Start configuration block for mesh partitioning"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a partitioning ... end block, used to
    specify the configuration for mesh partitioning. Keywords allowed
    in a partitioning ... end block: )" + std::string("\'")
    + algorithm::string() + "\'.";
  }
};
using partitioning = keyword< partitioning_info, TAOCPP_PEGTL_STRING("partitioning") >;

struct move_info {
  static std::string name() { return "move"; }
  static std::string shortDescription() { return
    "Start configuration block configuring surface movement"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a move ... end block, used to
    configure surface movement for ALE simulations. Keywords allowed
    in a move ... end block: )" + std::string("\'")
    + sideset::string() + "\'.";
  }
};
using move = keyword< move_info, TAOCPP_PEGTL_STRING("move") >;

struct amr_uniform_info {
  using code = Code< u >;
  static std::string name() { return "uniform refine"; }
  static std::string shortDescription() { return
    "Select uniform initial mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to select uniform initial mesh refinement.)"; }
};
using amr_uniform = keyword< amr_uniform_info, TAOCPP_PEGTL_STRING("uniform") >;

struct amr_uniform_deref_info {
  using code = Code< d >;
  static std::string name() { return "uniform derefine"; }
  static std::string shortDescription() { return
    "Select uniform initial mesh de-refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to select uniform initial mesh de-refinement.)"; }
};
using amr_uniform_deref =
  keyword< amr_uniform_deref_info, TAOCPP_PEGTL_STRING("uniform_derefine") >;

struct amr_initial_cond_info {
  using code = Code< i >;
  static std::string name() { return "initial conditions"; }
  static std::string shortDescription() { return
    "Select initial-conditions-based initial mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to select initial-conditions-based initial mesh
       refinement.)"; }
};
using amr_initial_cond =
  keyword< amr_initial_cond_info, TAOCPP_PEGTL_STRING("ic") >;

struct amr_edgelist_info {
  using code = Code< e >;
  static std::string name() { return "edge list"; }
  static std::string shortDescription() { return
    "Configure edge-node pairs for initial refinement"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a list of edges that are explicitly
    tagged for initial refinement during setup in inciter. The keyword
    introduces an edgelist ... end block within an amr ... end block and must
    contain a list of integer pairs, i.e., the number of ids must be even,
    denoting the end-points of the nodes (=edge) which should be tagged for
    refinement.)"; }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 0;
    static std::string description() { return "two ints"; }
  };
};
using amr_edgelist =
  keyword< amr_edgelist_info, TAOCPP_PEGTL_STRING("edgelist") >;

struct amr_coords_info {
  using code = Code< c >;
  static std::string name() { return "coordinates"; }
  static std::string shortDescription() { return
    "Configure initial refinement using coordinate planes"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure entire volumes on a given side of a
    plane in 3D space. The keyword introduces an coords ... end block within
    an amr ... end block and must contain the either or multiple of the
    following keywords: x- <real>, x+ <real>, y- <real>, y+ <real>, z- <real>,
    z+ <real>. All edges of the input mesh will be tagged for refinement whose
    end-points lie less than (-) or larger than (+) the real number given.
    Example: 'x- 0.5' refines all edges whose end-point coordinates are less
    than 0.5. Multiple specifications are understood by combining with a logical
    AND. That is: 'x- 0.5 y+ 0.3' refines all edges whose end-point x
    coordinates are less than 0.5 AND y coordinates are larger than 0.3.)"; }
};
using amr_coords =
  keyword< amr_coords_info, TAOCPP_PEGTL_STRING("coords") >;

struct amr_initial_info {
  static std::string name() { return "Initial refinement typelist"; }
  static std::string shortDescription() { return
    "Configure initial mesh refinement (before time stepping)"; }
  static std::string longDescription() { return
    R"(This keyword is used to add to a list of initial mesh refinement types
    that happens before t = 0. Example: initial uniform initial ic inital
    uniform, which yiedls an initial uniform refinement, followed by a
    refinement based on the numerical error computed based on the initial
    conditions, followed by another step of unfirom refinement.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + amr_uniform::string() + "\' | \'"
                  + amr_uniform_deref::string()  + "\' | \'"
                  + amr_initial_cond::string() + "\' | \'"
                  + amr_edgelist::string() + "\' | \'"
                  + amr_coords::string() + '\'';
    }
  };
};
using amr_initial = keyword< amr_initial_info, TAOCPP_PEGTL_STRING("initial") >;

struct amr_refvar_info {
  static std::string name() { return "refinement variable(s)"; }
  static std::string shortDescription() { return
    "Configure dependent variables used for adaptive mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to configured a list of dependent variables that
    trigger adaptive mesh refinement based on estimating their numerical error.
    These refinement variables are used for both initial (i.e., before time
    stepping) mesh refinement as well as during time stepping. Only previously
    (i.e., earlier in the input file) selected dependent variables can be
    configured as refinement variables. Dependent variables are required to be
    defined in all equation system configuration blocks, e.g., transport ...
    end, by using the 'depvar' keyword. Example: transport depvar c end amr
    refvar c end end. Selecting a particular scalar component in a system is
    done by appending the equation number to the refvar: Example: transport
    depvar q ncomp 3 end amr refvar q1 q2 end end, which configures two
    refinement variables: the first and third scalar component of the previously
    configured transport equation system.)"; }
  struct expect {
    static std::string description() { return "strings"; }
  };
};
using amr_refvar = keyword< amr_refvar_info, TAOCPP_PEGTL_STRING("refvar") >;

struct amr_xminus_info {
  static std::string name() { return "initial refinement: x-"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates lower than an x-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the x coordinate of a plane perpendicular to
    coordinate x in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'x- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'x- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_xminus =
  keyword< amr_xminus_info, TAOCPP_PEGTL_STRING("x-") >;

struct amr_xplus_info {
  static std::string name() { return "initial refinement: x+"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates larger than an x-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the x coordinate of a plane perpendicular
    to coordinate x in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'x+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'x+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_xplus =
  keyword< amr_xplus_info, TAOCPP_PEGTL_STRING("x+") >;

struct amr_yminus_info {
  static std::string name() { return "initial refinement: y-"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates lower than an y-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the y coordinate of a plane perpendicular to
    coordinate y in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'y- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'y- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_yminus =
  keyword< amr_yminus_info, TAOCPP_PEGTL_STRING("y-") >;

struct amr_yplus_info {
  static std::string name() { return "initial refinement: y+"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates larger than an y-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the y coordinate of a plane perpendicular
    to coordinate y in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'y+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'y+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_yplus =
  keyword< amr_yplus_info, TAOCPP_PEGTL_STRING("y+") >;

struct amr_zminus_info {
  static std::string name() { return "initial refinement: z-"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates lower than an z-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the z coordinate of a plane perpendicular to
    coordinate z in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'z- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'z- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_zminus =
  keyword< amr_zminus_info, TAOCPP_PEGTL_STRING("z-") >;

struct amr_zplus_info {
  static std::string name() { return "initial refinement: z+"; }
  static std::string shortDescription() { return "Configure initial refinement "
    "for coordinates larger than an z-normal plane"; }
  static std::string longDescription() { return
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the z coordinate of a plane perpendicular
    to coordinate z in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'z+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'z+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using amr_zplus =
  keyword< amr_zplus_info, TAOCPP_PEGTL_STRING("z+") >;

struct amr_jump_info {
  static std::string name() { return "jump"; }
  static std::string shortDescription() { return
    "Error estimation based on the jump in the solution normalized by solution";
  }
  static std::string longDescription() { return
    R"(This keyword is used to select the jump-based error indicator for
    solution-adaptive mesh refinement. The error is estimated by computing the
    magnitude of the jump in the solution value normalized by the solution
    value.)"; }
};
using amr_jump =
  keyword< amr_jump_info, TAOCPP_PEGTL_STRING("jump") >;

struct amr_hessian_info {
  static std::string name() { return "Hessian"; }
  static std::string shortDescription() { return
    "Error estimation based on the Hessian normalized by solution value"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Hessian-based error indicator for
    solution-adaptive mesh refinement. The error is estimated by computing the
    Hessian (2nd derivative matrix) of the solution normalized by sum of the
    absolute values of the gradients at edges-end points.)"; }
};
using amr_hessian = keyword< amr_hessian_info, TAOCPP_PEGTL_STRING("hessian") >;

struct amr_error_info {
  static std::string name() { return "Error estimator"; }
  static std::string shortDescription() { return
    "Configure the error type for solution-adaptive mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the algorithm used to estimate the error
    for solution-adaptive mesh refinement.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + amr_jump::string() + "\' | \'"
                  + amr_hessian::string() + '\'';
    }
  };
};
using amr_error = keyword< amr_error_info, TAOCPP_PEGTL_STRING("error") >;

struct amr_t0ref_info {
  static std::string name() { return "Mesh refinement at t<0"; }
  static std::string shortDescription() { return
    "Enable mesh refinement at t<0"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable initial mesh refinement, which can be
    configured to perform multiple levels of mesh refinement based on various
    refinement criteria and configuration settings.)";
  }
  struct expect {
    using type = bool;
    static std::string choices() { return "true | false"; }
    static std::string description() { return "string"; }
  };
};
using amr_t0ref = keyword< amr_t0ref_info, TAOCPP_PEGTL_STRING("t0ref") >;

struct amr_dtref_info {
  static std::string name() { return "Mesh refinement at t>0"; }
  static std::string shortDescription() { return
    "Enable mesh refinement at t>0"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable soution-adaptive mesh refinement during "
    "time stepping.)";
  }
  struct expect {
    using type = bool;
    static std::string choices() { return "true | false"; }
    static std::string description() { return "string"; }
  };
};
using amr_dtref = keyword< amr_dtref_info, TAOCPP_PEGTL_STRING("dtref") >;

struct amr_dtref_uniform_info {
  static std::string name() { return "Uniform-only mesh refinement at t>0"; }
  static std::string shortDescription() { return
    "Enable mesh refinement at t>0 but only perform uniform refinement"; }
  static std::string longDescription() { return R"(This keyword is used to force
    uniform-only soution-adaptive mesh refinement during time stepping.)";
  }
  struct expect {
    using type = bool;
    static std::string choices() { return "true | false"; }
    static std::string description() { return "string"; }
  };
};
using amr_dtref_uniform =
  keyword< amr_dtref_uniform_info, TAOCPP_PEGTL_STRING("dtref_uniform") >;

struct amr_dtfreq_info {
  static std::string name() { return "Mesh refinement frequency"; }
  static std::string shortDescription() { return
    "Set mesh refinement frequency during time stepping"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure the frequency of mesh refinement
    during time stepping. The default is 3, which means that mesh refinement
    will be performed every 3rd time step.)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max();
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using amr_dtfreq = keyword< amr_dtfreq_info, TAOCPP_PEGTL_STRING("dtfreq") >;

struct amr_maxlevels_info {
  static std::string name() { return "Maximum mesh refinement levels"; }
  static std::string shortDescription() { return
    "Set maximum allowed mesh refinement levels"; }
  static std::string longDescription() { return
    R"(This keyword is used to configure the maximum allowed mesh refinement
    levels. The default is 2.)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max();
    static std::string description() { return "uint"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using amr_maxlevels = keyword< amr_maxlevels_info,
  TAOCPP_PEGTL_STRING("maxlevels") >;

struct amr_tolref_info {
  static std::string name() { return "refine tolerance"; }
  static std::string shortDescription() { return "Configure refine tolerance"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the tolerance used to tag an edge for
    refinement if the relative error exceeds this value.)"; }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static constexpr type upper = 1.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using amr_tolref =
  keyword< amr_tolref_info, TAOCPP_PEGTL_STRING("tol_refine") >;

struct amr_tolderef_info {
  static std::string name() { return "derefine tolerance"; }
  static std::string shortDescription() {
    return "Configure derefine tolerance"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the tolerance used to tag an edge for
    derefinement if the relative error is below this value.)"; }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static constexpr type upper = 1.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using amr_tolderef =
  keyword< amr_tolderef_info, TAOCPP_PEGTL_STRING("tol_derefine") >;

struct amr_info {
  static std::string name() { return "AMR"; }
  static std::string shortDescription() { return
    "Start configuration block configuring adaptive mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the amr ... end block, used to
    configure adaptive mesh refinement. Keywords allowed
    in this block: )" + std::string("\'")
    + amr_t0ref::string() + "\' | \'"
    + amr_dtref::string() + "\' | \'"
    + amr_dtref_uniform::string() + "\' | \'"
    + amr_dtfreq::string() + "\' | \'"
    + amr_maxlevels::string() + "\' | \'"
    + amr_initial::string() + "\' | \'"
    + amr_refvar::string() + "\' | \'"
    + amr_tolref::string() + "\' | \'"
    + amr_tolderef::string() + "\' | \'"
    + amr_error::string() + "\' | \'"
    + amr_coords::string() + "\' | \'"
    + amr_edgelist::string() + "\'.";
  }
};
using amr = keyword< amr_info, TAOCPP_PEGTL_STRING("amr") >;

struct filename_info {
  static std::string name() { return "filename"; }
  static std::string shortDescription() { return "Set filename"; }
  static std::string longDescription() { return
    R"(Set filename, e.g., mesh filename for solver coupling.)";
  }
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using filename = keyword< filename_info, TAOCPP_PEGTL_STRING("filename") >;

struct location_info {
  static std::string name() { return "location"; }
  static std::string shortDescription() { return "Configure location"; }
  static std::string longDescription() { return
    R"(Configure location of a mesh relative to another, e.g., for solver
       coupling.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using location = keyword< location_info, TAOCPP_PEGTL_STRING("location") >;

struct orientation_info {
  static std::string name() { return "orientation"; }
  static std::string shortDescription() { return "Configure orientation"; }
  static std::string longDescription() { return
    R"(Configure orientation of a mesh relative to another, e.g., for solver
       coupling.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using orientation =
  keyword< orientation_info, TAOCPP_PEGTL_STRING("orientation") >;

struct mesh_info {
  static std::string name() { return "Mesh specification block"; }
  static std::string shortDescription() { return
    "Start configuration block assigning a mesh to a solver"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a mesh ... end block, used to
    assign and configure a mesh to a solver.)";
  }
};
using mesh = keyword< mesh_info, TAOCPP_PEGTL_STRING("mesh") >;

struct reference_info {
  static std::string name() { return "Mesh transformation"; }
  static std::string shortDescription() { return
    "Specify mesh transformation relative to a mesh of another solver"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a solver, given with a dependent
       variable, configured upstream in the input file, whose mesh is used as a
       reference to which the mesh being configured is transformed relative
       to.)";
  }
  struct expect {
    using type = char;
    static std::string description() { return "character"; }
  };
};
using reference = keyword< reference_info, TAOCPP_PEGTL_STRING("reference") >;

// This will go away once all the keywords below are documented
struct undefined_info {
  static std::string name() { return "undef"; }
  static std::string shortDescription() { return "undefined"; }
  static std::string longDescription() { return "Undefined."; }
};

} // kw::

#endif // Keywords_h
