// *****************************************************************************
/*!
  \file      src/Control/InciterConfig.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Lua parser for Inciter's control file
  \see       https://github.com/edubart/minilua
  \see       https://www.codingwiththomas.com/blog/a-lua-c-api-cheat-sheet
*/
// *****************************************************************************

#include "Compiler.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include "minilua.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#include "InciterConfig.hpp"

#include "NoWarning/charm++.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "TaggedTuple.hpp"
#include "PrintTaggedTupleDeep.hpp"
#include "Writer.hpp"

namespace inciter {
namespace ctr {

static constexpr auto largeint = std::numeric_limits< int64_t >::max();
static constexpr auto largeuint = std::numeric_limits< uint64_t >::max();
static constexpr auto largereal = std::numeric_limits< double >::max();

void
Config::cmdline( int argc, char** argv )
// *****************************************************************************
//! Contructor: parse inciter command line
//! \param[in] argc Number of arguments to executable
//! \param[in] argv Arguments to executable
// *****************************************************************************
{
  if (!argc) return;

  // Defaults
  if (!tk::git_commit().empty()) get< tag::commit >() = tk::git_commit();
  get< tag::output >() = "out";
  get< tag::diag >() = "diag";
  get< tag::checkpoint >() = "restart";
  get< tag::lbfreq >() = 1;
  get< tag::rsfreq >() = 1000;

  if (argc == 1) {
    help( argv );
    CkExit( EXIT_FAILURE );
  }
  tk::Print print;

  // Process command line arguments
  int c;
  while ((c = getopt( argc, argv, "bc:d:fh?i:l:no:qr:s:u:v" )) != -1) {
    switch (c) {
      case '?':
      case 'h':
      default:
        help( argv );
        CkExit();
        break;
      case 'b':
        get< tag::benchmark >() = true;
        break;
      case 'c':
        get< tag::control >() = optarg;
        break;
      case 'd':
        get< tag::diag >() = optarg;
        break;
      case 'f':
        get< tag::feedback >() = true;
        break;
      case 'i':
        get< tag::input >() = optarg;
        break;
      case 'l':
        get< tag::lbfreq >() = std::stoul( optarg );
        break;
      case 'n':
        get< tag::nonblocking >() = true;
        break;
      case 'o':
        get< tag::output >() = optarg;
        break;
      case 'r':
        get< tag::rsfreq >() = std::stoul( optarg );
        break;
      case 'q':
        get< tag::quiescence >() = true;
        break;
      case 'u':
        get< tag::virt >() = std::stod( optarg );
        break;
      case 'v':
        print << '\n';
        print.version( tk::inciter_executable(), tk::git_commit() );
        CkExit();
        break;
    }
  }

  if (optind != argc) {
    print << "\nA non-option was supplied";
    help( argv );
    CkExit( EXIT_FAILURE );
  }
  ErrChk( not get< tag::input >().empty(),
          "Mandatory input mesh file not specified. Use -i <filename>." );
  ErrChk( not get< tag::control >().empty(),
          "Mandatory control file not specified. Use -c <filename>." );
}

void
Config::help( char** argv )
// *****************************************************************************
// Echo help on command line arguments
//! \param[in] argv Arguments to executable
// *****************************************************************************
{
  tk::Print() <<
    "\nUsage: " << argv[0] << " -i <in.exo> -c <config.q> [OPTION]...\n"
    "\n"
    "  -h, -?        Print out this help\n"
    "  -b            Benchmark mode, "
                     "default: " << get< tag::benchmark >() << "\n" <<
    "  -c <config.q> Specify control file\n"
    "  -d <diag>     Specify diagnostics file, "
                     "default: " << get< tag::diag >() << "\n" <<
    "  -f            Extra feedback, "
                     "default: " << get< tag::feedback >() << "\n" <<
    "  -i <in.exo>   Specify input mesh file\n"
    "  -l <int>      Load balancing frequency, "
                     "default: " << get< tag::lbfreq >() << "\n" <<
    "  -n            Non-blocking migration, "
                     "default: " << get< tag::nonblocking >() << "\n" <<
    "  -o <outfile>  Base-filename for field output, "
                     "default: " << get< tag::output >() << "\n" <<
    "  -r <int>      Checkpoint frequency, "
                     "default: " << get< tag::rsfreq >() << "\n" <<
    "  -q            Enable quiescence detection, "
                     "default: " << get< tag::quiescence >() << "\n" <<
    "  -u <real>     Virtualization, "
                     "default: " << get< tag::virt >() << "\n" <<
    "  -v            Print revision information\n"
    "\n";
}

[[maybe_unused]] static void
dumpstack( lua_State *L )
// *****************************************************************************
// Dump lua stack for debugging
//! \param[in] L Lua state
// *****************************************************************************
{
  int top=lua_gettop(L);
  for (int i=1; i <= top; i++) {
    printf("%d\t%s\t", i, luaL_typename(L,i));
    switch (lua_type(L, i)) {
      case LUA_TNUMBER:
        printf("%g\n",lua_tonumber(L,i));
        break;
      case LUA_TSTRING:
        printf("%s\n",lua_tostring(L,i));
        break;
      case LUA_TBOOLEAN:
        printf("%s\n", (lua_toboolean(L, i) ? "true" : "false"));
        break;
      case LUA_TNIL:
        printf("%s\n", "nil");
        break;
      default:
        printf("%p\n",lua_topointer(L,i));
        break;
    }
  }
}

static int64_t
sigint( lua_State* L,
        const char* name,
        int64_t def = largeint,
        bool global = false )
// *****************************************************************************
// Parse integer from global scope
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  int64_t a = def;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) {
      lua_getfield( L, -1, name );
    } else {
      return a;
    }
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_isinteger( L, -1 ), std::string(name) + " must be an integer" );
    a = lua_tointeger( L, -1 );
  }

  lua_pop( L, 1 );

  return a;
}

static uint64_t
unsigint( lua_State* L,
          const char* name,
          uint64_t def = largeuint,
          bool global = false )
// *****************************************************************************
// Parse unsigned integer from global scope
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  uint64_t a = def;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) {
      lua_getfield( L, -1, name );
    } else {
      return a;
    }
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_isinteger( L, -1 ), std::string(name) + " must be an integer" );
    a = static_cast< uint64_t >( lua_tointeger( L, -1 ) );
  }

  lua_pop( L, 1 );

  return a;
}

static bool
boolean( lua_State* L, const char* name, bool def = false, bool global = false )
// *****************************************************************************
// Parse boolean from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  bool a = def;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) {
      lua_getfield( L, -1, name );
    } else {
      return a;
    }
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_isboolean( L, -1 ), std::string(name) + " must be a boolean" );
    a = lua_toboolean( L, -1 );
  }

  lua_pop( L, 1 );

  return a;
}

static double
real( lua_State* L,
      const char* name,
      double def = largereal,
      bool global = false )
// *****************************************************************************
// Parse real from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  double a = def;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, name ); else return a;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_isnumber( L, -1 ), std::string(name) + " must be a number" );
    a = lua_tonumber( L, -1 );
  }

  lua_pop( L, 1 );

  return a;
}

static std::string
string( lua_State* L,
        const char* name,
        const char* def = "default",
        bool global = false )
// *****************************************************************************
// Parse string from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  std::string a = def;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, name ); else return a;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_isstring( L, -1 ), std::string(name) + " must be a string" );
    a = lua_tostring( L, -1 );
  }

  lua_pop( L, 1 );

  return a;
}

static std::vector< double >
vector( lua_State* L,
        const char* name,
        double def = largereal,
        bool global = false )
// *****************************************************************************
// Parse vector table from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] def Default if does not exist
//! \param[in] global True to parse from global scope, false from table on stack
//! \return Vector components parsed
// *****************************************************************************
{
  std::vector< double > v( 3, def );

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, name ); else return v;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "vector must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_isnumber( L, -1 ), "vector components must be numbers" );
      v[ static_cast<std::size_t>(i-1) ] = lua_tonumber( L, -1 );
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );

  return v;
}

static std::vector< std::string >
stringlist( lua_State* L, const char* name, bool global = false )
// *****************************************************************************
// Parse string list table from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] global True to parse from global scope, false from table on stack
//! \return List of strings parsed
// *****************************************************************************
{
  std::vector< std::string > v;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, name ); else return v;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "stringlist must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_isstring( L, -1 ), "stringlist components must be strings" );
      v.push_back( lua_tostring( L, -1 ) );
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );

  return v;
}

static std::vector< uint64_t >
unsigints( lua_State* L, const char* name, bool global = false )
// *****************************************************************************
// Parse table of unsigned integers from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] global True to parse from global scope, false from table on stack
//! \return List of unsigned integers parsed
// *****************************************************************************
{
  std::vector< uint64_t > v;

  if (global) {
    lua_getglobal( L, name );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, name ); else return v;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "unsigints must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_isinteger( L, -1 ), "unsigints components must be numbers" );
      v.push_back( static_cast< uint64_t >( lua_tointeger( L, -1 ) ) );
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );

  return v;
}

static std::vector< int >
sideset( lua_State* L, bool global = false )
// *****************************************************************************
// Parse sideset table from global scope or table
//! \param[in] L Lua state
//! \param[in] global True to parse from global scope, false from table on stack
//! \return Vector of side set ids parsed
// *****************************************************************************
{
  std::vector< int > v;

  if (global) {
    lua_getglobal( L, "sideset" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "sideset" ); else return v;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "sideset must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_isinteger( L, -1 ), "sideset id must be a number" );
      int a = static_cast< int >( lua_tointeger( L, -1 ) );
      ErrChk( a >= 0, "sideset id must be non-negative" );
      v.push_back( a );
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );

  return v;
}

static std::vector< std::vector< double > >
range( lua_State* L, bool global = false )
// *****************************************************************************
// Parse range(s) table from global scope or table
//! \param[in] L Lua state
//! \param[in] global True to parse from global scope, false from table on stack
//! \return Vector of vectors of range(s)
// *****************************************************************************
{
  std::vector< std::vector< double > > v;

  if (global) {
    lua_getglobal( L, "range" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "range" ); else return v;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "range must be a table" );
    int64_t n = luaL_len( L, -1 );
    v.emplace_back();
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      if (lua_isnumber( L, -1 )) {
        v.back().push_back( lua_tonumber( L, -1 ) );
      } else {
        ErrChk( lua_istable( L, -1 ), "non-number range must be a table" );
        if (i>1) v.emplace_back();
        int64_t m = luaL_len( L, -1 );
        for (int64_t j=1; j<=m; ++j) {
          lua_geti( L, -1, j );
          ErrChk( lua_isnumber( L, -1 ), "vector components must be numbers" );
          v.back().push_back( lua_tonumber( L, -1 ) );
          lua_pop( L, 1 );
        }
      }
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );

  return v;
}

static void
fieldout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse fieldout table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "fieldout" );

  cfg.get< tag::fieldout_iter >() = unsigint( L, "iter" );
  cfg.get< tag::fieldout_time >() = real( L, "time" );
  cfg.get< tag::fieldout_range >() = range( L );
  cfg.get< tag::fieldout >() = sideset( L );

  lua_pop( L, 1 );
}

static void
histout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse histout table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "histout" );

  cfg.get< tag::histout_iter >() = unsigint( L, "iter" );
  cfg.get< tag::histout_time >() = real( L, "time" );
  cfg.get< tag::histout_range >() = range( L );
  cfg.get< tag::histout_precision >() = sigint( L, "precision", 8 );
  cfg.get< tag::histout_format >() = string( L, "format" );

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "points" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "histout points must be a table" );
      auto& r = cfg.get< tag::histout >();
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_istable( L, -1 ), "histout point must be a table" );
        r.emplace_back();
        auto& p = r.back();
        int64_t m = luaL_len( L, -1 );
        for (int64_t j=1; j<=m; ++j) {
          lua_geti( L, -1, j );
          ErrChk( lua_isnumber( L, -1 ), "point coordinate must be a number" );
          p.push_back( lua_tonumber( L, -1 ) );
          lua_pop( L, 1 );
        }
        lua_pop( L, 1 );
      }
    }
    lua_pop( L, 1 );
  }

  lua_pop( L, 1 );
}

static void
integout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse integout table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "integout" );

  cfg.get< tag::integout_iter >() = unsigint( L, "iter" );
  cfg.get< tag::integout_time >() = real( L, "time" );
  cfg.get< tag::integout_range >() = range( L );
  cfg.get< tag::integout_precision >() = sigint( L, "precision", 8 );
  cfg.get< tag::integout_format >() = string( L, "format" );
  cfg.get< tag::integout >() = sideset( L );

  lua_pop( L, 1 );
}

static void
diag( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse diag table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "diag" );

  cfg.get< tag::diag_iter >() = unsigint( L, "iter", 1 );
  cfg.get< tag::diag_precision >() = sigint( L, "precision", 8 );
  cfg.get< tag::diag_format >() = string( L, "format" );

  lua_pop( L, 1 );
}

static void
bc_dir( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_dir table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "bc_dir" );

  if (lua_istable( L, -1 )) {
    auto& s = cfg.get< tag::bc_dir >();
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_istable( L, -1 ), "bc_dir must be a table" );
      s.emplace_back();
      auto& b = s.back();
      int64_t m = luaL_len( L, -1 );
      for (int64_t j=1; j<=m; ++j) {
        lua_geti( L, -1, j );
        ErrChk( lua_isinteger( L, -1 ), "bc_dir entry must be an integer" );
        b.push_back( static_cast< int >( lua_tointeger( L, -1 ) ) );
        lua_pop( L, 1 );
      }
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );
}

static void
bc_sym( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_sym table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "bc_sym" );

  cfg.get< tag::bc_sym >() = sideset( L );

  lua_pop( L, 1 );
}

static void
bc_far( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_far table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "bc_far" );

  cfg.get< tag::bc_far >() = sideset( L );
  cfg.get< tag::bc_far_velocity >() = vector( L, "velocity" );
  cfg.get< tag::bc_far_density >() = real( L, "density" );
  cfg.get< tag::bc_far_pressure >() = real( L, "pressure" );

  lua_pop( L, 1 );
}

static void
bc_pre( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_pre table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "bc_pre" );

  if (lua_istable( L, -1 )) {
    auto bc_pre_set = [&]( const char* name ) {
      lua_getfield( L, -1, name );
      if (!lua_isnil( L, -1 )) {
        ErrChk( lua_istable( L, -1 ), "bc_pre must be a table" );
        cfg.get< tag::bc_pre >().push_back( sideset( L ) );
        cfg.get< tag::bc_pre_density >().push_back( real(L,"density") );
        cfg.get< tag::bc_pre_pressure >().push_back( real(L,"pressure") );
      }
      lua_pop( L, 1 );
    };

    bc_pre_set( "inlet" );
    bc_pre_set( "outlet" );
  }

  lua_pop( L, 1 );
}

static void
ic( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse ic table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "ic" );

  cfg.get< tag::ic_velocity >() = vector( L, "velocity" );
  cfg.get< tag::ic_density >() = real( L, "density" );
  cfg.get< tag::ic_pressure >() = real( L, "pressure" );
  cfg.get< tag::ic_energy >() = real( L, "energy" );
  cfg.get< tag::ic_temperature >() = real( L, "temperature" );

  auto box_extent = [&]( const char* axis, auto& v ) {
    lua_getfield( L, -1, axis );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "ic box extents must be a table" );
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_isnumber( L, -1 ), "ic extent must be a number" );
        v.push_back( lua_tonumber( L, -1 ) );
        lua_pop( L, 1 );
      }
    }
    lua_pop( L, 1 );
  };

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "boxes" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "ic boxes must be a table" );
      auto& boxes = cfg.get< tag::ic >();
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_istable( L, -1 ), "ic box must be a table" );
        boxes.emplace_back();
        auto& box = boxes.back();
        box_extent( "x", box.get< tag::x >() );
        box_extent( "y", box.get< tag::y >() );
        box_extent( "z", box.get< tag::z >() );
        box.get< tag::ic_velocity >() = vector( L, "velocity" );
        box.get< tag::ic_density >() = real( L, "density" );
        box.get< tag::ic_pressure >() = real( L, "pressure" );
        box.get< tag::ic_energy >() = real( L, "energy" );
        box.get< tag::ic_temperature >() = real( L, "temperature" );
        lua_pop( L, 1 );
      }
    }
    lua_pop( L, 1 );
  }

  lua_pop( L, 1 );
}

static void
mat( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse mat table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "mat" );

  cfg.get< tag::mat_spec_heat_ratio >() = real( L, "spec_heat_ratio" );
  cfg.get< tag::mat_spec_heat_const_vol >() = real( L, "spec_heat_const_vol" );
  cfg.get< tag::mat_spec_gas_const >() = real(L, "spec_gas_const", 287.052874);
  cfg.get< tag::mat_heat_conductivity >() = real( L, "heat_conductivity" );

  lua_pop( L, 1 );
}

static void
problem( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse problem table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "problem" );

  cfg.get< tag::problem >() = string( L, "name", "userdef" );
  cfg.get< tag::problem_alpha >() = real( L, "alpha" );
  cfg.get< tag::problem_kappa >() = real( L, "kappa" );
  cfg.get< tag::problem_r0 >() = real( L, "r0" );
  cfg.get< tag::problem_p0 >() = real( L, "p0" );
  cfg.get< tag::problem_ce >() = real( L, "ce" );
  cfg.get< tag::problem_beta >() = vector( L, "beta" );

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "src" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "problem source must be a table" );
      auto& s = cfg.get< tag::problem_src >();
      s.get< tag::location >() = vector( L, "location" );
      s.get< tag::radius >() = real( L, "radius" );
      s.get< tag::release_time >() = real( L, "release_time" );
    }
    lua_pop( L, 1 );
  }

  const auto& solver = cfg.get< tag::solver >();
  const auto& problem = cfg.get< tag::problem >();

  auto& n = cfg.get< tag::problem_ncomp >();
  n = 5;

  if (problem == "slot_cyl" || problem == "point_src") ++n;

  if (solver == "chocg") {
    //if (problem.find("poisson") != std::string::npos)
    n = 1;
  }

  lua_pop( L, 1 );
}

static void
href( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse href table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "href" );

  cfg.get< tag::href_t0 >() = boolean( L, "t0" );
  cfg.get< tag::href_dt >() = boolean( L, "dt" );
  cfg.get< tag::href_dtfreq >() = unsigint( L, "dtfreq", 5 );
  cfg.get< tag::href_init >() = stringlist( L, "init" );
  cfg.get< tag::href_refvar >() = unsigints( L, "refvar" );
  cfg.get< tag::href_error >() = string( L, "error", "jump" );
  cfg.get< tag::href_maxlevels >() = unsigint( L, "maxlevels", 2 );

  lua_pop( L, 1 );
}

static void
deactivate( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse deactivate table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "deactivate" );

  cfg.get< tag::deafreq >() = unsigint( L, "freq", 0 );
  cfg.get< tag::deatol >() = real( L, "tol", 1.0e-3 );
  cfg.get< tag::deadif >() = real( L, "dif", 0.0 );
  cfg.get< tag::deasys >() = unsigints( L, "sys" );
  cfg.get< tag::deatime >() = real( L, "time", 0.0 );

  cfg.get< tag::deactivate >() = cfg.get< tag::deafreq >();  // on if freq > 0

  lua_pop( L, 1 );
}

static void
pressure( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse pressure table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "pressure" );

  cfg.get< tag::pre_iter >() = unsigint( L, "iter", 10 );
  cfg.get< tag::pre_tol >() = real( L, "tol", 1.0e-3 );
  cfg.get< tag::pre_verbose >() = unsigint( L, "verbose", 0 );

  lua_pop( L, 1 );
}

static void
lb( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse lb (load balancing configuration) table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "lb" );

  cfg.get< tag::lbtime >() = real( L, "time", 0.0 );

  lua_pop( L, 1 );
}

void
Config::control()
// *****************************************************************************
// Parse control file
// *****************************************************************************
{
  const auto& controlfile = get< tag::control >();

  tk::Print print;
  print.section( "Control file: " + controlfile );

  lua_State* L = luaL_newstate();
  luaL_openlibs( L );

  std::string err;

  if (luaL_dofile( L, controlfile.c_str() ) == LUA_OK) {

    get< tag::nstep >() = unsigint( L, "nstep", largeuint, true );
    get< tag::term >() = real( L, "term", largereal, true );
    get< tag::ttyi >() = unsigint( L, "ttyi", 1, true );
    get< tag::cfl >() = real( L, "cfl", 0.0, true );
    get< tag::dt >() = real( L, "dt", 0.0, true );
    get< tag::turkel >() = real( L, "turkel", 0.5, true );
    get< tag::velinf >() = vector( L, "velinf", 1.0, true );
    get< tag::t0 >() = real( L, "t0", 0.0, true );
    get< tag::reorder >() = boolean( L, "reorder", false, true );
    get< tag::flux >() = string( L, "flux", "rusanov", true );
    get< tag::steady >() = boolean( L, "steady", false, true );
    get< tag::residual >() = real( L, "residual", 0.0, true );
    get< tag::rescomp >() = unsigint( L, "rescomp", 1, true );
    get< tag::part >() = string( L, "part", "rcb", true );
    get< tag::zoltan_params >() = stringlist( L, "zoltan_params", true );
    get< tag::solver >() = string( L, "solver", "riecg", true );
    get< tag::stab2 >() = boolean( L, "stab2", false, true );
    get< tag::stab2coef >() = real( L, "stab2coef", 0.2, true );
    get< tag::fct >() = boolean( L, "fct", true, true );
    get< tag::fctdif >() = real( L, "fctdif", 1.0, true );
    get< tag::fctclip >() = boolean( L, "fctclip", false, true );
    get< tag::fctsys >() = unsigints( L, "fctsys", true );
    get< tag::fctfreeze >() = real( L, "fctfreeze", -largereal, true );
    get< tag::freezeflow >() = real( L, "freezeflow", 1.0, true );
    get< tag::freezetime >() = real( L, "freezetime", 0.0, true );

    ic( L, *this );
    bc_dir( L, *this );
    bc_sym( L, *this );
    bc_far( L, *this );
    bc_pre( L, *this );
    problem( L, *this );
    mat( L, *this );
    fieldout( L, *this );
    histout( L, *this );
    integout( L, *this );
    diag( L, *this );
    href( L, *this );
    deactivate( L, *this );
    pressure( L, *this );
    lb( L, *this );

    print << "Solver: " << get< tag::solver >() << '\n';

  } else {
    err = lua_tostring( L, -1 );
  }

  lua_close( L );

  ErrChk( err.empty(), "Lua error during control file parsing: " + err );

  // Output command line object to file
  auto logfilename = tk::inciter_executable() + "_control.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), *this );
  print << "Control data saved in " << logfilename << '\n';
}

} // ctr::
} // inciter::
