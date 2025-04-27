// *****************************************************************************
/*!
  \file      src/Control/InciterConfig.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
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

extern int g_nrestart;

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
        if (!g_nrestart) get< tag::input >().push_back( optarg );
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
        if (!g_nrestart) get< tag::virt >().push_back( std::stod( optarg ) );
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

  // Augment virtualization parameters if necessary
  auto& vir = get< tag::virt>();
  auto inps = get< tag::input >().size();
  if (inps > vir.size()) vir.resize( inps, 0.0 );

  // Basic error handling
  ErrChk( not get< tag::input >().empty(),
          "Mandatory input mesh file not specified. Use -i <filename>." );
  ErrChk( get< tag::input >().size() <= 2,
          "The maximum number of meshes for coupled problems is 2. If you "
          "need more, put them into the same mesh file." );
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
    "  -i <in.exo>   Specify an input mesh file. Use it another time to "
                     "specify a second mesh file for coupled problems.\n"
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
    "  -u <real>     Virtualization, default: 0.0\n" <<
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
intergrid_( lua_State* L, Config& cfg,
            std::vector< std::vector< int > >& sets,
            std::vector< std::vector< uint64_t > >& lays,
            std::vector< std::string >& syms )
// *****************************************************************************
// Parse integrid_* table from table for multple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
//! \param[in,out] sets State to push back intergrid setids to (outer vec: mesh)
//! \param[in,out] lays State to push back intergrid layers to (outer vec: mesh)
//! \param[in,out] syms State to push back holes-ymmetries to (vec: mesh)
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  sets.resize( nf );
  lays.resize( nf );
  syms.resize( nf );
  std::string basename = "intergrid_";

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    if (lua_istable(L, -1)) lua_getfield( L, -1, name.c_str() ); else return;

    if (!lua_isnil( L, -1 )) {
      sets[k] = sideset( L );
      lays[k] = unsigints( L, "layers" );
      syms[k] = string( L, "sym", "" );
    }

    lua_pop( L, 1 );

  }
}

static void
overset( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse overset table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "overset" );

  auto& tf = cfg.get< tag::overset>();
  intergrid_( L, cfg,
              tf.get< tag::intergrid_ >(),
              tf.get< tag::layers_ >(),
              tf.get< tag::sym_ >() );

  lua_pop( L, 1 );
}

static void
fieldout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse fieldout table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "fieldout" );

  auto& tf = cfg.get< tag::fieldout >();
  tf.get< tag::sidesets >() = sideset( L );
  tf.get< tag::iter >() = unsigint( L, "iter" );
  tf.get< tag::time >() = real( L, "time" );
  tf.get< tag::range >() = range( L );

  lua_pop( L, 1 );
}

static void
fieldout_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse fieldout_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "fieldout_";
  auto& tf = cfg.get< tag::fieldout_ >();
  tf.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& tfk = tf[k];
    tfk.get< tag::sidesets >() = sideset( L );
    tfk.get< tag::iter >() = unsigint( L, "iter" );
    tfk.get< tag::time >() = real( L, "time" );
    tfk.get< tag::range >() = range( L );

    lua_pop( L, 1 );

  }
}

static void
histout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse histout table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "histout" );

  auto& th = cfg.get< tag::histout >();
  th.get< tag::iter >() = unsigint( L, "iter" );
  th.get< tag::time >() = real( L, "time" );
  th.get< tag::range >() = range( L );
  th.get< tag::precision >() = sigint( L, "precision", 8 );
  th.get< tag::format >() = string( L, "format" );

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "points" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "histout points must be a table" );
      auto& r = th.get< tag::points >();
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
histout_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse histout_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "histout_";
  auto& th = cfg.get< tag::histout_ >();
  th.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& thk = th[k];

    thk.get< tag::iter >() = unsigint( L, "iter" );
    thk.get< tag::time >() = real( L, "time" );
    thk.get< tag::range >() = range( L );
    thk.get< tag::precision >() = sigint( L, "precision", 8 );
    thk.get< tag::format >() = string( L, "format" );

    if (lua_istable( L, -1 )) {
      lua_getfield( L, -1, "points" );
      if (!lua_isnil( L, -1 )) {
        ErrChk( lua_istable( L, -1 ), "histout points must be a table" );
        auto& r = thk.get< tag::points >();
        int64_t n = luaL_len( L, -1 );
        for (int64_t i=1; i<=n; ++i) {
          lua_geti( L, -1, i );
          ErrChk( lua_istable( L, -1 ), "histout point must be a table" );
          r.emplace_back();
          auto& p = r.back();
          int64_t m = luaL_len( L, -1 );
          for (int64_t j=1; j<=m; ++j) {
            lua_geti( L, -1, j );
            ErrChk( lua_isnumber(L,-1), "point coordinate must be a number" );
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
}

static void
integout( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse integout table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "integout" );

  auto& ti = cfg.get< tag::integout >();
  ti.get< tag::sidesets >() = sideset( L );
  ti.get< tag::integrals >() = stringlist( L, "integrals" );
  ti.get< tag::iter >() = unsigint( L, "iter" );
  ti.get< tag::time >() = real( L, "time" );
  ti.get< tag::range >() = range( L );
  ti.get< tag::precision >() = sigint( L, "precision", 8 );
  ti.get< tag::format >() = string( L, "format" );

  lua_pop( L, 1 );
}

static void
integout_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse integout_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "integout_";
  auto& ti = cfg.get< tag::integout_ >();
  ti.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& tik = ti[k];

    tik.get< tag::sidesets >() = sideset( L );
    tik.get< tag::integrals >() = stringlist( L, "integrals" );
    tik.get< tag::iter >() = unsigint( L, "iter" );
    tik.get< tag::time >() = real( L, "time" );
    tik.get< tag::range >() = range( L );
    tik.get< tag::precision >() = sigint( L, "precision", 8 );
    tik.get< tag::format >() = string( L, "format" );

    lua_pop( L, 1 );

  }
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
bc_dir( lua_State* L,
        std::vector< std::vector< int > >& mask,
        bool global = false )
// *****************************************************************************
// Parse bc_dir table
//! \param[in,out] L Lua state
//! \param[in] global True to parse from global scope, false from table on stack
//! \param[in,out] mask Config state to store Dirichlet BC setids and mask
// *****************************************************************************
{
  if (global) {
    lua_getglobal( L, "bc_dir" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "bc_dir" ); else return;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "bc_dir must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_istable( L, -1 ), "bc_dir table entry must be a table" );
      mask.emplace_back();
      auto& b = mask.back();
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
bc_dirval( lua_State* L,
           std::vector< std::vector< double > >& val,
           bool global = false )
// *****************************************************************************
// Parse bc_dirval table
//! \param[in,out] L Lua state
//! \param[in] global True to parse from global scope, false from table on stack
//! \param[in,out] val Config state to store Dirichlet BC setids and values
// *****************************************************************************
{
  if (global) {
    lua_getglobal( L, "bc_dirval" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "bc_dirval" ); else return;
  }

  if (!lua_isnil( L, -1 )) {
    ErrChk( lua_istable( L, -1 ), "bc_dirval must be a table" );
    int64_t n = luaL_len( L, -1 );
    for (int64_t i=1; i<=n; ++i) {
      lua_geti( L, -1, i );
      ErrChk( lua_istable( L, -1 ), "bc_dirval table entry must be a table" );
      val.emplace_back();
      auto& b = val.back();
      int64_t m = luaL_len( L, -1 );
      for (int64_t j=1; j<=m; ++j) {
        lua_geti( L, -1, j );
        ErrChk( lua_isnumber( L, -1 ), "bc_dirval entry must be an real" );
        b.push_back( static_cast< double >( lua_tonumber( L, -1 ) ) );
        lua_pop( L, 1 );
      }
      lua_pop( L, 1 );
    }
  }

  lua_pop( L, 1 );
}

static void
bc_dir_( lua_State* L,
         std::vector< std::vector< std::vector< int > > >& mask,
         std::size_t nf,
         bool global = false )
// *****************************************************************************
// Parse bc_dir_* table from global scope or table for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] mask State to push back Dirichlet BC setids and mask to
//!                (outer vec: mesh)
//! \param[in] nf Number of mesh files specified on command line
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  if (nf == 1) return;

  mask.resize( nf );
  std::string basename = "bc_dir_";

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);

    if (global) {
      lua_getglobal( L, name.c_str() );
    } else {
      if (lua_istable(L, -1)) lua_getfield( L, -1, name.c_str() ); else return;
    }

    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), name + " must be a table" );
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_istable( L, -1 ), name + " table entry must be a table" );
        mask[k].emplace_back();
        auto& b = mask[k].back();
        int64_t m = luaL_len( L, -1 );
        for (int64_t j=1; j<=m; ++j) {
          lua_geti( L, -1, j );
          ErrChk( lua_isinteger( L, -1 ), name + " entry must be an integer" );
          b.push_back( static_cast< int >( lua_tointeger( L, -1 ) ) );
          lua_pop( L, 1 );
        }
        lua_pop( L, 1 );
      }
    }

    lua_pop( L, 1 );

  }
}

static void
bc_dirval_( lua_State* L,
            std::vector< std::vector< std::vector< double > > >& val,
            std::size_t nf,
            bool global = false )
// *****************************************************************************
// Parse bc_dirval table
//! \param[in,out] L Lua state
//! \param[in,out] val Config state to store Dirichlet BC setids and values
//! \param[in] nf Number of mesh files specified on command line
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  if (nf == 1) return;

  val.resize( nf );
  std::string basename = "bc_dirval_";

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);

    if (global) {
      lua_getglobal( L, name.c_str() );
    } else {
      if (lua_istable(L, -1)) lua_getfield( L, -1, name.c_str() ); else return;
    }

    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "bc_dirval must be a table" );
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_istable( L, -1 ), "bc_dirval table entry must be a table" );
        val[k].emplace_back();
        auto& b = val[k].back();
        int64_t m = luaL_len( L, -1 );
        for (int64_t j=1; j<=m; ++j) {
          lua_geti( L, -1, j );
          ErrChk( lua_isnumber( L, -1 ), "bc_dirval entry must be an real" );
          b.push_back( static_cast< double >( lua_tonumber( L, -1 ) ) );
          lua_pop( L, 1 );
        }
        lua_pop( L, 1 );
      }
    }

    lua_pop( L, 1 );

  }
}

static void
bc_sym( lua_State* L, std::vector< int >& s, bool global = false )
// *****************************************************************************
// Parse bc_sym table from global scope or table
//! \param[in,out] L Lua state
//! \param[in] global True to parse from global scope, false from table on stack
//! \param[in,out] s Config state
// *****************************************************************************
{
  if (global) {
    lua_getglobal( L, "bc_sym" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "bc_sym" ); else return;
  }

  if (!lua_isnil( L, -1 )) {
    s = sideset( L );
  }

  lua_pop( L, 1 );
}

static void
bc_sym_( lua_State* L,
         std::vector< std::vector< int > >& s,
         std::size_t nf,
         bool global = false )
// *****************************************************************************
// Parse bc_sym table from global scope or table for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] s State to push back symmetry BC to (outer vec: mesh)
//! \param[in] nf Number of mesh files specified on command line
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  if (nf == 1) return;

  s.resize( nf );
  std::string basename = "bc_sym_";

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);

    if (global) {
      lua_getglobal( L, name.c_str() );
    } else {
      if (lua_istable(L, -1)) lua_getfield( L, -1, name.c_str() ); else return;
    }

    if (!lua_isnil( L, -1 )) {
      s[k] = sideset( L );
    }

    lua_pop( L, 1 );

  }
}

static void
bc_noslip( lua_State* L, std::vector< int >& s, bool global = false )
// *****************************************************************************
// Parse bc_noslip table from global scope or table
//! \param[in,out] L Lua state
//! \param[in,out] s State to push back no-slip BC setids to
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  if (global) {
    lua_getglobal( L, "bc_noslip" );
  } else {
    if (lua_istable( L, -1 )) lua_getfield( L, -1, "bc_noslip" ); else return;
  }

  if (!lua_isnil( L, -1 )) {
    s = sideset( L );
  }

  lua_pop( L, 1 );
}

static void
bc_noslip_( lua_State* L,
            std::vector< std::vector< int > >& s,
            std::size_t nf,
            bool global = false )
// *****************************************************************************
// Parse bc_noslip_* table from global scope or table for multple meshes
//! \param[in,out] L Lua state
//! \param[in,out] s State to push back no-slip BC setids to (outer vec: mesh)
//! \param[in] nf Number of mesh files specified on command line
//! \param[in] global True to parse from global scope, false from table on stack
// *****************************************************************************
{
  if (nf == 1) return;

  s.resize( nf );
  std::string basename = "bc_noslip_";

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);

    if (global) {
      lua_getglobal( L, name.c_str() );
    } else {
      if (lua_istable(L, -1)) lua_getfield( L, -1, name.c_str() ); else return;
    }

    if (!lua_isnil( L, -1 )) {
      s[k] = sideset( L );
    }

    lua_pop( L, 1 );

  }
}

static void
bc_far( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_far table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "bc_far" );

  auto& tf = cfg.get< tag::bc_far >();
  tf.get< tag::sidesets >() = sideset( L );
  tf.get< tag::velocity >() = vector( L, "velocity" );
  tf.get< tag::density >() = real( L, "density" );
  tf.get< tag::pressure >() = real( L, "pressure" );

  lua_pop( L, 1 );
}

static void
bc_far_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_far_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "bc_far_";
  auto& tf = cfg.get< tag::bc_far_ >();
  tf.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& tfk = tf[k];
    tfk.get< tag::sidesets >() = sideset( L );
    tfk.get< tag::velocity >() = vector( L, "velocity" );
    tfk.get< tag::density >() = real( L, "density" );
    tfk.get< tag::pressure >() = real( L, "pressure" );

    lua_pop( L, 1 );

  }
}

static void
bc_pre( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_pre table from global scope
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
        auto& tb = cfg.get< tag::bc_pre >();
        tb.get< tag::sidesets >().push_back( sideset( L ) );
        tb.get< tag::density >().push_back( real(L,"density") );
        tb.get< tag::pressure >().push_back( real(L,"pressure") );
      }
      lua_pop( L, 1 );
    };

    bc_pre_set( "inlet" );
    bc_pre_set( "outlet" );
  }

  lua_pop( L, 1 );
}

static void
bc_pre_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse bc_pre_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "bc_pre_";
  auto& tb = cfg.get< tag::bc_pre_ >();
  tb.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );
    auto& tbk = tb[k];

    if (lua_istable( L, -1 )) {
      auto bc_pre_set = [&]( const char* n ) {
        lua_getfield( L, -1, n );
        if (!lua_isnil( L, -1 )) {
          ErrChk( lua_istable( L, -1 ), "bc_pre must be a table" );
          tbk.get< tag::sidesets >().push_back( sideset( L ) );
          tbk.get< tag::density >().push_back( real(L,"density") );
          tbk.get< tag::pressure >().push_back( real(L,"pressure") );
        }
        lua_pop( L, 1 );
      };

      bc_pre_set( "inlet" );
      bc_pre_set( "outlet" );
    }

    lua_pop( L, 1 );

  }
}

static void
ic( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse ic table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "ic" );

  auto& tic = cfg.get< tag::ic >();
  tic.get< tag::velocity >() = vector( L, "velocity" );
  tic.get< tag::density >() = real( L, "density" );
  tic.get< tag::pressure >() = real( L, "pressure" );
  tic.get< tag::energy >() = real( L, "energy" );
  tic.get< tag::temperature >() = real( L, "temperature" );

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
      auto& boxes = tic.get< tag::boxes >();
      int64_t n = luaL_len( L, -1 );
      for (int64_t i=1; i<=n; ++i) {
        lua_geti( L, -1, i );
        ErrChk( lua_istable( L, -1 ), "ic box must be a table" );
        boxes.emplace_back();
        auto& box = boxes.back();
        box_extent( "x", box.get< tag::box_x >() );
        box_extent( "y", box.get< tag::box_y >() );
        box_extent( "z", box.get< tag::box_z >() );
        box.get< tag::box_velocity >() = vector( L, "velocity" );
        box.get< tag::box_density >() = real( L, "density" );
        box.get< tag::box_pressure >() = real( L, "pressure" );
        box.get< tag::box_energy >() = real( L, "energy" );
        box.get< tag::box_temperature >() = real( L, "temperature" );
        lua_pop( L, 1 );
      }
    }
    lua_pop( L, 1 );
  }

  lua_pop( L, 1 );
}

static void
ic_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse ic table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "ic_";
  auto& tic = cfg.get< tag::ic_ >();
  tic.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& tick = tic[k];
    tick.get< tag::velocity >() = vector( L, "velocity" );
    tick.get< tag::density >() = real( L, "density" );
    tick.get< tag::pressure >() = real( L, "pressure" );
    tick.get< tag::energy >() = real( L, "energy" );
    tick.get< tag::temperature >() = real( L, "temperature" );

    auto box_extent = [&]( const char* axis, auto& v ) {
      lua_getfield( L, -1, axis );
      if (!lua_isnil( L, -1 )) {
        ErrChk( lua_istable( L, -1 ), name + " box extents must be a table" );
        int64_t n = luaL_len( L, -1 );
        for (int64_t i=1; i<=n; ++i) {
          lua_geti( L, -1, i );
          ErrChk( lua_isnumber( L, -1 ), name + " extent must be a number" );
          v.push_back( lua_tonumber( L, -1 ) );
          lua_pop( L, 1 );
        }
      }
      lua_pop( L, 1 );
    };

    if (lua_istable( L, -1 )) {
      lua_getfield( L, -1, "boxes" );
      if (!lua_isnil( L, -1 )) {
        ErrChk( lua_istable( L, -1 ), name + " boxes must be a table" );
        auto& boxes = tick.get< tag::boxes >();
        int64_t n = luaL_len( L, -1 );
        for (int64_t i=1; i<=n; ++i) {
          lua_geti( L, -1, i );
          ErrChk( lua_istable( L, -1 ), name + " box must be a table" );
          boxes.emplace_back();
          auto& box = boxes.back();
          box_extent( "x", box.get< tag::box_x >() );
          box_extent( "y", box.get< tag::box_y >() );
          box_extent( "z", box.get< tag::box_z >() );
          box.get< tag::box_velocity >() = vector( L, "velocity" );
          box.get< tag::box_density >() = real( L, "density" );
          box.get< tag::box_pressure >() = real( L, "pressure" );
          box.get< tag::box_energy >() = real( L, "energy" );
          box.get< tag::box_temperature >() = real( L, "temperature" );
          lua_pop( L, 1 );
        }
      }
      lua_pop( L, 1 );
    }

    lua_pop( L, 1 );

  }
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

  cfg.get< tag::mat_spec_heat_ratio >() = real( L, "spec_heat_ratio", 1.4 );
  cfg.get< tag::mat_spec_heat_const_vol >() = real( L, "spec_heat_const_vol" );
  cfg.get< tag::mat_spec_gas_const >() = real(L, "spec_gas_const", 287.052874);
  cfg.get< tag::mat_heat_conductivity >() = real( L, "heat_conductivity" );
  cfg.get< tag::mat_dyn_viscosity >() = real( L, "dyn_viscosity", 0.0 );
  cfg.get< tag::mat_dyn_diffusivity >() = real( L, "dyn_diffusivity", 0.0 );

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
  if ( problem == "slot_cyl" or
       problem == "point_src" or
       problem == "sheardiff" )
  {
     ++n;
  }

  if (solver == "chocg") {
    n -= 2;
  }
  else if (solver == "lohcg") {
    n -= 1;
  }

  lua_pop( L, 1 );
}

static void
href( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse href table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "href" );

  auto& ht = cfg.get< tag::href >();
  ht.get< tag::t0 >() = boolean( L, "t0" );
  ht.get< tag::dt >() = boolean( L, "dt" );
  ht.get< tag::dtfreq >() = unsigint( L, "dtfreq", 5 );
  ht.get< tag::init >() = stringlist( L, "init" );
  ht.get< tag::refvar >() = unsigints( L, "refvar" );
  ht.get< tag::error >() = string( L, "error", "jump" );
  ht.get< tag::maxlevels >() = unsigint( L, "maxlevels", 2 );

  lua_pop( L, 1 );
}


static void
href_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse href table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "href_";
  auto& ht = cfg.get< tag::href_ >();
  ht.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& htk = ht[k];
    htk.get< tag::t0 >() = boolean( L, "t0" );
    htk.get< tag::dt >() = boolean( L, "dt" );
    htk.get< tag::dtfreq >() = unsigint( L, "dtfreq", 5 );
    htk.get< tag::init >() = stringlist( L, "init" );
    htk.get< tag::refvar >() = unsigints( L, "refvar" );
    htk.get< tag::error >() = string( L, "error", "jump" );
    htk.get< tag::maxlevels >() = unsigint( L, "maxlevels", 6 );

    lua_pop( L, 1 );

  }
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
// Parse pressure table from global scope
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "pressure" );

  auto& tp = cfg.get< tag::pressure >();
  tp.get< tag::iter >() = unsigint( L, "iter", 10 );
  tp.get< tag::tol >() = real( L, "tol", 1.0e-3 );
  tp.get< tag::verbose >() = unsigint( L, "verbose", 0 );
  tp.get< tag::hydrostat >() = unsigint( L, "hydrostat" );
  tp.get< tag::pc >() = string( L, "pc", "none" );
  bc_dir( L, tp.get< tag::bc_dir >() );
  bc_dirval( L, tp.get< tag::bc_dirval >() );
  bc_sym( L, tp.get< tag::bc_sym >() );

  lua_pop( L, 1 );
}

static void
pressure_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse pressure_* table from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "pressure_";
  auto& tp = cfg.get< tag::pressure_ >();
  tp.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    auto& tpk = tp[k];
    tpk.get< tag::iter >() = unsigint( L, "iter", 10 );
    tpk.get< tag::tol >() = real( L, "tol", 1.0e-3 );
    tpk.get< tag::verbose >() = unsigint( L, "verbose", 0 );
    tpk.get< tag::hydrostat >() = unsigint( L, "hydrostat" );
    tpk.get< tag::pc >() = string( L, "pc", "none" );
    bc_dir( L, tpk.get< tag::bc_dir >() );
    bc_dirval( L, tpk.get< tag::bc_dirval >() );
    bc_sym( L, tpk.get< tag::bc_sym >() );

    lua_pop( L, 1 );

  }
}

static void
part_( lua_State* L, Config& cfg, const char* def = "default" )
// *****************************************************************************
// Parse part_* field from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
//! \param[in] def Default if does not exist
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "part_";
  auto& vp = cfg.get< tag::part_ >();
  if (vp.size() != nf ) vp.resize( nf, def );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );

    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_isstring( L, -1 ), std::string(name) + " must be a string" );
      vp[k] = lua_tostring( L, -1 );
    }

    lua_pop( L, 1 );

  }
}

static void
zoltan_params_( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse zoltan_params_* field from global scope for multiple meshes
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  auto nf = cfg.get< tag::input >().size();
  if (nf == 1) return;

  std::string basename = "zoltan_params_";
  auto& vl = cfg.get< tag::zoltan_params_ >();
  vl.resize( nf );

  for (std::size_t k=0; k<nf; ++k) {

    std::string name = basename + std::to_string(k+1);
    lua_getglobal( L, name.c_str() );
    auto& v = vl[k];

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

  }
}

static void
momentum( lua_State* L, Config& cfg )
// *****************************************************************************
// Parse momentum table
//! \param[in,out] L Lua state
//! \param[in,out] cfg Config state
// *****************************************************************************
{
  lua_getglobal( L, "momentum" );

  cfg.get< tag::mom_iter >() = unsigint( L, "iter", 10 );
  cfg.get< tag::mom_tol >() = real( L, "tol", 1.0e-3 );
  cfg.get< tag::mom_verbose >() = unsigint( L, "verbose", 0 );
  cfg.get< tag::mom_pc >() = string( L, "pc", "none" );

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
    get< tag::rk >() = unsigint( L, "rk", 1, true );
    get< tag::theta >() = real( L, "theta", 0.0, true );
    get< tag::dt >() = real( L, "dt", 0.0, true );
    get< tag::turkel >() = real( L, "turkel", 0.5, true );
    get< tag::soundspeed >() = real( L, "soundspeed", 1.0, true );
    get< tag::velinf >() = vector( L, "velinf", 1.0, true );
    get< tag::t0 >() = real( L, "t0", 0.0, true );
    get< tag::reorder >() = boolean( L, "reorder", false, true );
    get< tag::flux >() = string( L, "flux", "rusanov", true );
    get< tag::steady >() = boolean( L, "steady", false, true );
    get< tag::residual >() = real( L, "residual", 0.0, true );
    get< tag::rescomp >() = unsigint( L, "rescomp", 1, true );

    get< tag::part >() = string( L, "part", "rcb", true );
    part_( L, *this, "rcb" );

    get< tag::zoltan_params >() = stringlist( L, "zoltan_params", true );
    zoltan_params_( L, *this );

    get< tag::solver >() = string( L, "solver", "riecg", true );
    get< tag::stab >() = boolean( L, "stab", true, true );
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
    ic_( L, *this );

    bc_dir( L, get< tag::bc_dir >(), true );
    bc_dir_( L, get< tag::bc_dir_ >(), get< tag::input >().size(), true );
    bc_dirval( L, get< tag::bc_dirval >(), true );
    bc_dirval_( L, get< tag::bc_dirval_ >(), get< tag::input >().size(), true );

    bc_sym( L, get< tag::bc_sym >(), true );
    bc_sym_( L, get< tag::bc_sym_ >(), get< tag::input >().size(), true );

    bc_noslip( L, get< tag::bc_noslip >(), true );
    bc_noslip_( L, get< tag::bc_noslip_ >(), get< tag::input >().size(), true );

    bc_far( L, *this );
    bc_far_( L, *this );

    bc_pre( L, *this );
    bc_pre_( L, *this );

    problem( L, *this );
    mat( L, *this );
    overset( L, *this );
    fieldout( L, *this );
    fieldout_( L, *this );
    histout( L, *this );
    histout_( L, *this );
    integout( L, *this );
    integout_( L, *this );
    diag( L, *this );
    href( L, *this );
    href_( L, *this );
    deactivate( L, *this );
    pressure( L, *this );
    pressure_( L, *this );
    momentum( L, *this );
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
