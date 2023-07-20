// *****************************************************************************
/*!
  \file      src/Control/LuaParser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Lua parser
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

#include "LuaParser.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

static constexpr auto largeint = std::numeric_limits< int64_t >::max();
static constexpr auto largeuint = std::numeric_limits< uint64_t >::max();
static constexpr auto largereal = std::numeric_limits< double >::max();

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
vector( lua_State* L, const char* name, bool global = false )
// *****************************************************************************
// Parse vector table from global scope or table
//! \param[in] L Lua state
//! \param[in] name Label to parse
//! \param[in] global True to parse from global scope, false from table on stack
//! \return Vector components parsed
// *****************************************************************************
{
  std::vector< double > v( 3, largereal );

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
fieldout( lua_State* L )
// *****************************************************************************
// Parse fieldout table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "fieldout" );

  g_inputdeck.get< tag::fieldout_iter >() = unsigint( L, "iter" );
  g_inputdeck.get< tag::fieldout_time >() = real( L, "time" );
  g_inputdeck.get< tag::fieldout_range >() = range( L );
  g_inputdeck.get< tag::fieldout >() = sideset( L );

  lua_pop( L, 1 );
}

static void
histout( lua_State* L )
// *****************************************************************************
// Parse histout table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "histout" );

  g_inputdeck.get< tag::histout_iter >() = unsigint( L, "iter" );
  g_inputdeck.get< tag::histout_time >() = real( L, "time" );
  g_inputdeck.get< tag::histout_range >() = range( L );
  g_inputdeck.get< tag::histout_precision >() = sigint( L, "precision", 8 );
  g_inputdeck.get< tag::histout_format >() = string( L, "format" );

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "points" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "histout points must be a table" );
      auto& r = g_inputdeck.get< tag::histout >();
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
integout( lua_State* L )
// *****************************************************************************
// Parse integout table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "integout" );

  g_inputdeck.get< tag::integout_iter >() = unsigint( L, "iter" );
  g_inputdeck.get< tag::integout_time >() = real( L, "time" );
  g_inputdeck.get< tag::integout_range >() = range( L );
  g_inputdeck.get< tag::integout_precision >() = sigint( L, "precision", 8 );
  g_inputdeck.get< tag::integout_format >() = string( L, "format" );
  g_inputdeck.get< tag::integout >() = sideset( L );

  lua_pop( L, 1 );
}

static void
diag( lua_State* L )
// *****************************************************************************
// Parse diag table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "diag" );

  g_inputdeck.get< tag::diag_iter >() = unsigint( L, "iter", 1 );
  g_inputdeck.get< tag::diag_precision >() = sigint( L, "precision", 8 );
  g_inputdeck.get< tag::diag_format >() = string( L, "format" );

  lua_pop( L, 1 );
}

static void
bc_dir( lua_State* L )
// *****************************************************************************
// Parse bc_dir table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "bc_dir" );

  if (lua_istable( L, -1 )) {
    auto& s = g_inputdeck.get< tag::bc_dir >();
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
bc_sym( lua_State* L )
// *****************************************************************************
// Parse bc_sym table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "bc_sym" );

  g_inputdeck.get< tag::bc_sym >() = sideset( L );

  lua_pop( L, 1 );
}

static void
bc_far( lua_State* L )
// *****************************************************************************
// Parse bc_far table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "bc_far" );

  g_inputdeck.get< tag::bc_far >() = sideset( L );
  g_inputdeck.get< tag::bc_far_velocity >() = vector( L, "velocity" );
  g_inputdeck.get< tag::bc_far_density >() = real( L, "density" );
  g_inputdeck.get< tag::bc_far_pressure >() = real( L, "pressure" );

  lua_pop( L, 1 );
}

static void
bc_pre( lua_State* L )
// *****************************************************************************
// Parse bc_pre table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "bc_pre" );

  if (lua_istable( L, -1 )) {
    auto bc_pre_set = [&]( const char* name ) {
      lua_getfield( L, -1, name );
      if (!lua_isnil( L, -1 )) {
        ErrChk( lua_istable( L, -1 ), "bc_pre must be a table" );
        g_inputdeck.get< tag::bc_pre >().push_back( sideset( L ) );
        g_inputdeck.get< tag::bc_pre_density >().push_back( real(L,"density") );
        g_inputdeck.get< tag::bc_pre_pressure >().push_back( real(L,"pressure") );
      }
      lua_pop( L, 1 );
    };

    bc_pre_set( "inlet" );
    bc_pre_set( "outlet" );
  }

  lua_pop( L, 1 );
}

static void
ic( lua_State* L )
// *****************************************************************************
// Parse ic table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "ic" );

  g_inputdeck.get< tag::ic_velocity >() = vector( L, "velocity" );
  g_inputdeck.get< tag::ic_density >() = real( L, "density" );
  g_inputdeck.get< tag::ic_pressure >() = real( L, "pressure" );
  g_inputdeck.get< tag::ic_energy >() = real( L, "energy" );
  g_inputdeck.get< tag::ic_temperature >() = real( L, "temperature" );

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
      auto& boxes = g_inputdeck.get< tag::ic >();
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
mat( lua_State* L )
// *****************************************************************************
// Parse mat table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "mat" );

  g_inputdeck.get< tag::mat_spec_heat_ratio >() =
    real( L, "spec_heat_ratio" );

  g_inputdeck.get< tag::mat_spec_heat_const_vol >() =
    real( L, "spec_heat_const_vol" );

  g_inputdeck.get< tag::mat_heat_conductivity >() =
    real( L, "heat_conductivity" );

  lua_pop( L, 1 );
}

static void
problem( lua_State* L )
// *****************************************************************************
// Parse problem table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "problem" );

  g_inputdeck.get< tag::problem >() = string( L, "name", "userdef" );
  g_inputdeck.get< tag::problem_alpha >() = real( L, "alpha" );
  g_inputdeck.get< tag::problem_kappa >() = real( L, "kappa" );
  g_inputdeck.get< tag::problem_r0 >() = real( L, "r0" );
  g_inputdeck.get< tag::problem_p0 >() = real( L, "p0" );
  g_inputdeck.get< tag::problem_ce >() = real( L, "ce" );
  g_inputdeck.get< tag::problem_beta >() = vector( L, "beta" );

  if (lua_istable( L, -1 )) {
    lua_getfield( L, -1, "source" );
    if (!lua_isnil( L, -1 )) {
      ErrChk( lua_istable( L, -1 ), "problem source must be a table" );
      auto& s = g_inputdeck.get< tag::problem_source >();
      s.get< tag::location >() = vector( L, "location" );
      s.get< tag::radius >() = real( L, "radius" );
      s.get< tag::release_time >() = real( L, "release_time" );
    }
    lua_pop( L, 1 );
  }

  const auto& problem = g_inputdeck.get< tag::problem >();
  auto& n = g_inputdeck.get< tag::problem_ncomp >();
  n = 5;
  if (problem == "slot_cyl" || problem == "point_src") ++n;

  lua_pop( L, 1 );
}

static void
href( lua_State* L )
// *****************************************************************************
// Parse href table
//! \param[in] L Lua state
// *****************************************************************************
{
  lua_getglobal( L, "href" );

  g_inputdeck.get< tag::href_t0 >() = boolean( L, "t0" );
  g_inputdeck.get< tag::href_dt >() = boolean( L, "dt" );
  g_inputdeck.get< tag::href_dtfreq >() = unsigint( L, "dtfreq", 5 );
  g_inputdeck.get< tag::href_init >() = stringlist( L, "init" );
  g_inputdeck.get< tag::href_refvar >() = unsigints( L, "refvar" );
  g_inputdeck.get< tag::href_error >() = string( L, "error", "jump" );
  g_inputdeck.get< tag::href_maxlevels >() = unsigint( L, "maxlevels", 2 );

  lua_pop( L, 1 );
}

std::string
parseLua( const char* inputfile )
// *****************************************************************************
// Parse lua input file
//! \param[in] inputfile Input file to parse
// *****************************************************************************
{
  lua_State* L = luaL_newstate();
  luaL_openlibs( L );

  std::string err;

  if (luaL_dofile( L, inputfile ) == LUA_OK) {

    g_inputdeck.get< tag::nstep >() = unsigint( L, "nstep", largeuint, true );
    g_inputdeck.get< tag::term >() = real( L, "term", largereal, true );
    g_inputdeck.get< tag::ttyi >() = unsigint( L, "ttyi", 1, true );
    g_inputdeck.get< tag::cfl >() = real( L, "cfl", 0.0, true );
    g_inputdeck.get< tag::dt >() = real( L, "dt", 0.0, true );
    g_inputdeck.get< tag::t0 >() = real( L, "t0", 0.0, true );
    g_inputdeck.get< tag::reorder >() = boolean( L, "reorder", false, true );
    g_inputdeck.get< tag::steady >() = boolean( L, "steady", false, true );
    g_inputdeck.get< tag::residual >() = real( L, "residual", 1.0e-14, true );
    g_inputdeck.get< tag::rescomp >() = unsigint( L, "rescomp", 1, true );
    g_inputdeck.get< tag::part >() = string( L, "part", "rcb", true );

    ic( L );
    bc_dir( L );
    bc_sym( L );
    bc_far( L );
    bc_pre( L );
    problem( L );
    mat( L );
    fieldout( L );
    histout( L );
    integout( L );
    diag( L );
    href( L );

  } else {
    err = lua_tostring( L, -1 );
  }

  lua_close( L );

  return err;
}

} // ::inciter
