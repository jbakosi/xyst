// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     NodeDiagnostics class for collecting nodal diagnostics
  \details   NodeDiagnostics class for collecting nodal diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************

#include "Diagnostics.hpp"
#include "NodeDiagnostics.hpp"
#include "DiagReducer.hpp"
#include "Discretization.hpp"
#include "Problems.hpp"

namespace inciter {

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::NodeDiagnostics;

void
NodeDiagnostics::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details This routine is supposed to be called from a Charm++ initnode
//!   routine. Since the runtime system executes initnode routines exactly once
//!   on every logical node early on in the Charm++ init sequence, they must be
//!   static as they are called without an object. See also: Section
//!   "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  DiagMerger = CkReduction::addReducer( diagnostics::mergeDiag );
}

bool
NodeDiagnostics::rhocompute( Discretization& d,
                             const tk::Fields& u,
                             const tk::Fields& un,
                             uint64_t diag_iter ) const
// *****************************************************************************
//  Compute diagnostics for density-based solvers
//! \param[in] d Discretization proxy to read from
//! \param[in] u Current solution vector
//! \param[in] un Previous solution vector
//! \param[in] diag_iter Diagnostics output frequency
//! \return True if diagnostics have been computed
//! \details Diagnostics are defined as some norm, e.g., L2 norm, of a quantity,
//!   computed in mesh nodes, A, as ||A||_2 = sqrt[ sum_i(A_i)^2 V_i ],
//!   where the sum is taken over all mesh nodes and V_i is the nodal volume.
//!   We send multiple sets of quantities to the host for aggregation across
//!   the whole mesh. The final aggregated solution will end up in
//!   Transporter::diagnostics(). Aggregation of the partially computed
//!   diagnostics is done via potentially different policies for each field.
//! \see inciter::mergeDiag(), src/Inciter/Diagnostics.hpp
// *****************************************************************************
{
  using namespace diagnostics;

  // Only compute diagnostics if needed in this time step
  if ( (d.It()+1) % diag_iter ) return false;

  auto ncomp = u.nprop();

  // Diagnostics vector (of vectors) during aggregation. See
  // Inciter/Diagnostics.h.
  std::vector< std::vector< tk::real > >
    diag( NUMDIAG, std::vector< tk::real >( ncomp, 0.0 ) );

  const auto& v = d.V();  // nodal volumes without contributions from others

  // query function to evaluate analytic solution (if defined)
  auto sol = problems::SOL();

  // Evaluate analytic solution (if defined)
  auto an = u;
  if (sol) {
    const auto& coord = d.Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (std::size_t i=0; i<u.nunk(); ++i) {
      auto s = sol( x[i], y[i], z[i], d.T()+d.Dt() );
      s[1] /= s[0];
      s[2] /= s[0];
      s[3] /= s[0];
      s[4] = s[4] / s[0] - 0.5*(s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
    }
  }

  // Put in norms sweeping our mesh chunk
  for (std::size_t i=0; i<u.nunk(); ++i) {
    // Compute sum for L2 norm of the numerical solution
    for (std::size_t c=0; c<ncomp; ++c)
      diag[L2SOL][c] += u(i,c) * u(i,c) * v[i];
    // Compute sum for L2 norm of the residual
    for (std::size_t c=0; c<ncomp; ++c)
      diag[L2RES][c] += (u(i,c)-un(i,c)) * (u(i,c)-un(i,c)) * v[i];
    // Compute sum for the total energy over the entire domain (first entry)
    diag[TOTALEN][0] += u(i,4) * v[i];
    // Compute sum for L2 norm of the numerical-analytic solution
    if (sol) {
      auto nu = u[i];
      nu[1] /= nu[0];
      nu[2] /= nu[0];
      nu[3] /= nu[0];
      nu[4] = nu[4] / nu[0] - 0.5*(nu[1]*nu[1] + nu[2]*nu[2] + nu[3]*nu[3]);
      for (std::size_t c=0; c<5; ++c) {
        auto du = nu[c] - an(i,c);
        diag[L2ERR][c] += du * du * v[i];
        diag[L1ERR][c] += std::abs( du ) * v[i];
      }
      for (std::size_t c=5; c<ncomp; ++c) {
        auto du = u(i,c) - an(i,c);
        diag[L2ERR][c] += du * du * v[i];
        diag[L1ERR][c] += std::abs( du ) * v[i];
      }
    }
  }

  // Append diagnostics vector with metadata on the current time step
  // ITER:: Current iteration count (only the first entry is used)
  // TIME: Current physical time (only the first entry is used)
  // DT: Current physical time step size (only the first entry is used)
  diag[ITER][0] = static_cast< tk::real >( d.It()+1 );
  diag[TIME][0] = d.T() + d.Dt();
  diag[DT][0] = d.Dt();

  // Contribute to diagnostics
  auto stream = serialize( d.MeshId(), ncomp, diag );
  d.contribute( stream.first, stream.second.get(), DiagMerger,
    CkCallback(CkIndex_Transporter::rhodiagnostics(nullptr), d.Tr()) );

  return true;        // diagnostics have been computed
}

bool
NodeDiagnostics::precompute( Discretization& d,
                             const tk::Fields& u,
                             const tk::Fields& un,
                             const std::vector< tk::real >& p,
                             const std::vector< tk::real >& dp,
                             uint64_t diag_iter ) const
// *****************************************************************************
//  Compute diagnostics for pressure-based solvers
//! \param[in] d Discretization proxy to read from
//! \param[in] u Current solution vector
//! \param[in] un Previous solution vector
//! \param[in] p Current pressure solution
//! \param[in] dp Recent pressure solution increment
//! \param[in] diag_iter Diagnostics output frequency
//! \return True if diagnostics have been computed
//! \details Diagnostics are defined as some norm, e.g., L2 norm, of a quantity,
//!   computed in mesh nodes, A, as ||A||_2 = sqrt[ sum_i(A_i)^2 V_i ],
//!   where the sum is taken over all mesh nodes and V_i is the nodal volume.
//!   We send multiple sets of quantities to the host for aggregation across
//!   the whole mesh. The final aggregated solution will end up in
//!   Transporter::diagnostics(). Aggregation of the partially computed
//!   diagnostics is done via potentially different policies for each field.
//! \see inciter::mergeDiag(), src/Inciter/Diagnostics.hpp
// *****************************************************************************
{
  using namespace diagnostics;

  // Only compute diagnostics if needed in this time step
  if ( (d.It()+1) % diag_iter ) return false;

  Assert( p.size() == u.nunk(), "Size mismatch" );
  Assert( p.size() == dp.size(), "Size mismatch" );
  Assert( u.nunk() == un.nunk(), "Size mismatch" );

  auto ncomp = u.nprop();

  const auto& v = d.V();  // nodal volumes without contributions from others

  // query function to evaluate analytic solution (if defined)
  auto pressure_sol = problems::PRESSURE_SOL();

  // Evaluate analytic solution (if defined)
  auto ap = p;
  if (pressure_sol) {
    ncomp = 0;
    const auto& coord = d.Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (std::size_t i=0; i<p.size(); ++i) {
      ap[i] = pressure_sol( x[i], y[i], z[i] );
    }
  }

  // query function to evaluate analytic solution (if defined)
  auto sol = problems::SOL();

  // Evaluate analytic solution (if defined)
  auto an = u;
  if (sol) {
    const auto& coord = d.Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (std::size_t i=0; i<u.nunk(); ++i) {
      auto s = sol( x[i], y[i], z[i], d.T()+d.Dt() );
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
    }
  }

  // Diagnostics vector (of vectors) during aggregation. See
  // Inciter/Diagnostics.h.
  std::vector< std::vector< tk::real > >
    diag( NUMDIAG, std::vector< tk::real >( ncomp+1, 0.0 ) );

  // Put in norms sweeping our mesh chunk
  for (std::size_t i=0; i<u.nunk(); ++i) {
    // Compute sum for L2 norm of the numerical solution
    diag[L2SOL][0] += p[i] * p[i] * v[i];
    for (std::size_t c=0; c<ncomp; ++c) {
      diag[L2SOL][c+1] += u(i,c) * u(i,c) * v[i];
    }
    // Compute sum for L2 norm of the residual
    diag[L2RES][0] += dp[i] * dp[i] * v[i];
    for (std::size_t c=0; c<ncomp; ++c) {
      diag[L2RES][c+1] += (u(i,c)-un(i,c)) * (u(i,c)-un(i,c)) * v[i];
    }
    // Compute sum for the total energy over the entire domain
    diag[TOTALEN][0] += 0.0 * v[i];
    // Compute sum for norms of the numerical-analytic pressure solution
    if (pressure_sol) {
      auto pd = p[i] - ap[i];
      diag[L2ERR][0] += pd * pd * v[i];
      diag[L1ERR][0] += std::abs( pd ) * v[i];
    }
    // Compute sum for norms of the numerical-analytic adv/diff solution
    if (sol) {
      for (std::size_t c=0; c<ncomp; ++c) {
        auto du = u(i,c) - an(i,c);
        diag[L2ERR][c+1] += du * du * v[i];
        diag[L1ERR][c+1] += std::abs( du ) * v[i];
      }
    }
  }

  // Append diagnostics vector with metadata on the current time step
  // ITER:: Current iteration count (only the first entry is used)
  // TIME: Current physical time (only the first entry is used)
  // DT: Current physical time step size (only the first entry is used)
  diag[ITER][0] = static_cast< tk::real >( d.It() );
  diag[TIME][0] = d.T();
  diag[DT][0] = d.Dt();

  // Contribute to diagnostics
  auto stream = serialize( d.MeshId(), ncomp+1, diag );
  d.contribute( stream.first, stream.second.get(), DiagMerger,
    CkCallback(CkIndex_Transporter::prediagnostics(nullptr), d.Tr()) );

  return true;        // diagnostics have been computed
}

bool
NodeDiagnostics::accompute( Discretization& d,
                            const tk::Fields& u,
                            const tk::Fields& un,
                            uint64_t diag_iter ) const
// *****************************************************************************
//  Compute diagnostics for artificial compressibility solvers
//! \param[in] d Discretization proxy to read from
//! \param[in] u Current solution vector
//! \param[in] un Previous solution vector
//! \param[in] diag_iter Diagnostics output frequency
//! \return True if diagnostics have been computed
//! \details Diagnostics are defined as some norm, e.g., L2 norm, of a quantity,
//!   computed in mesh nodes, A, as ||A||_2 = sqrt[ sum_i(A_i)^2 V_i ],
//!   where the sum is taken over all mesh nodes and V_i is the nodal volume.
//!   We send multiple sets of quantities to the host for aggregation across
//!   the whole mesh. The final aggregated solution will end up in
//!   Transporter::diagnostics(). Aggregation of the partially computed
//!   diagnostics is done via potentially different policies for each field.
//! \see inciter::mergeDiag(), src/Inciter/Diagnostics.hpp
// *****************************************************************************
{
  using namespace diagnostics;

  // Only compute diagnostics if needed in this time step
  if ( (d.It()+1) % diag_iter ) return false;

  auto ncomp = u.nprop();

  // Diagnostics vector (of vectors) during aggregation. See
  // Inciter/Diagnostics.h.
  std::vector< std::vector< tk::real > >
    diag( NUMDIAG, std::vector< tk::real >( ncomp, 0.0 ) );

  const auto& v = d.V();  // nodal volumes without contributions from others

  // query function to evaluate analytic solution (if defined)
  auto sol = problems::SOL();

  // Evaluate analytic solution (if defined)
  auto an = u;
  if (sol) {
    const auto& coord = d.Coord();
    const auto& x = coord[0];
    const auto& y = coord[1];
    const auto& z = coord[2];
    for (std::size_t i=0; i<u.nunk(); ++i) {
      auto s = sol( x[i], y[i], z[i], d.T()+d.Dt() );
      for (std::size_t c=0; c<s.size(); ++c) an(i,c) = s[c];
    }
  }

  // Put in norms sweeping our mesh chunk
  for (std::size_t i=0; i<u.nunk(); ++i) {
    // Compute sum for L2 norm of the numerical solution
    for (std::size_t c=0; c<ncomp; ++c)
      diag[L2SOL][c] += u(i,c) * u(i,c) * v[i];
    // Compute sum for L2 norm of the residual
    for (std::size_t c=0; c<ncomp; ++c)
      diag[L2RES][c] += (u(i,c)-un(i,c)) * (u(i,c)-un(i,c)) * v[i];
    // Compute sum for the total energy over the entire domain
    diag[TOTALEN][0] += 0.0 * v[i];
    // Compute sum for norms of the numerical-analytic adv/diff solution
    if (sol) {
      for (std::size_t c=1; c<ncomp; ++c) {
        auto du = u(i,c) - an(i,c);
        diag[L2ERR][c] += du * du * v[i];
        diag[L1ERR][c] += std::abs( du ) * v[i];
      }
    }
  }

  // Append diagnostics vector with metadata on the current time step
  // ITER:: Current iteration count (only the first entry is used)
  // TIME: Current physical time (only the first entry is used)
  // DT: Current physical time step size (only the first entry is used)
  diag[ITER][0] = static_cast< tk::real >( d.It() );
  diag[TIME][0] = d.T();
  diag[DT][0] = d.Dt();

  // Contribute to diagnostics
  auto stream = serialize( d.MeshId(), ncomp, diag );
  d.contribute( stream.first, stream.second.get(), DiagMerger,
    CkCallback(CkIndex_Transporter::acdiagnostics(nullptr), d.Tr()) );

  return true;        // diagnostics have been computed
}
