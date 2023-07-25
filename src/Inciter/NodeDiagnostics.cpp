// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
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
NodeDiagnostics::compute( Discretization& d,
                          const tk::Fields& u,
                          const tk::Fields& un,
                          uint64_t diag_iter ) const
// *****************************************************************************
//  Compute diagnostics, e.g., residuals, norms of errors, etc.
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
  // Optionally collect diagnostics and send for aggregation across all workers

  using namespace diagnostics;

  if ( !((d.It()+1) % diag_iter) ) {     // if remainder, don't dump

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
        for (std::size_t c=0; c<s.size(); ++c) an(i,c,0) = s[c];
      }
    }

    // Put in norms sweeping our mesh chunk
    for (std::size_t i=0; i<u.nunk(); ++i) {
      // Compute sum for L2 norm of the numerical solution
      for (std::size_t c=0; c<ncomp; ++c)
        diag[L2SOL][c] += u(i,c,0) * u(i,c,0) * v[i];
      // Compute sum for L2 norm of the residual
      for (std::size_t c=0; c<ncomp; ++c)
        diag[L2RES][c] += (u(i,c,0)-un(i,c,0)) * (u(i,c,0)-un(i,c,0)) * v[i];
      // Compute sum for the total energy over the entire domain (first entry)
      diag[TOTALSOL][0] += u(i,4,0) * v[i];
      // Compute sum for L2 norm of the numerical-analytic solution
      if (sol) {
        auto nu = u[i];
        nu[1] /= nu[0];
        nu[2] /= nu[0];
        nu[3] /= nu[0];
        nu[4] = nu[4] / nu[0] - 0.5*(nu[1]*nu[1] + nu[2]*nu[2] + nu[3]*nu[3]);
        for (std::size_t c=0; c<5; ++c) {
          auto du = nu[c] - an(i,c,0);
          diag[L2ERR][c] += du * du * v[i];
          diag[L1ERR][c] += std::abs( du ) * v[i];
        }
        for (std::size_t c=5; c<ncomp; ++c) {
          auto du = u(i,c,0) - an(i,c,0);
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
    auto stream = serialize( d.MeshId(), diag );
    d.contribute( stream.first, stream.second.get(), DiagMerger,
      CkCallback(CkIndex_Transporter::diagnostics(nullptr), d.Tr()) );

    return true;        // diagnostics have been computed
  }

  return false;         // diagnostics have not been computed
}
