// *****************************************************************************
/*!
  \file      src/Inciter/Diagnostics.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Common data for collecting diagnostics
  \details   Common data for collecting (node-, elem-, etc) diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************
#pragma once

#include <cstddef>

namespace inciter {
namespace diagnostics {

//! Number of entries in diagnostics vector (of vectors)
constexpr std::size_t NUMDIAG = 8;

//! Diagnostics labels
enum Diag { L2SOL=0,    //!< L2 norm of numerical solution
            L2RES,      //!< L2 norm of the residual
            TOTALEN,    //!< Total energy over entire domain
            L2ERR,      //!< L2 norm of numerical-analytic solution
            L1ERR,      //!< L1 norm of numerical-analytic solution
            ITER,       //!< Iteration count
            TIME,       //!< Physical time
            DT };       //!< Time step size

} // diagnostics::
} // inciter::
