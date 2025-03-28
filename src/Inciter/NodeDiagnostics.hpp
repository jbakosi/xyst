// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     NodeDiagnostics class for collecting diagnostics
  \details   NodeDiagnostics class for collecting diagnostics, e.g., residuals,
    and various norms of errors while solving partial differential equations.
*/
// *****************************************************************************
#pragma once

#include "Discretization.hpp"
#include "PUPUtil.hpp"

namespace inciter {

//! NodeDiagnostics class used to compute diagnostics while integrating PDEs
class NodeDiagnostics {

  public:
    //! Configure Charm++ custom reduction types initiated from this class
    static void registerReducers();

    //! Compute diagnostics for density-based solvers
    bool rhocompute( Discretization& d,
                     const tk::Fields& u,
                     const tk::Fields& un,
                     uint64_t diag_iter ) const;

    //! Compute diagnostics for pressure-based solvers
    bool precompute( Discretization& d,
                     const tk::Fields& u,
                     const tk::Fields& un,
                     const std::vector< tk::real >& p,
                     const std::vector< tk::real >& dp,
                     uint64_t diag_iter ) const;

    //! Compute diagnostics for artificial compressibility solvers
    bool accompute( Discretization& d,
                    const tk::Fields& u,
                    const tk::Fields& un,
                    uint64_t diag_iter ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    void pup( PUP::er & ) {}
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Diagnostics object reference
    friend void operator|( PUP::er& p, NodeDiagnostics& d ) { d.pup(p); }
    //@}
};

} // inciter::
