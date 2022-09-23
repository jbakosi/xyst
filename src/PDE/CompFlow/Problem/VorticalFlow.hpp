// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/VorticalFlow.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the single-material compressible flow
    equations
  \details   This file defines a Problem policy class for the single-material
    compressible flow equations, defined under PDE/CompFlow/. See
    PDE/CompFlow/Problem.h for general requirements on Problem policy classes
    for CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowProblemVorticalFlow_h
#define CompFlowProblemVorticalFlow_h

#include <string>
#include <unordered_set>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! CompFlow system of PDEs problem: vortical flow
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemVorticalFlow {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::compflow;

  public:
    //! Initialize numerical solution
    static tk::InitializeFn::result_type
    initialize( ncomp_t system, ncomp_t, tk::real x, tk::real y,
                tk::real z, tk::real );

    //! Evaluate analytical solution at (x,y,z) for all components
    static tk::InitializeFn::result_type
    analyticSolution( ncomp_t system, ncomp_t, tk::real x, tk::real y,
                      tk::real z, tk::real );

    //! Compute and return source term for vortical flow manufactured solution
    //! \param[in] system Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in,out] r Density source
    //! \param[in,out] ru X momentum source
    //! \param[in,out] rv Y momentum source
    //! \param[in,out] rw Z momentum source
    //! \param[in,out] re Specific total energy source
    //! \note The function signature must follow tk::SrcFn
    static tk::CompFlowSrcFn::result_type
    src( ncomp_t system, tk::real x, tk::real y, tk::real z, tk::real,
         tk::real& r, tk::real& ru, tk::real& rv, tk::real& rw, tk::real& re )
    {
      using tag::param; using tag::compflow;

      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
      const auto& b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
      // ratio of specific heats
      auto g = g_inputdeck.get< param, compflow, tag::gamma >()[ system ][0];
      // evaluate solution at x,y,z
      auto s = initialize( system, 5, x, y, z, 0.0 );

      // density source
      r = 0.0;
      // momentum source
      ru = a*s[1]/s[0] - b*s[2]/s[0];
      rv = b*s[1]/s[0] + a*s[2]/s[0];
      rw = 0.0;
      // energy source
      re = (ru*s[1] + rv*s[2])/s[0] + 8.0*a*a*a*z*z/(g-1.0);
    }

    //! Return analytic field names to be output to file
    std::vector< std::string > analyticFieldNames( ncomp_t ) const;

    //! Return names of integral variables to be output to diagnostics file
    std::vector< std::string > names( ncomp_t ) const;

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::VORTICAL_FLOW; }
};

} // inciter::

#endif // CompFlowProblemVorticalFlow_h
