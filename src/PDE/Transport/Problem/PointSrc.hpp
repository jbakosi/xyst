// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/PointSrc.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for scalar transport equations
  \details   This file declares a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.h
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.h for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************
#ifndef TransportProblemPointSrc_h
#define TransportProblemPointSrc_h

#include <vector>
#include <array>

#include "Types.hpp"
#include "Fields.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

/*! Transport PDE problem: diffusion of a shear layer
    \details This class sets up a dispersion from point source.
    \see Bakosi, Franzese, Boybeyi, Joint PDF modeling of
       turbulent flow and dispersion in an urban street canyon, Boundary-Layer
        Meteorol., 131, 2, 2009. http://dx.doi.org/10.1007/s10546-009-9370-x
*/
class TransportProblemPointSrc {
  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::transport;

  public:
    //! Initialize numerical solution
    static std::vector< tk::real >
    initialize( ncomp_t system, ncomp_t ncomp, tk::real x, tk::real y,
                tk::real z, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t, ncomp_t, tk::real, tk::real, tk::real, tk::real )
    { return {}; }

    //! Add source
    static void src( ncomp_t system, ncomp_t offset, tk::real t,
                     const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::Fields& U );

    //! Do error checking on PDE parameters
    void errchk( ncomp_t system, ncomp_t ) const;

    //! Assign prescribed shear velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t, ncomp_t ncomp, tk::real, tk::real, tk::real,
                        tk::real )
    {
      return std::vector< std::array<tk::real,3> >( ncomp, {0,0,0} );
    }

    //! Return true if velocity is prescribed by this problem
    static bool prescribedVel() noexcept { return false; }

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::POINT_SRC; }
};

} // inciter::

#endif // TransportProblemPointSrc_h
