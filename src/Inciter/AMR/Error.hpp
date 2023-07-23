// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"
#include "AMR/edge.hpp"

namespace AMR {

//! Class for computing error estimates for mesh refinement
class Error {

  public:
    //! Compute error estimate for a scalar quantity
    tk::real scalar( const tk::Fields& u,
                     const edge_t& edge,
                     uint64_t c,
                     const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const std::pair< std::vector< std::size_t >,
                                      std::vector< std::size_t > >& esup,
                     const std::string& err ) const;

  private:
    //! Estimate error for scalar quantity on edge based on jump in solution
    tk::real
    error_jump( const tk::Fields& u,
                const edge_t& edge,
                uint64_t c ) const;

    //! Estimate error for scalar quantity on edge based on Hessian of solution
    tk::real
    error_hessian( const tk::Fields& u,
                   const edge_t& edge,
                   uint64_t c,
                   const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& esup ) const;
};

} // AMR::
