// *****************************************************************************
/*!
  \file      src/Inciter/Integrals.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Common data for collecting integrals
  \details   Common data for collecting integrals.
*/
// *****************************************************************************
#pragma once

namespace inciter {
namespace integrals {

//! Number of entries in integrals vector (of vectors)
const std::size_t NUMINT = 7;

//! Integral labels
enum Ints { ITER=0              //!< Iteration count
          , TIME                //!< Physical time
          , DT                  //!< Time step size
          , MASS_FLOW_RATE      //!< Mass flow rate
          , FORCE_X             //!< Force in x direction
          , FORCE_Y             //!< Force in y direction
          , FORCE_Z             //!< Force in z direction
          };

} // integrals::
} // inciter::
