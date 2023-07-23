// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************
#ifndef MeshConvDriver_h
#define MeshConvDriver_h

#include <iosfwd>

#include "MeshConvCmdLine.hpp"

//! Mesh converter declarations and definitions
namespace meshconv {

//! Mesh converter driver used polymorphically with tk::Driver
class MeshConvDriver {

  public:
    void convert( const std::string& input,
                  const std::string& output,
                  bool reorder ) const;
};

} // meshconv::

#endif // MeshConvDriver_h
