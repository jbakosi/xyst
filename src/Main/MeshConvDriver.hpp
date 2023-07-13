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

#include "MeshConv/CmdLine/CmdLine.hpp"

//! Mesh converter declarations and definitions
namespace meshconv {

//! Mesh converter driver used polymorphically with tk::Driver
class MeshConvDriver {

  public:
    //! Constructor
    explicit MeshConvDriver( const ctr::CmdLine& cmdline );

    //! Execute
    void execute( int sig ) const;

  private:
    const bool m_reorder;               //!< Whether to also reorder mesh nodes
    std::string m_input;                //!< Input file name
    std::string m_output;               //!< Output file name
};

} // meshconv::

#endif // MeshConvDriver_h
