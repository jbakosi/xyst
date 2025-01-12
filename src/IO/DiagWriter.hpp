// *****************************************************************************
/*!
  \file      src/IO/DiagWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Text diagnostics writer definition
  \details   This file feines the ASCII diagnostics writer class that
     facilitates outputing diagnostics to text files.
*/
// *****************************************************************************
#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "Types.hpp"
#include "Writer.hpp"

namespace tk {

//! \brief DiagWriter : tk::Writer
//! \details ASCII diagnostics writer class that facilitates outputing
//!   diagnostics to text files.
class DiagWriter : public tk::Writer {

  public:
    //! Constructor
    explicit DiagWriter( const std::string& filename,
                         const std::string& format = "default",
                         std::streamsize precision = std::cout.precision(),
                         std::ios_base::openmode mode = std::ios_base::out );

    //! Write out diagnostics file header
    void header( const std::vector< std::string >& name ) const;

    //! Write diagnostics file
    std::size_t write( uint64_t it,
                       tk::real t,
                       tk::real dt,
                       const std::vector< tk::real >& diagnostics );

  private:
    int m_precision;     //!< Floating-point precision in digits
    int m_width;         //!< Floating-point number width
};

} // tk::
