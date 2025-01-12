// *****************************************************************************
/*!
  \file      src/IO/PDFWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     PDF writer class declaration
  \details   This file declares a PDF writer class that facilitates outputing
    probability density functions (PDFs) into files in various formats using
    various configurations.
*/
// *****************************************************************************
#pragma once

#include <string>
#include <iostream>

#include "Writer.hpp"
#include "UniPDF.hpp"
#include "StatCtr.hpp"

namespace tk {

//! PDFWriter : Writer
class PDFWriter : public tk::Writer {

  public:
    //! Constructor
    explicit PDFWriter(
      const std::string& filename,
      const std::string& format = "default",
      std::streamsize precision = std::cout.precision() );

    //! Write univariate PDF to text file
    void writeTxt( const UniPDF& pdf, const tk::ctr::PDFInfo& info ) const;

  private:
    //! Assert the number of sample space dimensions given
    template< std::size_t size, class Container >
    void assertSampleSpaceDimensions( [[maybe_unused]] const Container& c )
    const {
      Assert( c.size() == size,
              "Number of sample space variables must equal " +
              std::to_string( size ) + " in PDF writer." );
    }

    //! Assert the number of sample space extents given
    template< std::size_t size, class Container >
    void assertSampleSpaceExtents( const Container& c ) const {
      if (!c.empty())
        Assert( c.size() == size*2,
                "PDF user-specified sample space extents must be defined by " +
                std::to_string( size*2 ) +" real numbers: minx, maxx, ..." );
    }

    //! Query extents and other metadata of univariate PDF sample space
    void extents( const UniPDF& pdf,
                  const std::vector< tk::real >& uext,
                  std::size_t& nbi,
                  tk::real& min,
                  tk::real& max,
                  tk::real& binsize,
                  std::array< long, 2*UniPDF::dim >& ext,
                  std::vector< tk::real >& outpdf ) const;
};

} // tk::
