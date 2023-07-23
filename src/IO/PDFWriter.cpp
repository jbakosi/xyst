// *****************************************************************************
/*!
  \file      src/IO/PDFWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Univariate PDF writer
  \brief     PDF writer class definition
  \details   This file defines a PDF writer class that facilitates outputing
    probability density functions (PDFs) into files in various formats using
    various configurations.
*/
// *****************************************************************************

#include <iomanip>

#include "NoWarning/exodusII.hpp"

#include "PDFWriter.hpp"
#include "Exception.hpp"

using tk::PDFWriter;

PDFWriter::PDFWriter( const std::string& filename,
                      const std::string& format,
                      std::streamsize precision ) :
  Writer( filename )
// *****************************************************************************
//  Constructor
//! \param[in] filename Output filename to which output the PDF
//! \param[in] format Configure floating-point output format for ASCII output
//! \param[in] precision Configure precision for floating-point ASCII output
// *****************************************************************************
{
  // Set floating-point format for output file stream
  if (format == "default")
    {} //m_outFile << std::defaultfloat;   GCC does not yet support this
  else if (format == "fixed")
    m_outFile << std::fixed;
  else if (format == "scientific")
    m_outFile << std::scientific;
  else Throw( "Text floating-point format not recognized." );

  // Set numeric precision for output file stream if the input makes sense
  if (precision > 0 && precision < std::numeric_limits< tk::real >::digits10+2)
    m_outFile << std::setprecision( static_cast<int>(precision) );
}

void
PDFWriter::writeTxt( const UniPDF& pdf, const tk::ctr::PDFInfo& info ) const
// *****************************************************************************
//  Write out standardized univariate PDF to file
//! \param[in] pdf Univariate PDF
//! \param[in] info PDF metadata
// *****************************************************************************
{
  const auto& name = info.name;
  const auto& uext = info.exts;
  const auto& vars = info.vars;
  const auto& it = info.it;
  const auto& time = info.time;

  assertSampleSpaceDimensions< 1 >( vars );
  assertSampleSpaceExtents< 1 >( uext );

  // Query and optionally override number of bins and minimum of sample space if
  // user-specified extents were given and copy probabilities from pdf to an
  // array for output
  std::size_t nbi;
  tk::real min, max;
  std::vector< tk::real > outpdf;
  tk::real binsize;
  std::array< long, 2*UniPDF::dim > ext;
  extents( pdf, uext, nbi, min, max, binsize, ext, outpdf );

  // Output header
  m_outFile << "# vim: filetype=sh:\n#\n"
            << "# Univariate PDF: " << name << '(' << vars[0] << ')' << '\n'
            << "# -----------------------------------------------\n"
            << "# Numeric precision: " << m_outFile.precision() << '\n'
            << "# Bin size: " << binsize << '\n'
            << "# Number of bins estimated: " << ext[1] - ext[0] + 1
            << '\n'
            << "# Number of bins output: " << nbi << '\n'
            << "# Sample space extent: [" << min << " : " << max << "]\n"
            << "# Integral: " << pdf.integral() << "\n"
            << "# Iteration: " << it << "\n"
            << "# Physical time: " << time << "\n#\n"
            << "# Example step-by-step visualization with gnuplot\n"
            << "# -----------------------------------------------\n"
            << "# gnuplot> set grid\n"
            << "# gnuplot> unset key\n"
            << "# gnuplot> set xlabel \"" << vars[0] << "\"\n"
            << "# gnuplot> set ylabel \"" << name << "(" << vars[0] << ")\"\n"
            << "# gnuplot> plot ";
  if (!uext.empty()) m_outFile << "[" << uext[0] << ':' << uext[1] << "] ";
  m_outFile << "\"" << m_filename << "\" with points\n#\n"
            << "# Gnuplot one-liner for quick copy-paste\n"
            << "# -----------------------------------------------\n"
            << "# set grid; unset key; set xlabel \"" << vars[0]
            << "\"; set ylabel \"" << name << "(" << vars[0]
            << ")\"; plot";
  if (!uext.empty()) m_outFile << " [" << uext[0] << ':' << uext[1] << "]";
  m_outFile << " \"" << m_filename << "\" w p\n#\n"
            << "# Data columns: " << vars[0] << ", " << name << "(" << vars[0]
            << ")\n# -----------------------------------------------\n";

  // If no user-specified sample space extents, output pdf map directly
  if (uext.empty()) {
    for (const auto& p : pdf.map())
      m_outFile << binsize * static_cast<tk::real>(p.first) << '\t'
                << static_cast<tk::real>(p.second) / binsize /
                   static_cast<tk::real>(pdf.nsample())
                << std::endl;
  } else { // If user-specified sample space extents, output outpdf array
    std::size_t bin = 0;
    for (const auto& p : outpdf)
      m_outFile << binsize * static_cast<tk::real>(bin++) + uext[0] << '\t'
                << p << std::endl;
  }
}
void
PDFWriter::extents( const UniPDF& pdf,
                    const std::vector< tk::real >& uext,
                    std::size_t& nbi,
                    tk::real& min,
                    tk::real& max,
                    tk::real& binsize,
                    std::array< long, 2*UniPDF::dim >& ext,
                    std::vector< tk::real >& outpdf ) const
// *****************************************************************************
//  Query extents and other metadata of univariate PDF sample space
//! \details Query and optionally override number of bins and minimum of sample
//!    space if user-specified extents were given and copy probabilities from
//!    pdf to an array for output for plotting univariate PDF.
//! \param[in] pdf Univariate PDF object
//! \param[in] uext User-specified extents of sample space
//! \param[inout] nbi Number of bins
//! \param[inout] min Minimum value of sample space
//! \param[inout] max Maximum value of sample space
//! \param[inout] binsize Bin size
//! \param[inout] ext Extents of sample space
//! \param[inout] outpdf PDF ready to be written out to file
// *****************************************************************************
{
  assertSampleSpaceExtents< 1 >( uext );

  // Query bin size and extents of sample space from PDF
  binsize = pdf.binsize();
  ext = pdf.extents();

  // Compute number of bins of sample space (min bins: 1)
  Assert( ext[1] >= ext[0], "Wrong extents in PDFWriter::extents" );
  nbi = static_cast< std::size_t >( ext[1] - ext[0] + 1 );

  // Compute minimum and maximum of sample space
  min = binsize * static_cast< tk::real >( ext[0] );
  max = binsize * static_cast< tk::real >( ext[1] );

  // Override number of bins and minimum if user-specified extents were given,
  // and copy probabilities from pdf to an array for output
  if (!uext.empty()) {
    // Override number of bins by that based on user-specified extents
    Assert( uext[1] >= uext[0],
            "Wrong user-defined extents in PDFWriter::extents" );
    nbi = static_cast< std::size_t >(
            std::lround( (uext[1] - uext[0]) / binsize ) );
    // Override extents
    min = uext[0];
    max = uext[1];

    // Size output pdf to user-requested dimensions to overridden nbi and
    // initialize output probabilities to zero
    outpdf = std::vector< tk::real >( nbi, 0.0 );

    // Fill requested region of pdf to be output from computed pdf
    for (const auto& p : pdf.map()) {
      // Compute (i.e., shift) bin indices relative to user-requested extents
      const auto bin = p.first - std::lround( uext[0] / binsize );
      // Only copy probability value if shifted bin indices fall within
      // user-requested extents (lower inclusive, upper exclusive)
      if (bin >= 0 && bin < std::lround( (uext[1] - uext[0]) / binsize )) {
        Assert( static_cast<std::size_t>(bin) < nbi,
                "Bin overflow in user-specified-extent-based bin "
                "calculation of univariate PDF extents." );
        // Copy normalized probability to output pdf
        outpdf[ static_cast<std::size_t>(bin) ] =
          p.second / binsize / static_cast<tk::real>(pdf.nsample());
      }
    }
  }
}
