// *****************************************************************************
/*!
  \file      src/Statistics/PDFReducer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging PDFs across PEs
  \details   Custom Charm++ reducer for merging PDFs across PEs.
*/
// *****************************************************************************
#ifndef PDFReducer_h
#define PDFReducer_h

#include <tuple>
#include <vector>

#include "NoWarning/charm++.hpp"

#include "UniPDF.hpp"

namespace tk {

//! Serialize univariate PDF to raw memory stream
std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid, const std::vector< tk::UniPDF >& u );

//! \brief Charm++ custom reducer for merging a univariate PDF during reduction
//!    across PEs
CkReductionMsg*
mergeUniPDFs( int nmsg, CkReductionMsg **msgs );

} // tk::

#endif // PDFReducer_h
