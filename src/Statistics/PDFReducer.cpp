// *****************************************************************************
/*!
  \file      src/Statistics/PDFReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging PDFs across PEs
  \details   Custom Charm++ reducer for merging PDFs across PEs.
*/
// *****************************************************************************

#include "PDFReducer.hpp"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid, const std::vector< tk::UniPDF >& u )
// *****************************************************************************
// Serialize univariate PDFs to raw memory stream
//! \param[in] meshid Mesh ID
//! \param[in] u Univariate PDFs
//! \return Pair of the length and the raw stream containing the serialized PDFs
// *****************************************************************************
{
  // Prepare for serializing PDF to a raw binary stream, compute size
  PUP::sizer sizer;
  // cppcheck-suppress uninitvar
  sizer | meshid;
  sizer | const_cast< std::vector< tk::UniPDF >& >( u );

  // Create raw character stream to store the serialized PDF
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize PDF, the message will contain a univariate PDF
  PUP::toMem packer( flatData.get() );
  // cppcheck-suppress uninitvar
  packer | meshid;
  packer | const_cast< std::vector< tk::UniPDF >& >( u );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeUniPDFs( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging a univariate PDFs during reduction across
// PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized PDF
//! \return Aggregated PDF built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized univariate PDFs
  std::size_t meshid;
  std::vector< tk::UniPDF > updf;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize PDFs from raw stream
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | updf;

  for (int m=1; m<nmsg; ++m) {
    // Unpack PDF
    std::size_t mid;
    std::vector< tk::UniPDF > u;
    PUP::fromMem curCreator( msgs[m]->getData() );
    // cppcheck-suppress uninitvar
    curCreator | mid;
    curCreator | u;
    // Merge PDFs
    // cppcheck-suppress uninitvar
    meshid = mid;
    std::size_t i = 0;
    // cppcheck-suppress knownEmptyContainer
    // cppcheck-suppress containerOutOfBounds
    for (const auto& p : u) updf[i++].addPDF( p );
  }

  // Serialize vector of merged PDF to raw stream
  // cppcheck-suppress uninitvar
  auto stream = tk::serialize( meshid, updf );

  // Forward serialized PDFs
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
