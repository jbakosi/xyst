// *****************************************************************************
/*!
  \file      src/Inciter/PartsReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging mesh part assignments across PEs
*/
// *****************************************************************************

#include "PartsReducer.hpp"
#include "Exception.hpp"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize(
  const std::unordered_map< std::size_t, std::size_t >& d )
// *****************************************************************************
// Serialize parts to raw memory stream
//! \param[in] d Mesh graph to aggregate
//! \return Pair of the length and the raw stream containing the serialized data
// *****************************************************************************
{
  // Prepare for serializing parts to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_map< std::size_t, std::size_t >& >( d );

  // Create raw character stream to store the serialized parts
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain graph
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_map< std::size_t, std::size_t >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeParts( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging mesh graphs during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing serialized data
//! \return Aggregated parts built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized mesh graph
  std::unordered_map< std::size_t, std::size_t > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize from raw stream
  // cppcheck-suppress uninitvar
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack graph
    std::unordered_map< std::size_t, std::size_t > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    // cppcheck-suppress uninitvar
    curCreator | w;
    // Aggregate
    for (const auto& [g,p] : w) {
      Assert( v.find(g) == end(v), "Conflicting partition assignments" );
      v[g] = p;
    }
  }

  // Serialize concatenated parts to raw stream
  auto stream = tk::serialize( v );

  // Forward serialized parts
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
