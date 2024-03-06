// *****************************************************************************
/*!
  \file      src/Inciter/GraphReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging mesh graphs across PEs
*/
// *****************************************************************************

#include "GraphReducer.hpp"

namespace tk {

std::pair< int, std::unique_ptr<char[]> >
serialize(
  const std::unordered_map< std::size_t, std::vector< std::size_t > >& d )
// *****************************************************************************
// Serialize diagnostics to raw memory stream
//! \param[in] d Mesh graph to aggregate
//! \return Pair of the length and the raw stream containing the serialized data
// *****************************************************************************
{
  // Prepare for serializing diagnostics to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_map< std::size_t,
                        std::vector< std::size_t > >& >( d );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain graph
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_map< std::size_t,
                        std::vector< std::size_t > >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeGraph( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging mesh graphs during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing serialized data
//! \return Aggregated diagnostics built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized mesh graph
  std::unordered_map< std::size_t, std::vector< std::size_t > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize from raw stream
  // cppcheck-suppress uninitvar
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack graph
    std::unordered_map< std::size_t, std::vector< std::size_t > > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    // cppcheck-suppress uninitvar
    curCreator | w;
    // Aggregate
    for (const auto& [g,n] : w) {
      auto& sup = v[g];
      if (sup.empty()) {
        sup.insert( end(sup), begin(n), end(n) );
      } else {
        // lower compute node id owns
        sup[0] = std::min( sup[0], n[0] );
        sup.insert( end(sup), n.begin()+1, n.end() );
      }
    }
  }

  // Serialize concatenated diagnostics vector to raw stream
  auto stream = tk::serialize( v );

  // Forward serialized diagnostics
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // tk::
