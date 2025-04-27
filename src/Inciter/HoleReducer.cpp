// *****************************************************************************
/*!
  \file      src/Inciter/HoleReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for aggregating hole data across PEs
  \details   Custom Charm++ reducer for aggregating hole data across PEs.
*/
// *****************************************************************************

#include <stddef.h>
#include <type_traits>
#include <memory>

#include "HoleReducer.hpp"
#include "Exception.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

std::pair< int, std::unique_ptr<char[]> >
serialize( const std::unordered_map< std::size_t, std::vector< tk::real > >& d )
// *****************************************************************************
// Serialize hole surface data to raw memory stream
//! \param[in] d Hole data structure
//! \return Pair of the length and the raw stream containing serialized data
// *****************************************************************************
{
  // Prepare for serializing hole data to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::unordered_map< std::size_t,
                                          std::vector< tk::real > >& >( d );

  // Create raw character stream to store the serialized data
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain a vector
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::unordered_map< std::size_t,
                                           std::vector< tk::real > >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeHole( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging hole data during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized
//!   hole surface data
//! \return Aggregated hole surface data built for further aggregation if needed
// *****************************************************************************
{
  std::unordered_map< std::size_t, std::vector< tk::real > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack partial hole data
    std::unordered_map< std::size_t, std::vector< tk::real > > inholes;
    PUP::fromMem curCreator( msgs[m]->getData() );
    curCreator | inholes;
    // Aggregate hole data
    for (auto&& [hid,data] : inholes) {
      std::move( begin(data), end(data), std::back_inserter(v[hid]) );
     }
  }

  // Serialize concatenated hole surface data to raw stream
  auto stream = serialize( v );

  // Forward serialized hole surface data
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // inciter::
