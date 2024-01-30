// *****************************************************************************
/*!
  \file      src/Inciter/IntegralReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging integrals across PEs
  \details   Custom Charm++ reducer for merging integrals across PEs.
*/
// *****************************************************************************

#include <memory>

#include "Integrals.hpp"
#include "IntegralReducer.hpp"
#include "Exception.hpp"

namespace inciter {
namespace integrals {

std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid, const std::vector< std::map<int,tk::real> >& d )
// *****************************************************************************
// Serialize integrals to raw memory stream
//! \param[in] meshid Mesh ID
//! \param[in] d Integral contributions
//! \return Pair of the length and the raw stream containing the serialized
//!   vectors
// *****************************************************************************
{
  // Prepare for serializing integrals to a raw binary stream, compute size
  PUP::sizer sizer;
  // cppcheck-suppress uninitvar
  sizer | meshid;
  sizer | const_cast< std::vector< std::map< int, tk::real > >& >( d );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize integrals
  PUP::toMem packer( flatData.get() );
  // cppcheck-suppress uninitvar
  packer | meshid;
  packer | const_cast< std::vector< std::map< int, tk::real > >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeIntegrals( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging integrals during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized
//!   integrals
//! \return Aggregated integrals built for further aggregation if needed
// *****************************************************************************
{
  // Will store deserialized integrals
  std::size_t meshid;
  std::vector< std::map< int, tk::real > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack vector
    std::size_t mid;
    std::vector< std::map< int, tk::real > > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    // cppcheck-suppress uninitvar
    curCreator | mid;
    curCreator | w;
    // Aggregate integrals
    // cppcheck-suppress uninitvar
    // cppcheck-suppress unreadVariable
    meshid = mid;
    Assert( v.size() == w.size(), "Size mismatch during integrals aggregation");
    Assert( v.size() == NUMINT, "Size mismatch during integrals aggregation" );
    // Aggregate applying integrals aggregation policy
    // Copy ITER, TIME, DT
    // cppcheck-suppress containerOutOfBounds
    v[ITER] = w[ITER];
    // cppcheck-suppress containerOutOfBounds
    v[TIME] = w[TIME];
    // cppcheck-suppress containerOutOfBounds
    v[DT] = w[DT];
    // Sum integrals
    // cppcheck-suppress containerOutOfBounds
    for (const auto& [s,d] : w[MASS_FLOW_RATE]) v[MASS_FLOW_RATE][s] += d;
  }

  // Serialize concatenated diagnostics vector to raw stream
  // cppcheck-suppress uninitvar
  auto stream = serialize( meshid, v );

  // Forward serialized diagnostics
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // integrals::
} // inciter::
