// *****************************************************************************
/*!
  \file      src/Inciter/DiagReducer.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Custom Charm++ reducer for merging diagnostics across PEs
  \details   Custom Charm++ reducer for merging diagnostics across PEs.
*/
// *****************************************************************************

#include <stddef.h>
#include <type_traits>
#include <memory>

#include "DiagReducer.hpp"
#include "Diagnostics.hpp"
#include "Exception.hpp"

namespace inciter {
namespace diagnostics {

std::pair< int, std::unique_ptr<char[]> >
serialize( std::size_t meshid, const std::vector< std::vector< tk::real > >& d )
// *****************************************************************************
// Serialize diagnostics to raw memory stream
//! \param[in] meshid Mesh ID
//! \param[in] d Diagnostics vector of vectors (of eq components)
//! \return Pair of the length and the raw stream containing the serialized
//!   vectors
// *****************************************************************************
{
  // Prepare for serializing diagnostics to a raw binary stream, compute size
  PUP::sizer sizer;
  sizer | meshid;
  sizer | const_cast< std::vector< std::vector< tk::real > >& >( d );

  // Create raw character stream to store the serialized vectors
  std::unique_ptr<char[]> flatData = std::make_unique<char[]>( sizer.size() );

  // Serialize vector, each message will contain a vector
  PUP::toMem packer( flatData.get() );
  packer | meshid;
  packer | const_cast< std::vector< std::vector< tk::real > >& >( d );

  // Return size of and raw stream
  return { sizer.size(), std::move(flatData) };
}

CkReductionMsg*
mergeDiag( int nmsg, CkReductionMsg **msgs )
// *****************************************************************************
// Charm++ custom reducer for merging diagnostics during reduction across PEs
//! \param[in] nmsg Number of messages in msgs
//! \param[in] msgs Charm++ reduction message containing the serialized
//!   diagnostics
//! \return Aggregated diagnostics built for further aggregation if needed
// *****************************************************************************
{
  using namespace diagnostics;

  // Will store deserialized diagnostics vector of vectors
  std::size_t meshid;
  std::vector< std::vector< tk::real > > v;

  // Create PUP deserializer based on message passed in
  PUP::fromMem creator( msgs[0]->getData() );

  // Deserialize vector from raw stream
  // cppcheck-suppress uninitvar
  creator | meshid;
  creator | v;

  for (int m=1; m<nmsg; ++m) {
    // Unpack vector
    std::size_t mid;
    std::vector< std::vector< tk::real > > w;
    PUP::fromMem curCreator( msgs[m]->getData() );
    // cppcheck-suppress uninitvar
    curCreator | mid;
    curCreator | w;
    // Aggregate diagnostics vector
    // cppcheck-suppress uninitvar
    // cppcheck-suppress unreadVariable
    meshid = mid;
    Assert( v.size() == w.size(),
            "Size mismatch during diagnostics aggregation" );
    Assert( v.size() == NUMDIAG,
            "Size mismatch during diagnostics aggregation" );
    // cppcheck-suppress unsignedLessThanZero
    for (std::size_t i=0; i<v.size(); ++i)
      Assert( v[i].size() == w[i].size(),
              "Size mismatch during diagnostics aggregation" );
    // Apply diagnostics aggregation policy
    // Sum for L2 normal of the numerical solution for all scalar components
    for (std::size_t i=0; i<v[L2SOL].size(); ++i) v[L2SOL][i] += w[L2SOL][i];
    // Sum for the L2 norm of the residual of all components
    for (std::size_t i=0; i<v[L2RES].size(); ++i) v[L2RES][i] += w[L2RES][i];
    // Sum of the total energy over the entire domain
    v[TOTALEN][0] += w[TOTALEN][0];
    // Sum for the L2 norm of the numerical - analytical solution for all comps
    for (std::size_t i=0; i<v[L2ERR].size(); ++i) v[L2ERR][i] += w[L2ERR][i];
    // Sum for the L1 norm of the numerical - analytical solution for all comps
    for (std::size_t i=0; i<v[L1ERR].size(); ++i) v[L1ERR][i] += w[L1ERR][i];
    // Copy ITER, TIME, DT
    for (std::size_t i=0; i<v[ITER].size(); ++i) v[ITER][i] = w[ITER][i];
    for (std::size_t i=0; i<v[TIME].size(); ++i) v[TIME][i] = w[TIME][i];
    for (std::size_t i=0; i<v[DT].size(); ++i) v[DT][i] = w[DT][i];
  }

  // Serialize concatenated diagnostics vector to raw stream
  auto stream = serialize( meshid, v );

  // Forward serialized diagnostics
  return CkReductionMsg::buildNew( stream.first, stream.second.get() );
}

} // diagnostics::
} // inciter::
