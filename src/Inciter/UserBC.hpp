// *****************************************************************************
/*!
  \file      src/Inciter/UserBC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Setup boundary conditions based on user input
  \details   Setup boundary conditions based on user input.
*/
// *****************************************************************************
#ifndef UserBC_h
#define UserBC_h

#include "CartesianProduct.hpp"

namespace inciter {

//! Function object for querying the side set ids the user configured
//! \details Used to query and collect the side set ids the user has
//!   configured for all PDE types querying all BC types. Used on a
//!   Carteisan product of 2 type lists: PDE types and BC types.
struct UserBC {
  const ctr::InputDeck& inputdeck;
  std::unordered_set< int >& userbc;
  explicit UserBC( const ctr::InputDeck& i, std::unordered_set< int >& u )
    : inputdeck(i), userbc(u) {}
  template< typename U > void operator()( brigand::type_<U> ) {
    using tag::param;
    using eq = typename brigand::front< U >;
    using bc = typename brigand::back< U >;
    for (const auto& s : inputdeck.get< param, eq, tag::bc, bc >())
      for (const auto& i : s) userbc.insert( std::stoi(i) );
  }
};

//! Function object for querying the side set ids for time dependent BCs
//! \details Used to query and collect the side set ids the user has
//!   configured for all PDE types querying time dependent BCs. Used on
//!   PDE type list.
struct UserTimedepBC {
  const ctr::InputDeck& inputdeck;
  std::unordered_set< int >& userbc;
  explicit UserTimedepBC( const ctr::InputDeck& i,
    std::unordered_set< int >& u )
    : inputdeck(i), userbc(u) {}
  template< typename eq > void operator()( brigand::type_<eq> ) {
    using tag::param;
    for (const auto& sys : inputdeck.get< param, eq, tag::bctimedep >()) {
      for (const auto& b : sys) {
        for (auto i : b.template get< tag::sideset >())
          userbc.insert( std::stoi(i) );
      }
    }
  }
};

} // inciter::

#endif // UserBC_h
