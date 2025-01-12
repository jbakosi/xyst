// *****************************************************************************
/*!
  \file      src/Base/Callback.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple types used for passing Charm++ callbacks
  \details   Tagged tuple types used for passing Charm++ callbacks.
*/
// *****************************************************************************
#pragma once

#include "NoWarning/charm++.hpp"

#include "TaggedTuple.hpp"

namespace tag {
struct load;
struct partitioned;
struct distributed;
struct refinserted;
struct queried;
struct responded;
struct compatibility;
struct bndint;
struct matched;
struct refined;
struct discinserted;
struct workinserted;
} // tag::

namespace tk {

using PartitionerCallback =
  tk::TaggedTuple< brigand::list<
      tag::queried,        CkCallback
    , tag::responded,      CkCallback
    , tag::load,           CkCallback
    , tag::partitioned,    CkCallback
    , tag::distributed,    CkCallback
    , tag::refinserted,    CkCallback
  > >;

using RefinerCallback =
  tk::TaggedTuple< brigand::list<
      tag::queried,        CkCallback
    , tag::responded,      CkCallback
    , tag::compatibility,  CkCallback
    , tag::bndint,         CkCallback
    , tag::matched,        CkCallback
    , tag::refined,        CkCallback
  > >;

using SorterCallback =
  tk::TaggedTuple< brigand::list<
      tag::queried,        CkCallback
    , tag::responded,      CkCallback
    , tag::discinserted,   CkCallback
    , tag::workinserted,   CkCallback
  > >;

} // tk::
