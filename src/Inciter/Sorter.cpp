// *****************************************************************************
/*!
  \file      src/Inciter/Sorter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh sorter for global distributed mesh reordering
  \see       Sorter.h for more info.
*/
// *****************************************************************************

#include <vector>
#include <algorithm>

#include "Sorter.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"
#include "InciterConfig.hpp"

namespace inciter {

extern ctr::Config g_cfg;

} // inciter::

using inciter::Sorter;

Sorter::Sorter( std::size_t meshid,
                const CProxy_Transporter& transporter,
                const tk::CProxy_MeshWriter& meshwriter,
                const tk::SorterCallback& cbs,
                const CProxy_Discretization& discretization,
                const CProxy_RieCG& riecg,
                const CProxy_LaxCG& laxcg,
                const CProxy_ZalCG& zalcg,
                const CProxy_KozCG& kozcg,
                CkCallback reorderRefiner,
                const std::vector< std::size_t >& ginpoel,
                const tk::UnsMesh::CoordMap& coordmap,
                const tk::UnsMesh::Chunk& el,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::vector< std::size_t >& triinpoel,
                const std::map< int, std::vector< std::size_t > >& bnode,
                int nchare ) :
  m_meshid( meshid ),
  m_host( transporter ),
  m_meshwriter( meshwriter ),
  m_cbs( cbs ),
  m_discretization( discretization ),
  m_riecg( riecg ),
  m_laxcg( laxcg ),
  m_zalcg( zalcg ),
  m_kozcg( kozcg ),
  m_reorderRefiner( reorderRefiner ),
  m_ginpoel( ginpoel ),
  m_coordmap( coordmap ),
  m_el( el ),
  m_nbnd( 0 ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_bnode( bnode ),
  m_nchare( nchare ),
  m_nodeset( begin(ginpoel), end(ginpoel) ),
  m_noffset( 0 ),
  m_nodech(),
  m_chnode(),
  m_nodeCommMap(),
  m_reordcomm(),
  m_start( 0 ),
  m_newnodes(),
  m_newcoordmap(),
  m_reqnodes(),
  m_lower( 0 ),
  m_upper( 0 )
// *****************************************************************************
//  Constructor: prepare owned mesh node IDs for reordering
//! \param[in] meshid Mesh ID
//! \param[in] transporter Transporter (host) Charm++ proxy
//! \param[in] meshwriter Mesh writer Charm++ proxy
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] discretization Discretization Charm++ proxy
//! \param[in] riecg RieCG Charm++ proxy
//! \param[in] laxcg RieCG Charm++ proxy
//! \param[in] zalcg ZalCG Charm++ proxy
//! \param[in] kozcg KozCG Charm++ proxy
//! \param[in] reorderRefiner Callback to use to send reordered mesh to Refiner
//! \param[in] ginpoel Mesh connectivity (this chare) using global node IDs
//! \param[in] coordmap Mesh node coordinates (this chare) for global node IDs
//! \param[in] el Elements of the mesh chunk we operate on
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-faces
//! \param[in] bnode Node ids mapped to side set ids
//! \param[in] nchare Total number of Charm++ worker chares
// *****************************************************************************
{
  // Ensure boundary face ids will not index out of face connectivity
  Assert( std::all_of( begin(m_bface), end(m_bface),
            [&](const auto& s)
            { return std::all_of( begin(s.second), end(s.second),
                       [&](auto f){ return f*3+2 < m_triinpoel.size(); } ); } ),
          "Boundary face data structures inconsistent" );
}

void
Sorter::setup( std::size_t npoin )
// *****************************************************************************
// Setup chare mesh boundary node communication map
//! \param[in] npoin Total number of mesh points in mesh. Note that the number
//!   of mesh points does not have to be exactly the total number of points in
//!   the mesh. It can be a larger number, but not less. This is only used here
//!   to assign nodes to workers that will assign ids to mesh nodes during node
//!   reordering.
// *****************************************************************************
{
  // Compute the number of nodes (chunksize) a chare will build a node
  // communication map for. We compute two values of chunksize: one for when
  // the global node ids are abounded between [0...npoin-1], inclusive, and
  // another one for when the global node ids are assigned by a hash algorithm
  // during initial mesh refinement. In the latter case, the maximum
  // representable value of a std::size_t is assumed to be the large global node
  // id and is used to compute the chunksize. To compute the bin id, we attempt
  // to use the first chunksize first: if it gives a chare id that is
  // (strictly) lower than the number of chares, that's good. If not, we compute
  // the bin id based on the second chunksize, which almost always will give a
  // bin id strictly lower than the number of chares, except if the global node
  // id assigned by the hash algorithm in Refiner hits the maximum
  // representable number in std::size_t. If that is the case, we just assign
  // that node to the last chare.
  auto N = static_cast< std::size_t >( m_nchare );
  std::array< std::size_t, 2 > chunksize{{
     npoin / N, std::numeric_limits< std::size_t >::max() / N }};

  // Find chare-boundary nodes of our mesh chunk. This algorithm collects the
  // global mesh node ids on the chare boundary. A node is on a chare boundary
  // if it belongs to a face of a tetrahedron that has no neighbor tet at a
  // face. The nodes are categorized to bins that will be sent to different
  // chares to build point-to-point communication maps across all chares. The
  // binning is determined by the global node id divided by the chunksizes. See
  // discussion above on how we use two chunksizes for global node ids assigned
  // by the hash algorithm in Refiner (if initial mesh refinement has been
  // done).
  std::map< int, std::unordered_set< std::size_t > > chbnd;
  auto el = tk::global2local( m_ginpoel );      // generate local mesh data
  const auto& inpoel = std::get< 0 >( el );     // local connectivity
  auto esup = tk::genEsup( inpoel, 4 );         // elements surrounding points
  auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)
      if (esuel[mark+f] == -1)
        for (std::size_t n=0; n<3; ++n) {
          auto g = m_ginpoel[ mark+tk::lpofa[f][n] ];
          auto bin = g / chunksize[0];
          if (bin >= N) bin = g / chunksize[1];
          if (bin >= N) bin = N - 1;
          Assert( bin < N, "Will index out of number of chares" );
          chbnd[ static_cast< int >( bin ) ].insert( g );
        }
  }

  // Send boundary data in bins to chares that will compute communication maps
  // for the data in the bin. These bins form a distributed table.  Note that
  // we only send data to those chares that have data to work on. The receiving
  // sides do not know in advance if they receive messages or not.  Completion
  // is detected by having the receiver respond back and counting the responses
  // on the sender side, i.e., this chare.
  m_nbnd = chbnd.size();
  if (m_nbnd == 0) {
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::queried >() );
  } else {
    for (const auto& [ targetchare, bnd ] : chbnd) {
      thisProxy[ targetchare ].query( thisIndex, bnd );
    }
  }
}

void
Sorter::query( int fromch, const std::unordered_set< std::size_t >& bnd  )
// *****************************************************************************
// Incoming query for a list of mesh nodes for which this chare compiles node
// communication maps
//! \param[in] fromch Sender chare ID
//! \param[in] bnd Chare-boundary data from another chare
// *****************************************************************************
{
  // Store incoming nodes in node->chare and its inverse, chare->node, maps
  for (auto n : bnd) m_nodech[ n ].push_back( fromch );
  m_chnode[ fromch ].insert( begin(bnd), end(bnd) );

  // Report back to chare message received from
  thisProxy[ fromch ].recvquery();
}

void
Sorter::recvquery()
// *****************************************************************************
// Receive receipt of boundary node lists to query
// *****************************************************************************
{
  if (--m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::queried >() );
}

void
Sorter::response()
// *****************************************************************************
//  Respond to boundary node list queries
// *****************************************************************************
{
  std::unordered_map< int, std::map< int, std::unordered_set< std::size_t > > >
    exp;

  // Compute node communication map to be sent back to chares
  for (const auto& [ neighborchare, bndnodes ] : m_chnode) {
    auto& nc = exp[ neighborchare ];
    for (auto n : bndnodes)
      for (auto d : tk::cref_find(m_nodech,n))
        if (d != neighborchare)
          nc[d].insert( n );
  }

  // Send communication maps to chares that issued a query to us. Communication
  // maps were computed above for those chares that queried this map from us.
  // This data form a distributed table and we only work on a chunk of it. Note
  // that we only send data back to those chares that have queried us. The
  // receiving sides do not know in advance if the receive messages or not.
  // Completion is detected by having the receiver respond back and counting
  // the responses on the sender side, i.e., this chare.
  m_nbnd = exp.size();
  if (m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::responded >() );
  else
    for (const auto& [ targetchare, maps ] : exp)
      thisProxy[ targetchare ].bnd( thisIndex, maps );
}

void
Sorter::bnd( int fromch,
  const std::map< int, std::unordered_set< std::size_t > >& nodeCommMap )
// *****************************************************************************
// Receive boundary node communication maps for our mesh chunk
//! \param[in] fromch Sender chare ID
//! \param[in] nodeCommMap Communication map assembled by chare fromch
// *****************************************************************************
{
  for (const auto& [ neighborchare, map ] : nodeCommMap) {
    m_nodeCommMap[ neighborchare ].insert( begin(map), end(map) );
  }

  // Report back to chare message received from
  thisProxy[ fromch ].recvbnd();
}

void
Sorter::recvbnd()
// *****************************************************************************
// Receive receipt of boundary node communication map
// *****************************************************************************
{
  if (--m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::responded >() );
}

void
Sorter::start()
// *****************************************************************************
//  Start reordering (if enabled)
// *****************************************************************************
{
  if (g_cfg.get< tag::feedback >()) m_host.chcomm();

  tk::destroy( m_nodech );
  tk::destroy( m_chnode );

  if (g_cfg.get< tag::reorder >())
    mask();   // continue with mesh node reordering if requested (or required)
  else
    createDiscWorkers();  // skip mesh node reordering
}

void
Sorter::mask()
// *****************************************************************************
//  Start preparing for mesh node reordering in parallel
// *****************************************************************************
{
  // Compute asymmetric communcation map that will be used for reordering. This
  // communication map is asymmetric because it associates global mesh node IDs
  // to chares only with lower IDs than thisIndex. That is because this chare
  // will need to receive new (reorderd) node IDs only from chares with lower
  // IDs than thisIndex during node reordering. Since it only stores data for
  // lower chare IDs, it is asymmetric. Note that because of this algorithm the
  // type of m_nodeCommMap is an ordered map, because of the std::none_of()
  // algorithm needs to look at ALL chares this chare potentially communicates
  // nodes with that have lower chare IDs that thisIndex. Since the map is
  // ordered, it can walk through from the beginning of m_nodeCommMap until the
  // outer loop variable c, which is the chare ID the outer loop works on in a
  // given cycle.
  for (auto c=m_nodeCommMap.cbegin(); c!=m_nodeCommMap.cend(); ++c) {
    if (thisIndex > c->first) {
      auto& n = m_reordcomm[ c->first ];
      for (auto j : c->second) {
        if (std::none_of( m_nodeCommMap.cbegin(), c,
             [j]( const auto& s ) {
               return s.second.find(j) != end(s.second); } ))
        {
          n.insert(j);
        }
      }
      if (n.empty()) m_reordcomm.erase( c->first );
    }
  }

  // Count up total number of nodes this chare will need to receive
  auto nrecv = tk::sumvalsize( m_reordcomm );

  if ( g_cfg.get< tag::feedback >() ) m_host.chmask();

  // Compute number of mesh node IDs we will assign IDs to
  auto nuniq = m_nodeset.size() - nrecv;

  // Start computing offsets for node reordering
  thisProxy.offset( thisIndex, nuniq );
}

void
Sorter::offset( int c, std::size_t u )
// *****************************************************************************
//  Receive number of uniquely assigned global mesh node IDs from chares with
//  lower IDs than thisIndex
//! \param[in] c Chare ID
//! \param[in] u Number of mesh node IDs chare c will assign IDs to
//! \details This function computes the offset each chare will need to start
//!   assigning its new node IDs from. The offset for a chare is the
//!   offset for the previous chare plus the number of node IDs the previous
//!   chare (uniquely) assigns new IDs for minus the number of node IDs the
//!   previous chare receives from others (lower chares). This is computed here
//!   in a parallel/distributed fashion by each chare sending its number of node
//!   IDs (that it uniquely assigns) to all chares. Note that each chare would
//!   only need to send this information to chares with higher IDs, but instead
//!   this function is called in a broadcast fashion, because that is more
//!   efficient than individual calls to only chares with higher IDs. Therefore
//!   when computing the offsets, we only count the lower chares. When this is
//!   done, we have the precise asymmetric communication map as well as the
//!   start offset on all chares and so we can start the distributed global mesh
//!   node ID reordering.
// *****************************************************************************
{
  if (c < thisIndex) m_start += u;
  if (++m_noffset == m_nchare) reorder();
}

void
Sorter::reorder()
// *****************************************************************************
//  Reorder global mesh node IDs
// *****************************************************************************
{
  // Activate SDAG waits for arriving requests from other chares requesting new
  // node IDs for node IDs we assign new IDs to during reordering; and for
  // computing/receiving lower and upper bounds of global node IDs our chare's
  // linear system will operate on after reordering.
  thisProxy[ thisIndex ].wait4prep();

  // Send out request for new global node IDs for nodes we do not reorder
  for (const auto& [ targetchare, nodes ] : m_reordcomm)
    thisProxy[ targetchare ].request( thisIndex, nodes );

  // Lambda to decide if node is assigned a new ID by this chare. If node is not
  // found in the asymmetric communication map, it is owned, i.e., this chare
  // assigns its new id.
  auto ownnode = [ this ]( std::size_t p ) {
    return std::all_of( m_reordcomm.cbegin(), m_reordcomm.cend(),
                        [&](const auto& s)
                        { return s.second.find(p) == s.second.cend(); } );
  };

  // Reorder our chunk of the mesh node IDs. Looping through all of our node
  // IDs, we test if we are to assign a new ID to a node ID, and if so, we
  // assign a new ID, i.e., reorder, by constructing a map associating new to
  // old IDs (m_newnodes). We also count up the reordered nodes, which serves as
  // the new node id. We also store the node coordinates associated to the new
  // node ID.
  for (auto p : m_nodeset)
    if (ownnode(p)) {
      m_newnodes[ p ] = m_start;        // assign new node ID (reorder)
      m_newcoordmap.emplace( m_start, tk::cref_find(m_coordmap,p) );
      ++m_start;
    }

  // Trigger SDAG wait indicating that reordering our node IDs are complete
  reorderowned_complete();

  // If all our nodes have new IDs assigned, reordering complete on this chare
  if (m_newnodes.size() == m_nodeset.size()) finish();
}

void
Sorter::request( int c, const std::unordered_set< std::size_t >& nd )
// *****************************************************************************
//  Request new global node IDs for old node IDs
//! \param[in] c Chare request coming from and to which we send new IDs to
//! \param[in] nd Set of old node IDs whose new IDs are requested
// *****************************************************************************
{
  // Queue up requesting chare and node IDs
  m_reqnodes.push_back( { c, nd } );
  // Trigger SDAG wait signaling that node IDs have been requested from us
  nodes_requested_complete();
}

void
// cppcheck-suppress unusedFunction
Sorter::prepare()
// *****************************************************************************
//  Find new node IDs for old ones and return them to the requestor(s)
// *****************************************************************************
{
  // Find and return new node IDs to sender
  for (const auto& [ requestorchare, nodes ] : m_reqnodes) {
    std::unordered_map< std::size_t,
      std::tuple< std::size_t, tk::UnsMesh::Coord > > n;
    for (auto p : nodes) {
      auto newid = tk::cref_find( m_newnodes, p );
      n.emplace( p,
        std::make_tuple( newid, tk::cref_find(m_newcoordmap,newid) ) );
    }
    thisProxy[ requestorchare ].neworder( n );
  }

  tk::destroy( m_reqnodes ); // Clear queue of requests just fulfilled

  // Re-enable SDAG wait for preparing new node requests
  thisProxy[ thisIndex ].wait4prep();

  // Re-enable trigger signaling that reordering of owned node IDs are
  // complete right away
  reorderowned_complete();
}

void
Sorter::neworder( const std::unordered_map< std::size_t,
                        std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes )
// *****************************************************************************
//  Receive new (reordered) global node IDs
//! \param[in] nodes Map associating new to old node IDs
// *****************************************************************************
{
  // Store new node IDs associated to old ones, and node coordinates associated
  // to new node IDs.
  for (const auto& [ oldid, newnodes ] : nodes) {
    auto newid = std::get< 0 >( newnodes );
    m_newnodes[ oldid ] = newid;
    m_newcoordmap.emplace( newid, std::get< 1 >( newnodes ) );
  }

  // If all our nodes have new IDs assigned, reorder complete on this PE
  if (m_newnodes.size() == m_nodeset.size()) finish();
}

void
Sorter::finish()
// *****************************************************************************
//  Compute final result of reordering
//! \details Reordering is now complete on this chare. We now remap all mesh
//!   data to reflect the new ordering.
// *****************************************************************************
{
  // Update elem connectivity with the reordered node IDs
  tk::remap( m_ginpoel, m_newnodes );

  // Update node coordinate map with the reordered IDs
  m_coordmap = m_newcoordmap;

  // Update mesh chunk data structure held in our state with new node order
  m_el = tk::global2local( m_ginpoel );

  // Update symmetric chare-node communication map with the reordered IDs
  for (auto& [ neighborchare, map ] : m_nodeCommMap) {
    std::unordered_set< std::size_t > n;
    for (auto p : map) n.insert( tk::cref_find( m_newnodes, p ) );
    map = std::move( n );
  }

  // Update boundary face-node connectivity with the reordered node IDs
  tk::remap( m_triinpoel, m_newnodes );

  // Update boundary node lists with the reordered node IDs
  for (auto& [ setid, nodes ] : m_bnode) tk::remap( nodes, m_newnodes );

  // Update mesh in Refiner after reordering
  m_reorderRefiner.send();

  // Progress report to host
  if ( g_cfg.get< tag::feedback >() ) m_host.chreordered();

  createDiscWorkers();
}

void
Sorter::mesh( std::vector< std::size_t >& ginpoel,
              tk::UnsMesh::CoordMap& coordmap,
              std::vector< std::size_t >& triinpoel,
              std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
// Update mesh data we hold for whoever calls this function
//! \param[in,out] ginpoel Mesh connectivity using global IDs
//! \param[in,out] coordmap Map of mesh node coordinates
//! \param[in,out] triinpoel Boundary face-node connectivity
//! \param[in] bnode Node lists of side sets
// *****************************************************************************
{
  ginpoel = m_ginpoel;
  coordmap = m_coordmap;
  triinpoel = m_triinpoel;
  bnode = m_bnode;
}

void
Sorter::createDiscWorkers()
// *****************************************************************************
//  Create Discretization chare array elements on this PE
//! \details We create chare array elements by calling the insert() member
//!   function, which allows specifying the PE on which the array element is
//!   created. and we send each chare array element the chunk of mesh it will
//!   operate on.
// *****************************************************************************
{
  // Create worker array element using Charm++ dynamic chare array element
  // insertion.

  m_discretization[ thisIndex ].insert( m_meshid, m_host,
    m_meshwriter, m_coordmap, m_el, m_nodeCommMap, m_nchare );

  contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
              m_cbs.get< tag::discinserted >() );
}

void
Sorter::createWorkers()
// *****************************************************************************
//  Create worker chare array element
// *****************************************************************************
{
  // Make sure (bound) base is already created and accessible
  Assert( m_discretization[ thisIndex ].ckLocal() != nullptr,
          "About to pass nullptr" );

  // Create worker array element using Charm++ dynamic chare array element
  // insertion.
  const auto& solver = g_cfg.get< tag::solver >();
  if (solver == "riecg") {
    m_riecg[ thisIndex ].insert( m_discretization, m_bface, m_bnode,
                                 m_triinpoel );
  }
  else if (solver == "laxcg") {
    m_laxcg[ thisIndex ].insert( m_discretization, m_bface, m_bnode,
                                 m_triinpoel );
  }
  else if (solver == "zalcg") {
    m_zalcg[ thisIndex ].insert( m_discretization, m_bface, m_bnode,
                                 m_triinpoel );
  }
  else if (solver == "kozcg") {
    m_kozcg[ thisIndex ].insert( m_discretization, m_bface, m_bnode,
                                 m_triinpoel );
  }
  else {
    Throw( "Unknown solver: " + solver );
  }

  if ( g_cfg.get< tag::feedback >() ) m_host.chcreated();

  contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
              m_cbs.get< tag::workinserted >() );

  // Free up some memory
  tk::destroy( m_ginpoel );
  tk::destroy( m_coordmap );
  tk::destroy( m_bface );
  tk::destroy( m_triinpoel );
  tk::destroy( m_bnode );
  tk::destroy( m_nodeset );
  tk::destroy( m_nodech );
  tk::destroy( m_chnode );
  tk::destroy( m_nodeCommMap );
  tk::destroy( m_reordcomm );
  tk::destroy( m_newnodes );
  tk::destroy( m_reqnodes );
}

#include "NoWarning/sorter.def.h"
