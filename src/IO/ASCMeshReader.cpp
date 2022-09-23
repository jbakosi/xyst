// *****************************************************************************
/*!
  \file      src/IO/ASCMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ASC mesh reader class definition
  \details   ASC mesh reader class definition. Mesh reader facilitating reading
             a mesh from a simple text file used by Jacob Waltz's Chicoma code.
*/
// *****************************************************************************

#include <array>
#include <istream>
#include <string>
#include <vector>
#include <cstddef>

#include "Types.hpp"
#include "Exception.hpp"
#include "UnsMesh.hpp"
#include "Reorder.hpp"
#include "ASCMeshReader.hpp"

using tk::ASCMeshReader;

void
ASCMeshReader::readHeader()
// *****************************************************************************
//  Read ASC mesh header
// *****************************************************************************
{
  std::string s;

  // ndim
  m_inFile >> s;
  Assert( s == "*ndim", "Invalid keyword, expected: '*ndim'" );
  int ndim;
  m_inFile >> ndim;
  Assert( ndim == 3, "Only 3D meshes are supported" );

  // numNodeSets (throw away)
  m_inFile >> s;
  Assert( s == "*numNodeSets", "Invalid keyword, expected: '*numNodeSets'" );
  int nn;
  m_inFile >> nn;

  // numSideSets (throw away)
  m_inFile >> s;
  Assert( s == "*numSideSets", "Invalid keyword, expected: '*numSideSets'" );
  int ns;
  m_inFile >> ns;
}

void
ASCMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read ASC mesh
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  // Read header
  readHeader();
  // Read nodes
  readNodes( mesh );
  // Read elements
  readElements( mesh );
}

void
ASCMeshReader::readNodes( UnsMesh& mesh )
// *****************************************************************************
//  Read nodes
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  std::string s;

  m_inFile >> s;
  Assert( s == "*nodes", "Invalid keyword, expected: '*nodes'" );
  int nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0,
          "Number of nodes must be greater than zero in file " + m_filename  );

  // Read in node coordinates: x-coord y-coord z-coord, ignore node IDs, assume
  // sorted
  for (int i=0; i<nnode; ++i) {
    int n;
    tk::real x, y, z;
    m_inFile >> n >> x >> y >> z;
    mesh.x().push_back( x );
    mesh.y().push_back( y );
    mesh.z().push_back( z );
  }
}

void
ASCMeshReader::readElements( UnsMesh& mesh )
// *****************************************************************************
//  Read element connectivity
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  std::string s;

  m_inFile >> s;
  Assert( s == "*cells", "Invalid keyword, expected: '*cells'" );

  int nel;
  m_inFile >> nel;
  ErrChk( nel > 0,
          "Number of cells must be greater than zero in file " + m_filename  );

  // Read in tetrahedra element tags and connectivity
  for (int i=0; i<nel; ++i) {
    int a, b, c;
    std::array< std::size_t, 4 > n;
    // ignore cell id, a, b
    m_inFile >> a >> b >> c >> n[3] >> n[0] >> n[1] >> n[2];
    mesh.tetinpoel().push_back( n[0] );
    mesh.tetinpoel().push_back( n[1] );
    // switch nodes 2 and 3 to enforce positive volume
    mesh.tetinpoel().push_back( n[3] );
    mesh.tetinpoel().push_back( n[2] );
  }

  // Shift node IDs to start from zero
  shiftToZero( mesh.tetinpoel() );
}
