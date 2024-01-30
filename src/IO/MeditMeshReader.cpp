// *****************************************************************************
/*!
  \file      src/IO/MeditMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Medit mesh reader class definition
  \details   Medit mesh reader class definition. Only supports tetrahedra.
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
#include "DerivedData.hpp"
#include "MeditMeshReader.hpp"

using tk::MeditMeshReader;

void
MeditMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Read Medit mesh
//! \param[in] mesh Unstructured mesh object
// *****************************************************************************
{
  std::size_t ntet = 0;
  // cppcheck-suppress unreadVariable
  std::size_t triid = 0;

  while (!m_inFile.eof()) {
    std::string s;
    m_inFile >> s;

    if (s == "Triangles") {

      int nel;
      m_inFile >> nel;
      //std::cout << "ntri: " << nel << '\n';
      ErrChk( nel > 0, "Number of triangles (surface elements) must be greater "
                       "than zero in file " + m_filename );
      for (int i=0; i<nel; ++i) {
        int tag;
        std::array< std::size_t, 3 > n;
        m_inFile >> n[0] >> n[1] >> n[2] >> tag;
        auto& triinpoel = mesh.triinpoel();
        triinpoel.push_back( n[0] );
        triinpoel.push_back( n[1] );
        triinpoel.push_back( n[2] );
        mesh.bface()[ tag ].push_back( triid++ );
        mesh.faceid()[ tag ].push_back( 0 );
      }

    } else if (s == "Tetrahedra") {

      int nel;
      m_inFile >> nel;
      //std::cout << "ntet: " << nel << '\n';
      ErrChk( nel > 0, "Number of tetrahedra (volume elements) must be greater "
                       "than zero in file " + m_filename );

      // Read in tetrahedra element tags and connectivity
      for (int i=0; i<nel; ++i) {
        int tag;
        std::array< std::size_t, 4 > n;
        m_inFile >> n[0] >> n[1] >> n[2] >> n[3] >> tag;
        mesh.tetinpoel().push_back( n[0] );
        mesh.tetinpoel().push_back( n[1] );
        mesh.tetinpoel().push_back( n[2] );
        mesh.tetinpoel().push_back( n[3] );
        ++ntet;
      }

    } else if (s == "Vertices") {

      int nnode;
      m_inFile >> nnode;
      //std::cout << "nnode: " << nnode << '\n';
      ErrChk( nnode > 0,
          "Number of nodes must be greater than zero in file " + m_filename  );

      // Read in node coordinates: x-coord y-coord z-coord
      int tag;
      for (int i=0; i<nnode; ++i) {
        tk::real x, y, z;
        m_inFile >> x >> y >> z >> tag;
        mesh.x().push_back( x );
        mesh.y().push_back( y );
        mesh.z().push_back( z );
      }
    }
  }

  // adjust boundary element ids in exo, since tets are written first
  for (auto& [sid,triangles] : mesh.bface()) {
    for (auto& tid : triangles) {
      // cppcheck-suppress useStlAlgorithm
      tid += ntet;
    }
  }

  // Shift node IDs to start from zero
  shiftToZero( mesh.triinpoel() );

  // Shift node IDs to start from zero
  shiftToZero( mesh.tetinpoel() );
}
