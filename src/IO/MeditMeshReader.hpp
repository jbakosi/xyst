// *****************************************************************************
/*!
  \file      src/IO/MeditMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Medit mesh reader class declaration
  \details   Medit mesh reader class declaration. Only supports tetrahedra.
*/
// *****************************************************************************
#ifndef MeditMeshReader_h
#define MeditMeshReader_h

#include <iosfwd>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief MeditMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html
class MeditMeshReader : public Reader {

  public:
    //! Constructor
    explicit MeditMeshReader( const std::string& filename ) :
      Reader( filename ) {}

    //! Read Medit mesh
    void readMesh( UnsMesh& mesh );
};

} // tk::

#endif // MeditMeshReader_h
