// *****************************************************************************
/*!
  \file      src/IO/GmshMeshIO.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Gmsh mesh reader and writer related types
  \details   Gmsh mesh reader and writer related types.
*/
// *****************************************************************************
#ifndef GmshMeshIO_h
#define GmshMeshIO_h

namespace tk {

//! Identifiers of supported Gmsh elements
enum GmshElemType { TRI = 2,
                    TET = 4,
                    PNT = 15 };

//! Supported Gmsh mesh file types
enum class GmshFileType { UNDEFINED = -1,
                          ASCII = 0,
                          BINARY = 1 };

} // tk::

#endif // GmshMeshIO_h
