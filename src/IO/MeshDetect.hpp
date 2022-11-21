// *****************************************************************************
/*!
  \file      src/IO/MeshDetect.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unstructured mesh file format detection
  \details   Unstructured mesh file format detection functions.
*/
// *****************************************************************************
#ifndef MeshDetect_h
#define MeshDetect_h

#include <iosfwd>

namespace tk {

//! Supported mesh readers
enum class MeshReaderType : uint8_t { GMSH = 0
                                    , NETGEN
                                    , EXODUSII
                                    , HYPER
                                    , ASC
                                    , UGRID
                                    , RDGFLO
                                    };

//! Supported mesh writers
enum class MeshWriterType : uint8_t { GMSH = 0
                                    , NETGEN
                                    , EXODUSII
                                    };

//! Detect input mesh file type
MeshReaderType
detectInput( const std::string& filename );

//! Determine output mesh file type
MeshWriterType
pickOutput( const std::string& filename );

} // tk::

#endif // MeshDetect_h
