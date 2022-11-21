// *****************************************************************************
/*!
  \file      src/IO/MeshDetect.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Unstructured mesh file format detection
  \details   Unstructured mesh file format detection functions.
*/
// *****************************************************************************

#include <string>
#include <fstream>
#include <stdexcept>

#include "MeshDetect.hpp"
#include "Exception.hpp"
#include "Reader.hpp"

namespace tk {

MeshReaderType
detectInput( const std::string& filename )
// *****************************************************************************
//  Detect input mesh file type
//! \param[in] filename File to open and detect its type
//! \return enum specifying the mesh reader type
// *****************************************************************************
{
  std::ifstream inFile;

  // Check if file exists, throw exception if it does not
  inFile.open( filename, std::ifstream::in );
  ErrChk( inFile.good(), "Failed to open file: " + filename );

  inFile.get();
  ErrChk( !inFile.fail(), "Failed to read from file: " + filename );

  // Close it
  inFile.close();
  ErrChk( !inFile.fail(), "Failed to close file: " + filename );

  // Get first three letters from input file
  std::string s( Reader( filename ).firstline().substr(0,4) );

  if ( s.find("$Me") != std::string::npos ) {
    return MeshReaderType::GMSH;
  } else if ( s.find("CDF") != std::string::npos ||
              s.find("HDF") != std::string::npos ) {
    return MeshReaderType::EXODUSII;
  } else if ( s.find("<?x") != std::string::npos ) {
    return MeshReaderType::HYPER;
  } else if ( s.find("*nd") != std::string::npos ) {
    return MeshReaderType::ASC;
  } else if ( s.find("   ") != std::string::npos ) {
    return MeshReaderType::UGRID;
  } else if ( s.find(" npo") != std::string::npos ) {
    return MeshReaderType::RDGFLO;
  } else {
    try {
      // cppcheck-suppress ignoredReturnValue
      std::stoi(s);    // try to convert to an integer
    } catch ( const std::invalid_argument& ) {
      Throw( "Input mesh file type could not be determined from header: " +
             filename );
    }
    // could also catch std::out_of_range, the other exception potentially
    // thrown by std::stoi(), but a three-digit integer will always fit into int

    // if we got here, the above string-to-integer conversion succeeded
    return MeshReaderType::NETGEN;
  }
}

MeshWriterType
pickOutput( const std::string& filename )
// *****************************************************************************
//  Determine output mesh file type
//! \param[in] filename Filename to pick its type based on extension given
//! \return enum specifying the mesh writer type
// *****************************************************************************
{
  // Get extension of input file name
  std::string fn = filename;
  std::string ext( fn.substr(fn.find_last_of(".") + 1) );

  if ( ext == "msh" ) {
    return MeshWriterType::GMSH;
  } else if ( ext == "exo" || ext == "h5" ) {
    return MeshWriterType::EXODUSII;
  } else if ( ext == "mesh" ) {
    return MeshWriterType::NETGEN;
  } else {
    Throw( "Output mesh file type could not be determined from extension of "
           "filename '" + filename + "'; valid extensions are: "
           "'msh' for Gmsh, 'exo' or 'h5' for ExodusII, 'mesh' for Netgen's "
           "neutral" );
  }
}

} // tk::
