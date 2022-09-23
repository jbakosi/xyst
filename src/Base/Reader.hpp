// *****************************************************************************
/*!
  \file      src/Base/Reader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reader base class declaration
  \details   Reader base class declaration. Reader base serves as a base class
    for various file readers. It does generic low-level I/O, e.g., opening and
    closing a file, and associated error handling.
*/
// *****************************************************************************
#ifndef Reader_h
#define Reader_h

#include <fstream>
#include <vector>
#include <cstddef>

namespace tk {

//! Reader base serves as a base class for various file readers. It does generic
//! low-level I/O, e.g., opening and closing a file, and associated error
//! handling.
class Reader {

  public:
    //! Constructor: Acquire file handle
    explicit Reader( const std::string& filename,
                     std::ios_base::openmode mode = std::ifstream::in );

    //! Destructor: Release file handle
    virtual ~Reader() noexcept;

    //! Unformatted read
    //! \param[in] data Buffer to read to
    //! \param[in] count Number of characters to read
    void read( char* data, std::streamsize count )
    { m_inFile.read( data, count ); }

    //! Return first line (for detection of file type based on header)
    std::string firstline();

    //! Read file and return a string for each line
    std::vector< std::string > lines();

    //! Read a given line from file
    std::string line( std::size_t lineNum );

  protected:
    const std::string m_filename;            //!< File name
    std::ifstream m_inFile;                  //!< File input stream
};

} // tk::

#endif // Reader_h
