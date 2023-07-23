// *****************************************************************************
/*!
  \file      src/Base/Print.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     General purpose pretty printer functionality
  \details   This file contains general purpose printer functions. Using the
    functions defined here provides formatting, and a consistent look with
    simple client-side code.
*/
// *****************************************************************************
#ifndef Print_h
#define Print_h

#include <iostream>
#include <cmath>

#include "Timer.hpp"
#include "Exception.hpp"
#include "PrintUtil.hpp"

namespace tk {

//! Pretty printer base. Contains general purpose printer functions. Using the
//! functions defined here provides formatting, and a consistent look with
//! simple client-side code.
class Print {

  public:
    //! Constructor
    explicit Print() : m_stream( std::cout ) {}

    //! Operator << for printing any type to the verbose stream.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] t Reference to an arbitrary object of type T. T must define
    //! operator<< for std::ostream-compatible streams.
    //! \return The internal stream buffer of the stream
    template< typename T >
    friend const Print& operator<<( const Print& os, const T& t ) {
      os.m_stream << t; return os;
    }

    //! Formatted print of section title
    //! \param[in] t Section title to be printed
    void section( const std::string& t ) const {
      std::string underline( t.size(), '-' );
      m_stream << '\n' << t.c_str() << '\n' << underline.c_str() << '\n';
    }

    //! Formatted print of item: name : value
    //! \param[in] name Item name to be printed
    //! \param[in] value Item value to be printed
    template< typename T >
    void item( const std::string& name, const T& value ) const {
      m_stream << name.c_str() << ": ";
      if constexpr( std::is_same_v< T, std::string > ) {
        m_stream << value.c_str();
      } else {
        m_stream << value;
      }
      m_stream << '\n';
    }

    //! Formatted print of item: h:m:s.
    //! \param[in] name Item name to be printed
    //! \param[in] watch Watch (in hours, minutes, seconds) to be printed as
    //!   item value
    void item( const std::string& name, const tk::Timer::Watch& watch ) const {
      m_stream << name.c_str() << ": "
               << watch.hrs.count() << ':'
               << watch.min.count() << ':'
               << watch.sec.count() << '\n';
    }

    //! Formatted print of a performance statistic (an item of a list)
    //! \param[in] name Performance statistic name to be printed
    //! \param[in] value Performance statistic value
    void perfitem( const std::string& name, tk::real value ) const {
      m_stream << name.c_str() << " : " << value << '\n';
    }

    //! Formatted print of elapsed times
    //! \param[in] t Title of section containing a list of elapsed times
    //! \param[in] clock std::vector of strings (clock names) and associated
    //!   timers which could be in various formats as long as there is a
    //!   corresponding item() overload that can apply operator << for outputing
    //!   their value to an output stream. Examples of allowed ClockFormats are:
    //!   tk::Timer::Watch, which is a struct containing a timestamp in h:m:s
    //!   format, and the return value of Timer::dsec(), which is a tk::real.
    template< class ClockFormat >
    void time( const std::string& t,
               const std::vector<
                 std::pair< std::string, ClockFormat > >& clock ) const
    {
      section( t );
      for (const auto& c : clock) item( c.first, c.second );
      m_stream << '\n';
    }

    //! Echo formatted print of a diagnostics message within a progress section
    //! \param[in] labels Label parts of diagnostics message
    //! \param[in] values Value parts of diagnostics message
    //! \param[in] precr If true start with a CR/LF, if false end with it
    //! \note The number of labels and values must equal.
    void diag( const std::vector< std::string >& labels,
               const std::vector< std::string >& values,
               bool precr = true ) const
    {
      Assert( labels.size() == values.size(), "Size mismatch" );
      if (!labels.empty()) {
        m_stream << (precr ? "\n" : "") << labels[0] << ": " << values[0];
        for (std::size_t i=1; i<labels.size(); ++i) {
          m_stream << ", " << labels[i] << ": " << values[i];
        }
        m_stream << (precr ? " " : "\n");
      }
    }

    //! Start formatted print of a diagnostics message
    //! Start formatted print of a diagnostics message
    //! \param[in] msg First part of message to print as a diagnostics message
    void diagstart( const std::string& msg ) const {
      m_stream << msg.c_str() << ' ';
    }

    //! Finish formatted print of a diagnostics message
    //! \param[in] msg Last part of message to print as a diagnostics message
    void diagend( const std::string& msg ) const {
      m_stream << msg.c_str() << '\n';
    }

    //! Echo formatted print of a progress message
    //! \param[in] prefix Strings to output prefixing the progress report
    //! \param[in] done Array of integers indicating how many have been done
    //! \param[in] max Array of integers indicating how many to be done
    //! \param[in] progress_size Size of previous progress report (to overwrite)
    //! \details All input arrays are the same size. The prefix strings
    //!   are optional, i.e., they can be empty strings. The function generates
    //!   an output to the stream configured in the following fashion:
    //!   pre1[done1/max1], pre2[done2/max2], ..., e.g., r:[1/3], b[2/8].
    //!   Whenever this function is called, a number of backspaces are put into
    //!   the stream so that the new progress report string overwrites the old
    //!   one. In order to backtrack the correct amount, the length of the old
    //!   progress report is stored (by whatever object holds us) and passed in
    //!   by reference in progress_size, which is overwritten here once it has
    //!   been used for backtracking. Therefore, for restarting a new series of
    //!   progress reports, this variable must be zeroed. Also, it is best to
    //!   not to interleave multiple tasks, because even if a different
    //!   progress_size is kept for each, there is no regard as to which line we
    //!   output to in the stream. In other words, multiple task outputs will
    //!   be intermingled, leading to confusing output.
    template< std::size_t N >
    void progress( const std::array< std::string, N >& prefix,
                   const std::array< int, N >& done,
                   const std::array< int, N >& max,
                   std::size_t& progress_size ) const
    {
      // lambda to determine the number of digits in an integer
      auto numdig = []( int i ) -> std::size_t {
        return i > 0 ?
          static_cast< std::size_t >( std::log10(static_cast<double>(i)) ) + 1
          : 1; };
      // Backspace so that new progress can overwrite old one
      m_stream << std::string( progress_size, '\b' ).c_str();
      std::stringstream ss;
      auto ip = prefix.cbegin();
      auto id = done.cbegin();
      auto im = max.cbegin();
      progress_size = 0;
      while (ip != prefix.cend()) {
        // Compute new length of progress string
        progress_size += 4 + ip->size() + numdig(*id) + numdig(*im);
        // Construct and output new progress string to stream
        ss << *ip << ":[" << *id << '/' << *im << ']';
        ++ip; ++id; ++im;
        // if next subprogress is not the last one, put in a comma
        if (ip != prefix.cend()) {
          ss << ", ";
          progress_size += 2;
        } else {
          ss << ' ';
          ++progress_size;
        }
      }
      m_stream << ss.str().c_str();
    }

    //! Print version information
    //! \param[in] executable Name of executable to output version for
    //! \param[in] git_commit Git commit sha1 to output
    void version( const std::string& executable,
                  const std::string& git_commit ) const {
      m_stream << "\nXyst::" << executable.c_str()
               << ", revision " << git_commit.c_str() << "\n\n";
    }

    //! Print mandatory arguments information
    //! \param[in] args Mandaatory-arguments infor to output
    void mandatory( const std::string& args ) const {
      m_stream << "\n>>> ERROR: " << args.c_str() << '\n';
    }

    //! Print example usage information
    //! \param[in] example Example command line to output
    //! \param[in] msg Message to output after example
    void usage( const std::string& example, const std::string& msg ) const {
      m_stream << "\nUsage: " << example.c_str() << '\n'
               << msg.c_str() << ". See also -h." << "\n\n";
    }

    //! Print lower and upper bounds for a keyword if defined
    template< typename Info >
    void bounds( const Info& info ) const {
      if (info.lower)
        m_stream << splitLines( *info.lower, "  ", "Lower bound: " ).c_str();
      if (info.upper)
        m_stream << splitLines( *info.upper, "  ", "Upper bound: " ).c_str();
    }

    //! Print unit tests header (with legend)
    //! \param[in] t Section title
    //! \param[in] group String attempting to match unit test groups
    void unithead( const std::string& t, const std::string& group ) const {
      section( t );
      m_stream << "Groups: " + (group.empty() ? "all" : group) +
                  " (use -g str to match groups)\n" +
                  "Legend: [done/failed] group/test: result\n\n";
    }

    //! Print one-liner info for test
    //! \details Columns:
    //!   [done/failed]
    //!   - done: number of tests completed so far
    //!   - failed: number of failed tests so far
    //!   name of the test group
    //!   name of the test
    //!   result (with additional info if failed)
    //!   Assumed fields for status:
    //!   - status[0]: test group name
    //!   - status[1]: test name
    //!   - status[2]: result (tut::test_result::result_type as string)
    //!   - status[3]: exception message for failed test
    //!   - status[4]: exception type id for failed test
    void test( std::size_t ncomplete,
               std::size_t nfail,
               const std::vector< std::string >& status )
    {
      if (status[2] != "8") {             // if not dummy
        std::stringstream ss;
        ss << "[" << ncomplete << "/" << nfail << "] " << status[0] << ":"
           << status[1];
        m_stream << ss.str() << ": "
                 << result( status[2], status[3], status[4] ) << '\n';
      }
    }

    //! Print Inciter header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag.
    void headerInciter() const {
       m_stream << R"(
____  ___                __    __   .___              .__  __                
\   \/  /___.__. _______/  |_  \ \  |   | ____   ____ |__|/  |_  ___________ 
 \     /<   |  |/  ___/\   __\  \ \ |   |/    \_/ ___\|  \   __\/ __ \_  __ \
 /     \ \___  |\___ \  |  |    / / |   |   |  \  \___|  ||  | \  ___/|  | \/
/___/\  \/ ____/____  > |__|   /_/  |___|___|  /\___  >__||__|  \___  >__|   
      \_/\/         \/                       \/     \/              \/)"
      << '\n';
    }

    //! Print UnitTest header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag.
    void headerUnitTest() const {
       m_stream << R"(
____  ___                __    __    ____ ___      .__  __ ___________              __   
\   \/  /___.__. _______/  |_  \ \  |    |   \____ |__|/  |\__    ___/___   _______/  |_ 
 \     /<   |  |/  ___/\   __\  \ \ |    |   /    \|  \   __\|    |_/ __ \ /  ___/\   __\
 /     \ \___  |\___ \  |  |    / / |    |  /   |  \  ||  |  |    |\  ___/ \___ \  |  |  
/___/\  \/ ____/____  > |__|   /_/  |______/|___|  /__||__|  |____| \___  >____  > |__|  
      \_/\/         \/                           \/                     \/     \/)"
      << '\n';
    }

    //! Print MeshConv header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag.
    void headerMeshConv() const {
      m_stream << R"(
____  ___                __    __      _____                .__    _________                     
\   \/  /___.__. _______/  |_  \ \    /     \   ____   _____|  |__ \_   ___ \  ____   _______  __
 \     /<   |  |/  ___/\   __\  \ \  /  \ /  \_/ __ \ /  ___/  |  \/    \  \/ /  _ \ /    \  \/ /
 /     \ \___  |\___ \  |  |    / / /    Y    \  ___/ \___ \|   Y  \     \___(  <_> )   |  \   / 
/___/\  \/ ____/____  > |__|   /_/  \____|__  /\___  >____  >___|  /\______  /\____/|___|  /\_/  
      \_/\/         \/                      \/     \/     \/     \/        \/            \/)"
      << '\n';
    }

  private:
    std::ostream& m_stream;     //!< Output stream

    //! Return human-readable test result based on result code
    //! \param[in] code Result code
    //! \param[in] msg Message to append
    //! \param[in] ex Expection message to attach to exceptions cases
    std::string result( const std::string& code,
                        const std::string& msg,
                        const std::string& ex ) const
    {
      if (code == "0") return "ok";
      else if (code == "1") return "fail: " + msg;
      else if (code == "2") return "except: " + msg + ex;
      else if (code == "3") return "warning: " + msg;
      else if (code == "4") return "terminate: " + msg;
      else if (code == "5") return "ex_ctor: " + msg + ex;
      else if (code == "6") return "rethrown: " + msg + ex;
      else if (code == "7") return "skipped: " + msg;
      else if (code == "8") return "dummy";
      else Throw( "No such unit test result code found" );
    }
};

} // tk::

#endif // Print_h
