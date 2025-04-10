// *****************************************************************************
/*!
  \file      src/Base/Timer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Timer declaration
  \details   Timer declaration. Timer is a simple class to do timing various
    parts of the code in a portable way. The functionality is intended to be
    very minimal and simple, but still convenient to use, with as little state
    as possible. For an example client code, see walker::Main in Main/Walker.C.
*/
// *****************************************************************************
#ifndef Timer_h
#define Timer_h

#include <cstdint>
#include <chrono>

#include "NoWarning/pup.hpp"      // for er
#include "Types.hpp"

namespace tk {

//! \brief Simple class to do timing
//! \details Timing use at various parts of the code in a portable way. The
//!   functionality is intended to be very minimal and simple, but still
//!   convenient to use, with as little state as possible.
class Timer {

  public:
    // Shorthands for seconds duration, hours, minutes, seconds
    using Dsec = std::chrono::duration< real >;
    using hours = std::chrono::hours;
    using minutes = std::chrono::minutes;
    using seconds = std::chrono::seconds;
    // Shorthand for clock, setting clock type
    using clock = std::chrono::high_resolution_clock;

    //! Watch stores time in hours:minutes:seconds
    struct Watch {
      hours hrs;
      minutes min;
      seconds sec;
      //! Zero constructor. Zeros hours, minutes, and seconds
      explicit Watch() :
        hrs( std::chrono::duration_cast< hours >( clock::duration::zero()) ),
        min( std::chrono::duration_cast< minutes >( clock::duration::zero() ) ),
        sec( std::chrono::duration_cast< seconds >( clock::duration::zero() ) )
      {}
      //! Fill constructor. Initialize hours, minutes, and seconds given
      explicit Watch( hours&& h, minutes&& m, seconds&& s ) :
        hrs( std::move(h) ), min( std::move(m) ), sec( std::move(s) ) {}
    };

    //! Constructor: initialize clock to current time stamp
    explicit Timer() : m_start( clock::now() ), m_prev( m_start ) {}

    //! Zero timer
    void zero() { m_start = clock::now(); m_prev = m_start; }

    //! Query time in second since the constructor call
    //! \return Time elapsed between start and stop as a real number
    real dsec() const {
      return std::chrono::duration_cast< Dsec >(clock::now() - m_start).count();
    }

    //! Query time in second since the constructor call
    Watch hms() const;

    //! Estimate time for accomplishment
    void eta( real term, real time, uint64_t nstep, uint64_t it,
              tk::real res0, tk::real res, tk::real rest,
              Watch& elapsedWatch, Watch& estimatedWatch );

    /** @name Pack/Unpack: Serialize Timer object for Charm++ */
    ///@{
    //! Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p( reinterpret_cast<char*>(&m_start), sizeof(clock::time_point) );
      p( reinterpret_cast<char*>(&m_prev), sizeof(clock::time_point) );
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] t Timer object reference
    friend void operator|( PUP::er& p, Timer& t ) { t.pup(p); } 
    ///@}

  private:
    clock::time_point m_start;  //!< Time stamp at Timer() or zero()
    clock::time_point m_prev;   //!< Time stamp at previous update in eta()
};

//! Convert existing time stamp as a real to Watch (global scope)
Timer::Watch
hms( tk::real stamp );

} // tk::

#endif // Timer_h
