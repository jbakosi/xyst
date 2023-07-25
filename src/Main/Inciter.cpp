// *****************************************************************************
/*!
  \file      src/Main/Inciter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <unordered_map>
#include <vector>
#include <iostream>

#include "XystBuildConfig.hpp"
#ifdef XYST_AMPI
  #include "NoWarning/mpi.hpp"
#endif

#include "Types.hpp"
#include "Init.hpp"
#include "Timer.hpp"
#include "Exception.hpp"
#include "ProcessException.hpp"
#include "InciterConfig.hpp"
#include "LBSwitch.hpp"

#include "NoWarning/inciter.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! Load balancer switch group proxy
tk::CProxy_LBSwitch LBSwitchProxy;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! Inciter declarations and definitions
namespace inciter {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! Configuration data structure, containing all input data
//! \details This object is in global scope, it contains all of user input, and
//!   thus it is made available to all PEs for convenience reasons. The runtime
//!   system distributes it to all PEs during initialization. Once distributed,
//!   the object does not change.
ctr::Config g_cfg;

//! Number of times restarted counter
int g_nrestart;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // inciter::

#ifdef XYST_AMPI
//! Main function for Charm++'s AMPI
//! \note If this is not defined, Charm++ does not wait for CkExit().
int main( int, char** ) {
  MPI_Init( nullptr, nullptr );
  return 0;
}
#endif

//! \brief Charm++ main chare for the shock hydrodynamics executable, inciter.
//! \details In inciter the Charm++ runtime system is initialized only after the
//!   mesh has been read in, partitioned, and the necessary data structures,
//!   e.g., communication maps, have been generated. This delayed initialization
//!   of the Charm++ runtime system is required since the mesh partitioning is
//!   done by Zoltan, an MPI library. Note that this Charm++ main chare object
//!   should not be in a namespace.
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details Inciter's main chare constructor is the entry point of the
    //!   Charm++ portion of inciter, called by the Charm++ runtime system. The
    //!   constructor does basic initialization steps, prints out some useful
    //!   information to screen (in verbose mode), and instantiates a driver.
    //!   Since Charm++ is fully asynchronous, the constructor usually spawns
    //!   asynchronous objects and immediately exits. Thus in the body of the
    //!   main chare constructor we fire up an 'execute' chare, which then calls
    //!   back to Main::execute(). Finishing the main chare constructor the
    //!   Charm++ runtime system then starts the network-migration of all
    //!   global-scope data (if any). The execute chare calling back to
    //!   Main::execute() signals the end of the migration of the global-scope
    //!   data. Then we are ready to execute the driver. Since inciter is
    //!   parallel and asynchronous, its driver fires up additional Charm++
    //!   chare objects which then call back to Main::finalize() at some point
    //!   in the future when all work has been finished. finalize() then exits
    //!   by calling Charm++'s CkExit(), shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    explicit Main( CkArgMsg* msg )
    try :
      m_timer(1)
    {
      tk::setSignalHandlers();
      using inciter::g_cfg;
      // Parse command line
      g_cfg.cmdline( msg->argc, msg->argv );
      tk::echoHeader( tk::HeaderType::INCITER );
      tk::echoBuildEnv( msg->argv[0] );
      tk::echoRunEnv(msg->argc, msg->argv, g_cfg.get< tag::quiescence >());
      delete msg;
      mainProxy = thisProxy;
      if (g_cfg.get< tag::quiescence >()) {
        CkStartQD( CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
      }
      m_timer.emplace_back();
      // Parse control file
      inciter::g_cfg.control();
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
    } catch (...) { tk::processExceptionCharm(); }

    //! Migrate constructor: returning from a checkpoint
    explicit Main( CkMigrateMessage* msg )
    try :
      CBase_Main( msg ), m_timer(1)
    {
      tk::setSignalHandlers();
      using inciter::g_cfg;
      // Parse command line after restart
      g_cfg.cmdline( reinterpret_cast<CkArgMsg*>(msg)->argc,
                     reinterpret_cast<CkArgMsg*>(msg)->argv );
      tk::echoHeader( tk::HeaderType::INCITER );
      tk::echoBuildEnv( reinterpret_cast<CkArgMsg*>(msg)->argv[0] );
      tk::echoRunEnv( reinterpret_cast<CkArgMsg*>(msg)->argc,
                      reinterpret_cast<CkArgMsg*>(msg)->argv,
                      g_cfg.get< tag::quiescence >() );
      // increase number of restarts (available for Transporter on PE 0)
      ++inciter::g_nrestart;
      mainProxy = thisProxy;
      if (g_cfg.get< tag::quiescence >()) {
        CkStartQD( CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
      }
      m_timer.emplace_back();
      // Parse control file after restart
      g_cfg.control();
    } catch (...) { tk::processExceptionCharm(); }

    //! Execute driver created and initialized by constructor
    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        // Instantiate Transporter chare on PE 0 which drives time-integration
        inciter::CProxy_Transporter::ckNew( 0 );
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Towards normal exit but collect chare state first (if any)
    void finalize() {
      tk::finalize( m_timer, m_timestamp );
    }

    //! Entry method triggered when quiescence is detected
    [[noreturn]] void quiescence() { Throw( "Quiescence detected" ); }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ mainchare, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_timer;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] m Mainchare object reference
    friend void operator|( PUP::er& p, Main& m ) { m.pup(p); }
    //@}

  private:
    std::vector< tk::Timer > m_timer;           //!< Timers
    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
class execute : public CBase_execute {
  public:
    //! Constructor
    explicit execute() { mainProxy.execute(); }
    //! Migrate constructor
    explicit execute( CkMigrateMessage* m ) : CBase_execute( m ) {}
};

#include "NoWarning/inciter.def.h"
