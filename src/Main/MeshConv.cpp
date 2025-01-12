// *****************************************************************************
/*!
  \file      src/Main/MeshConv.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh file converter Charm++ main chare
  \details   Mesh file converter Charm++ main chare. This file contains the
    definition of the Charm++ main chare, equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <vector>
#include <utility>
#include <iostream>

#include "Timer.hpp"
#include "Types.hpp"
#include "XystConfig.hpp"
#include "Init.hpp"
#include "MeshConvDriver.hpp"
#include "MeshConvConfig.hpp"
#include "ProcessException.hpp"

#include "NoWarning/charm.hpp"
#include "NoWarning/meshconv.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! \brief Charm++ main chare for the mesh converter executable, meshconv.
//! \details Note that this object should not be in a namespace.
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details MeshConv's main chare constructor is the entry point of the
    //!   program, called by the Charm++ runtime system. The constructor does
    //!   basic initialization steps, e.g., parser the command-line, prints out
    //!   some useful information to screen (in verbose mode), and instantiates
    //!   a driver. Since Charm++ is fully asynchronous, the constructor
    //!   usually spawns asynchronous objects and immediately exits. Thus in the
    //!   body of the main chare constructor we fire up an 'execute' chare,
    //!   which then calls back to Main::execute(). Finishing the main chare
    //!   constructor the Charm++ runtime system then starts the
    //!   network-migration of all global-scope data (if any). The execute chare
    //!   calling back to Main::execute() signals the end of the migration of
    //!   the global-scope data. Then we are ready to execute the driver which
    //!   calls back to Main::finalize() when it finished. Then finalize() exits
    //!   by calling Charm++'s CkExit(), shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    explicit Main( CkArgMsg* msg )
    try :
      m_timer(1)
    {
      tk::setSignalHandlers();
      // Parse command line
      m_cfg.cmdline( msg->argc, msg->argv );
      tk::echoHeader( tk::HeaderType::MESHCONV );
      tk::echoBuildEnv( msg->argv[0] );
      tk::echoRunEnv(msg->argc, msg->argv, m_cfg.get< tag::quiescence >());
      delete msg;
      mainProxy = thisProxy;
      if (m_cfg.get< tag::quiescence >()) {
        CkStartQD( CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
      }
      m_timer.emplace_back();
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        meshconv::MeshConvDriver().convert( m_cfg.get< tag::input >(),
          m_cfg.get< tag::output >(), m_cfg.get< tag::reorder >() );
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Towards normal exit but collect chare state first (if any)
    void finalize() {
      tk::finalize( m_timer, m_timestamp );
    }

    //! Add a time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp ) {
      try {
        m_timestamp.emplace_back( label, tk::hms( stamp ) );
      } catch (...) { tk::processExceptionCharm(); }
    }
    //! Add multiple time stamps contributing to final timers output
    void timestamp( const std::vector< std::pair< std::string, tk::real > >& s )
    { for (const auto& t : s) timestamp( t.first, t.second ); }

    //! Entry method triggered when quiescence is detected
    [[noreturn]] void quiescence() { Throw( "Quiescence detected" ); }

  private:
    meshconv::ctr::Config m_cfg;                //!< Config parsed from cmdline
    std::vector< tk::Timer > m_timer;           //!< Timers
    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
class execute : public CBase_execute {
  public: explicit execute() { mainProxy.execute(); }
};

#include "NoWarning/meshconv.def.h"
