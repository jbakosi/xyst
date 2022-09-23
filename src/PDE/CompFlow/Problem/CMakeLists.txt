add_library(CompFlowProblem
            BoxInitialization.cpp
            FieldOutput.cpp
            NLEnergyGrowth.cpp
            RayleighTaylor.cpp
            SodShocktube.cpp
            RotatedSodShocktube.cpp
            SedovBlastwave.cpp
            TaylorGreen.cpp
            UserDefined.cpp
            VorticalFlow.cpp
)

target_include_directories(CompFlowProblem PUBLIC
                           ${PROJECT_SOURCE_DIR}/..
                           ${PROJECT_SOURCE_DIR}/../Base
                           ${PROJECT_SOURCE_DIR}/../Control
                           ${PROJECT_SOURCE_DIR}/../Control/Inciter
                           ${PROJECT_SOURCE_DIR}/../Mesh
                           ${PROJECT_SOURCE_DIR}/../PDE
                           ${PROJECT_SOURCE_DIR}/../PDE/CompFlow
                           ${PROJECT_SOURCE_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS})

set_target_properties(CompFlowProblem
                      PROPERTIES LIBRARY_OUTPUT_NAME xyst_compflowproblem)

INSTALL(TARGETS CompFlowProblem
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
