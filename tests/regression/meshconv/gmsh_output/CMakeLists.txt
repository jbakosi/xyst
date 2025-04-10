# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(asc2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse.msh
                    POSTPROCESS_PROG ${GMSH_EXECUTABLE}
                    POSTPROCESS_PROG_ARGS -check cube_coarse.msh
                                  COMMAND grep -v binary
                                  COMMAND grep -v dense
                                  COMMAND grep -vi Reading
                                  COMMAND grep -v current
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE cube_coarse.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)

add_regression_test(exo2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.exo
                    ARGS -i box_24.exo -o box_24.msh
                    POSTPROCESS_PROG ${GMSH_EXECUTABLE}
                    POSTPROCESS_PROG_ARGS -check box_24.msh
                                  COMMAND grep -v binary
                                  COMMAND grep -v dense
                                  COMMAND grep -vi Reading
                                  COMMAND grep -v current
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_box_24.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)

add_regression_test(gmshtxt2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh
                    ARGS -i box_24.txt.msh -o box_24.msh
                    POSTPROCESS_PROG ${GMSH_EXECUTABLE}
                    POSTPROCESS_PROG_ARGS -check box_24.msh
                                  COMMAND grep -v binary
                                  COMMAND grep -v dense
                                  COMMAND grep -vi Reading
                                  COMMAND grep -v current
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_box_24.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)

add_regression_test(netgen2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES torus.mesh
                    ARGS -i torus.mesh -o torus.msh
                    POSTPROCESS_PROG ${GMSH_EXECUTABLE}
                    POSTPROCESS_PROG_ARGS -check torus.msh
                                  COMMAND grep -v binary
                                  COMMAND grep -v dense
                                  COMMAND grep -vi Reading
                                  COMMAND grep -v current
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_torus.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)
