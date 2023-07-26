# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_pipe ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES pipe.q rectangle_01_1.5k.exo
                    ARGS -c pipe.q -i rectangle_01_1.5k.exo
                    BIN_BASELINE pipe.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                                  pipe_int.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS riecg)

# Parallel

add_regression_test(compflow_euler_pipe ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES pipe.q rectangle_01_1.5k.exo
                    ARGS -c pipe.q -i rectangle_01_1.5k.exo
                    BIN_BASELINE pipe_pe2.std.exo.0
                                 pipe_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                                  pipe_int.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_pipe_int_timeint ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES pipe_int_timeint.q rectangle_01_1.5k.exo
                    ARGS -c pipe_int_timeint.q -i rectangle_01_1.5k.exo
                    TEXT_BASELINE diag.std
                                  pipe_int_timeint.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_pipe_int_timerange ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES pipe_int_timerange.q rectangle_01_1.5k.exo
                    ARGS -c pipe_int_timerange.q -i rectangle_01_1.5k.exo
                    TEXT_BASELINE diag.std
                                  pipe_int_timerange.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS riecg)

# Parallel + virtualization + migration

add_regression_test(compflow_euler_pipe_u0.5_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES pipe.q rectangle_01_1.5k.exo
                    ARGS -c pipe.q -i rectangle_01_1.5k.exo -u 0.5 -l 3
                         +balancer RecBipartLB +LBDebug 1
                    TEXT_BASELINE diag.std
                                  pipe_int.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS riecg)
