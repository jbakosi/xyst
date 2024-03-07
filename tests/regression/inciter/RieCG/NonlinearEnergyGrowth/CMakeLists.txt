# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(riecg_nleg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES nleg.q unitcube_1k.exo
                    ARGS -c nleg.q -i unitcube_1k.exo
                    BIN_BASELINE nleg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

# Parallel + no virtualization

add_regression_test(riecg_nleg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg.q unitcube_1k.exo
                    ARGS -c nleg.q -i unitcube_1k.exo
                    BIN_BASELINE nleg_pe4.std.exo.0
                                 nleg_pe4.std.exo.1
                                 nleg_pe4.std.exo.2
                                 nleg_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

# Parallel + virtualization

add_regression_test(riecg_nleg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg.q unitcube_1k.exo
                    ARGS -c nleg.q -i unitcube_1k.exo -u 0.5
                    BIN_BASELINE nleg_pe4_u0.5.std.exo.0
                                 nleg_pe4_u0.5.std.exo.1
                                 nleg_pe4_u0.5.std.exo.2
                                 nleg_pe4_u0.5.std.exo.3
                                 nleg_pe4_u0.5.std.exo.4
                                 nleg_pe4_u0.5.std.exo.5
                                 nleg_pe4_u0.5.std.exo.6
                                 nleg_pe4_u0.5.std.exo.7
                    BIN_RESULT out.e-s.0.8.0
                               out.e-s.0.8.1
                               out.e-s.0.8.2
                               out.e-s.0.8.3
                               out.e-s.0.8.4
                               out.e-s.0.8.5
                               out.e-s.0.8.6
                               out.e-s.0.8.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)
