# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(chocg_ldc ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ldc.q unitcube_1k.exo
                    ARGS -c ldc.q -i unitcube_1k.exo
                    BIN_BASELINE ldc.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_ldc.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)
