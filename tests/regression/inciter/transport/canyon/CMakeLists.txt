# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel

add_regression_test(transport_canyon ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES canyon.q canyon.exo
                    ARGS -c canyon.q -i canyon.exo -v
                    BIN_BASELINE canyon_pe2.std.exo.0
                                 canyon_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(transport_canyon_farfield ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES canyon_farfield.q canyon.exo
                    ARGS -c canyon_farfield.q -i canyon.exo -v
                    BIN_BASELINE canyon_farfield_pe2.std.exo.0
                                 canyon_farfield_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_farfield.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(transport_canyon_histwrongsize ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES canyon_wrongsize.q canyon.exo
                    ARGS -c canyon_wrongsize.q -i canyon.exo -v
                    EXTRA_PASS_REGEXP "size of the parameter vector is incorrect"
                    SKIP_RETURN_CODE true
                    LABELS riecg)
