# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(chocg_poisson_const ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES poisson_const.q unitcube_01_1k.exo
                    ARGS -c poisson_const.q -i unitcube_01_1k.exo
                    BIN_BASELINE poisson_const.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_poisson_const.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

# Parallel

add_regression_test(chocg_poisson_const
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poisson_const.q unitcube_01_1k.exo
                    ARGS -c poisson_const.q -i unitcube_01_1k.exo
                    BIN_BASELINE poisson_const_pe2.std.exo.0
                                 poisson_const_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poisson_const.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_poisson_sine
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poisson_sine.q unitcube_01_1k.exo
                    ARGS -c poisson_sine.q -i unitcube_01_1k.exo
                    BIN_BASELINE poisson_sine_pe2.std.exo.0
                                 poisson_sine_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poisson_sine.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_poisson_sine3
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poisson_sine3.q unitcube_01_1k.exo
                    ARGS -c poisson_sine3.q -i unitcube_01_1k.exo
                    BIN_BASELINE poisson_sine3_pe2.std.exo.0
                                 poisson_sine3_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poisson_sine3.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_poisson_neumann
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poisson_neumann.q unitcube_0pidiv4_1k.exo
                    ARGS -c poisson_neumann.q -i unitcube_0pidiv4_1k.exo
                    BIN_BASELINE poisson_neumann_pe2.std.exo.0
                                 poisson_neumann_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poisson_neumann.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

# Parallel + virtualization

add_regression_test(chocg_poisson_const
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 4
                    INPUTFILES poisson_const_phg.q unitcube_01_1k.exo
                    ARGS -c poisson_const_phg.q -i unitcube_01_1k.exo -u 0.8
                    TEXT_BASELINE diag_poisson_const.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)
