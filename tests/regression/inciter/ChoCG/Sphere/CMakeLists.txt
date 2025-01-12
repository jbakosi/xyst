# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(chocg_elemblock_unsupported ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS -i sphere_bar.exo -c sphere.q
                    INPUTFILES sphere_bar.exo sphere.q
                    EXTRA_PASS_REGEXP "CALL TRACE"
                    SKIP_RETURN_CODE true
                    LABELS chocg)

add_regression_test(chocg_inviscid_sphere ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES inviscid_sphere.q sphere2_5K.exo
                    ARGS -c inviscid_sphere.q -i sphere2_5K.exo
                    BIN_BASELINE inviscid_sphere.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_inviscid_sphere.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_viscous_sphere ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sphere_chocg_viscous_test.q sphere2_5K.exo
                    ARGS -c sphere_chocg_viscous_test.q -i sphere2_5K.exo
                    BIN_BASELINE sphere_chocg_viscous_test.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_sphere_chocg_viscous_test.std
                                  inte_sphere_chocg_viscous_test.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS chocg)

# Parallel + virtualization

add_regression_test(chocg_inviscid_sphere_u0.6
                    ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES inviscid_sphere.q sphere2_5K.exo
                    ARGS -c inviscid_sphere.q -i sphere2_5K.exo -u 0.6
                    TEXT_BASELINE diag_inviscid_sphere.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)
