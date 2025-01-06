# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(lohcg_viscous_sphere ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sphere_lohcg_viscous_test.q sphere2_5K.exo
                    ARGS -c sphere_lohcg_viscous_test.q -i sphere2_5K.exo
                    BIN_BASELINE sphere_lohcg_viscous_test.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_sphere_lohcg_viscous_test.std
                                  inte_sphere_lohcg_viscous_test.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS lohcg)

# Parallel + virtualization

add_regression_test(lohcg_viscous_sphere_u0.6 ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES sphere_lohcg_viscous_test.q sphere2_5K.exo
                    ARGS -c sphere_lohcg_viscous_test.q -i sphere2_5K.exo -u 0.6
                    TEXT_BASELINE diag_sphere_lohcg_viscous_test.std
                                  inte_sphere_lohcg_viscous_test.std
                    TEXT_RESULT diag
                                out.int
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        int.ndiff.cfg
                    LABELS lohcg)
