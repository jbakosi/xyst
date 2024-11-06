# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(lohcg_poiseuille_damp2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES poiseuille_damp2.q poiseuille1tetz.exo
                    ARGS -c poiseuille_damp2.q -i poiseuille1tetz.exo
                    BIN_BASELINE poiseuille_damp2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_poiseuille_damp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg)

add_regression_test(lohcg_poiseuille_damp4 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES poiseuille_damp4.q poiseuille1tetz.exo
                    ARGS -c poiseuille_damp4.q -i poiseuille1tetz.exo
                    BIN_BASELINE poiseuille_damp4.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_poiseuille_damp4.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg)

add_regression_test(lohcg_poiseuille_rk4 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES poiseuille_rk4.q poiseuille1tetz.exo
                    ARGS -c poiseuille_rk4.q -i poiseuille1tetz.exo
                    BIN_BASELINE poiseuille_rk4.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_poiseuille_rk4.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg)

# Parallel

add_regression_test(lohcg_poiseuille
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poiseuille_damp2.q poiseuille1tetz.exo
                    ARGS -c poiseuille_damp2.q -i poiseuille1tetz.exo
                    BIN_BASELINE poiseuille_pe2.std.exo.0
                                 poiseuille_pe2.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poiseuille_damp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg)

# Parallel + virtualization

add_regression_test(lohcg_poiseuille_u0.8
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poiseuille_damp2.q poiseuille1tetz.exo
                    ARGS -c poiseuille_damp2.q -i poiseuille1tetz.exo -u 0.8
                    BIN_BASELINE poiseuille_pe2_u0.8.std.exo.0
                                 poiseuille_pe2_u0.8.std.exo.1
                                 poiseuille_pe2_u0.8.std.exo.2
                                 poiseuille_pe2_u0.8.std.exo.3
                                 poiseuille_pe2_u0.8.std.exo.4
                                 poiseuille_pe2_u0.8.std.exo.5
                                 poiseuille_pe2_u0.8.std.exo.6
                                 poiseuille_pe2_u0.8.std.exo.7
                                 poiseuille_pe2_u0.8.std.exo.8
                                 poiseuille_pe2_u0.8.std.exo.9
                    BIN_RESULT out.e-s.0.10.0
                               out.e-s.0.10.1
                               out.e-s.0.10.2
                               out.e-s.0.10.3
                               out.e-s.0.10.4
                               out.e-s.0.10.5
                               out.e-s.0.10.6
                               out.e-s.0.10.7
                               out.e-s.0.10.8
                               out.e-s.0.10.9
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poiseuille_damp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg)

# Parallel + virtualization + migration

add_regression_test(lohcg_poiseuille_u0.8_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES poiseuille_damp2.q poiseuille1tetz.exo
                    ARGS -c poiseuille_damp2.q -i poiseuille1tetz.exo -u 0.8
                         +balancer RecBipartLB +LBDebug 1 -l 4
                    BIN_BASELINE poiseuille_pe2_u0.8.std.exo.0
                                 poiseuille_pe2_u0.8.std.exo.1
                                 poiseuille_pe2_u0.8.std.exo.2
                                 poiseuille_pe2_u0.8.std.exo.3
                                 poiseuille_pe2_u0.8.std.exo.4
                                 poiseuille_pe2_u0.8.std.exo.5
                                 poiseuille_pe2_u0.8.std.exo.6
                                 poiseuille_pe2_u0.8.std.exo.7
                                 poiseuille_pe2_u0.8.std.exo.8
                                 poiseuille_pe2_u0.8.std.exo.9
                    BIN_RESULT out.e-s.0.10.0
                               out.e-s.0.10.1
                               out.e-s.0.10.2
                               out.e-s.0.10.3
                               out.e-s.0.10.4
                               out.e-s.0.10.5
                               out.e-s.0.10.6
                               out.e-s.0.10.7
                               out.e-s.0.10.8
                               out.e-s.0.10.9
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_poiseuille_damp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS lohcg migration)
