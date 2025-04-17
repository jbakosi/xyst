# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(chocg_slotcyl ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl.q -i unitsquare_01_3.6k.exo
                    BIN_BASELINE slot_cyl.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_slotcyl_damp4 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_damp4.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_damp4.q -i unitsquare_01_3.6k.exo
                    BIN_BASELINE slot_cyl_damp4.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_damp4.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)

add_regression_test(chocg_slotcyl_damp4_freeze ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_damp4_freeze.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_damp4_freeze.q -i unitsquare_01_3.6k.exo
                    BIN_BASELINE slot_cyl_damp4_freeze.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_damp4_freeze.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg)


# Parallel + virtualization

add_regression_test(chocg_slotcyl_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitsquare_01_3.6k.exo
                    ARGS -f -c slot_cyl.q -i unitsquare_01_3.6k.exo -u 0.5
                    BIN_BASELINE slot_cyl_pe4_u0.5.std.exo.0
                                 slot_cyl_pe4_u0.5.std.exo.1
                                 slot_cyl_pe4_u0.5.std.exo.2
                                 slot_cyl_pe4_u0.5.std.exo.3
                                 slot_cyl_pe4_u0.5.std.exo.4
                                 slot_cyl_pe4_u0.5.std.exo.5
                                 slot_cyl_pe4_u0.5.std.exo.6
                                 slot_cyl_pe4_u0.5.std.exo.7
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
                    LABELS chocg)

# Parallel + virtualization + migration

add_regression_test(chocg_slotcyl_u0.9_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl.q -i unitsquare_01_3.6k.exo -u 0.9
                         -l 10 +balancer RecBipartLB +LBDebug 1 +cs
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS chocg migration)
