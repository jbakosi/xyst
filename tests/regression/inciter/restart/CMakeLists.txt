# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial, checkpoint on 1 PE, restart on 1 PE

add_regression_test(riecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_riecg.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_riecg_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE vortical_flow_riecg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg restart)

add_regression_test(riecg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_riecg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_riecg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE vortical_flow_riecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    CHECKPOINT riecg_checkpoint
                    LABELS riecg restart)

add_regression_test(kozcg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_kozcg.q unitcube_1k.exo
                    ARGS -c vortical_flow_kozcg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_kozcg_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_kozcg.cfg
                    TEXT_BASELINE vortical_flow_kozcg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.kozcg.cfg
                    LABELS kozcg restart)

add_regression_test(kozcg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_kozcg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_kozcg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_kozcg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_kozcg.cfg
                    TEXT_BASELINE vortical_flow_kozcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.kozcg.cfg
                    CHECKPOINT kozcg_checkpoint
                    LABELS kozcg restart)

add_regression_test(zalcg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_zalcg.q unitcube_1k.exo
                    ARGS -c vortical_flow_zalcg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_zalcg_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_zalcg.cfg
                    TEXT_BASELINE vortical_flow_zalcg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.zalcg.cfg
                    LABELS zalcg restart)

add_regression_test(zalcg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_zalcg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_zalcg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_zalcg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_zalcg.cfg
                    TEXT_BASELINE vortical_flow_zalcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.zalcg.cfg
                    CHECKPOINT zalcg_checkpoint
                    LABELS zalcg restart)

add_regression_test(chocg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ldc_chocg.q unitcube_1k.exo
                    ARGS -c ldc_chocg.q -i unitcube_1k.exo
                    BIN_BASELINE ldc_chocg_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_chocg.cfg
                    TEXT_BASELINE ldc_chocg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.chocg.cfg
                    LABELS chocg restart)

add_regression_test(chocg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ldc_chocg_restart.q unitcube_1k.exo
                    ARGS -c ldc_chocg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE ldc_chocg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_chocg.cfg
                    TEXT_BASELINE ldc_chocg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.chocg.cfg
                    CHECKPOINT chocg_checkpoint
                    LABELS chocg restart)

# Parallel, checkpoint on 4 PEs, restart on 4 PEs

add_regression_test(riecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_riecg.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_riecg_checkpoint.std.exo.0
                                 vortical_flow_riecg_checkpoint.std.exo.1
                                 vortical_flow_riecg_checkpoint.std.exo.2
                                 vortical_flow_riecg_checkpoint.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_riecg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg restart)

add_regression_test(riecg_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_riecg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_riecg.std.exo.0
                                 vortical_flow_riecg.std.exo.1
                                 vortical_flow_riecg.std.exo.2
                                 vortical_flow_riecg.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_riecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    CHECKPOINT riecg_checkpoint
                    LABELS riecg restart)

add_regression_test(kozcg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_kozcg.q unitcube_1k.exo
                    ARGS -c vortical_flow_kozcg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_kozcg_checkpoint.std.exo.0
                                 vortical_flow_kozcg_checkpoint.std.exo.1
                                 vortical_flow_kozcg_checkpoint.std.exo.2
                                 vortical_flow_kozcg_checkpoint.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_kozcg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_kozcg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.kozcg.cfg
                    LABELS kozcg restart)

add_regression_test(kozcg_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_kozcg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_kozcg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_kozcg.std.exo.0
                                 vortical_flow_kozcg.std.exo.1
                                 vortical_flow_kozcg.std.exo.2
                                 vortical_flow_kozcg.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_kozcg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_kozcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.kozcg.cfg
                    CHECKPOINT kozcg_checkpoint
                    LABELS kozcg restart)

add_regression_test(zalcg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_zalcg.q unitcube_1k.exo
                    ARGS -c vortical_flow_zalcg.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_zalcg_checkpoint.std.exo.0
                                 vortical_flow_zalcg_checkpoint.std.exo.1
                                 vortical_flow_zalcg_checkpoint.std.exo.2
                                 vortical_flow_zalcg_checkpoint.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_zalcg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_zalcg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.zalcg.cfg
                    LABELS zalcg restart)

add_regression_test(zalcg_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_zalcg_restart.q unitcube_1k.exo
                    ARGS -c vortical_flow_zalcg_restart.q -i unitcube_1k.exo
                         +restart restart
                    BIN_BASELINE vortical_flow_zalcg.std.exo.0
                                 vortical_flow_zalcg.std.exo.1
                                 vortical_flow_zalcg.std.exo.2
                                 vortical_flow_zalcg.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_zalcg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE vortical_flow_zalcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.zalcg.cfg
                    CHECKPOINT zalcg_checkpoint
                    LABELS zalcg restart)
