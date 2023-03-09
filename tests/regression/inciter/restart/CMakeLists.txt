# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial, checkpoint on 1 PE, restart on 1 PE

add_regression_test(riecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_riecg.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg.q -i unitcube_1k.exo -v
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
                    ARGS -c vortical_flow_riecg_restart.q -i unitcube_1k.exo -v
                         +restart restart
                    BIN_BASELINE vortical_flow_riecg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE vortical_flow_riecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    CHECKPOINT riecg_checkpoint
                    LABELS riecg restart)

# Parallel, checkpoint on 4 PEs, restart on 4 PEs

add_regression_test(riecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_riecg.q unitcube_1k.exo
                    ARGS -c vortical_flow_riecg.q -i unitcube_1k.exo -v
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
                    ARGS -c vortical_flow_riecg_restart.q -i unitcube_1k.exo -v
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
