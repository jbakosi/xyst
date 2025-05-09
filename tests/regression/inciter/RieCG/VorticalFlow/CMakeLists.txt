# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(riecg_vorticalflow ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow.q unitcube_1k.exo
                    ARGS -c vortical_flow.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(riecg_vorticalflow_stab2
                    ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_stab2.q unitcube_1k.exo
                    ARGS -c vortical_flow_stab2.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_stab2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_stab2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(riecg_vorticalflow_steady
                    ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_steady.q unitcube_1k.exo
                    ARGS -c vortical_flow_steady.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_steady.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_steady.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

# Parallel + no virtualization

add_regression_test(riecg_vorticalflow_steady
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES vortical_flow_steady.q unitcube_1k.exo
                    ARGS -c vortical_flow_steady.q -i unitcube_1k.exo
                    BIN_BASELINE vortical_flow_steady_pe4.std.exo.0
                                 vortical_flow_steady_pe4.std.exo.1
                                 vortical_flow_steady_pe4.std.exo.2
                                 vortical_flow_steady_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_steady.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

# Parallel + virtualization

add_regression_test(riecg_vorticalflow_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow.q unitcube_1k.exo
                    ARGS -f -c vortical_flow.q -i unitcube_1k.exo -u 0.5
                    BIN_BASELINE vortical_flow_pe4_u0.5.std.exo.0
                                 vortical_flow_pe4_u0.5.std.exo.1
                                 vortical_flow_pe4_u0.5.std.exo.2
                                 vortical_flow_pe4_u0.5.std.exo.3
                                 vortical_flow_pe4_u0.5.std.exo.4
                                 vortical_flow_pe4_u0.5.std.exo.5
                                 vortical_flow_pe4_u0.5.std.exo.6
                                 vortical_flow_pe4_u0.5.std.exo.7
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

add_regression_test(riecg_vorticalflow_hllc_u0.5
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES vortical_flow_hllc.q unitcube_1k.exo
                    ARGS -f -c vortical_flow_hllc.q -i unitcube_1k.exo -u 0.5
                    BIN_BASELINE vortical_flow_hllc_pe4_u0.5.std.exo.0
                                 vortical_flow_hllc_pe4_u0.5.std.exo.1
                                 vortical_flow_hllc_pe4_u0.5.std.exo.2
                                 vortical_flow_hllc_pe4_u0.5.std.exo.3
                                 vortical_flow_hllc_pe4_u0.5.std.exo.4
                                 vortical_flow_hllc_pe4_u0.5.std.exo.5
                                 vortical_flow_hllc_pe4_u0.5.std.exo.6
                                 vortical_flow_hllc_pe4_u0.5.std.exo.7
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
                    TEXT_BASELINE diag_hllc.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(riecg_vorticalflow_hllc_stab2_u0.5
                    ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_hllc_stab2.q unitcube_1k.exo
                    ARGS -f -c vortical_flow_hllc_stab2.q -i unitcube_1k.exo
                         -u 0.5
                    BIN_BASELINE vortical_flow_hllc_stab2_pe4_u0.5.std.exo.0
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.1
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.2
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.3
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.4
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.5
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.6
                                 vortical_flow_hllc_stab2_pe4_u0.5.std.exo.7
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
                    TEXT_BASELINE diag_hllc_stab2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)


add_regression_test(riecg_vorticalflow_reord_u0.5
                    ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_reord.q unitcube_1k.exo
                    ARGS -c vortical_flow_reord.q -i unitcube_1k.exo -u 0.5
                         -f
                    BIN_BASELINE vortical_flow_pe4_u0.5.std.exo.0
                                 vortical_flow_pe4_u0.5.std.exo.1
                                 vortical_flow_pe4_u0.5.std.exo.2
                                 vortical_flow_pe4_u0.5.std.exo.3
                                 vortical_flow_pe4_u0.5.std.exo.4
                                 vortical_flow_pe4_u0.5.std.exo.5
                                 vortical_flow_pe4_u0.5.std.exo.6
                                 vortical_flow_pe4_u0.5.std.exo.7
                    BIN_RESULT out.e-s.0.8.0
                               out.e-s.0.8.1
                               out.e-s.0.8.2
                               out.e-s.0.8.3
                               out.e-s.0.8.4
                               out.e-s.0.8.5
                               out.e-s.0.8.6
                               out.e-s.0.8.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_reord.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

# Parallel + virtualization + migration

add_regression_test(riecg_vorticalflow_u0.9_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES vortical_flow.q unitcube_1k.exo
                    ARGS -c vortical_flow.q -i unitcube_1k.exo -u 0.9
                         -l 10 +balancer RecBipartLB +LBDebug 1 +cs
                    BIN_BASELINE vortical_flow_pe4_u0.9.std.exo.0
                                 vortical_flow_pe4_u0.9.std.exo.1
                                 vortical_flow_pe4_u0.9.std.exo.2
                                 vortical_flow_pe4_u0.9.std.exo.3
                                 vortical_flow_pe4_u0.9.std.exo.4
                                 vortical_flow_pe4_u0.9.std.exo.5
                                 vortical_flow_pe4_u0.9.std.exo.6
                                 vortical_flow_pe4_u0.9.std.exo.7
                                 vortical_flow_pe4_u0.9.std.exo.8
                                 vortical_flow_pe4_u0.9.std.exo.9
                                 vortical_flow_pe4_u0.9.std.exo.10
                                 vortical_flow_pe4_u0.9.std.exo.11
                                 vortical_flow_pe4_u0.9.std.exo.12
                                 vortical_flow_pe4_u0.9.std.exo.13
                                 vortical_flow_pe4_u0.9.std.exo.14
                                 vortical_flow_pe4_u0.9.std.exo.15
                                 vortical_flow_pe4_u0.9.std.exo.16
                                 vortical_flow_pe4_u0.9.std.exo.17
                                 vortical_flow_pe4_u0.9.std.exo.18
                                 vortical_flow_pe4_u0.9.std.exo.19
                                 vortical_flow_pe4_u0.9.std.exo.20
                                 vortical_flow_pe4_u0.9.std.exo.21
                                 vortical_flow_pe4_u0.9.std.exo.22
                                 vortical_flow_pe4_u0.9.std.exo.23
                                 vortical_flow_pe4_u0.9.std.exo.24
                                 vortical_flow_pe4_u0.9.std.exo.25
                                 vortical_flow_pe4_u0.9.std.exo.26
                                 vortical_flow_pe4_u0.9.std.exo.27
                                 vortical_flow_pe4_u0.9.std.exo.28
                                 vortical_flow_pe4_u0.9.std.exo.29
                                 vortical_flow_pe4_u0.9.std.exo.30
                                 vortical_flow_pe4_u0.9.std.exo.31
                                 vortical_flow_pe4_u0.9.std.exo.32
                                 vortical_flow_pe4_u0.9.std.exo.33
                                 vortical_flow_pe4_u0.9.std.exo.34
                                 vortical_flow_pe4_u0.9.std.exo.35
                                 vortical_flow_pe4_u0.9.std.exo.36
                                 vortical_flow_pe4_u0.9.std.exo.37
                    BIN_RESULT out.e-s.0.38.0
                               out.e-s.0.38.1
                               out.e-s.0.38.2
                               out.e-s.0.38.3
                               out.e-s.0.38.4
                               out.e-s.0.38.5
                               out.e-s.0.38.6
                               out.e-s.0.38.7
                               out.e-s.0.38.8
                               out.e-s.0.38.9
                               out.e-s.0.38.10
                               out.e-s.0.38.11
                               out.e-s.0.38.12
                               out.e-s.0.38.13
                               out.e-s.0.38.14
                               out.e-s.0.38.15
                               out.e-s.0.38.16
                               out.e-s.0.38.17
                               out.e-s.0.38.18
                               out.e-s.0.38.19
                               out.e-s.0.38.20
                               out.e-s.0.38.21
                               out.e-s.0.38.22
                               out.e-s.0.38.23
                               out.e-s.0.38.24
                               out.e-s.0.38.25
                               out.e-s.0.38.26
                               out.e-s.0.38.27
                               out.e-s.0.38.28
                               out.e-s.0.38.29
                               out.e-s.0.38.30
                               out.e-s.0.38.31
                               out.e-s.0.38.32
                               out.e-s.0.38.33
                               out.e-s.0.38.34
                               out.e-s.0.38.35
                               out.e-s.0.38.36
                               out.e-s.0.38.37
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg migration)
