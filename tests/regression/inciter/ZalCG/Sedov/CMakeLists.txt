# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_sedov_zalcg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sedov.q sedov_coarse.exo
                    ARGS -c sedov.q -i sedov_coarse.exo
                    BIN_BASELINE sedov.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS zalcg)

# Parallel + virtualization + migration

add_regression_test(compflow_euler_sedov_zalcg_u0.8_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES sedov.q sedov_coarse.exo
                    ARGS -c sedov.q -i sedov_coarse.exo -u 0.8
                         +balancer RecBipartLB +LBDebug 1 -l 5
                    BIN_BASELINE sedov_pe4_u0.8.std.exo.0
                                 sedov_pe4_u0.8.std.exo.1
                                 sedov_pe4_u0.8.std.exo.2
                                 sedov_pe4_u0.8.std.exo.3
                                 sedov_pe4_u0.8.std.exo.4
                                 sedov_pe4_u0.8.std.exo.5
                                 sedov_pe4_u0.8.std.exo.6
                                 sedov_pe4_u0.8.std.exo.7
                                 sedov_pe4_u0.8.std.exo.8
                                 sedov_pe4_u0.8.std.exo.9
                                 sedov_pe4_u0.8.std.exo.10
                                 sedov_pe4_u0.8.std.exo.11
                                 sedov_pe4_u0.8.std.exo.12
                                 sedov_pe4_u0.8.std.exo.13
                                 sedov_pe4_u0.8.std.exo.14
                                 sedov_pe4_u0.8.std.exo.15
                                 sedov_pe4_u0.8.std.exo.16
                                 sedov_pe4_u0.8.std.exo.17
                                 sedov_pe4_u0.8.std.exo.18
                    BIN_RESULT out.e-s.0.19.0
                               out.e-s.0.19.1
                               out.e-s.0.19.2
                               out.e-s.0.19.3
                               out.e-s.0.19.4
                               out.e-s.0.19.5
                               out.e-s.0.19.6
                               out.e-s.0.19.7
                               out.e-s.0.19.8
                               out.e-s.0.19.9
                               out.e-s.0.19.10
                               out.e-s.0.19.11
                               out.e-s.0.19.12
                               out.e-s.0.19.13
                               out.e-s.0.19.14
                               out.e-s.0.19.15
                               out.e-s.0.19.16
                               out.e-s.0.19.17
                               out.e-s.0.19.18
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS zalcg migration)
