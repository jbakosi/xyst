# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_bump_zalcg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES bump.q bump.exo
                    ARGS -c bump.q -i bump.exo
                    BIN_BASELINE bump.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS zalcg)

# Parallel + virtualization + migration

add_regression_test(compflow_euler_bump_zalcg_u0.8_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES bump.q bump.exo
                    ARGS -c bump.q -i bump.exo -u 0.8
                         +balancer RecBipartLB +LBDebug 1 -l 5
                    BIN_BASELINE bump_pe4_u0.8.std.exo.0
                                 bump_pe4_u0.8.std.exo.1
                                 bump_pe4_u0.8.std.exo.2
                                 bump_pe4_u0.8.std.exo.3
                                 bump_pe4_u0.8.std.exo.4
                                 bump_pe4_u0.8.std.exo.5
                                 bump_pe4_u0.8.std.exo.6
                                 bump_pe4_u0.8.std.exo.7
                                 bump_pe4_u0.8.std.exo.8
                                 bump_pe4_u0.8.std.exo.9
                                 bump_pe4_u0.8.std.exo.10
                                 bump_pe4_u0.8.std.exo.11
                                 bump_pe4_u0.8.std.exo.12
                                 bump_pe4_u0.8.std.exo.13
                                 bump_pe4_u0.8.std.exo.14
                                 bump_pe4_u0.8.std.exo.15
                                 bump_pe4_u0.8.std.exo.16
                                 bump_pe4_u0.8.std.exo.17
                                 bump_pe4_u0.8.std.exo.18
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
