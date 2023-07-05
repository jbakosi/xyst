# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_sod ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod.q rectangle_01_1.5k.exo
                    ARGS -c sod.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                                  sod_hist.p1.std
                                  sod_hist.p2.std
                    TEXT_RESULT diag
                                out.hist.p1
                                out.hist.p2
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_sod_hist_timeint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_hist_timeint.q rectangle_01_1.5k.exo
                    ARGS -c sod_hist_timeint.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_hist_timeint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                                  sod_hist_timeint.p1.std
                                  sod_hist_timeint.p2.std
                    TEXT_RESULT diag
                                out.hist.p1
                                out.hist.p2
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_sod_surf ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_surf.q rectangle_01_1.5k.exo
                    ARGS -c sod_surf.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod.std.exo
                                 sod-surf.2.std.exo
                                 sod-surf.4.std.exo
                                 sod-surf.5.std.exo
                                 sod-surf.6.std.exo
                    BIN_RESULT out.e-s.0.1.0
                               out-surf.2.e-s.0.1.0
                               out-surf.4.e-s.0.1.0
                               out-surf.5.e-s.0.1.0
                               out-surf.6.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_sod_surf_timeint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_surf_timeint.q rectangle_01_1.5k.exo
                    ARGS -c sod_surf_timeint.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_timeint.std.exo
                                 sod-surf_timeint.2.std.exo
                                 sod-surf_timeint.4.std.exo
                                 sod-surf_timeint.5.std.exo
                                 sod-surf_timeint.6.std.exo
                    BIN_RESULT out.e-s.0.1.0
                               out-surf.2.e-s.0.1.0
                               out-surf.4.e-s.0.1.0
                               out-surf.5.e-s.0.1.0
                               out-surf.6.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_sod_userdef ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_userdef.q rectangle_01_1.5k.exo
                    ARGS -c sod_userdef.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                                  sod_hist.p1.std
                                  sod_hist.p2.std
                    TEXT_RESULT diag
                                out.hist.p1
                                out.hist.p2
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                    LABELS riecg)

add_regression_test(nodepvar ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_nodepvar.q rectangle_01_1.5k.exo
                    ARGS -c sod_nodepvar.q -i rectangle_01_1.5k.exo -v
                    EXTRA_PASS_REGEXP "Dependent variable not specified"
                    SKIP_RETURN_CODE true)

add_regression_test(badrange ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_badrange.q rectangle_01_1.5k.exo
                    ARGS -c sod_badrange.q -i rectangle_01_1.5k.exo -v
                    EXTRA_PASS_REGEXP "time range"
                    SKIP_RETURN_CODE true)

add_regression_test(badrange_dt ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_badrange_dt.q rectangle_01_1.5k.exo
                    ARGS -c sod_badrange_dt.q -i rectangle_01_1.5k.exo -v
                    EXTRA_PASS_REGEXP "time range"
                    SKIP_RETURN_CODE true)

# Parallel

add_regression_test(compflow_euler_sod_hist_range ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_hist_range.q rectangle_01_1.5k.exo
                    ARGS -c sod_hist_range.q -i rectangle_01_1.5k.exo -v
                    TEXT_BASELINE diag_hist_range.std
                                  sod_hist_range.p1.std
                                  sod_hist_range.p2.std
                                  sod_hist_range.p3.std
                    TEXT_RESULT diag
                                out.hist.p1
                                out.hist.p2
                                out.hist.p3
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                    LABELS riecg)

add_regression_test(compflow_euler_sod_parsurf ${INCITER_EXECUTABLE}
                    NUMPES 2 PPN 1
                    INPUTFILES sod_surf_rcb.q rectangle_01_1.5k.exo
                    ARGS -c sod_surf_rcb.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod-surf.1.std.exo.0
                                 sod-surf.1.std.exo.1
                                 sod-surf.5.std.exo.0
                                 sod-surf.5.std.exo.1
                    BIN_RESULT out-surf.1.e-s.0.2.0
                               out-surf.1.e-s.0.2.1
                               out-surf.5.e-s.0.2.0
                               out-surf.5.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg)

# Parallel + virtualization + migration

add_regression_test(compflow_euler_sod_u0.8_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES sod.q rectangle_01_1.5k.exo
                    ARGS -c sod.q -i rectangle_01_1.5k.exo -v -u 0.8
                         +balancer RecBipartLB +LBDebug 1
                    BIN_BASELINE sod_pe4_u0.8.std.exo.0
                                 sod_pe4_u0.8.std.exo.1
                                 sod_pe4_u0.8.std.exo.2
                                 sod_pe4_u0.8.std.exo.3
                                 sod_pe4_u0.8.std.exo.4
                                 sod_pe4_u0.8.std.exo.5
                                 sod_pe4_u0.8.std.exo.6
                                 sod_pe4_u0.8.std.exo.7
                                 sod_pe4_u0.8.std.exo.8
                                 sod_pe4_u0.8.std.exo.9
                                 sod_pe4_u0.8.std.exo.10
                                 sod_pe4_u0.8.std.exo.11
                                 sod_pe4_u0.8.std.exo.12
                                 sod_pe4_u0.8.std.exo.13
                                 sod_pe4_u0.8.std.exo.14
                                 sod_pe4_u0.8.std.exo.15
                                 sod_pe4_u0.8.std.exo.16
                                 sod_pe4_u0.8.std.exo.17
                                 sod_pe4_u0.8.std.exo.18
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
                                  sod_hist.p1.std
                                  sod_hist.p2.std
                    TEXT_RESULT diag
                                out.hist.p1
                                out.hist.p2
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                                        hist.ndiff.cfg
                                        hist.ndiff.cfg
                    LABELS riecg migration)
