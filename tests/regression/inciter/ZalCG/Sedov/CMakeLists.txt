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
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag.ndiff.cfg
                    LABELS zalcg migration)

add_regression_test(compflow_euler_sedov_zalcg_dea_u0.8_migr
                    ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES sedov_dea.q sedov_coarse.exo
                    ARGS -c sedov_dea.q -i sedov_coarse.exo -u 0.8
                         +balancer RecBipartLB +LBDebug 1 -l 5
                    TEXT_BASELINE diag_dea.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF diag_dea.ndiff.cfg
                    LABELS zalcg migration)
