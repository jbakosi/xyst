# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# See regression/inciter/mesh_refinement/README.txt for naming conventions

# Serial

add_regression_test(amr_t0ref_u_trans_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_riecg.q -i unitsquare_01_955.exo -v -f
                    BIN_BASELINE slot_cyl_riecg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg amr)

add_regression_test(amr_t0ref_u_trans_reord_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_reord_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_reord_riecg.q
                         -i unitsquare_01_955.exo -v -f
                    BIN_BASELINE slot_cyl_riecg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg amr)
