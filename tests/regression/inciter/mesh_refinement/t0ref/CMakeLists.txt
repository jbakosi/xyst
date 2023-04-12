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

add_regression_test(amr_t0ref_idiu_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ic_uniform_deref.q unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref.q -i unitsquare_01_141.exo -v
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.1.0
                                 ic_uniform_deref_t0ref.e-s.1.1.0
                                 ic_uniform_deref_t0ref.e-s.2.1.0
                                 ic_uniform_deref_t0ref.e-s.3.1.0
                                 ic_uniform_deref_t0ref.e-s.4.1.0
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               t0ref.e-s.4.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(amr_t0ref_uudd_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES uudd.q unitsquare_01_141.exo
                    ARGS -c uudd.q -i unitsquare_01_141.exo -v
                    BIN_BASELINE uudd_t0ref.e-s.0.1.0
                                 uudd_t0ref.e-s.1.1.0
                                 uudd_t0ref.e-s.2.1.0
                                 uudd_t0ref.e-s.3.1.0
                                 uudd_t0ref.e-s.4.1.0
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               t0ref.e-s.4.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(amr_t0ref_iid_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES iid.q unitsquare_01_141.exo
                    ARGS -c iid.q -i unitsquare_01_141.exo -v
                    BIN_BASELINE iid_t0ref.e-s.0.1.0
                                 iid_t0ref.e-s.1.1.0
                                 iid_t0ref.e-s.2.1.0
                                 iid_t0ref.e-s.3.1.0
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)
