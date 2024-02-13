# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# See regression/inciter/mesh_refinement/README.txt for naming conventions

# Serial

add_regression_test(href_t0ref_u_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_riecg.q -i unitsquare_01_955.exo -f
                    BIN_BASELINE slot_cyl_riecg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg amr)

add_regression_test(href_t0ref_u_reord_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_reord_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_reord_riecg.q
                         -i unitsquare_01_955.exo -f
                    BIN_BASELINE slot_cyl_riecg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg amr)

#add_regression_test(href_t0ref_cc_riecg ${INCITER_EXECUTABLE}
#                    NUMPES 1
#                    INPUTFILES slot_cyl_amr_coords.q unitcube_01_364.exo
#                    ARGS -c slot_cyl_amr_coords.q -i unitcube_01_364.exo
#                    BIN_BASELINE amr_init_coords.1.std.exo.0
#                    BIN_RESULT out.e-s.0.1.0
#                    BIN_DIFF_PROG_ARGS -m
#                    BIN_DIFF_PROG_CONF exodiff.cfg
#                    LABELS riecg amr)

add_regression_test(href_t0ref_idiu_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ic_uniform_deref.q unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref.q -i unitsquare_01_141.exo
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

add_regression_test(href_t0ref_uudd_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES uudd.q unitsquare_01_141.exo
                    ARGS -c uudd.q -i unitsquare_01_141.exo
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

add_regression_test(href_t0ref_iid_riecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES iid.q unitsquare_01_141.exo
                    ARGS -c iid.q -i unitsquare_01_141.exo
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
                    LABELS riecg amr)

# Parallel, no virtualization

add_regression_test(href_t0ref_u_riecg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_riecg.q -i unitsquare_01_955.exo -f
                    BIN_BASELINE slot_cyl_riecg_pe4_u0.0.std.exo.0
                                 slot_cyl_riecg_pe4_u0.0.std.exo.1
                                 slot_cyl_riecg_pe4_u0.0.std.exo.2
                                 slot_cyl_riecg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS riecg amr)

# See 90811bfa
#add_regression_test(href_t0ref_ii_riecg ${INCITER_EXECUTABLE}
#                    NUMPES 2
#                    INPUTFILES ic_ic.q unitsquare_01_141.exo
#                    ARGS -c ic_ic.q -i unitsquare_01_141.exo
#                    BIN_BASELINE ic_ic_t0ref.e-s.0.2.0
#                                 ic_ic_t0ref.e-s.0.2.1
#                                 ic_ic_t0ref.e-s.1.2.0
#                                 ic_ic_t0ref.e-s.1.2.1
#                                 ic_ic_t0ref.e-s.2.2.0
#                                 ic_ic_t0ref.e-s.2.2.1
#                                 amr_ic_ic.std.exo.0
#                                 amr_ic_ic.std.exo.1
#                    BIN_RESULT t0ref.e-s.0.2.0
#                               t0ref.e-s.0.2.1
#                               t0ref.e-s.1.2.0
#                               t0ref.e-s.1.2.1
#                               t0ref.e-s.2.2.0
#                               t0ref.e-s.2.2.1
#                               out.e-s.0.2.0
#                               out.e-s.0.2.1
#                    BIN_DIFF_PROG_ARGS -m
#                    BIN_DIFF_PROG_CONF exodiff.cfg
#                    LABELS riecg amr)

# See 90811bfa
#add_regression_test(href_t0ref_iu_riecg ${INCITER_EXECUTABLE}
#                    NUMPES 2
#                    INPUTFILES ic_uniform.q unitsquare_01_141.exo
#                    ARGS -c ic_uniform.q -i unitsquare_01_141.exo
#                    BIN_BASELINE amr_ic_uniform.std.exo.0
#                                 amr_ic_uniform.std.exo.1
#                    BIN_RESULT  out.e-s.0.2.0
#                                out.e-s.0.2.1
#                    BIN_DIFF_PROG_ARGS -m
#                    BIN_DIFF_PROG_CONF exodiff.cfg
#                    LABELS riecg amr)

add_regression_test(href_t0ref_idiu_riecg ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform_deref.q unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref.q -i unitsquare_01_141.exo
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.2.0
                                 ic_uniform_deref_t0ref.e-s.0.2.1
                                 ic_uniform_deref_t0ref.e-s.1.2.0
                                 ic_uniform_deref_t0ref.e-s.1.2.1
                                 ic_uniform_deref_t0ref.e-s.2.2.0
                                 ic_uniform_deref_t0ref.e-s.2.2.1
                                 ic_uniform_deref_t0ref.e-s.3.2.0
                                 ic_uniform_deref_t0ref.e-s.3.2.1
                                 ic_uniform_deref_t0ref.e-s.4.2.0
                                 ic_uniform_deref_t0ref.e-s.4.2.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(href_t0ref_idddiu_riecg ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform_deref2.q unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref2.q -i unitsquare_01_141.exo
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.2.0
                                 ic_uniform_deref_t0ref.e-s.0.2.1
                                 ic_uniform_deref_t0ref.e-s.1.2.0
                                 ic_uniform_deref_t0ref.e-s.1.2.1
                                 ic_uniform_deref_t0ref.e-s.2.2.0
                                 ic_uniform_deref_t0ref.e-s.2.2.1
                                 ic_uniform_deref_t0ref.e-s.3.2.0
                                 ic_uniform_deref_t0ref.e-s.3.2.1
                                 ic_uniform_deref_t0ref.e-s.4.2.0
                                 ic_uniform_deref_t0ref.e-s.4.2.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.4.2.0  # skip two no-op derefs
                               t0ref.e-s.4.2.1
                               t0ref.e-s.5.2.0
                               t0ref.e-s.5.2.1
                               t0ref.e-s.6.2.0
                               t0ref.e-s.6.2.1
                    BIN_DIFF_PROG_CONF exodiff_coords_but_ignore_time.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(href_t0ref_uudd_riecg ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES uudd.q unitsquare_01_141.exo
                    ARGS -c uudd.q -i unitsquare_01_141.exo
                    BIN_BASELINE uudd_t0ref.e-s.0.2.0
                                 uudd_t0ref.e-s.0.2.1
                                 uudd_t0ref.e-s.1.2.0
                                 uudd_t0ref.e-s.1.2.1
                                 uudd_t0ref.e-s.2.2.0
                                 uudd_t0ref.e-s.2.2.1
                                 uudd_t0ref.e-s.3.2.0
                                 uudd_t0ref.e-s.3.2.1
                                 uudd_t0ref.e-s.4.2.0
                                 uudd_t0ref.e-s.4.2.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(href_t0ref_uuuddd_riecg ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES uuuddd.q unitsquare_01_141.exo
                    ARGS -c uuuddd.q -i unitsquare_01_141.exo
                    BIN_BASELINE uuuddd_t0ref.e-s.0.2.0
                                 uuuddd_t0ref.e-s.0.2.1
                                 uuuddd_t0ref.e-s.1.2.0
                                 uuuddd_t0ref.e-s.1.2.1
                                 uuuddd_t0ref.e-s.2.2.0
                                 uuuddd_t0ref.e-s.2.2.1
                                 uuuddd_t0ref.e-s.3.2.0
                                 uuuddd_t0ref.e-s.3.2.1
                                 uuuddd_t0ref.e-s.4.2.0
                                 uuuddd_t0ref.e-s.4.2.1
                                 uuuddd_t0ref.e-s.5.2.0
                                 uuuddd_t0ref.e-s.5.2.1
                                 uuuddd_t0ref.e-s.6.2.0
                                 uuuddd_t0ref.e-s.6.2.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                               t0ref.e-s.5.2.0
                               t0ref.e-s.5.2.1
                               t0ref.e-s.6.2.0
                               t0ref.e-s.6.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

add_regression_test(href_t0ref_iid ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES iid.q unitsquare_01_141.exo
                    ARGS -c iid.q -i unitsquare_01_141.exo
                    BIN_BASELINE iid_t0ref.e-s.0.2.0
                                 iid_t0ref.e-s.1.2.0
                                 iid_t0ref.e-s.2.2.0
                                 iid_t0ref.e-s.3.2.0
                                 iid_t0ref.e-s.0.2.1
                                 iid_t0ref.e-s.1.2.1
                                 iid_t0ref.e-s.2.2.1
                                 iid_t0ref.e-s.3.2.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.1.2.0
                               t0ref.e-s.2.2.0
                               t0ref.e-s.3.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS riecg amr)

# Parallel, virtualization

add_regression_test(href_t0ref_u_iecg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr_riecg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_riecg.q -i unitsquare_01_955.exo
                         -u 0.5
                    BIN_BASELINE slot_cyl_riecg_pe4_u0.5.std.exo.0
                                 slot_cyl_riecg_pe4_u0.5.std.exo.1
                                 slot_cyl_riecg_pe4_u0.5.std.exo.2
                                 slot_cyl_riecg_pe4_u0.5.std.exo.3
                                 slot_cyl_riecg_pe4_u0.5.std.exo.4
                                 slot_cyl_riecg_pe4_u0.5.std.exo.5
                                 slot_cyl_riecg_pe4_u0.5.std.exo.6
                                 slot_cyl_riecg_pe4_u0.5.std.exo.7
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
                    LABELS riecg amr)
