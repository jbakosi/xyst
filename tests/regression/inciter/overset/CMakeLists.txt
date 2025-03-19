# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(lohcg_overset_linear ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES lohcg_overset_linear.q
                               half_viscous_sphere_bg_11k.exo
                               half_viscous_sphere_overset.exo
                    ARGS -c lohcg_overset_linear.q
                    ARGS -i half_viscous_sphere_bg_11k.exo
                         -i half_viscous_sphere_overset.exo
                    BIN_BASELINE lohcg_overset_linear.std.0.exo
                                 lohcg_overset_linear.std.1.exo
                    BIN_RESULT out.0.e-s.0.1.0
                               out.1.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -p
                    LABELS lohcg overset)

add_regression_test(lohcg_overset_intergrid ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES lohcg_overset_linear.q
                               half_viscous_sphere_bg_11k.exo
                               half_viscous_sphere_overset.exo
                    ARGS -c lohcg_overset_linear.q
                    ARGS -i half_viscous_sphere_bg_11k.exo
                         -i half_viscous_sphere_overset.exo
                    BIN_BASELINE lohcg_overset_intergrid.std.0.exo
                                 lohcg_overset_intergrid.std.1.exo
                    BIN_RESULT out.0.e-s.0.1.0
                               out.1.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_intergrid.cfg
                    LABELS lohcg overset)

# Parallel + virtualization

add_regression_test(lohcg_overset_linear ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES lohcg_overset_linear.q
                               half_viscous_sphere_bg_11k.exo
                               half_viscous_sphere_overset.exo
                    ARGS -c lohcg_overset_linear.q
                    ARGS -i half_viscous_sphere_bg_11k.exo
                         -i half_viscous_sphere_overset.exo
                         -u 0.3 -u 0.6
                    BIN_BASELINE lohcg_overset_linear.std.0.exo.0
                                 lohcg_overset_linear.std.0.exo.1
                                 lohcg_overset_linear.std.0.exo.2
                                 lohcg_overset_linear.std.0.exo.3
                                 lohcg_overset_linear.std.0.exo.4
                                 lohcg_overset_linear.std.1.exo.0
                                 lohcg_overset_linear.std.1.exo.1
                                 lohcg_overset_linear.std.1.exo.2
                                 lohcg_overset_linear.std.1.exo.3
                                 lohcg_overset_linear.std.1.exo.4
                                 lohcg_overset_linear.std.1.exo.5
                                 lohcg_overset_linear.std.1.exo.6
                                 lohcg_overset_linear.std.1.exo.7
                                 lohcg_overset_linear.std.1.exo.8
                    BIN_RESULT out.0.e-s.0.5.0
                               out.0.e-s.0.5.1
                               out.0.e-s.0.5.2
                               out.0.e-s.0.5.3
                               out.0.e-s.0.5.4
                               out.1.e-s.0.9.0
                               out.1.e-s.0.9.1
                               out.1.e-s.0.9.2
                               out.1.e-s.0.9.3
                               out.1.e-s.0.9.4
                               out.1.e-s.0.9.5
                               out.1.e-s.0.9.6
                               out.1.e-s.0.9.7
                               out.1.e-s.0.9.8
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS lohcg overset)
