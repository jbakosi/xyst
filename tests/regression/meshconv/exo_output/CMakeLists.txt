# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(asc2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse.exo
                    BIN_BASELINE cube_coarse.exo.std
                    BIN_RESULT cube_coarse.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(ugrid2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sphere_tetra.1.ugrid
                    ARGS -i sphere_tetra.1.ugrid -o sphere_tetra.exo
                    BIN_BASELINE sphere_tetra.exo.std
                    BIN_RESULT sphere_tetra.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(rdgflo2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES rdgflo.domn.bump
                    ARGS -i rdgflo.domn.bump -o bump.exo
                    BIN_BASELINE bump.exo.std
                    BIN_RESULT bump.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(netgen2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.mesh
                    ARGS -i box_24.mesh -o box_24.exo
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(medit2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube86.mesh
                    ARGS -i cube86.mesh -o cube86.exo
                    BIN_BASELINE cube86.exo.std
                    BIN_RESULT cube86.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(medittrionly2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube86_trionly.mesh
                    ARGS -i cube86_trionly.mesh -o cube86_trionly.exo
                    BIN_BASELINE cube86_trionly.exo.std
                    BIN_RESULT cube86_trionly.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshtxt2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh
                    ARGS -i box_24.txt.msh -o box_24.exo
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshbin2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.msh
                    ARGS -i box_24.msh -o box_24.exo
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(multiblockexo2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_5blocks.exo
                    ARGS -i shear_5blocks.exo -o shear.exo
                    BIN_BASELINE shear.exo.std
                    BIN_RESULT shear.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(stitchexo2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES out-surf.1.e-s.0.30.0
                               out-surf.1.e-s.0.30.1
                               out-surf.1.e-s.0.30.2
                               out-surf.1.e-s.0.30.3
                               out-surf.1.e-s.0.30.4
                               out-surf.1.e-s.0.30.5
                               out-surf.1.e-s.0.30.6
                               out-surf.1.e-s.0.30.7
                               out-surf.1.e-s.0.30.8
                               out-surf.1.e-s.0.30.9
                               out-surf.1.e-s.0.30.10
                               out-surf.1.e-s.0.30.11
                               out-surf.1.e-s.0.30.12
                               out-surf.1.e-s.0.30.13
                               out-surf.1.e-s.0.30.14
                               out-surf.1.e-s.0.30.15
                               out-surf.1.e-s.0.30.16
                               out-surf.1.e-s.0.30.17
                               out-surf.1.e-s.0.30.18
                               out-surf.1.e-s.0.30.19
                               out-surf.1.e-s.0.30.20
                               out-surf.1.e-s.0.30.21
                               out-surf.1.e-s.0.30.22
                               out-surf.1.e-s.0.30.23
                               out-surf.1.e-s.0.30.24
                               out-surf.1.e-s.0.30.25
                               out-surf.1.e-s.0.30.26
                               out-surf.1.e-s.0.30.27
                               out-surf.1.e-s.0.30.28
                               out-surf.1.e-s.0.30.29
                    ARGS -i out-surf.1.e-s.0.30.% -o out.exo
                    BIN_BASELINE out-surf.1.exo
                    BIN_RESULT out.exo
                    BIN_DIFF_PROG_CONF exodiff_stitch.cfg)

add_regression_test(asc2reordexo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse_reord.exo -r
                    BIN_BASELINE cube_coarse_reord.exo.std
                    BIN_RESULT cube_coarse_reord.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

if(EXCEPTIONS_WRITE_TO_CERR)
  add_regression_test(percent ${MESHCONV_EXECUTABLE}
                      NUMPES 1
                      ARGS -i out-surf.1.e-s.0.blah.% -o out.exo
                      EXTRA_PASS_REGEXP "percent sign"
                      SKIP_RETURN_CODE true)
endif()

add_regression_test(signal_unhandled ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse_reord.exo
                         --raise_signal 7
                    EXPECT_SIGNAL 7
                    EXCLUSIVE_PASS_REGEXP "handled signal: 7"
                    SKIP_RETURN_CODE true)

add_regression_test(signal_handled ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse_reord.exo
                         --raise_signal 6
                    EXPECT_SIGNAL 6
                    EXCLUSIVE_PASS_REGEXP "handled signal: 6"
                    SKIP_RETURN_CODE true)
