# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(chocg_elemblock_unsupported ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS -i sphere_bar.exo -c sphere.q
                    INPUTFILES sphere_bar.exo sphere.q
                    EXTRA_PASS_REGEXP "must only contain TRI"
                    SKIP_RETURN_CODE true
                    LABELS chocg)
