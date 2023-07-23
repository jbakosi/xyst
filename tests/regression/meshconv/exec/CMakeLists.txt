# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    ARGS -v
                    EXTRA_PASS_REGEXP "revision")

add_regression_test(help ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    ARGS -h
                    EXTRA_PASS_REGEXP "Usage")
