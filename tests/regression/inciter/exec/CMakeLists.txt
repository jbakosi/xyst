# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --version
                    EXTRA_PASS_REGEXP "revision")

add_regression_test(help ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --help
                    EXTRA_PASS_REGEXP "inciter Command-line Parameters")

add_regression_test(unrecognized ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --nosucharg
                    EXTRA_PASS_REGEXP "Unknown keyword"
                    SKIP_RETURN_CODE true)

add_regression_test(lbfreq ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --lbfreq 0
                    EXTRA_PASS_REGEXP "out of bounds"
                    SKIP_RETURN_CODE true)

add_regression_test(rsfreq ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --rsfreq 0
                    EXTRA_PASS_REGEXP "out of bounds"
                    SKIP_RETURN_CODE true)
