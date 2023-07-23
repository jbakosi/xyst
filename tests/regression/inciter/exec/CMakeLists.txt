# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS -v
                    EXTRA_PASS_REGEXP "revision")

add_regression_test(help ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS -h
                    EXTRA_PASS_REGEXP "Usage")

add_regression_test(unrecognized ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --nosucharg
                    EXTRA_PASS_REGEXP "invalid option"
                    SKIP_RETURN_CODE true)
