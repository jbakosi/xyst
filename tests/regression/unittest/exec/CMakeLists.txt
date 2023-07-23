# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS -v
                    EXTRA_PASS_REGEXP "revision")

add_regression_test(help ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS -h
                    EXTRA_PASS_REGEXP "Usage")

add_regression_test(singlegroup ${UNITTEST_EXECUTABLE}
                    NUMPES 2
                    ARGS -g DerivedData
                    EXTRA_PASS_REGEXP "Groups: DerivedData")

add_regression_test(nogroup ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS -g group_definitely_does_not_exist
                    EXTRA_PASS_REGEXP "No test groups match"
                    SKIP_RETURN_CODE true)
