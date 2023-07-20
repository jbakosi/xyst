# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --version
                    EXTRA_PASS_REGEXP "version")

add_regression_test(help ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --help
                    EXTRA_PASS_REGEXP "unittest Command-line Parameters")

add_regression_test(singlegroup ${UNITTEST_EXECUTABLE}
                    NUMPES 2
                    ARGS --group DerivedData
                    EXTRA_PASS_REGEXP "Groups: DerivedData")

add_regression_test(nogroup ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --group group_definitely_does_not_exist
                    EXTRA_PASS_REGEXP "No test groups to be executed")
