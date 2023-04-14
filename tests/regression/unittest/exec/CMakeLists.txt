# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --version
                    EXTRA_PASS_REGEXP "version"
                    LABELS unittest)

add_regression_test(help ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --help
                    EXTRA_PASS_REGEXP "unittest Command-line Parameters"
                    LABELS unittest)

add_regression_test(helpkw ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --helpkw group
                    EXTRA_PASS_REGEXP "unittest command-line keyword"
                    LABELS unittest)

add_regression_test(singlegroup ${UNITTEST_EXECUTABLE}
                    NUMPES 2
                    ARGS --group DerivedData -v
                    EXTRA_PASS_REGEXP "Groups: DerivedData"
                    LABELS unittest)

add_regression_test(nogroup ${UNITTEST_EXECUTABLE}
                    NUMPES 1
                    ARGS --group group_definitely_does_not_exist -v
                    EXTRA_PASS_REGEXP "No test groups to be executed"
                    LABELS unittest)
