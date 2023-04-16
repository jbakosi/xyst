# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --version
                    EXTRA_PASS_REGEXP "version"
                    LABELS inciter)

add_regression_test(help ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --help
                    EXTRA_PASS_REGEXP "inciter Command-line Parameters"
                    LABELS inciter)

add_regression_test(helpctr ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --helpctr
                    EXTRA_PASS_REGEXP "inciter Control File Keywords"
                    LABELS inciter)

add_regression_test(helpkw ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --helpkw compflow
                    EXTRA_PASS_REGEXP "inciter control file keyword"
                    LABELS inciter)

add_regression_test(unrecognized ${INCITER_EXECUTABLE}
                    NUMPES 1
                    ARGS --nosucharg
                    EXTRA_PASS_REGEXP "Unknown keyword"
                    SKIP_RETURN_CODE true
                    LABELS inciter)
