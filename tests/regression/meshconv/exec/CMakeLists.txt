# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(version ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    ARGS --version
                    EXTRA_PASS_REGEXP "version")

add_regression_test(help ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    ARGS --help
                    EXTRA_PASS_REGEXP "meshconv Command-line Parameters")

add_regression_test(helpkw ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    ARGS --helpkw reorder
                    EXTRA_PASS_REGEXP "meshconv command-line keyword")
