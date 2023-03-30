find_package(Git QUIET)

# Macro: RunGitCommand
# Description: short-hand macro for calling a git function. Outputs are the
#              "exit_code" and "output" variables. The "_permit_git_failure"
#              variable can locally override the exit code checking- use it
#              with caution.
macro(RunGitCommand)
    execute_process(COMMAND
        "${GIT_EXECUTABLE}" ${ARGV}
        WORKING_DIRECTORY "${_working_dir}"
        RESULT_VARIABLE exit_code
        OUTPUT_VARIABLE output
        ERROR_VARIABLE stderr
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT exit_code EQUAL 0 AND NOT _permit_git_failure)
        set(ENV{GIT_RETRIEVED_STATE} "false")

        # Issue 26: git info not properly set
        #
        # Check if we should fail if any of the exit codes are non-zero.
        # Most methods have a fall-back default value that's used in case of non-zero
        # exit codes. If you're feeling risky, disable this safety check and use
        # those default values.
        if(GIT_FAIL_IF_NONZERO_EXIT )
            string(REPLACE ";" " " args_with_spaces "${ARGV}")
            message(FATAL_ERROR "${stderr} (${GIT_EXECUTABLE} ${args_with_spaces})")
        endif()
    endif()
endmacro()
