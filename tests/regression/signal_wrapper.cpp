#include <cstdio>
#include <cstdlib>
#include <sys/wait.h>
#include <unistd.h>

int main( int, char** argv )
{
  int expected = atoi( argv[1] );
  printf( "SignalWrapper> Expecting signal: %d\n", expected );

  pid_t pid = fork();
  if (pid == -1) {
    // Signal_wrapper failed to fork
    return 1;
  }
  else if (pid) {
    // Parent - wait child and interpret its result
    int status = 0;
    wait( &status );
    printf( "SignalWrapper> Child process status: %d, WIFSIGNALED: %d, "
            "WTERMSIG: %d\n", status, WIFSIGNALED(status), WTERMSIG(status) );
    if (WIFSIGNALED(status)) {
      if (WTERMSIG(status) == expected) {
        printf( "SignalWrapper> Child process correctly not handled "
                "signal: %d\n", WTERMSIG(status) );
        return 0;
      } else {
        printf( "SignalWrapper> Child process returned unexpected "
                "signal: %d\n", WTERMSIG(status) );
        return EXIT_FAILURE;
      }
    } else {
      printf( "SignalWrapper> Child process handled signal: %d\n", expected );
      return EXIT_SUCCESS;
    }
  }
  else {
    // Child - execute wrapped command
    execvp( argv[2], argv + 2 );
    exit( EXIT_FAILURE );
  }
}
