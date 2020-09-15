#include <stdbool.h>

#define UNITY_MPI_TEST_PASS 0
#define UNITY_MPI_TEST_IGNORE 1
#define UNITY_MPI_TEST_FAIL 2


// getting the error messages from unity
#define UNITY_OUTPUT_CHAR(a) mpi_addon_putchar(a);
int mpi_addon_putchar(int a);
void mpi_addon_err_msg_start();
void mpi_addon_err_msg_stop();

// write a custom error message
int mpi_addon_write_err_msg(const char * s);

// communicate with other ranks and with the user
int communicate_test_result(int current_test_ignored, int current_test_failed);
void unity_mpi_collect_and_print_err_msg();


// init
void init_unity_mpi_addon();
void free_unity_mpi_addon();

