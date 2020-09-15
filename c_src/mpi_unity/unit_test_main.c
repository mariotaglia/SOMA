#include "unity_fixture.h"
#include <mpi.h>
#include <sys/time.h>
#include "run_unit_tests.h"
#include "err_handling.h"

#define RESET "\x1B[0m"
#define RED   "\x1B[31m"
#define GREEN "\x1B[32m"

static int final_test_evaluation(int ret_unity, long long duration);
static int print_final_test_results(int * ret_vals, int ret_val_size, long long duration);

static long long get_current_time_millis();

int main(int argc, const char ** argv)
{
	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    change_test_world_size (world_size);

	long long start_time = get_current_time_millis();
	int ret_unity = UnityMain(argc, argv, run_all_unit_tests);
    long long stop_time = get_current_time_millis();
    long long duration = stop_time - start_time;
	
	int test_ret = final_test_evaluation(ret_unity, duration);

	MPI_Finalize();

	// should I return failures?
    // Ends MPI ungracefully but makes it easier to check for errors from commandline
	(void) test_ret;
	return 0;
}

static int print_final_test_results(int * ret_vals, int ret_val_size, long long duration)
{

	printf("\n------------\n\n");
	// see what ranks have failures
	int *fail_ranks = (int*)malloc(ret_val_size*sizeof(int));
	int num_failing_ranks = 0;

	for (int i=0; i < ret_val_size; i++)
	{
		if (ret_vals[i] != 0)
		{
			fail_ranks[num_failing_ranks++] = i;
		}
	}

    printf("%lld ms\n", duration);

	if (num_failing_ranks == 0)
	{
		printf(GREEN "OK\n");
	}
	else if (num_failing_ranks == ret_val_size)
	{
		printf(RED "FAIL (on all ranks)\n");
	}
	else
	{
		printf(RED "FAIL (on rank(s) ");
		for (int i=0; i < num_failing_ranks; i++)
		{
			printf("%d,", fail_ranks[i]);
		}
		printf(")\n");
	}

	printf(RESET);

	for (int i=0; i < num_failing_ranks; i++)
	{
		printf("rank %d failed %d tests\n",
				fail_ranks[i], ret_vals[fail_ranks[i]]);
	}

	free(fail_ranks);

	return num_failing_ranks;
}


static int final_test_evaluation(int ret_unity, long long duration)
{

	int worldrank, worldsize;
	MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
	int *ret_unity_all = (int*)malloc(worldsize*sizeof(int));
	MPI_Gather(&ret_unity, 1, MPI_INT,
			ret_unity_all, 1, MPI_INT,
			0, MPI_COMM_WORLD);

	int ret;
	if (worldrank != 0)
		ret = 0;
	else
		ret = print_final_test_results(ret_unity_all, worldsize, duration);
	free(ret_unity_all);
	return ret;
}

static long long get_current_time_millis()
{
    struct timeval t;
    gettimeofday(&t, NULL);

    return (t.tv_sec) * 1000 + (t.tv_usec) / 1000 ;
}

