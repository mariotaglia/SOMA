#include <mpi.h>
#include "unity_mpi_addon.h"
#include "unity_fixture.h"
#include "err_handling.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define MALLOC_ERROR_CHECK(ptr) {\
	if (ptr == NULL)\
	{\
		fprintf(stderr, "malloc error on rank %d (%s:%d:%s)",\
				world_rank, __FILE__, __LINE__, __func__);\
		exit(1);\
	}\
}

enum listening_state { write_through, record_err_msg, ignore };
static enum listening_state current_state = ignore;
static enum listening_state default_state = ignore;

static char *err_msg = NULL;
static int err_len = 0;
static int err_max_len = 1000;

static int world_rank;

static int err_msg_append_char(int ch)
{
	if (err_len == err_max_len - 1)
	{
		err_max_len *= 2;
		err_msg = realloc(err_msg, err_max_len*sizeof(*err_msg));
		MALLOC_ERROR_CHECK(err_msg);
	}

	err_msg[err_len++] = (char)ch;

	return ch;
}

static void delete_terminating_null_character()
{
    if (err_len == 0)
        return;

    if (err_msg[err_len - 1] == '\0')
    {
        err_len --;
    }
}

static void clear_err_msg()
{
	err_len = 0;
	err_msg[0] = '\0';
}

void unity_mpi_collect_and_print_err_msg()
{
    int test_world_size;
    MPI_Comm_size(get_test_world(), & test_world_size);
	// gather sizes of error messages
	int size_on_rank[ test_world_size];
	MPI_Gather(&err_len, 1, MPI_INT,
			size_on_rank, 1, MPI_INT,
			0, get_test_world());

	if (world_rank == 0)
	{

		int total_size = 0;
		int dspls[ test_world_size];
		for (int i=0; i < test_world_size; i++)
		{
			dspls[i] = total_size;
			total_size += size_on_rank[i];
		}
		assert(total_size > 0);
		char * recv = malloc(total_size * sizeof(char));
		MALLOC_ERROR_CHECK(recv);

		MPI_Gatherv(err_msg, err_len, MPI_CHAR,
				recv, size_on_rank, dspls, MPI_CHAR,
				0, get_test_world());


		int fail_ranks[test_world_size];
		int num_fail_ranks = 0;
		for (int i=0; i <test_world_size; i++)
		{
			if (size_on_rank[i] > 0)
			{
				fail_ranks[num_fail_ranks++] = i;
			}
		}
		// someone must have an error message.
		assert(num_fail_ranks > 0);

		if (num_fail_ranks ==test_world_size)
		{
			printf("\nfailure on all %d ranks:\n", test_world_size);
		}
		else
		{
			printf("\nfailure on rank(s) (");
			for (int i=0; i<num_fail_ranks; i++)
			{
				printf("%d,", fail_ranks[i]);
			}
			printf(") (there are %d ranks)\n", test_world_size);
		}
		for (int i=0; i<num_fail_ranks; i++)
		{
			printf("rank %d: %s\n", fail_ranks[i], &recv[dspls[fail_ranks[i]]]);
		}

		free(recv);
		}
	else
	{
		MPI_Gatherv(err_msg, err_len, MPI_CHAR,
				NULL, NULL, NULL, MPI_CHAR,
				0, get_test_world());

	}

	clear_err_msg();

}



void mpi_addon_err_msg_start()
{
    delete_terminating_null_character();
	current_state = record_err_msg;
}

void mpi_addon_err_msg_stop()
{
	err_msg_append_char('\0');
	current_state = default_state;
}

int mpi_addon_putchar(int ch)
{
	if (current_state == ignore)
		return ch;
	
	else if (current_state == write_through)
        {
            int ret = putchar(ch);
            fflush(stdout);
            return ret;
        }

	else if (current_state == record_err_msg)
	{
		return err_msg_append_char(ch);
	}
	else
	{
		fprintf(stderr, "illegal state in mpi_addon_putchar:"
				" current_state is %d "
				"(rank %d, file %s, line %d, func %s)\n",
				current_state, world_rank, __FILE__, __LINE__, __func__);
		exit(1);
	}

	

}

int communicate_test_result(int current_test_ignored, int current_test_failed)
{
	// you cant fail a test you ignored!
	assert(!current_test_ignored || !current_test_failed);

	int my_result;

	if (current_test_ignored)
	{
		my_result = UNITY_MPI_TEST_IGNORE;
	}
	else if (current_test_failed)
	{
		my_result = UNITY_MPI_TEST_FAIL;
	}
	else // success
	{
		my_result = UNITY_MPI_TEST_PASS;
	}

	int result;
	MPI_Allreduce(&my_result, &result, 1,
		       	MPI_INT, MPI_MAX, get_test_world());

	// if this rank ignores a test, all other ranks should have ignored it, too.
	assert( my_result != UNITY_MPI_TEST_IGNORE || result == UNITY_MPI_TEST_IGNORE);

	return result;
}

void init_unity_mpi_addon()
{

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	err_msg = malloc(err_max_len * sizeof(*err_msg));
	MALLOC_ERROR_CHECK(err_msg);


	if (world_rank == 0)
		default_state = write_through;
	else
		default_state = ignore;

	current_state = default_state;
}

void free_unity_mpi_addon()
{
	err_len = 0;
	err_max_len = 0;
	free(err_msg);
	err_msg = NULL;
}


int mpi_addon_write_err_msg (const char *s)
{
    delete_terminating_null_character();

    while (*s)
        {
            mpi_addon_putchar(*s);
            s++;
        }

    mpi_addon_putchar('\0');

    return 0;
}
