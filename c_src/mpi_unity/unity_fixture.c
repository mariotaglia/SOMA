/* Copyright (c) 2010 James Grenning and Contributed to Unity Project
 * ==========================================
 *  Unity Project - A Test Framework for C
 *  Copyright (c) 2007 Mike Karlesky, Mark VanderVoord, Greg Williams
 *  [Released under MIT License. Please refer to license.txt for details]
 * ========================================== */

#include "unity_mpi_addon.h"
#include "unity_fixture.h"
#include "unity_internals.h"
#include <string.h>
#include <signal.h>
#include <mpi.h>
#include <stdio.h>
#include <stdbool.h>

#define TERMINAL_YELLOW   "\x1B[33m"
#define COLOR_RESET "\x1B[0m"

static MPI_Comm test_world;
static bool test_world_initialized = false;

MPI_Comm get_test_world()
{
    if (test_world_initialized)
        return test_world;

    fprintf(stderr, "ERROR: requested test world before first call to change_test_world_size\n");
    abort();
}

int get_test_world_rank()
{
    int rank;
    MPI_Comm_rank(get_test_world(), &rank);
    return rank;
}

int get_test_world_size()
{
    int size;
    MPI_Comm_size(get_test_world(), &size);
    return size;
}

int change_test_world_size(int newsize)
{
    if (newsize < 1)
    {
        fprintf(stderr, "ERROR in testing harness: size of new test-communicator must be"
                        "1 or greater, was %d, aborting\n", newsize);
        abort(); // this is not a recoverable mistake
    }
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (newsize > world_size)
    {
        if (rank==0)
            {
                fprintf(stderr, TERMINAL_YELLOW
                "\nERROR: the tests need %d ranks, but there are only %d available,"
                " some tests will be skipped. (use mpiexec -n <number> to get more ranks,"
                " use --use-hwthread-cpus or even --oversubscribe if necessary)\n"
                COLOR_RESET, newsize, world_size);
                fflush(stderr);
            }
        return -1;
    }
    int * sizes = (int *)malloc(sizeof(int) * world_size);
    MPI_Gather(&newsize, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int error = 0;
    if (rank == 0)
    {
        for (int i=0; i<world_size; i++)
        {
            if (sizes[i] != newsize)
            {
                if (!error)
                {
                    fprintf(stderr, "ERROR: all ranks on world must request the same newsize for test-world\n"
                                "rank 0 has received newsize==%d\n", newsize);
                    error = 1;
                }
                fprintf(stderr, "rank %d: newsize = %d\n", i, sizes[i]);
            }
        }
        if (error)
        {
            fprintf(stderr, "aborting.\n");
            abort();
        }
    }


    int ret;

    if (rank < newsize)
    {
        // this rank is part of the new test-communicator
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &test_world);
        ret = 1;
    }
    else
    {
        // not part of new tests. (MPI_UNDEFINED in split -> comm is MPI_COMM_NULL)
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &test_world);
        ret = 0;
    }

    free(sizes);
    test_world_initialized = true;

    return ret;
}

void signal_handler(int sig)
{
    char * sig_name;
    switch (sig)
        {
            case SIGFPE:
                sig_name = "SIGFPE";
                break;
            case SIGABRT:
                sig_name = "SIGABRT";
                break;
            case SIGSEGV:
                sig_name = "SIGSEGV";
                break;
            default:
                fprintf(stderr, "signal handler called with signal %d for which it was not designed. Aborting.", sig);
                exit(1);
        }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "world-rank %d received signal %d ( %s ) in test %s (line %lu, file %s)\n",
            rank,
            sig,
            sig_name,
            Unity.CurrentTestName,
            Unity.CurrentTestLineNumber,
            Unity.TestFile);
}
struct UNITY_FIXTURE_T UnityFixture;

/* If you decide to use the function pointer approach.
 * Build with -D UNITY_OUTPUT_CHAR=outputChar and include <stdio.h>
 * int (*outputChar)(int) = putchar; */

void setUp(void)    { /*does nothing*/ }
void tearDown(void) { /*does nothing*/ }

static void announceTestRun(unsigned int runNumber)
{
    UnityPrint("Unity test run ");
    UnityPrintNumberUnsigned(runNumber+1);
    UnityPrint(" of ");
    UnityPrintNumberUnsigned(UnityFixture.RepeatCount);
    UNITY_PRINT_EOL();
}

int UnityMain(int argc, const char* argv[], void (*runAllTests)(void))
{
    init_unity_mpi_addon();
    signal(SIGFPE, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);
    int result = UnityGetCommandLineOptions(argc, argv);
    unsigned int r;
    if (result != 0)
        return result;

    for (r = 0; r < UnityFixture.RepeatCount; r++)
    {
        UnityBegin(argv[0]);
        announceTestRun(r);
        runAllTests();
        if (!UnityFixture.Verbose) UNITY_PRINT_EOL();
        UnityEnd();
    }
    free_unity_mpi_addon();

    return (int)Unity.TestFailures;
}

static int selected(const char* filter, const char* name)
{
    if (filter == 0)
        return 1;
    return strstr(name, filter) ? 1 : 0;
}

static int testSelected(const char* test)
{
    return selected(UnityFixture.NameFilter, test);
}

static int groupSelected(const char* group)
{
    return selected(UnityFixture.GroupFilter, group);
}

void UnityTestRunner(unityfunction* setup,
                     unityfunction* testBody,
                     unityfunction* teardown,
                     const char* printableName,
                     const char* group,
                     const char* name,
                     const char* file,
                     unsigned int line)
{
    if (testSelected(name) && groupSelected(group))
    {
        Unity.TestFile = file;
        Unity.CurrentTestName = printableName;
        Unity.CurrentTestLineNumber = line;
        if (UnityFixture.Verbose)
        {
            UnityPrint(printableName);
        #ifndef UNITY_REPEAT_TEST_NAME
            Unity.CurrentTestName = NULL;
        #endif
        }
        else if (UnityFixture.Silent)
        {
            /* Do Nothing */
        }
        else
        {
            UNITY_OUTPUT_CHAR('.');
        }

        Unity.NumberOfTests++;
        UnityPointer_Init();

        UNITY_EXEC_TIME_START();

        if (TEST_PROTECT())
        {
            setup();
            testBody();
        }
        if (TEST_PROTECT())
        {
            teardown();
        }
        if (TEST_PROTECT())
        {
            UnityPointer_UndoAllSets();
        }
        UnityConcludeFixtureTest();
    }
}

void UnityIgnoreTest(const char* printableName, const char* group, const char* name)
{
    if (testSelected(name) && groupSelected(group))
    {
        Unity.NumberOfTests++;
        Unity.TestIgnores++;
        if (UnityFixture.Verbose)
        {
            UnityPrint(printableName);
            UNITY_PRINT_EOL();
        }
        else if (UnityFixture.Silent)
        {
            /* Do Nothing */
        }
        else
        {
            UNITY_OUTPUT_CHAR('!');
        }
    }
}

/*-------------------------------------------------------- */
/*Automatic pointer restoration functions */
struct PointerPair
{
    void** pointer;
    void* old_value;
};

static struct PointerPair pointer_store[UNITY_MAX_POINTERS];
static int pointer_index = 0;

void UnityPointer_Init(void)
{
    pointer_index = 0;
}

void UnityPointer_Set(void** pointer, void* newValue, UNITY_LINE_TYPE line)
{
    if (pointer_index >= UNITY_MAX_POINTERS)
    {
        UNITY_TEST_FAIL(line, "Too many pointers set");
    }
    else
    {
        pointer_store[pointer_index].pointer = pointer;
        pointer_store[pointer_index].old_value = *pointer;
        *pointer = newValue;
        pointer_index++;
    }
}

void UnityPointer_UndoAllSets(void)
{
    while (pointer_index > 0)
    {
        pointer_index--;
        *(pointer_store[pointer_index].pointer) =
            pointer_store[pointer_index].old_value;
    }
}

int UnityGetCommandLineOptions(int argc, const char* argv[])
{
    int i;
    UnityFixture.Verbose = 0;
    UnityFixture.Silent = 0;
    UnityFixture.GroupFilter = 0;
    UnityFixture.NameFilter = 0;
    UnityFixture.RepeatCount = 1;

    if (argc == 1)
        return 0;

    for (i = 1; i < argc; )
    {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            /* Usage */
            UnityPrint("Runs a series of unit tests.");
            UNITY_PRINT_EOL();
            UNITY_PRINT_EOL();
            UnityPrint("When no flag is specified, all tests are run.");
            UNITY_PRINT_EOL();
            UNITY_PRINT_EOL();
            UnityPrint("Optional flags:");
            UNITY_PRINT_EOL();
            UnityPrint("  -v          Verbose output: show all tests executed even if they pass");
            UNITY_PRINT_EOL();
            UnityPrint("  -s          Silent mode: minimal output showing only test failures");
            UNITY_PRINT_EOL();
            UnityPrint("  -g NAME     Only run tests in groups that contain the string NAME");
            UNITY_PRINT_EOL();
            UnityPrint("  -n NAME     Only run tests whose name contains the string NAME");
            UNITY_PRINT_EOL();
            UnityPrint("  -r NUMBER   Repeatedly run all tests NUMBER times");
            UNITY_PRINT_EOL();
            UnityPrint("  -h, --help  Display this help message");
            UNITY_PRINT_EOL();
            UNITY_PRINT_EOL();
#ifdef UNITY_CUSTOM_HELP_MSG
            /* User-defined help message, e.g. to point to project-specific documentation */
            UnityPrint(UNITY_CUSTOM_HELP_MSG);
            UNITY_PRINT_EOL();
#else
            /* Default help suffix if a custom one is not defined */
            UnityPrint("More information about Unity: https://www.throwtheswitch.org/unity");
            UNITY_PRINT_EOL();
#endif
            return 1;  /* Exit without running the tests */
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            UnityFixture.Verbose = 1;
            i++;
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            UnityFixture.Silent = 1;
            i++;
        }
        else if (strcmp(argv[i], "-g") == 0)
        {
            i++;
            if (i >= argc)
                return 1;
            UnityFixture.GroupFilter = argv[i];
            i++;
        }
        else if (strcmp(argv[i], "-n") == 0)
        {
            i++;
            if (i >= argc)
                return 1;
            UnityFixture.NameFilter = argv[i];
            i++;
        }
        else if (strcmp(argv[i], "-r") == 0)
        {
            UnityFixture.RepeatCount = 2;
            i++;
            if (i < argc)
            {
                if (*(argv[i]) >= '0' && *(argv[i]) <= '9')
                {
                    unsigned int digit = 0;
                    UnityFixture.RepeatCount = 0;
                    while (argv[i][digit] >= '0' && argv[i][digit] <= '9')
                    {
                        UnityFixture.RepeatCount *= 10;
                        UnityFixture.RepeatCount += (unsigned int)argv[i][digit++] - '0';
                    }
                    i++;
                }
            }
        }
        else
        {
            /* ignore unknown parameter */
            i++;
        }
    }
    return 0;
}

void UnityConcludeFixtureTest(void)
{

    int res = communicate_test_result(Unity.CurrentTestIgnored, Unity.CurrentTestFailed);
    if (res == UNITY_MPI_TEST_FAIL)
        unity_mpi_collect_and_print_err_msg();


    if (Unity.CurrentTestIgnored)
    {
        Unity.TestIgnores++;
        UNITY_PRINT_EOL();
    }
    else if (!Unity.CurrentTestFailed)
    {
        if (UnityFixture.Verbose)
        {
            UnityPrint(" ");
            UnityPrint(UnityStrPass);
            UNITY_EXEC_TIME_STOP();
            UNITY_PRINT_EXEC_TIME();
            UNITY_PRINT_EOL();
        }
    }
    else /* Unity.CurrentTestFailed */
    {
        Unity.TestFailures++;
        UNITY_PRINT_EOL();
    }

    Unity.CurrentTestFailed = 0;
    Unity.CurrentTestIgnored = 0;
}
