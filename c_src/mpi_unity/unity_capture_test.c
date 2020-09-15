#include "unity_fixture.h"
#include "unity_capture.h"
#include "unity_capture_internals.h"
#include "err_handling.h"

#include <stdbool.h>

#include <unistd.h>
#include <fcntl.h>

#define NO_BARRIER (unsigned)0
#define BARRIER_BEFORE (unsigned)1
#define BARRIER_AFTER (unsigned)2
static int fenced_system (const char *cmd, unsigned barr)
{
    int ret = 0;

    if ((barr | BARRIER_BEFORE) == barr)
        {
            MPI_Barrier (get_test_world ());
        }

    int rank;
    MPI_Comm_rank (get_test_world (), &rank);
    if (rank == 0)
        {
            ret = system (cmd);
        }

    if ((barr | BARRIER_AFTER) == barr)
        {
            MPI_Barrier (get_test_world ());
        }

    return ret;
}

#define ENSURE_EXISTS(file_name_lit) fenced_system("touch " file_name_lit, BARRIER_AFTER)
#define ENSURE_NOT_EXISTS(file_name_lit) fenced_system("touch " file_name_lit " && rm " file_name_lit, BARRIER_AFTER)

#define FILE_IS_EMPTY(file_name_lit) file_is_empty("find . -empty -name " file_name_lit)
static bool file_is_empty (const char *find_file)
{
    FILE *fp = popen (find_file, "r");
    TEST_ASSERT_NOT_NULL(fp);
    char output[1000];
    bool ret = (NULL == (fgets (output, sizeof (output), fp)));
    pclose (fp);
    return ret;
}

struct H5_error_handling_data {
    herr_t (*func)(long, void*);
    void *client_data;
};

static void silence_h5_errors (struct H5_error_handling_data *err_data)
{

    H5Eget_auto (H5E_DEFAULT, &err_data->func, &err_data->client_data);

    H5Eset_auto (H5E_DEFAULT, NULL, NULL);
}

static void restore_h5_errors (struct H5_error_handling_data err_data)
{
    H5Eset_auto (H5E_DEFAULT, err_data.func, err_data.client_data);
}

TEST_GROUP(string_handling);

TEST_SETUP(string_handling)
{
    UnityMalloc_StartTest ();
}

TEST_TEAR_DOWN(string_handling)
{
    UnityMalloc_EndTest ();
}

TEST(string_handling, simple_append_num_to_name)
{
    char *name = "x.h5";
    int num = 0;
    char *newname = append_number_to_filename (name, num);
    TEST_ASSERT_EQUAL_STRING("x_0.h5", newname);
    free (newname);
}

TEST(string_handling, append_num_to_complex_name)
{
    char *name = "x->vari_able.h5";
    int  num= 0;
    char *newname = append_number_to_filename (name, num);
    TEST_ASSERT_EQUAL_STRING("x->vari_able_0.h5", newname);
    free (newname);
}

TEST(string_handling, append_longer_num)
{
    char *name = "x.txt";
    int num = 112;
    char *newname = append_number_to_filename (name, num);
    TEST_ASSERT_EQUAL_STRING("x_112.txt", newname);
    free (newname);
}

#define H5FILE "test_file.h5"

TEST_GROUP(file_creation_w1_4);

TEST_SETUP(file_creation_w1_4)
{
    UnityMalloc_StartTest ();
}

TEST_TEAR_DOWN(file_creation_w1_4)
{
    UnityMalloc_EndTest ();
}

TEST(file_creation_w1_4, create_not_existing)
{

    ENSURE_NOT_EXISTS(H5FILE);
    capture_file_t file = capture_file_open (H5FILE, get_test_world (), FAIL_IF_EXISTS);
    TEST_ASSERT_TRUE(file.file_id > 0);
    capture_file_close (&file);
}

TEST(file_creation_w1_4, fail_to_create_existing_if_mode_is_fail_if_exists)
{
    ENSURE_EXISTS(H5FILE);

    struct H5_error_handling_data err_data;
    silence_h5_errors (&err_data);
    capture_file_t file = capture_file_open (H5FILE, get_test_world (), FAIL_IF_EXISTS);
    restore_h5_errors (err_data);
    TEST_ASSERT_TRUE(file.file_id < 0);

    capture_file_close (&file);
}

TEST(file_creation_w1_4, override_existing_if_mode_says_so)
{
    ENSURE_EXISTS(H5FILE);

    capture_file_t file = capture_file_open (H5FILE, get_test_world (), OVERRIDE_IF_EXISTS);
    TEST_ASSERT_TRUE(file.file_id > 0);

    TEST_ASSERT_TRUE(access(H5FILE, F_OK) != 1);
    TEST_ASSERT_FALSE(FILE_IS_EMPTY (H5FILE));

    capture_file_close (&file);
}

TEST(file_creation_w1_4, ignore_existing_and_do_nothing_if_mode_is_no_op_on_existing)
{
    ENSURE_EXISTS(H5FILE);

    capture_file_t file = capture_file_open (H5FILE, get_test_world (), NO_OP_IF_EXISTS);

    TEST_ASSERT_EQUAL(file.file_id, -1);

    capture_file_close (&file);
}

TEST(file_creation_w1_4, can_create_with_noop_if_exists)
{

    ENSURE_EXISTS(H5FILE);
    capture_file_t file = capture_file_open (H5FILE, get_test_world (), NO_OP_IF_EXISTS);

    TEST_ASSERT(file.act == IGNORE_CALLS);

    capture_file_close (&file);
}

TEST(file_creation_w1_4, can_create_with_noop_if_not_exists)
{
    ENSURE_NOT_EXISTS(H5FILE);
    capture_file_t file = capture_file_open (H5FILE, get_test_world (), NO_OP_IF_EXISTS);

    TEST_ASSERT(file.file_id > 0);
    TEST_ASSERT(file.act == FILE_WRITE);

    capture_file_close (&file);

}

TEST(file_creation_w1_4, cannot_create_with_readonly)
{
    ENSURE_NOT_EXISTS(H5FILE);

    struct H5_error_handling_data err_data;
    silence_h5_errors(&err_data);
    capture_file_t file = capture_file_open (H5FILE, get_test_world(), READ_ONLY);
    restore_h5_errors(err_data);
    TEST_ASSERT(file.file_id < 0);

}

TEST(file_creation_w1_4, has_proper_group_structure)
{
    ENSURE_NOT_EXISTS(H5FILE);
    capture_file_t file = capture_file_open (H5FILE, get_test_world (), FAIL_IF_EXISTS);
    TEST_ASSERT(file.file_id > 0);

    htri_t inputs_exist = H5Lexists (file.file_id, "inputs", H5P_DEFAULT);
    htri_t results_exist = H5Lexists (file.file_id, "results", H5P_DEFAULT);
    htri_t special_exist = H5Lexists (file.file_id, "special", H5P_DEFAULT);

    // no failures
    TEST_ASSERT_TRUE(inputs_exist >= 0 && results_exist >= 0 && special_exist >= 0);

    // test they exist
    TEST_ASSERT_TRUE(inputs_exist > 0);
    TEST_ASSERT_TRUE(results_exist > 0);
    TEST_ASSERT_TRUE(special_exist > 0);

    capture_file_close (&file);
}

TEST(file_creation_w1_4, evade_if_exists_simply_creates_file_if_non_existant)
{
    ENSURE_NOT_EXISTS(H5FILE);

    capture_file_t file = capture_file_open(H5FILE, get_test_world(), EVADE_IF_EXISTS);
    TEST_ASSERT_TRUE(file.file_id > 0);
    TEST_ASSERT_TRUE(!FILE_IS_EMPTY(H5FILE));
    capture_file_close(&file);
}

TEST(file_creation_w1_4, evade_if_exists_evades_one)
{
    ENSURE_EXISTS("file_to_evade.h5");
    ENSURE_NOT_EXISTS("file_to_evade_1.h5");

    capture_file_t file = capture_file_open("file_to_evade.h5", get_test_world(), EVADE_IF_EXISTS);
    TEST_ASSERT_TRUE(file.file_id > 0);
    TEST_ASSERT_TRUE(access("file_to_evade.h5", F_OK) != -1);
    TEST_ASSERT_TRUE(access("file_to_evade_1.h5", F_OK) != -1);
    TEST_ASSERT_TRUE(!FILE_IS_EMPTY("file_to_evade_1.h5"));

    capture_file_close(&file);
}



TEST_GROUP(read_and_write_single_w1_3_6);
static capture_file_t rws_file;


TEST_SETUP(read_and_write_single_w1_3_6)
{
    UnityMalloc_StartTest();
    ENSURE_NOT_EXISTS(H5FILE);
    rws_file = capture_file_open (H5FILE, get_test_world (), FAIL_IF_EXISTS);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);
}

TEST_TEAR_DOWN(read_and_write_single_w1_3_6)
{
    capture_file_close (&rws_file);
    UnityMalloc_EndTest();
}

TEST(read_and_write_single_w1_3_6, write_native_int)
{
    int x = 2;
    int err = SAVE_INT(x, rws_file, INPUTS);
    TEST_ASSERT_EQUAL_INT(0, err);
}

TEST(read_and_write_single_w1_3_6, write_and_read_native_int)
{
    int x = 1000 + rws_file.comm_rank;
    int err = SAVE_INT(x, rws_file, INPUTS);
    TEST_ASSERT_EQUAL_INT(0, err);
    capture_file_close (&rws_file);

    rws_file = capture_file_open (H5FILE, get_test_world (), READ_ONLY);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);


    int a = -1;
    err = read_int (rws_file, INPUTS, "x", &a, 1);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_INT(x, a);
}

TEST(read_and_write_single_w1_3_6, do_not_write_to_noop_file)
{
    int x = 2;
    capture_file_t file;
    file.file_id = -1;
    file.act = IGNORE_CALLS;

    int err = SAVE_INT(x, file, INPUTS);
    TEST_ASSERT_TRUE(err == 0);
}

TEST(read_and_write_single_w1_3_6, fail_to_read_with_wrong_sized_communicator)
{
    int x = 2;
    int err = SAVE_INT(x, rws_file, INPUTS);
    TEST_ASSERT_EQUAL_INT(0, err);
    capture_file_close(&rws_file);

    int color = (get_test_world_rank() < get_test_world_size()/2) ? 0 : MPI_UNDEFINED;
    MPI_Comm small_comm;
    MPI_Comm_split(get_test_world(), color, 0, &small_comm);

    if (small_comm != MPI_COMM_NULL)
    {
        //redirect stderr temporarily
        int backup, newoutput;
        fflush(stdout);
        backup = dup( STDOUT_FILENO);
        newoutput = open("/dev/null", O_WRONLY);
        dup2(newoutput, STDOUT_FILENO);
        close(newoutput);

        rws_file = capture_file_open(H5FILE, small_comm, READ_ONLY);

        //restore stderr
        fflush(stdout);
        dup2(backup, STDOUT_FILENO);
        close(backup);


        TEST_ASSERT_TRUE(rws_file.file_id < 0);

        MPI_Comm_free(&small_comm);
    }


}

TEST(read_and_write_single_w1_3_6, write_and_read_native_double)
{
    double x = 1 + (double) rws_file.comm_rank / 10;
    int err = SAVE_DOUBLE (x, rws_file, RESULTS) ;
    TEST_ASSERT_EQUAL_INT(0, err);
    capture_file_close (&rws_file);


    rws_file = capture_file_open(H5FILE, get_test_world (), READ_ONLY);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);

    double a = -1;
    err = read_double (rws_file, RESULTS, "x", &a, 1);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_DOUBLE(x, a);
}

TEST(read_and_write_single_w1_3_6, write_and_read_various_types)
{
    float float_in = 1 + (float ) rws_file.comm_rank / 10;
    uint16_t short_in = rws_file.comm_rank * 1000;
    uint64_t long_in = rws_file.comm_rank * 1e16;
    int err = SAVE_FLOAT (float_in, rws_file, RESULTS);
    err |= SAVE_UINT16(short_in, rws_file, INPUTS) ;
    err |= SAVE_UINT64(long_in, rws_file, INPUTS) ;
    TEST_ASSERT_EQUAL_INT(0, err);
    capture_file_close (&rws_file);


    rws_file = capture_file_open(H5FILE, get_test_world (), READ_ONLY);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);

    float float_out = -1;
    uint16_t short_out = -1;
    uint64_t long_out = -1;
    err = read_float (rws_file, RESULTS, "float_in", &float_out, 1);
    err |= read_uint16 (rws_file, INPUTS, "short_in", &short_out, 1);
    err |= read_uint64 (rws_file, INPUTS, "long_in", &long_out, 1);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_UINT16(short_in, short_out);
    TEST_ASSERT_EQUAL_UINT64(long_in, long_out);
    TEST_ASSERT_EQUAL_FLOAT(float_in, float_out);

}

struct example_struct
{
    int a;
    double b;
};

TEST(read_and_write_single_w1_3_6, rw_weird_characters)
{
    struct example_struct es;
    struct example_struct * esp = &es;
    esp->a = 3;
    esp->b = 2.0;
    int err = SAVE_DOUBLE(esp->b, rws_file, INPUTS);
    err |= SAVE_INT(esp->a, rws_file, INPUTS);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);
    TEST_ASSERT_EQUAL_INT(0, err);

    int ar;
    double br;
    err |= read_int(rws_file, INPUTS, "esp->a", &ar, 1);
    err |= read_double(rws_file, INPUTS, "esp->b", &br, 1);
    TEST_ASSERT_EQUAL_INT(0, err);

    TEST_ASSERT_EQUAL_DOUBLE(br, esp->b);
    TEST_ASSERT_EQUAL_INT(ar, esp->a);



}

TEST_GROUP(read_write_array_w1_4);

TEST_SETUP(read_write_array_w1_4)
{
    UnityMalloc_StartTest();
    ENSURE_NOT_EXISTS(H5FILE);
    rws_file = capture_file_open(H5FILE, get_test_world(), FAIL_IF_EXISTS);
    TEST_ASSERT_TRUE(rws_file.file_id > 0);
    TEST_ASSERT_TRUE(rws_file.act == FILE_WRITE);
}

TEST_TEAR_DOWN(read_write_array_w1_4)
{
    capture_file_close(&rws_file);
    UnityMalloc_EndTest();
}

TEST(read_write_array_w1_4, read_int)
{

    int size = 100;
    int * data = malloc(sizeof(int) *size);
    int * data_read = malloc(sizeof(int) * size);
    TEST_ASSERT_NOT_NULL(data);
    TEST_ASSERT_NOT_NULL(data_read);

    for (int i=0; i<size; i++)
    {
        data[i] = 10*i;
    }

    int err = SAVE_INT_ARR(data, size, rws_file, SPECIAL);
    TEST_ASSERT_EQUAL_INT(0, err);

    capture_file_close(&rws_file);

    rws_file = capture_file_open(H5FILE, get_test_world(), READ_ONLY);

    err = read_int(rws_file, SPECIAL, "data", data_read, size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_INT_ARRAY(data, data_read, size);

    free(data_read);
    free(data);
}

static void allocate_both(void ** a, void ** b, size_t size)
{
    *a = malloc(size);
    TEST_ASSERT_NOT_NULL(a);
    *b = malloc(size);
    TEST_ASSERT_NOT_NULL(b);
}

static void free_both(void *a, void *b)
{
    free(a);
    free(b);
}

#define INIT_DATA(data, data_read, size, type) {\
    allocate_both((void **) &data, (void **) &data_read, size * sizeof(type));\
    \
    for (int i=0; i<size; i++)\
    {\
        data[i] = (type) (100*i + get_test_world_rank()*10);\
    }\
}

static void reopen()
{
    capture_file_close(&rws_file);
    rws_file = capture_file_open(H5FILE, get_test_world(), READ_ONLY);
}


TEST(read_write_array_w1_4, read_uint16)
{
    uint16_t *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, uint16_t);

    int err = SAVE_UINT16_ARR(data, size, rws_file, SPECIAL);
    TEST_ASSERT_EQUAL_INT(0, err);

    reopen();

    err = read_uint16(rws_file, SPECIAL, "data", data_read, size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(data, data_read, size);

    free_both(data, data_read);
}


TEST(read_write_array_w1_4, read_double)
{
    double *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, double);

    int err = SAVE_DOUBLE_ARR(data, size, rws_file, SPECIAL);
    TEST_ASSERT_EQUAL_INT(0, err);

    reopen();

    err = read_double(rws_file, SPECIAL, "data", data_read, size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(data, data_read, size);

    free_both(data, data_read);
}

TEST(read_write_array_w1_4, read_float)
{
    float *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, float);

    int err = SAVE_FLOAT_ARR(data, size, rws_file, SPECIAL);
    TEST_ASSERT_EQUAL_INT(0, err);

    reopen();

    err = read_float(rws_file, SPECIAL, "data", data_read, size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(data, data_read, size);

    free_both(data, data_read);
}

TEST(read_write_array_w1_4, read_uint64)
{
    uint64_t *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, uint64_t);

    int err = SAVE_UINT64_ARR(data, size, rws_file, SPECIAL);
    TEST_ASSERT_EQUAL_INT(0, err);


    reopen();

    err = read_uint64(rws_file, SPECIAL, "data", data_read, size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_UINT64_ARRAY(data, data_read, size);

    free_both(data, data_read);
}

IGNORE_TEST(read_write_array_w1_4, write_unequal_amounts)
{
    DPRINT("write unequal amounts")
    uint16_t *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, uint16_t);

    int actual_size = (get_test_world_rank()+1)*25;

    int err = SAVE_UINT16_ARR(data, actual_size, rws_file, INPUTS);
    TEST_ASSERT_EQUAL_INT(0, err);
    DPRINT("wrote unequal amounts")

    reopen();

    fenced_system("h5dump " H5FILE, BARRIER_BEFORE | BARRIER_AFTER);

    err = read_uint16(rws_file, INPUTS, "data", data_read, actual_size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(data, data_read, actual_size);

    free_both(data, data_read);
}

IGNORE_TEST(read_write_array_w1_4, write_only_rank_0)
{
    uint16_t *data, *data_read;
    int size = 100;
    INIT_DATA(data, data_read, size, uint16_t);

    int actual_size = (get_test_world_rank() == 0)? size : 0;


    int err = SAVE_UINT16_ARR(data, actual_size, rws_file, INPUTS);
    TEST_ASSERT_EQUAL_INT(0, err);

    reopen();

    fenced_system("h5dump " H5FILE, BARRIER_BEFORE | BARRIER_AFTER);

    err = read_uint16(rws_file, INPUTS, "data", data_read, actual_size);
    TEST_ASSERT_EQUAL_INT(0, err);
    TEST_ASSERT_EQUAL_UINT16_ARRAY(data, data_read, actual_size);

    free_both(data, data_read);
}

TEST_GROUP(convenience_api_w1_3);

TEST_SETUP(convenience_api_w1_3)
{
    UnityMalloc_StartTest();
}

TEST_TEAR_DOWN(convenience_api_w1_3)
{
    UnityMalloc_EndTest();
}

TEST(convenience_api_w1_3, sample_use)
{
    ENSURE_NOT_EXISTS(H5FILE);

    // writing
    int x = 2;
    uint64_t large_number = 1e15;
    int arr[30] = {1,2,3,4,0};
    arr[29] = 100;
    file_open(H5FILE, get_test_world());
    SAVE(x);
    SAVE(large_number);

    change_group(SPECIAL);
    SAVE_ARR(arr, 30);

    file_close();

    // reading
    capture_file_t cf = capture_file_open(H5FILE, get_test_world(), READ_ONLY);

    int x_r;
    uint64_t large_number_r;
    int arr_r[30];

    read_int(cf, INPUTS, "x", &x_r, 1);
    read_uint64(cf, INPUTS, "large_number", &large_number_r, 1);
    read_int(cf, SPECIAL, "arr", arr_r, 30);

    TEST_ASSERT_EQUAL_INT_ARRAY(arr, arr_r, 30);
    TEST_ASSERT_EQUAL_INT(x, x_r);
    TEST_ASSERT_EQUAL_UINT64(large_number, large_number_r);
}