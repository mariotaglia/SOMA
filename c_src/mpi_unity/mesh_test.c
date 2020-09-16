#include "unity_fixture.h"
#include <memory.h>
#include <assert.h>
#include "err_handling.h"

#include "../mesh.h"
#include "../phase.h"
#include "../init.h"

static int init_for_acc_if_necessary(Phase * p)
{

#ifdef _OPENACC

    p->args.gpus_given = 1;
    p->args.gpus_arg = 4;

    p->args.only_gpu_given = 0;
    p->args.only_gpu_arg = 0;

    p->args.omp_threads_given = 0;
    p->args.omp_threads_arg = 0;

    return set_openacc_devices(p);
#else
    return 0;
#endif

}

static void phase_copyin_if_necessary(Phase * p)
{
    if (p->present_on_device)
        {
            fprintf(stderr, "WARNING: copyin of phase that is already on device!");
        }
#ifdef _OPENACC
#pragma acc enter data copyin(p[0:1])
#pragma acc enter data copyin(p->fields_unified[0:p->n_types*p->n_cells_local])
    if (p->args.N_domains_arg > 1)
        {
//#pragma acc enter data copyin(p->right_tmp_buffer[0:p->args.domain_buffer_arg * p->ny * p->nz])
//#pragma acc enter data copyin(p->left_tmp_buffer[0:p->args.domain_buffer_arg * p->ny * p->nz])
        }
#endif
    p->present_on_device = true;
}

static void phase_copyout_if_necessary(Phase * p)
{
    if (!p->present_on_device)
        {
            fprintf(stderr, "WARNING: copyout of phase that is not on device!");
        }
#ifdef _OPENACC
DPRINT("copyout fields_unified")
#pragma acc update host(p->fields_unified[0:p->n_types*p->n_cells_local])
DPRINT("copyout phase")
#pragma acc update host(p[0:1])
DPRINT("fields_unified copied")
//#pragma acc exit data delete(p->fields_unified)
//#pragma acc exit data delete(p[0:1])

DPRINT("delete successful")
    if (p->args.N_domains_arg > 1)
        {
//#pragma acc exit data copyout(p->right_tmp_buffer[0:p->args.domain_buffer_arg * p->ny * p->nz])
//#pragma acc exit data copyout(p->left_tmp_buffer[0:p->args.domain_buffer_arg * p->ny * p->nz])
        }
#endif
    p->present_on_device = false;
}

#define INIT_COMM_DENSITY_TEST \
    Phase p = init_communicate_density_fields_phase(n_domains, nx, ny, nz, ntypes, domain_buffer_arg);\
    uint16_t * expected = (uint16_t *) malloc(p.n_types * p.n_cells_local*sizeof(uint16_t));\
    TEST_ASSERT_NOT_NULL(expected);

#define LOOP_DENSITY_FIELD(body) \
for (unsigned type=0; type<p.n_types; type++)\
        for (int x= p.local_nx_low; x<p.local_nx_high; x++)\
            for (unsigned y=0; y<p.ny; y++)\
                for (unsigned z=0; z<p.nz; z++)\
                    {\
                        body\
                    }

#define FINALIZE_COMM_DENSITY_TEST \ 
    phase_copyin_if_necessary(&p);\
    DPRINT("calling comm");\
    communicate_density_fields(&p);\
    DPRINT("testing equality");\
    TEST_ASSERT_EQUAL_UINT16_ARRAY(expected, p.fields_unified, p.n_types*p.n_cells_local);\
    DPRINT("equality confirmed");\
    DPRINT("calling copyout");\
    phase_copyout_if_necessary(&p);\
    DPRINT("copyout successful");\
    free_communicate_density_field_phase(&p);\
    free(expected);


TEST_GROUP(communicate_density_fields_w2);

TEST_SETUP(communicate_density_fields_w2)
{
#if (ENABLE_DOMAIN_DECOMPOSITION != 1)
        TEST_IGNORE_MESSAGE("need to compile with domain decomp to test communicate density fields\n");
#endif
    UnityMalloc_StartTest();
}

TEST_TEAR_DOWN(communicate_density_fields_w2){
    UnityMalloc_EndTest();
}

static int get_my_domain(const Phase * p)
{
    int domain_size = get_test_world_size() / p->args.N_domains_arg;
    TEST_ASSERT_EQUAL_INT(p->info_MPI.domain_size, domain_size);
    return get_test_world_rank() / domain_size;
}

// initializes info_mpi for a domain decomposition using test_world as the sim_comm
static void init_info_mpi(Info_MPI * inf, int n_domains)
{
    int world_size = get_test_world_size();
    TEST_ASSERT_MESSAGE(world_size%n_domains == 0, "number of domains does not evenly divide number of ranks (there must be an error in the test code in " __FILE__ ")");

    inf->SOMA_comm_sim = get_test_world();
    inf->sim_size = get_test_world_size();
    inf->sim_rank = get_test_world_rank();

    inf->SOMA_comm_world = get_test_world();
    inf->world_size = get_test_world_size();
    inf->world_rank = get_test_world_rank();

    inf->domain_size = get_test_world_size() / n_domains;
    int my_domain = get_test_world_rank() / inf->domain_size;
    MPI_Comm_split(inf->SOMA_comm_sim, my_domain, 0, &(inf->SOMA_comm_domain));

    MPI_Comm_size(inf->SOMA_comm_domain, &inf->domain_size);
    MPI_Comm_rank(inf->SOMA_comm_domain, &inf->domain_rank);
}

// initializes metadata for fields_unified and mallocs the required arrays which are: fields_unified and left/right_temp_buffer
void init_field_metadata (Phase *p, int nx, int ny, int nz, int ntypes, int domain_buffer_arg, int n_domains)
{
    p->nx = nx;
    p->ny = ny;
    p->nz = nz;
    p->n_types = ntypes;
    p->n_cells = p->nx*p->ny*p->nz;

    TEST_ASSERT_EQUAL_INT(n_domains, (p->args.N_domains_arg));

    int my_domain = get_my_domain(p);

    p->args.domain_buffer_arg = domain_buffer_arg;

    p->local_nx_low = (p->nx / n_domains * my_domain) - domain_buffer_arg;
    p->local_nx_high = p->local_nx_low + p->nx / p->args.N_domains_arg + 2*p->args.domain_buffer_arg;
    p->n_cells_local = (p->local_nx_high - p->local_nx_low) * p->ny*p->nz;
    p->fields_unified = (uint16_t *) malloc(p->n_cells_local * p->n_types * sizeof(uint16_t));
    TEST_ASSERT_NOT_NULL(p->fields_unified);
    if (n_domains > 1)
        {
            p->left_tmp_buffer = (uint16_t *) malloc(p->args.domain_buffer_arg * p->ny * p->nz * sizeof(uint16_t));
            TEST_ASSERT_NOT_NULL(p->left_tmp_buffer);
            p->right_tmp_buffer = (uint16_t *) malloc(p->args.domain_buffer_arg * p->ny * p->nz * sizeof(uint16_t));
            TEST_ASSERT_NOT_NULL(p->right_tmp_buffer);
        }
    else
        {
            p->left_tmp_buffer = NULL;
            p->right_tmp_buffer = NULL;
        }
}

static Phase init_communicate_density_fields_phase(int n_domains, int nx, int ny, int nz, int ntypes, int domain_buffer_arg)
{

    Phase p;
    memset(&p, -1, sizeof(p));

    p.args.N_domains_arg = n_domains;
    init_info_mpi(&p.info_MPI, n_domains);

    int err = init_for_acc_if_necessary(&p);
    TEST_ASSERT_EQUAL_INT(0, err);

    init_field_metadata (&p, nx, ny, nz, ntypes, domain_buffer_arg, n_domains);

    // dont run mpi_divergence as it is hard to test and would require more parameters
    p.args.load_balance_arg = 0;
    p.present_on_device = false;

    return p;
}

static void free_communicate_density_field_phase(Phase * p)
{
    free(p->right_tmp_buffer); p->right_tmp_buffer = NULL;
    free(p->left_tmp_buffer); p->left_tmp_buffer = NULL;
    free(p->fields_unified); p->fields_unified = NULL;
}


TEST(communicate_density_fields_w2, 1domain_2ranks)
{

    int n_domains = 1;
    int nx , ny , nz;
    nx = ny = nz = 5;
    int ntypes = 2;
    int domain_buffer_arg = 0;

DPRINT("init_comm");
    INIT_COMM_DENSITY_TEST

DPRINT("loop");
    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;

            uint16_t val = y + z*10;
            p.fields_unified[index] = val + p.info_MPI.domain_rank;
            expected[index] = 2*val + 1;
        )
	DPRINT("finalize");

    FINALIZE_COMM_DENSITY_TEST
}

TEST(communicate_density_fields_w2, 2domains_1rankperdomain_3types)
{
    int n_domains = 2;
    int nx = 10, ny = 9, nz = 7;
    int ntypes = 3;
    int domain_buffer_arg = 2;

    INIT_COMM_DENSITY_TEST

    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;
            assert(index != UINT64_MAX);
            assert(index < p.n_types*p.n_cells_local);

            p.fields_unified[index] = 1;
            if (x < p.local_nx_low + p.args.domain_buffer_arg*2 || x >= (p.local_nx_high - domain_buffer_arg*2))
                expected[index] = 2;
            else
                expected[index] = 1;
        )

    FINALIZE_COMM_DENSITY_TEST

}


TEST_GROUP(communicate_density_fields_w8);

TEST_SETUP(communicate_density_fields_w8)
{
#if (ENABLE_DOMAIN_DECOMPOSITION != 1)
    TEST_IGNORE_MESSAGE("need to compile with domain decomp to test communicate density fields\n");
#endif
    UnityMalloc_StartTest();
}

TEST_TEAR_DOWN(communicate_density_fields_w8){
    UnityMalloc_EndTest();
}



TEST(communicate_density_fields_w8, 2domains_4ranksperdomain_2types)
{
    int n_domains = 2;
    int nx = 20, ny = 9, nz = 7;
    int ntypes = 3;
    int domain_buffer_arg = 3;

    INIT_COMM_DENSITY_TEST


    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;
            assert(index != UINT64_MAX);
            assert(index < p.n_types*p.n_cells_local);

            p.fields_unified[index] = z+y;

            if (x < p.local_nx_low + p.args.domain_buffer_arg*2 || x >= (p.local_nx_high - domain_buffer_arg*2))
                expected[index] = (z+y)*8;
            else
                expected[index] = (z+y)*4;

        )

    FINALIZE_COMM_DENSITY_TEST
}

TEST(communicate_density_fields_w8, 8domains)
{
    int n_domains = 8;
    int nx = 80, ny = 5, nz = 3;
    int ntypes = 1;
    int domain_buffer_arg = 1;

    INIT_COMM_DENSITY_TEST

    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;
            p.fields_unified[index] = get_my_domain(&p);

            if (x < p.local_nx_low + p.args.domain_buffer_arg*2)
                expected[index] = get_my_domain(&p) + ((get_my_domain(&p) -1 + n_domains)%n_domains);
            else if (x >= (p.local_nx_high - domain_buffer_arg*2))
                expected[index] = get_my_domain(&p) + ((get_my_domain(&p) + 1)%n_domains);
            else
                expected[index] = get_my_domain(&p);
        )

    FINALIZE_COMM_DENSITY_TEST
}

TEST(communicate_density_fields_w8, 4domains_2rankseach_typesaredifferent)
{
    int n_domains = 4;
    int nx = 20, ny = 10, nz = 7;
    int ntypes = 2;
    int domain_buffer_arg = 2;

    INIT_COMM_DENSITY_TEST

    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;

            uint16_t  val =  type + (x+p.nx)%p.nx;
            p.fields_unified[index] = val;

            if (x < p.local_nx_low + p.args.domain_buffer_arg*2 || x >= (p.local_nx_high - domain_buffer_arg*2))
                expected[index] = val*4;
            else
                expected[index] = val*2;
        )

    FINALIZE_COMM_DENSITY_TEST
}

TEST(communicate_density_fields_w8, ranks_in_domain_are_different)
{
    int n_domains = 2;
    int nx = 20, ny = 10, nz = 10;
    int ntypes = 2;
    int domain_buffer_arg = 2;

    INIT_COMM_DENSITY_TEST

    LOOP_DENSITY_FIELD
        (
            uint64_t index = cell_coordinate_to_index(&p, x, y, z) + type*p.n_cells_local;

            uint16_t xt = (x+p.nx)%p.nx;
            uint16_t yt = (y+p.ny)%p.ny;
            uint16_t zt = (z+p.nz)%p.nz;
            uint16_t val[4];
            val[0] = type;
            val[1] = 3 + (xt < 4);
            val[2] = zt * yt + xt;
            val[3] = get_my_domain(&p);

            p.fields_unified[index] = val[p.info_MPI.domain_rank];

            uint16_t first_three_sum = val[0] + val[1] + val[2];
            uint16_t my_domain = val[3];
            if (x < p.local_nx_low + p.args.domain_buffer_arg*2)
                expected[index] = 2*first_three_sum + my_domain + ((my_domain - 1 + n_domains)%n_domains);
            else if (x >= (p.local_nx_high - domain_buffer_arg*2))
                expected[index] = 2*first_three_sum + my_domain + ((my_domain+1)%n_domains);
            else
                expected[index] = first_three_sum + val[3];

        )

    FINALIZE_COMM_DENSITY_TEST
}




