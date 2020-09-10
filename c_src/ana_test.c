
#include "unity_mpi_addon.h"
#include "unity_fixture.h"
#include "unity_capture.h"
#include "soma_config.h"

#include "phase.h"
#include "ana.h"

#include "err_handling.h"

TEST_GROUP(calc_non_bonded_energy_w2);

TEST_SETUP(calc_non_bonded_energy_w2){
    UnityMalloc_StartTest();
}

TEST_TEAR_DOWN(calc_non_bonded_energy_w2) {
    UnityMalloc_EndTest();
}

IGNORE_TEST(calc_non_bonded_energy_w2, example)
{
    capture_file_t cf = capture_file_open("/home/julian/SOMA_FINAL_MASTER/SOMA/build/testing/calc_non_bonded_energy_one_domain_two_ranks_test_data.h5", get_test_world(), READ_ONLY);

    int read_err = 0;

    Phase p;
    p.n_types = 1;
    read_err |= read_unsigned_int(cf, INPUTS, "p->n_types", &(p.n_types), 1);
    read_err |= read_unsigned_int(cf, INPUTS, "p->nx", &(p.nx), 1);
    read_err |= read_unsigned_int(cf, INPUTS, "p->ny", &(p.ny), 1);
    read_err |= read_unsigned_int(cf, INPUTS, "p->nz", &(p.nz), 1);
    read_err |= read_uint64(cf, INPUTS, "p->n_cells_local", &(p.n_cells_local), 1);
    size_t size = (p.n_types*p.nx*p.ny*p.nz);
    p.fields_unified = malloc(size * sizeof(*p.fields_unified));
    p.omega_field_unified = malloc(size * sizeof(*p.omega_field_unified ));
    TEST_ASSERT_TRUE(p.fields_unified != NULL && p.omega_field_unified != NULL);
    read_err |= read_uint16(cf, INPUTS, "p->fields_unified", (p.fields_unified), size);
    read_err |= read_double(cf, INPUTS, "p->omega_field_unified", p.omega_field_unified, size);
    read_err |= read_int(cf, INPUTS, "p->info_MPI.domain_size", &p.info_MPI.domain_size, 1);
    read_err |= read_int(cf, INPUTS, "p->args.N_domains_arg", &(p.args.N_domains_arg), 1);
    read_err |= read_int(cf, INPUTS, "p->local_nx_high", &(p.local_nx_high), 1);
    read_err |= read_int(cf, INPUTS, "p->local_nx_low", &(p.local_nx_low), 1);

    p.info_MPI.SOMA_comm_sim = get_test_world();
    p.info_MPI.sim_rank = get_test_world_rank();
    p.info_MPI.domain_rank = get_test_world_rank();

    soma_scalar_t * non_bonded_energy = malloc(p.n_types * sizeof(soma_scalar_t));
    soma_scalar_t * nb_energy_result = malloc(p.n_types * sizeof(soma_scalar_t));
    TEST_ASSERT_TRUE(non_bonded_energy != NULL && nb_energy_result != NULL);

    calc_non_bonded_energy(&p, non_bonded_energy);


    read_err |= read_double(cf, RESULTS, "non_bonded_energy", nb_energy_result, p.n_types);
    TEST_ASSERT_EQUAL_INT(0, read_err);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(nb_energy_result, non_bonded_energy, p.n_types);

    free(non_bonded_energy);
    free(nb_energy_result);
    free(p.omega_field_unified);
    free(p.fields_unified);
}


