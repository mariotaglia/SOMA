/* Copyright (C) 2016-2019 Ludwig Schneider

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

//! \file test.c
//! \brief Implementation of test.h

#include "test.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "rng.h"
#include "io.h"
#include "init.h"
#include "mesh.h"
#include "mc.h"
#include "independent_sets.h"
#include "server.h"


//todo: bring this test back to life by changing write_config_hdf5 to something current
/*
int test_read_write_hdf5(const struct Phase *const p)
{
    int status;
    unsigned int syserror;
    if ((status = write_config_hdf5(p, "/tmp/p1.h5")) != 0)
        {
            fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__, __LINE__);
            return status;
        }
#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

    struct Phase phase2;
    struct Phase *const p2 = &phase2;
    p2->info_MPI = p->info_MPI;
    if ((status = read_config_hdf5 (p2, "/tmp/p1.h5", NULL)) != 0)
        {
            fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__, __LINE__);
            return status;
        }
    init_phase (p2, NULL, NULL);

    if ((status = write_config_hdf5(p2, "/tmp/p2.h5")) != 0)
        {
            fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__, __LINE__);
            return status;
        }
    free_phase(p2);

    //Maybe the system() approach is not the best. But I think for testing purposes this is fine.
    if ((status = system("h5diff /tmp/p1.h5 /tmp/p2.h5")) != 0)
        {
            fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__, __LINE__);
            return status;
        }
    //Clean up files
    syserror = system("rm -f /tmp/p1.h5");
    if (syserror != 0)
        fprintf(stderr, "ERROR: removing p1.h5 failed\n");
    syserror = system("rm -f /tmp/p2.h5");
    if (syserror != 0)
        fprintf(stderr, "ERROR: removing p2.h5 failed\n");

    if (p->info_MPI.world_rank == 0)
        printf("INFO: At t= %d read_write_hdf5 test passed\n", p->time);
    return 0;
}
*/

int test_particle_types(const struct Phase *const p)
{
    for (uint64_t i = 0; i < p->n_poly_type; i++)
        {
            const unsigned int type_offset = p->poly_type_offset[i];
            const unsigned int N = p->poly_arch[type_offset];
            for (unsigned int j = 0; j < N; j++)
                {
                    const uint32_t info_bl = p->poly_arch[type_offset + 1 + j];
                    const unsigned int type = get_particle_type(info_bl);
                    if (type >= p->n_types)
                        {
                            fprintf(stderr,
                                    "ERROR: min. 1 particle has an undefined type. "
                                    "polytype= %u bead= %u type= %u n_types= %u info_bl=%d \n", (unsigned int)i, j,
                                    type, p->n_types, info_bl);
                            return type;
                        }
                }
        }
    if (p->info_MPI.world_rank == 0)
        printf("INFO: At t= %d particle_type test test passed\n", p->time);
    return 0;
}

int test_area51_violation(const struct Phase *const p)
{
    if (p->area51 != NULL)
        {

            for (uint64_t i = 0; i < p->n_polymers; i++)
                {
                    const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
                    for (unsigned int j = 0; j < N; j++)
                        {
                            const uint64_t index = coord_to_index(p, p->polymers[i].beads[j].x,
                                                                  p->polymers[i].beads[j].y,
                                                                  p->polymers[i].beads[j].z);
                            if (index >= p->n_cells_local)
                                {
                                    fprintf(stderr, "ERROR: domain Error %d %d. Particle out of domain %s:%d\n", (int)i,
                                            j, __FILE__, __LINE__);
                                    int x, y, z;
                                    coord_to_cell_coordinate(p, p->polymers[i].beads[j].x, p->polymers[i].beads[j].y,
                                                             p->polymers[i].beads[j].z, &x, &y, &z);
                                    fprintf(stderr, "\t %d %f %f %f\t %d %d %d\t %d %d\t %f %f %f\n",
                                            p->info_MPI.world_rank, p->polymers[i].beads[j].x,
                                            p->polymers[i].beads[j].y, p->polymers[i].beads[j].z, x, y, z,
                                            p->local_nx_low, p->local_nx_high, p->polymers[i].rcm.x,
                                            p->polymers[i].rcm.y, p->polymers[i].rcm.z);
                                    return index;
                                }
                            if (p->area51[index] == 1)
                                {
                                    fprintf(stderr, "ERROR: particle %u %u is in a forbidden area.\n", (unsigned int)i,
                                            j);
                                    return index;
                                }

                        }
                }
        }
    if (p->info_MPI.world_rank == 0)
        printf("INFO: At t= %d area51 violation test passed\n", p->time);
    return 0;
}

int test_independet_sets(const struct Phase *const p)
{
    if (p->sets == NULL)
        return 0;

    int ret = 0;
    unsigned int divergence = 0;
    for (unsigned int poly_type = 0; poly_type < p->n_poly_type; poly_type++)
        {
            struct IndependetSets *const set = &(p->sets[poly_type]);
            unsigned int largest_set = 0;
            unsigned int smallest_set = UINT_MAX;
            //printf("PolyType %d:\n",poly_type);
            for (unsigned int set_id = 0; set_id < set->n_sets; set_id++)
                {
                    //printf("\tSetId: %d\n\t\t",set_id);
                    if (set->set_length[set_id] > largest_set)
                        largest_set = set->set_length[set_id];
                    if (set->set_length[set_id] < smallest_set)
                        smallest_set = set->set_length[set_id];
                    for (unsigned int i = 0; i < set->set_length[set_id]; i++)
                        {
                            const unsigned int pi = set->sets[set_id * set->max_member + i];
                            //printf(" %d ",pi);
                            unsigned int n_neigh = 0;
                            for (unsigned int j = 0; j < set->set_length[set_id]; j++)
                                {
                                    const unsigned int pj = set->sets[set_id * set->max_member + j];
                                    const int start =
                                        get_bondlist_offset(p->poly_arch[p->poly_type_offset[poly_type] + pj + 1]);
                                    if (start > 0)
                                        {
                                            unsigned int i = start;
                                            unsigned int end;
                                            do
                                                {
                                                    const uint32_t info = p->poly_arch[i++];
                                                    end = get_end(info);
                                                    const int offset = get_offset(info);
                                                    const unsigned n_id = pj + offset;
                                                    if (n_id == pi)
                                                        n_neigh++;
                                            } while (end == 0);
                                        }
                                }
                            if (p->info_MPI.domain_rank == 0 && n_neigh > 0)
                                fprintf(stderr, "ERROR: mono %d of poly_type %d has neighbors in its set %d.\n",
                                        pi, poly_type, set_id);
                            ret += n_neigh;
                        }
                    //printf("\n");
                }
            if (largest_set - smallest_set > divergence)
                divergence = largest_set - smallest_set;
        }

    if (p->info_MPI.world_rank == 0 && ret == 0)
        printf("INFO: Test independet sets passed with max divergence %d.\n", divergence);
    return ret;
}

int test_area51_exact(const struct Phase *const p)
{
    unsigned int violations = 0;
    if (p->area51 != NULL)
        {
            for (uint64_t i = 0; i < p->n_polymers; i++)
                {
                    const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
                    for (unsigned int j = 0; j < N - 1; j++)
                        {
                            const Monomer a = p->polymers[i].beads[j];
                            const Monomer b = p->polymers[i].beads[j + 1];
                            Monomer dx;
                            dx.x = b.x - a.x;
                            dx.y = b.y - a.y;
                            dx.z = b.z - a.z;
                            const bool ok = possible_move_area51(p, a.x, a.y, a.z, dx.x, dx.y, dx.z, false);
                            if (!ok)
                                {
                                    //printf("violate %d %d %f %f %f\t%f %f %f\n", i, j, a.x, a.y, a.z, b.x, b.y, b.z);
                                    violations += 1;
                                }
                        }
                }
        }
#if ( ENABLE_MPI == 1 )
    MPI_Allreduce(MPI_IN_PLACE, &violations, 1, MPI_UNSIGNED, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
    if (p->info_MPI.sim_rank == 0)
        {
            if (violations == 0)
                printf("INFO: At t= %d area51 exact test passed\n", p->time);
            else
                printf
                    ("WARNING: At t= %d area51 exact test **FAILED** with %d violations. If area51 contains 'steps' this is not necessarily a problem.\n",
                     p->time, violations);
        }
    return violations;
}

int test_chains_in_domain(struct Phase *const p)
{
    if (p->args.N_domains_arg == 1)
        return 0;
    update_polymer_rcm(p);
    update_self_phase(p, 0);
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    unsigned int violations = 0;
    for (unsigned int i = 0; i < p->n_polymers; i++)
        {
            const int fold = rint(p->polymers[i].rcm.x / p->Lx);
            soma_scalar_t x = p->polymers[i].rcm.x - fold * p->Lx;
            if (x < 0)
                x += p->Lx;
            const unsigned int target_domain = x / p->Lx * p->args.N_domains_arg;
            if (target_domain != my_domain)
                {
                    violations += 1;
                    printf("ERROR: world rank %d owns a chain outside of its domain %f \n",
                           p->info_MPI.world_rank, p->polymers[i].rcm.x);
                }
        }
#if ( ENABLE_MPI == 1 )
    MPI_Allreduce(MPI_IN_PLACE, &violations, 1, MPI_UNSIGNED, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
    if (p->info_MPI.sim_rank == 0)
        {
            if (violations == 0)
                printf("INFO: At t= %d sim %d chain domain test passed.\n", p->time,
                       p->info_MPI.world_rank / p->info_MPI.sim_size);
            else
                printf("WARNING: At t= %d sim %d **FAILED** the domain chain test with %d violations.\n", p->time,
                       p->info_MPI.world_rank / p->info_MPI.sim_rank, violations);
        }
    return violations;
}

void test_field_sending_consistency(const struct Phase *p, struct sim_rank_info * sim_inf, uint64_t min_cell, uint64_t max_cell, uint64_t field_size)
{

    MESSAGE_ASSERT(field_size == max_cell - min_cell + 1, "field_size inconsistent with min and max cell");

    uint64_t n_cells_sending, local_sum; // sending local info
    uint64_t n_cells_sent_total=0, global_sum=0, max_cell_left, min_cell_right; //receiving global/nonlocal info

    const int req_num = 6;
    MPI_Request reqs[req_num];
    for (int i=0; i<req_num; i++)
        reqs[i] = MPI_REQUEST_NULL;

    if (sim_inf->field_comm.comm != MPI_COMM_NULL)
        {
            n_cells_sending = field_size;
            local_sum = 0;
            for (uint64_t i=min_cell; i <= max_cell; i++)
                {
                    local_sum += i;
                }
            const int left_rank =
                (((sim_inf->my_domain - 1) +
                  p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
                p->info_MPI.domain_rank;
            const int right_rank =
                (((sim_inf->my_domain + 1) +
                  p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
                p->info_MPI.domain_rank;
            MPI_Isend(&min_cell, 1, MPI_UINT64_T, left_rank, 0, sim_inf->sim_comm.comm, &reqs[2]);
            MPI_Isend(&max_cell, 1, MPI_UINT64_T, right_rank, 0, sim_inf->sim_comm.comm, &reqs[3]);
            MPI_Irecv(&min_cell_right, 1, MPI_UINT64_T, right_rank, 0, sim_inf->sim_comm.comm, &reqs[4]);
            MPI_Irecv(&max_cell_left, 1, MPI_UINT64_T, left_rank, 0, sim_inf->sim_comm.comm, &reqs[5]);
        }
    else
        {
            n_cells_sending = 0;
            local_sum = 0;
        }
    uint64_t n_cells_expected = p->nx*p->ny*p->nz;
    MPI_Ireduce(&n_cells_sending, &n_cells_sent_total, 1, MPI_UINT64_T, MPI_SUM, 0, sim_inf->sim_comm.comm, &reqs[0]);
    MPI_Ireduce(&local_sum, &global_sum, 1, MPI_UINT64_T, MPI_SUM, 0, sim_inf->sim_comm.comm, &reqs[1]);

    if (sim_inf->sim_comm.rank == 0)
        {
            MESSAGE_ASSERT(n_cells_sending == n_cells_expected, "Index calculation must be very wrong, not sending the correct amount of cells");
            MESSAGE_ASSERT((n_cells_expected-1)*(n_cells_expected-2)/2 == local_sum, "Index calculation gives wrong checksum.");
        }       MPI_Waitall(6, reqs, MPI_STATUSES_IGNORE);

    if (sim_inf->my_domain == 0)
        MESSAGE_ASSERT(max_cell_left == n_cells_expected - 1, "Mising cells in domain on the very right");
    else
        MESSAGE_ASSERT(max_cell_left == min_cell -1, "Missing cells or overlap in domain decomp");

    if (sim_inf->my_domain == p->args.N_domains_arg - 1)
        MESSAGE_ASSERT(min_cell_right == 0, "Missing cells in domain on the very left");
    else
        MESSAGE_ASSERT(min_cell_right == max_cell + 1, "Missing cells or overlap in domain decomp");

}