/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg

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

#include <string.h>
#include <assert.h>

#include "io_old.h"
#include "soma_config.h"
#include "phase.h"
#include "soma_memory.h"
#include "mesh.h"

int read_old_config(struct Phase *p, char *const filename)
{
    FILE *inputfile;

    /* open input file */
    if ((inputfile = fopen(filename, "r")) == NULL)
        {
            fprintf(stderr, "Cannot open coord file %s\n", filename);
            return -1;
        }
    if (p->info_MPI.world_size != 1)
        {                       //Make sure this is only called for single core processes.
            fprintf(stderr, "ERROR: convert can only be called with a single MPI-core.\n");
            return -11;
        }

    int iread, NX, NY, NZ, S, HYDRO, LEBD, TCHECK, TWRITE, VCHECK, ZA, ZB, ZC;
    unsigned int NpolyA, NmonoA, NpolyB, NmonoB, NpolyC, NmonoC;
    soma_scalar_t chiN, deltaN, kappaN, stiff, hN1, gN1, fN1, phia, phib, phic,
        D0, LY, LZ, time, F, f, xiN, gdot, forceN, kNcphi, dtcphi, Neq, Nav, Ntdgl, kT, dummy, dt;
    char *Cdummy;
    char cdummy[1024];
    unsigned int i, j;

    Cdummy = &cdummy[0];

    /* read first line of coord.dat */
    iread = fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
                   "chiN = %g; deltaN = %g; kappaN = %g; stiff = %g; hN1 = %g; gN1 = %g; fN1 = %g; phia = %g; phib = %g; phic = %g; D0 = %g; LY = %g; LZ = %g; N = %d %d %d; S = %d;\n",
#else                           //SINGLE_PRECISION
                   "chiN = %lg; deltaN = %lg; kappaN = %lg; stiff = %lg; hN1 = %lg; gN1 = %lg; fN1 = %lg; phia = %lg; phib = %lg; phic = %lg; D0 = %lg; LY = %lg; LZ = %lg; N = %d %d %d; S = %d;\n",
#endif                          //SINGLE_PRECISION
                   &chiN, &deltaN, &kappaN, &stiff, &hN1, &gN1, &fN1, &phia,
                   &phib, &phic, &D0, &LY, &LZ, &NX, &NY, &NZ, &S);

    if (iread != 17)
        {
            fprintf(stderr, "Old io ERROR: %s:%d. (%d)\n", __FILE__, __LINE__, iread);
            return -1;
        }

    /* read second line of coord.dat */
    iread = fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
                   "time =  %g; F = %g ( %g ); f = %g; xiN = %g; dt*N/xi = %g; gdot = %g; forceN = %g; kNcphi = %g; dtcphi = %g; Neq = %g; Nav = %g; Ntdgl = %g; kT = %g; HYDRO = %d; LEBD = %d; TCHECK = %d; TWRITE = %d; VCHECK = %d;\n",
#else                           //SINGLE_PRECISION
                   "time =  %lg; F = %lg ( %lg ); f = %lg; xiN = %lg; dt*N/xi = %lg; gdot = %lg; forceN = %lg; kNcphi = %lg; dtcphi = %lg; Neq = %lg; Nav = %lg; Ntdgl = %lg; kT = %lg; HYDRO = %d; LEBD = %d; TCHECK = %d; TWRITE = %d; VCHECK = %d;\n",
#endif                          //SINGLE_PRECISION
                   &time, &F, &dummy, &f, &xiN, &dt, &gdot, &forceN, &kNcphi,
                   &dtcphi, &Neq, &Nav, &Ntdgl, &kT, &HYDRO, &LEBD, &TCHECK, &TWRITE, &VCHECK);
    if (iread != 19)
        {
            fprintf(stderr, "Old io ERROR: %s:%d. (%d)\n", __FILE__, __LINE__, iread);
            return -2;
        }

    /* read third line of coord.dat - free comment */
    Cdummy = fgets(Cdummy, 1024, inputfile);

    /* read fourth line of coord.dat */
    iread =
        fscanf(inputfile, "%u %u %d  %u %u %d  %u %u %d\n", &NpolyA,
               &NmonoA, &ZA, &NpolyB, &NmonoB, &ZB, &NpolyC, &NmonoC, &ZC);
    if (iread != 9)
        {
            fprintf(stderr, "Old io ERROR: %s:%d. (%d)\n", __FILE__, __LINE__, iread);
            return -3;
        }

    p->reference_Nbeads = NmonoA;

    /* use the values we need */

    p->nx = NX;
    p->ny = NY;
    p->nz = NZ;
    p->n_cells = NX * NY * NZ;
    p->Lx = D0;
    p->Ly = LY;
    p->Lz = LZ;
    p->time = time;

    p->n_polymers = NpolyA;
    p->n_polymers_storage = p->n_polymers;
    assert(p->info_MPI.world_size == 1);        //Make sure this is only called for single core processes.
    p->n_polymers_global = p->n_polymers;
    p->n_types = 2;

    /* Assuming same diffusion constant for all monomers for this input file, but value can vary! */
    p->A = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->A == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    p->A[0] = dt / p->reference_Nbeads;
    p->A[1] = dt / p->reference_Nbeads;

    /* allocate XN interaction matrix */
    p->xn = (soma_scalar_t *) malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    if (p->xn == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }

    /* set parameters */
    p->xn[0 * p->n_types + 0] = kappaN;
    p->xn[1 * p->n_types + 1] = kappaN;
    p->xn[0 * p->n_types + 1] = deltaN;
    p->xn[1 * p->n_types + 0] = deltaN;

    //Create the polymer architecture list.
    //The old formate supports only 1 architecture -- a linear chain.
    p->n_poly_type = 1;
    p->poly_type_offset = (int *)malloc(p->n_poly_type * sizeof(int));
    if (p->poly_type_offset == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    p->poly_type_offset[0] = 0;
    //Length = 1(N) + N(mono) + 4 (+1,-1,+1-1)
    p->poly_arch_length = 1 + NmonoA + 4;

    p->poly_arch = (uint32_t *) malloc(p->poly_arch_length * sizeof(uint32_t));
    if (p->poly_arch == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }

    //Set up the architecture
    p->poly_arch[0] = NmonoA;

    p->poly_arch[0 + 1] = get_info_bl(NmonoA + 4, 0);
    for (unsigned int i = 1; i < NmonoA - 1; i++)
        {
            p->poly_arch[i + 1] = get_info_bl(NmonoA + 3, 0);
        }
    p->poly_arch[NmonoA] = get_info_bl(NmonoA + 1, 0);

    //Last monomer
    int end, offset, bond_type, info;
    bond_type = HARMONIC;
    end = 1;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[NmonoA + 1] = info;
    //First monomer
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[NmonoA + 2] = info;
    //Middle monomers
    end = 0;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[NmonoA + 3] = info;

    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    assert(NmonoA + 4 == p->poly_arch_length - 1);
    p->poly_arch[NmonoA + 4] = info;

    /* allocate space for polymers */
    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if (p->polymers == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }

    unsigned int counter = 0;
    for (i = 0; i < p->n_polymers; i++)
        {
            p->polymers[i].type = 0;
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            counter += N;
        }
    init_soma_memory(&(p->ph.beads), counter, sizeof(Monomer));
    init_soma_memory(&(p->ph.msd_beads), counter, sizeof(Monomer));
#if (ENABLE_MONOTYPE_CONVERSIONS == 1)
    init_soma_memory(&(p->ph.monomer_types), counter, sizeof(uint32_t));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    /* read in chains and set up particle data */
    for (i = 0; i < p->n_polymers; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            p->polymers[i].bead_offset = get_new_soma_memory_offset(&(p->ph.beads), N);
            if (p->polymers[i].bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            p->polymers[i].msd_bead_offset = get_new_soma_memory_offset(&(p->ph.msd_beads), N);
            if (p->polymers[i].msd_bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
#if (ENABLE_MONOTYPE_CONVERSIONS == 1)
            p->polymers[i].monomer_type_offset = get_new_soma_memory_offset(&(p->ph.monomer_types), N);
            if (p->polymers[i].monomer_type_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

            Monomer *beads = p->ph.beads.ptr;
            beads += p->polymers[i].bead_offset;
#if (ENABLE_MONOTYPE_CONVERSIONS == 1)
            uint32_t *monomer_types = (uint32_t *) p->ph.monomer_types.ptr;
            monomer_types += p->polymers[i].monomer_type_offset;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

            /* read the monomers of this chain */
            for (j = 0; j < N; j++)
                {
                    unsigned int tmp;
                    iread = fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
                                   "%g %g %g %u\n",
#else                           //SINGLE_PRECISION
                                   "%lg %lg %lg %u\n",
#endif                          //SINGLE_PRECISION
                                   &(beads[j].x), &(beads[j].y), &(beads[j].z), &tmp);
                    if (iread != 4)
                        {
                            fprintf(stderr, "Old io ERROR: %s:%d. (%d)\n", __FILE__, __LINE__, iread);
                            return -6;
                        }

                    tmp -= 1;   //First type in conf is 1. First type in program is 0
#if (ENABLE_MONOTYPE_CONVERSIONS == 1)
                    monomer_types[j] = tmp;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

                    //! \warning Assumed is that all polymers have the same
                    //! type. This includes, that the particle type
                    //! distribution is for all polymers identical. (The last
                    //! specified in the input format.) If this is not the 
                    //! case choose option ENABLE_MONOTYPE_CONVERSIONS.
                    uint32_t *const type_info = &(p->poly_arch[p->poly_type_offset[p->polymers[i].type] + j + 1]);
                    const int offset_bl = get_bondlist_offset(*type_info);
                    // Overwrite the type info every single step.
                    *type_info = get_info_bl(offset_bl, tmp);
                }
        }
    p->bead_data_read = true;
    p->mt_data_read = true;
    fclose(inputfile);

    p->area51 = NULL;
    p->external_field_unified = NULL;
    p->umbrella_field = NULL;
    p->hamiltonian = SCMF0;
    p->k_umbrella = (soma_scalar_t * const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->k_umbrella == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    memset(p->k_umbrella, 0, p->n_types * sizeof(soma_scalar_t));

    p->harmonic_normb_variable_scale = 1;
    p->cm_a = NULL;             // Deactivate CM movement with the old file format.
    p->pc.deltaMC = 0;
    p->pc.array = NULL;
    p->pc.input_type = NULL;
    p->pc.output_type = NULL;
    p->pc.reaction_end = NULL;
    p->pc.len_reactions = 0;

    p->field_scaling_type = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->field_scaling_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] = 1;

    return 0;
}

int read_old_geometry(struct Phase *p, const char *filename)
{
    if (p->args.N_domains_arg != 1)
        {
            fprintf(stderr, "ERROR: read_old_geometry works only with a single domain.\n");
            return -1;
        }

    /* check if also an old style geometry needs to be included */
    FILE *geofile = fopen(filename, "rb");
    if (geofile != NULL)
        {
            p->external_field_unified = (soma_scalar_t *) malloc(p->n_types * p->n_cells * sizeof(soma_scalar_t *));
            if (p->external_field_unified == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }

            p->area51 = (uint8_t *) malloc(p->n_cells * sizeof(uint8_t));
            if (p->area51 == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
            int x, y, z;
            soma_scalar_t wall, field, dielectric, compute_efield, fixed_potential;;

            // in the original definition external fields are symmetric to the two species
            for (unsigned int k = 0; k < p->n_cells; k++)
                {
#if (SINGLE_PRECISION == 1)
                    const int iread = fscanf(geofile, "%d %d %d %g %g %g %g %g", &x, &y, &z, &wall, &field, &dielectric,
                                             &compute_efield, &fixed_potential);
#else                           //SINGLE_PRECISION
                    const int iread =
                        fscanf(geofile, "%d %d %d %lg %lg %lg %lg %lg", &x, &y, &z, &wall, &field, &dielectric,
                               &compute_efield, &fixed_potential);
#endif                          //SINGLE_PRECISION
                    if (iread != 8)
                        {
                            fprintf(stderr, "Geometry io ERROR: %s:%d. (%d)\n", __FILE__, __LINE__, iread);
                            return -1;
                        }
                    if (wall > 0)
                        {
                            p->area51[k] = 0;
                        }
                    else
                        {
                            p->area51[k] = 1;
                        }
                    p->external_field_unified[cell_to_index_unified(p, k, 0)] = field;
                    p->external_field_unified[cell_to_index_unified(p, k, 1)] = -field;
                }
            fclose(geofile);
        }
    else
        {
            fprintf(stderr, "WARNING: Passed geometry file %s unreadable. Geometry will be ignored.\n", filename);
            p->area51 = NULL;
            p->external_field_unified = NULL;
        }
    return 0;
}

int read_beads0(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    //Distribute the polymers to the different cores.
    uint64_t n_polymers = p->n_polymers_global / p->info_MPI.sim_size;
    if ((unsigned int)p->info_MPI.sim_rank < p->n_polymers_global % p->info_MPI.sim_size)
        n_polymers += 1;
    p->n_polymers = n_polymers;
    p->n_polymers_storage = p->n_polymers;
    uint64_t n_polymer_offset = 0;
#if ( ENABLE_MPI == 1 )
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) & (p->n_polymers), &n_polymer_offset, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;
#endif                          //ENABLE_MPI

    if (p->info_MPI.sim_rank == p->info_MPI.sim_size - 1)
        assert(p->n_polymers + n_polymer_offset == p->n_polymers_global);

    //Polymers (only the local copies).
    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if (p->polymers == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    //Read in polymer related data, only the locally needed info
    unsigned int *const poly_type = (unsigned int *const)malloc(p->n_polymers_storage * sizeof(unsigned int));
    if (poly_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hsize_t hsize_polymers = p->n_polymers;
    hsize_t hsize_polymers_offset = n_polymer_offset;
    hid_t poly_type_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);
    hid_t poly_type_dataset = H5Dopen2(file_id, "/poly_type", H5P_DEFAULT);
    hid_t poly_type_dataspace = H5Dget_space(poly_type_dataset);
    H5Sselect_hyperslab(poly_type_dataspace, H5S_SELECT_SET, &(hsize_polymers_offset), NULL, &(hsize_polymers), NULL);

    if ((status =
         H5Dread(poly_type_dataset, H5T_NATIVE_UINT, poly_type_memspace, poly_type_dataspace, plist_id, poly_type)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(poly_type_memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(poly_type_dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(poly_type_dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            p->polymers[i].type = poly_type[i];
            p->polymers[i].tag = n_polymer_offset + i;  //The old file format doesn't contain tags, so we create some new ones.
        }
    free(poly_type);

    //Allocate space for monomers
    uint64_t counter = 0;
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            counter += N;
        }
    init_soma_memory(&(p->ph.beads), counter, sizeof(Monomer));
    init_soma_memory(&(p->ph.msd_beads), counter, sizeof(Monomer));
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            p->polymers[i].bead_offset = get_new_soma_memory_offset(&(p->ph.beads), N);
            if (p->polymers[i].bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            p->polymers[i].msd_bead_offset = get_new_soma_memory_offset(&(p->ph.msd_beads), N);
            if (p->polymers[i].msd_bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
        }

    if (H5Lexists(file_id, "/beads", H5P_DEFAULT) > 0)
        {
            //Get the data for the polymers
            //get the dimension
            unsigned int max_n_beads = 0;
            for (unsigned int i = 0; i < p->n_poly_type; i++)
                {
                    const unsigned int N = p->poly_arch[p->poly_type_offset[i]];
                    if (max_n_beads < N)
                        max_n_beads = N;
                }

            Monomer *const monomer_data = (Monomer * const)malloc(p->n_polymers * max_n_beads * sizeof(Monomer));
            if (monomer_data == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }

            hsize_t hsize_beads_memspace[2] = { p->n_polymers, max_n_beads };
            hid_t beads_memspace = H5Screate_simple(2, hsize_beads_memspace, NULL);
            hid_t beads_dataset = H5Dopen2(file_id, "/beads", H5P_DEFAULT);
            hid_t beads_dataspace = H5Dget_space(beads_dataset);
            hsize_t hsize_beads_offset[2] = { n_polymer_offset, 0 };
            H5Sselect_hyperslab(beads_dataspace, H5S_SELECT_SET, hsize_beads_offset, NULL, hsize_beads_memspace, NULL);

            hid_t monomer_memtype = get_monomer_memtype();

            if ((status =
                 H5Dread(beads_dataset, monomer_memtype, beads_memspace, beads_dataspace, plist_id, monomer_data)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            for (uint64_t i = 0; i < p->n_polymers; i++)
                {
                    const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
                    for (unsigned int j = 0; j < N; j++)
                        {
                            const Monomer data = monomer_data[i * max_n_beads + j];
                            *(((Monomer *) p->ph.beads.ptr) + p->polymers[i].bead_offset + j) = data;
                            *(((Monomer *) p->ph.msd_beads.ptr) + p->polymers[i].msd_bead_offset + j) = data;
                        }
                }
            free(monomer_data);

            if ((status = H5Tclose(monomer_memtype)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(beads_dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(beads_memspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Dclose(beads_dataset)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }
            p->bead_data_read = true;
        }
    else
        p->bead_data_read = false;
    return status;
}
