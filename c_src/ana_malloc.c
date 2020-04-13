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


#include "ana_malloc.h"
#include "phase.h"
#include <stdbool.h>
#include "ana.h"
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <stdio.h>

_Static_assert(sizeof(uint64_t) == sizeof(double), "double needs to be 8 bytes for now.");
_Static_assert(sizeof(uint64_t) == sizeof(soma_scalar_t), "for now, only builds with somascalar_t being double is supported");


union elem{
    uint64_t ival;
    double fval;
};

//"private" methods not to be included in the header
void init_sizes(uint64_t n_poly_type , uint64_t n_types, unsigned int q_size_dynamical, unsigned int q_size_static);
int rg_op(union elem *in, union elem *io);
int dvar_op(union elem *in, union elem *io);
int static_sf_op(union elem *in, union elem *io);
int dynamical_sf_op(union elem *in, union elem *io);
int re_op(union elem *in, union elem *io) ;
int anisotropy_op(union elem *in, union elem *io) ;
int msd_op(union elem *in, union elem *io);
int acc_ratio_op(union elem *in, union elem *io);
int bonded_energy_op(union elem *in, union elem *io);
int non_bonded_energy_op(union elem *in, union elem *io);
void sum_doubles(const double *in, double *io, size_t len);
void sum_ints(const uint64_t *in, uint64_t *io, size_t len);
void ana_reduce_func(void *inraw, void *inoutraw, int *len, MPI_Datatype *type);
void * get_pointer_from_offset(int offset, char * name);
void append_block(size_t *offset, int *block, size_t block_size, ssize_t *out_offset,
                    int *block_counts, MPI_Datatype *types, MPI_Aint *displs, MPI_Datatype type);

// cause  a compile-time-error if expression is false
// explanation: https://stackoverflow.com/questions/9229601/what-is-in-c-code
#define BUILD_BUG_OR_ZERO(e) (sizeof(struct { int:-!!(e); }));

#define BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE(a, b) (((a)<0) == ((b)<0))
//casting to void == ignore the value.
// This suppresses unused variable warning when intentionally not using x.
#define UNUSED(x) (void)(x)



// the send and receive buffers (operation is in-place)
// is NULL iff it doesn`t point to a valid buffer (ana_free resets it to null)
union elem *sendbuffer = NULL;
union elem *recvbuffer = NULL;

// MPI-Datatype of the buffers. Gets created and commited with ana_malloc and freed with ana_free.
MPI_Datatype buff_mpi_dtype;

// gets created once during the first call to ana_malloc, does not change afterwards
// use flag to achieve this
// operator is always assumed to be associative and commutative
MPI_Op ana_op;
bool op_is_valid = false;

// for each observable that needs to reduce something, here we have the
// offsets, so that the values can be found on the stack-array starting at stack+offset
// -1 means the value is not being used (the corresponding observable will not be calculated
ssize_t re_vals = -1; size_t re_vals_size;
ssize_t re_counter = -1; size_t re_counter_size;

ssize_t dvar_val = -1;

ssize_t rg_vals = -1; size_t  rg_vals_size;
ssize_t rg_counter = -1; size_t rg_counter_size;

ssize_t anisotropy_vals = -1; size_t anisotropy_vals_size;
ssize_t anisotropy_counter = -1; size_t anisotropy_counter_size;

ssize_t msd_vals = -1; size_t msd_vals_size;
ssize_t msd_counter = -1; size_t msd_counter_size;

// for acc_ratio
ssize_t accepts = -1;
ssize_t moves = -1;

ssize_t non_bonded_energy_val = -1; size_t non_bonded_energy_val_size;

ssize_t bonded_energy_val = -1; size_t bonded_energy_val_size;

ssize_t dynamical_sf_val = -1; size_t dynamical_sf_val_size;

ssize_t static_sf_val = -1; size_t static_sf_val_size;




void ana_reduce_func(void *inraw, void *inoutraw, int *len, MPI_Datatype *type){

    assert(*len==1); // in the context of this module, length must always be one.
    assert(*type == buff_mpi_dtype); //this operation is not defined on any other datatype.
    UNUSED(len); UNUSED(type); // make the compiler happy

    union elem *in = inraw;
    union elem *io = inoutraw;

    //Re
    assert(BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE((re_vals), (re_counter)));
    if (re_vals >= 0) {
        re_op(in, io);
    }

    //dvar
    if (dvar_val >= 0){
        dvar_op(in, io);
    }

    //Rg
    assert(BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE((rg_vals), (rg_counter)));
    if (rg_vals >= 0){
        rg_op(in, io);
    }

    //anisotropy
    assert(BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE((anisotropy_vals),(anisotropy_counter)));
    if (anisotropy_vals >= 0){
        anisotropy_op(in, io);
    }


    //msd
    assert( BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE((msd_vals), (msd_counter)) );
    if (msd_vals >= 0){
        msd_op(in, io);
    }

    //accratio
    assert(BOTH_NEGATIVE_OR_BOTH_NONNEGATIVE((moves), (accepts)));
    if (moves >= 0){
        acc_ratio_op(in, io);
    }

    // nonbonded energy
    if (non_bonded_energy_val >= 0){
        non_bonded_energy_op(in, io);
    }

    //bonded energy
    if (bonded_energy_val >= 0){
        bonded_energy_op(in, io);
    }

    // dynamical sf
    if (dynamical_sf_val >= 0){
        dynamical_sf_op(in, io);
    }

    // static sf
    if (static_sf_val >= 0){
        static_sf_op(in, io);
    }


}

void sum_ints(const uint64_t *in, uint64_t *io, size_t len){
    for (size_t i=0; i<len; i++){
        io[i] += in[i];
    }
}

void sum_doubles(const double *in, double *io, size_t len){
    for (size_t i = 0; i < len; i++) {
        io[i] += in[i];
    }
}

int non_bonded_energy_op(union elem *in, union elem *io){
    double * res_in = &(in[non_bonded_energy_val].fval);
    double * res_io = &(io[non_bonded_energy_val].fval);
    sum_doubles(res_in, res_io, non_bonded_energy_val_size);
    return MPI_SUCCESS;
}

int bonded_energy_op(union elem *in, union elem *io){
    double * res_in = &(in[bonded_energy_val].fval);
    double * res_io = &(io[bonded_energy_val].fval);
    sum_doubles(res_in, res_io, bonded_energy_val_size);
    return MPI_SUCCESS;
}

int acc_ratio_op(union elem *in, union elem *io){
    uint64_t * acc_in = &(in[accepts].ival);
    uint64_t * acc_io = &(io[accepts].ival);
    uint64_t * mov_in = &(in[moves].ival);
    uint64_t * mov_io = &(io[moves].ival);
    *mov_io += *mov_in;
    *acc_io += *acc_in;
    return MPI_SUCCESS;
}

int msd_op(union elem *in, union elem *io){
    double * in_res = &(in[msd_vals].fval);
    double * io_res = &(io[msd_vals].fval);
    uint64_t * in_counter = &(in[msd_counter].ival);
    uint64_t * io_counter = &(in[msd_counter].ival);
    sum_doubles(in_res, io_res, msd_vals_size);
    sum_ints(in_counter, io_counter, msd_counter_size);
    return MPI_SUCCESS;
}

int anisotropy_op(union elem *in, union elem *io) {
    double * in_res = &(in[anisotropy_vals].fval);
    double * io_res = &(io[anisotropy_vals].fval);
    uint64_t * in_counter = &(in[anisotropy_counter].ival);
    uint64_t * io_counter = &(in[anisotropy_counter].ival);

    sum_doubles(in_res, io_res, anisotropy_vals_size);
    sum_ints(in_counter, io_counter, anisotropy_counter_size);

    return MPI_SUCCESS;
}


int re_op(union elem *in, union elem *io) {
    double * in_res = &(in[re_vals].fval);
    double * io_res = &(io[re_vals].fval);
    uint64_t * in_counter = &(in[re_counter].ival);
    uint64_t * io_counter = &(in[re_counter].ival);

    sum_doubles(in_res, io_res, re_vals_size);
    sum_ints(in_counter, io_counter, re_counter_size);

    return MPI_SUCCESS;
}

int dynamical_sf_op(union elem *in, union elem *io){
    double * in_res = &(in[dynamical_sf_val].fval);
    double * io_res = &(io[dynamical_sf_val].fval);
    sum_doubles(in_res, io_res, dynamical_sf_val_size);
    return MPI_SUCCESS;
}

int static_sf_op(union elem *in, union elem *io){
    double * in_res = &(in[static_sf_val].fval);
    double * io_res = &(io[static_sf_val].fval);
    sum_doubles(in_res, io_res, static_sf_val_size);
    return MPI_SUCCESS;
}

int dvar_op(union elem *in, union elem *io){
    double * in_dvar = &(in[dvar_val].fval);
    double * io_dvar = &(io[dvar_val].fval);
    *io_dvar += *in_dvar;
    return MPI_SUCCESS;
}

int rg_op(union elem *in, union elem *io){
    double * in_res = &(in[rg_vals].fval);
    double * io_res = &(io[rg_vals].fval);
    uint64_t * in_counter = &(in[rg_counter].ival);
    uint64_t * io_counter = &(in[rg_counter].ival);

    sum_doubles(in_res, io_res, rg_vals_size);
    sum_ints(in_counter, io_counter, rg_counter_size);

    return MPI_SUCCESS;
}


// receives n_poly_type, the number of polymer architectures
// n_types, the number of monomer types
// also accesses global constant NUMBER_SOMA_BOND_TYPES
void init_sizes(uint64_t n_poly_type , uint64_t n_types, unsigned int q_size_dynamical, unsigned int q_size_static){

    //Re
    re_vals_size = 4*n_poly_type;
    re_counter_size = n_poly_type;

    //Rg
    rg_vals_size = 4*n_poly_type;
    rg_counter_size = n_poly_type;

    //anisotropy
    anisotropy_vals_size = 6*n_poly_type;
    anisotropy_counter_size = n_poly_type;

    //msd
    msd_vals_size = 8*n_poly_type;
    msd_counter_size = 2*n_poly_type;

    //non_bonded_energy
    non_bonded_energy_val_size = n_types;

    //bonded energies
    bonded_energy_val_size = NUMBER_SOMA_BOND_TYPES;

    //structure factors
    dynamical_sf_val_size = q_size_dynamical*n_poly_type*n_types*n_types;
    static_sf_val_size = q_size_static*n_poly_type*n_types*n_types;
}

int ana_malloc(uint64_t n_poly_type, uint64_t n_types, uint64_t q_size_dynamical, uint64_t q_size_static, const bool *needToDo){


    init_sizes(n_poly_type, n_types, q_size_dynamical, q_size_static);

    size_t offset = 0; // indicates current position in the buffers (size in units of 8byte-doubles
    int block = 0; // how many arrays are in the buffers so far

    const int num_obs = 13; //todo: merge this with SOMA_NUM_OBS from ana_info.c

    // information for creating the datatype
    // we need some maximumm length, so we choose 2 blocks per observable (but this is higher than needed)
    int block_counts[num_obs*2];
    MPI_Datatype types[num_obs*2];
    MPI_Aint displs[num_obs*2];

    // datatype will get created in global variable ana_mpi_dtype

    if (needToDo[Re]){

        append_block(&offset, &block, re_vals_size, &re_vals,
                block_counts, types, displs, MPI_DOUBLE);

        append_block(&offset, &block, re_counter_size, &re_counter,
                block_counts, types, displs, MPI_INT64_T);

    }

    if (needToDo[dvar]){
        append_block(&offset, &block, 1, &dvar_val,
                     block_counts, types, displs, MPI_DOUBLE);
	fflush(stdin);
    }

    if (needToDo[Rg]){

        append_block(&offset, &block, rg_vals_size, &rg_vals,
                     block_counts, types, displs, MPI_DOUBLE);

        append_block(&offset, &block, rg_counter_size, &rg_counter,
                     block_counts, types, displs, MPI_INT64_T);
    }

    if (needToDo[b_anisotropy]) {

        append_block(&offset, &block, anisotropy_vals_size, &anisotropy_vals,
                     block_counts, types, displs, MPI_DOUBLE);

        append_block(&offset, &block, anisotropy_counter_size, &anisotropy_counter,
                     block_counts, types, displs, MPI_INT64_T);
    }

    if (needToDo[MSD]){
        append_block(&offset, &block, msd_vals_size, &msd_vals,
                     block_counts, types, displs, MPI_DOUBLE);

        append_block(&offset, &block, msd_counter_size, &msd_counter,
                     block_counts, types, displs, MPI_INT64_T);
    }

    if (needToDo[acc_ratio]){
        append_block(&offset, &block, 1, &accepts,
                block_counts, types, displs, MPI_INT64_T);
        append_block(&offset, &block, 1, &moves,
                block_counts, types, displs, MPI_INT64_T);
    }

    if (needToDo[non_bonded_energy]){
        append_block(&offset, &block, non_bonded_energy_val_size, &non_bonded_energy_val,
                block_counts, types, displs, MPI_DOUBLE);
    }

    if (needToDo[bonded_energy]){
        append_block(&offset, &block, bonded_energy_val_size, &bonded_energy_val,
                     block_counts, types, displs, MPI_DOUBLE);
    }

    if (needToDo[dynamical_structure]){
        append_block(&offset, &block, dynamical_sf_val_size, &dynamical_sf_val,
                     block_counts, types, displs, MPI_DOUBLE);
    }

    if (needToDo[static_structure]){
        append_block(&offset, &block, static_sf_val_size, &static_sf_val,
                     block_counts, types, displs, MPI_DOUBLE);
    }
	
    printf("Type-create-arg: %d, %d, %ld, %d\n", block, block_counts[0], displs[0], types[0]==MPI_DOUBLE);
    MPI_Type_create_struct(block, block_counts, displs, types, &buff_mpi_dtype);
    MPI_Type_commit(&buff_mpi_dtype);

    if (!op_is_valid){
        // operator is assumed to be associative and commutative
        MPI_Op_create(ana_reduce_func, 1, &ana_op);
        op_is_valid = true;
    }

    sendbuffer = malloc(offset * sizeof(double));
    if (sendbuffer == NULL){
	    
    	fprintf(stderr, "ERROR: malloc failed. "
                        "(line %d of file %s)" , __LINE__, __FILE__);
	return -1;
	 
    }
    recvbuffer = malloc(offset * sizeof(double));
    if (recvbuffer == NULL){
	    
    	fprintf(stderr, "ERROR: malloc failed. "
                        "(line %d of file %s)" , __LINE__, __FILE__);
	free(sendbuffer);
	return -1;
	 
    }

    return 0;

}

int ana_reduce(int root, MPI_Comm comm){
    int rank;
    MPI_Comm_rank(comm, &rank);
    int ret;
    if (rank == root){
        ret = MPI_Reduce(MPI_IN_PLACE, sendbuffer, 1, buff_mpi_dtype, ana_op, root, comm);
    }
    else{
        ret = MPI_Reduce(sendbuffer, recvbuffer, 1, buff_mpi_dtype, ana_op, root, comm);
    }
    return ret;
}

void ana_free(){
    //free buffers, the datatypeset the global variables for the offset back to -1
    free(sendbuffer);
    free(recvbuffer);
    MPI_Type_free(&buff_mpi_dtype);
    sendbuffer = NULL;
    recvbuffer = NULL;
    re_vals = -1; re_counter = -1;
    dvar_val = -1;
    rg_vals = -1; rg_counter = -1;
    anisotropy_vals = -1; anisotropy_counter = -1;
    msd_vals = -1; msd_counter = -1;
    accepts = -1; moves = -1;
    non_bonded_energy_val = -1; bonded_energy_val = -1;
    dynamical_sf_val = -1; static_sf_val = -1;
}

// creates a block of size blocksize and mpi_type type on the global buffers
// and write the offset corresponding to the starting address of the block to out_offset
// increments offset and block to recognize the additional block
// write the information into block_counts, types and displs so that an MPI-struct may
// later be created using these three arrays
void append_block(size_t *offset, int *block, size_t block_size, ssize_t *out_offset,
        int *block_counts, MPI_Datatype *types, MPI_Aint *displs, MPI_Datatype type){

    *out_offset = *offset;
    displs[*block] = *offset* sizeof(double); // this currently relies on all datatypes in the buffer being the same length, later we will need an additional variable to track this.
    block_counts[*block] = block_size;
    types[*block] = type;

    *offset += block_size;
    (*block)++;

}

void * get_pointer_from_offset(int offset, char * name){
    if (sendbuffer == NULL){
        fprintf(stderr, "ERROR: invalid memory request "
                        "(line %d of file %s) for the space %s : "
                        "buffer is not allocated.\n", __LINE__, __FILE__, name);
        return NULL;
    }
    if (offset < 0){
        fprintf(stderr, "ERROR: invalid memory request "
                        "(line %d of file %s) for the space %s : "
                        "buffer is allocated, but observable "
                        "is not included in it.\n", __LINE__, __FILE__, name);
        return NULL;
    }

    return &sendbuffer[offset];
}

double * get_re_vals(){
    return get_pointer_from_offset(re_vals, "Re (values)");
}

uint64_t * get_re_counter(){
    return get_pointer_from_offset(re_counter, "Re (counter)");
}

double * get_dvar_val(){
    return get_pointer_from_offset(dvar_val, "dvar");
}

double * get_rg_vals(){
    return get_pointer_from_offset(rg_vals, "Rg (values)");
}

uint64_t * get_rg_counter(){
    return get_pointer_from_offset(rg_counter, "Rg (counter)");
}

double * get_anisotropy_vals(){
    return get_pointer_from_offset(anisotropy_vals, "anisotropy (values)");
}

uint64_t * get_anisotropy_counter(){
    return get_pointer_from_offset(anisotropy_counter, "anisotropy (counter)");
}

double * get_msd_vals(){
    return get_pointer_from_offset(msd_vals, "msd(values)");
}

uint64_t * get_msd_counter(){
    return get_pointer_from_offset(msd_counter, "msd (counter)");
}

uint64_t * get_accepts(){
    return get_pointer_from_offset(accepts, "accepts (for accratio)");
}

uint64_t * get_moves(){
    return get_pointer_from_offset(moves , "moves (for accratio)");
}

double * get_non_bonded_energy_val(){
    return get_pointer_from_offset(non_bonded_energy_val, "non-bonded energy");
}

double * get_bonded_energy_val(){
    return get_pointer_from_offset(bonded_energy_val, "bonded energy");
}

double * get_dynamical_sf_val(){
    return get_pointer_from_offset(dynamical_sf_val, "dynamical structure factor");
}

double * get_static_sf_val(){
    return get_pointer_from_offset(static_sf_val, "static structure factor");
}
