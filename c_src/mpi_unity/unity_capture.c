#include "unity_capture.h"
#include "unity_capture_internals.h"

#include "unity_fixture.h"

#include <mpi.h>
#include <hdf5.h>

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <unistd.h>

#include "err_handling.h"

static int create_groups_if_not_existing(hid_t file_id)
{
    char *groups[] = {"inputs", "results", "special"};
    for (unsigned i=0; i<3; i++)
        {
            htri_t exists = H5Lexists(file_id, groups[i], H5P_DEFAULT);

            if (exists < 0) return -1; // failure

            if (exists == 0) // does not exist
                {
                    hid_t group = H5Gcreate(file_id, groups[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    H5Gclose(group);
                }
        }
    return 0;
}

char * append_number_to_filename(const char * oldname, int num)
{
    char filename[10000];
    char suffix[100];
    sscanf(oldname, "%9999[^.].%99s", filename, suffix);

    size_t maxlen = 10000;
    char * newname = malloc(sizeof(char) * maxlen);
    if (newname == NULL)
    {
        fprintf(stderr, "Malloc error (file %s func %s line %d)\n",
            __FILE__, __func__, __LINE__);
        return NULL;
    }

    size_t written = snprintf(newname,maxlen-1, "%s_%d.%s", filename, num, suffix);
    if (written > maxlen)
    {
        fprintf(stderr, "ERROR writing filename (too long?) (file %s func %s line %d)\n", __FILE__, __func__, __LINE__);
        return NULL;
    }

    return newname;
}


static herr_t add_comm_size_attribute(hid_t file_id, MPI_Comm comm)
{
    hid_t attr = H5Acreate(file_id, "comm_size",H5T_NATIVE_INT, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    int size;
    MPI_Comm_size(comm, &size);
    herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &size);
    H5Aclose(attr);
    return status;
}

static int add_metadata_and_structure(hid_t file_id, MPI_Comm comm)
{
    herr_t status = add_comm_size_attribute(file_id, comm);
    if (status)
    {
        fprintf(stderr, "Error code %d, cannot add attribute (file %s, func %s, line %d)\n",
            status, __FILE__, __func__, __LINE__);
        return status;
    }
    int err = create_groups_if_not_existing(file_id);
    if (err)
    {
        fprintf(stderr, "ERROR: creating groups failed (error %d, file %s func %s line %d)\n",
                err, __FILE__, __func__, __LINE__);
        abort();
    }

    return 0;
}

static int get_file_comm_size(hid_t file_id)
{
    int size;
    hid_t attr = H5Aopen(file_id, "comm_size", H5P_DEFAULT);
    int err = H5Aread(attr, H5T_NATIVE_INT, &size);
    if (err)
    {
        fprintf(stderr, "ERROR: reading comm_size attribute from hdf5-file failed (error %d, file %s func %s line %d)\n",
                err, __FILE__, __func__, __LINE__);
        return -1;
    }
    H5Aclose(attr);
    return size;
}

static int get_comm_size(MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

static bool file_exists(const char * name)
{
    return (access(name, F_OK) != -1);
}

static int evade_name(const char * oldname, char ** newname)
{
    *newname = malloc(10000);
    if (*newname == NULL)
    {
        fprintf(stderr, "Malloc error (line %d file %s func %s)\n", __LINE__, __FILE__, __func__);
        return -1;
    }
    if (!file_exists(oldname))
    {
        strcpy(*newname, oldname);
        return 0;
    }

    int i = 1;
    do{

        free(*newname);
        *newname = append_number_to_filename (oldname, i);
        i++;
    } while (file_exists(*newname));

    return 0;
}

capture_file_t capture_file_open(const char * name, MPI_Comm comm_in, enum open_mode mode)
{
    capture_file_t cf;
    // values of capture file to be initialized
    MPI_Comm_dup(comm_in, &cf.comm);
    MPI_Comm_rank(cf.comm, &cf.comm_rank);
    MPI_Comm_size(cf.comm, &cf.comm_size);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, cf.comm, MPI_INFO_NULL);


    if (mode == EVADE_IF_EXISTS)
    {
        char *evasion_filename;
        evade_name(name, &evasion_filename);
        cf.file_id = H5Fcreate(evasion_filename, H5F_ACC_EXCL, H5P_DEFAULT, fapl);
        assert(cf.file_id > 0); // creating file should never fail, as we selected a unique name.
        add_metadata_and_structure(cf.file_id, comm_in);
        cf.act = FILE_WRITE;
        free(evasion_filename);
    }
    else if (mode == FAIL_IF_EXISTS)
    {
        cf.file_id = H5Fcreate(name, H5F_ACC_EXCL, H5P_DEFAULT, fapl);
        if (cf.file_id < 0) return cf;
        add_metadata_and_structure(cf.file_id, comm_in);
        cf.act = FILE_WRITE;
    }
    else if (mode == OVERRIDE_IF_EXISTS)
    {
        cf.file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        if (cf.file_id < 0) return cf;
        add_metadata_and_structure(cf.file_id, comm_in);
        cf.act = FILE_WRITE;
    }
    else if (mode == NO_OP_IF_EXISTS)
    {
        bool exists = file_exists(name);
        if (!exists)
        {
            // does not exist
            cf.file_id = H5Fcreate(name, H5F_ACC_EXCL, H5P_DEFAULT, fapl);
            add_metadata_and_structure(cf.file_id, comm_in);
            cf.act = FILE_WRITE;
        }
        else
        {
            cf.file_id = -1;
            cf.act = IGNORE_CALLS;
        }
    }
    else if (mode == READ_ONLY)
    {
        cf.file_id = H5Fopen(name, H5F_ACC_RDONLY, fapl);
        cf.act = FILE_READ;

        if (cf.file_id > 0)
        {
            int this_comm_size = get_comm_size(comm_in);
            int file_comm_size = get_file_comm_size(cf.file_id);

            if (this_comm_size != file_comm_size)
            {
                fprintf(stdout, "ERROR: attempting to open capture file with wrong number of ranks, file expects %d ranks, but was opened with %d ranks (bad arguments to function %s in file %s)\n", file_comm_size, this_comm_size, __func__, __FILE__);

                H5Fclose(cf.file_id);
                cf.file_id = -1;
            }
        }

    }
    else
    {
        fprintf(stderr, "ERROR: unknown mode %d, can not open capture_file (file %s func %s line %d)\n",
            mode, __FILE__, __func__, __LINE__);
        abort();
    }

    return cf;

}

void capture_file_close(capture_file_t * cf)
{
    if (cf->file_id > 0)
    {
        H5Fclose(cf->file_id);
        cf->file_id = -1;
    }
}

static char * group_name_from_enum(enum groups group)
{
    switch (group)
        {
            case INPUTS:
                return  "inputs";
            case RESULTS:
                return "results";
            case SPECIAL:
                return "special";
            default:
                fprintf(stderr, "ERROR: unkown group %d (file %s, function %s, line %d)",
                        group, __FILE__, __func__, __LINE__);
                return NULL;
        }
}

static int get_xfer_params(const capture_file_t * file, hid_t dset, hid_t *memspace, hid_t *filespace, hid_t *dxpl, size_t size)
{
    herr_t status;

    hsize_t offset[] = {file->comm_rank, 0};
    hsize_t count[] = {1, 1};
    hsize_t blocks[] = {1, size};

    *filespace = H5Dget_space(dset);
    status = H5Sselect_hyperslab(*filespace, H5S_SELECT_SET, offset, NULL, count, blocks);
    if (status != 0)
    {
        fprintf(stderr, "ERROR code %d: failed to select hyperslab,  (func %s, line %d, file %s)",
            status, __func__, __LINE__, __FILE__);
        return status;
    }
    *dxpl= H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(*dxpl, H5FD_MPIO_COLLECTIVE);

    hsize_t arr[] = {1, size};
    *memspace = H5Screate_simple(2, arr, arr);

    return 0;
}

int save_double (const double *val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
    return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_DOUBLE);
}

int save_float (const float *val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
    return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_FLOAT);
}

int save_uint16 (const uint16_t *val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
    return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_UINT16);
}

int save_uint64 (const uint64_t *val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
    return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_UINT64);
}

int save_int (const int *val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
    return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_INT);
}

int save_unsigned_int (const unsigned int * val, size_t size, const char *name, capture_file_t capture_file, enum groups group)
{
   return save_values ((void*)val, size, name, capture_file, group, H5T_NATIVE_UINT);
}

int save_values (void *val, size_t size, const char *name, capture_file_t capture_file, enum groups group, hid_t h5type)
{

    if (capture_file.act == FILE_READ)
    {
        fprintf(stderr, "ERROR: cannot write to file that is only open for reading.\n");
        return -1;
    }
    else if (capture_file.act == IGNORE_CALLS)
    {
        return 0;
    }
    else if (capture_file.act != FILE_WRITE)
    {
        fprintf(stderr, "ERROR: unknown capture_file action %d", capture_file.act);
        return -1;
    }
    // action of file is file_write, we are good to go

    char * group_name = group_name_from_enum(group);
    assert(group_name);

    hid_t group_id = H5Gopen(capture_file.file_id, group_name, H5P_DEFAULT);

    // want to use compact for access speed, but it does not seem to allow parallel writes
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(dcpl, H5D_CONTIGUOUS);

    hsize_t dims[] = {capture_file.comm_size, size};
    hid_t dspace = H5Screate_simple(2, dims, dims);

    hid_t dset = H5Dcreate(group_id, name, h5type, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Sclose(dspace);

    hid_t memspace, filespace, dxpl;
    int err = get_xfer_params(&capture_file, dset, &memspace, &filespace, &dxpl, size);
    if (err != 0)
    {
        fprintf(stderr, "failed to set xfer parameters\n");
        return err;
    }

    herr_t status = H5Dwrite(dset, h5type, memspace, filespace,
		      dxpl, val);
    if (status < 0)
    {
        fprintf(stderr, "ERROR code %d, H5Dwrite failed (file %s func %s line %d)\n",
            status, __FILE__, __func__, __LINE__);
        return status;
    }

    H5Pclose(dxpl);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
    H5Gclose(group_id);

    return status;
}

int read_int (capture_file_t cf, enum groups group, const char *name, int *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_INT);
}

int read_double (capture_file_t cf, enum groups group, const char *name, double *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_DOUBLE);
}
int read_uint16 (capture_file_t cf, enum groups group, const char *name, uint16_t *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_UINT16);
}
int read_uint64 (capture_file_t cf, enum groups group, const char *name, uint64_t *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_UINT64);
}
int read_float (capture_file_t cf, enum groups group, const char *name, float *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_FLOAT);
}

int read_unsigned_int (capture_file_t cf, enum groups group, const char *name, unsigned int *value, size_t size)
{
    return read_values (cf, group, name, value, size, H5T_NATIVE_UINT);
}
int read_values (capture_file_t cf, enum groups group, const char *name, void *value, size_t size, hid_t h5type)
{
    hid_t group_id = H5Gopen(cf.file_id, group_name_from_enum(group), H5P_DEFAULT);
    htri_t dset_exists = H5Lexists(group_id, name, H5P_DEFAULT);
    if (dset_exists < 0)
    {
        fprintf(stderr, "ERROR code %d checking dset existence (file %s, func %s, line %d)\n",
            dset_exists, __FILE__, __func__, __LINE__);
        return -1;
    }
    else if (dset_exists == 0)
    {
        fprintf(stderr, "ERROR: tried to read non-existing dataset \"%s\" (file %s func %s line %d)\n",
            name, __FILE__, __func__, __LINE__);
        return -1;
    }

    hid_t dset = H5Dopen(group_id, name, H5P_DEFAULT);

    hid_t memspace, filespace, dxpl;
    get_xfer_params(&cf, dset, &memspace, &filespace, &dxpl, size);

    herr_t status = H5Dread (dset, h5type, memspace, filespace,
                      dxpl, value);
    if (status < 0)
    {
        fprintf(stderr, "ERROR code %d H5Dread failed (file %s, func %s, line %d)",
            status, __FILE__, __func__, __LINE__);
        return -1;
    }

    H5Pclose(dxpl);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
    H5Gclose(group_id);

    return 0;
}

// convenience api
static capture_file_t current_file;
static enum groups current_group;

capture_file_t get_current_file()
{
    return current_file;
}

enum groups get_current_group()
{
    return current_group;
}

void change_group(enum groups group)
{
    switch (group)
        {
            case INPUTS:
            case RESULTS:
            case SPECIAL:
                current_group = group;
                break;
            default:
                fprintf(stderr,"ERROR: unknown group %d\n",group);
        }
}

void file_open (const char *name_suggestion, MPI_Comm comm)
{
    current_file = capture_file_open(name_suggestion, comm, EVADE_IF_EXISTS);
    assert(current_file.file_id > 0);
    current_group = INPUTS;
}
void file_close()
{
    capture_file_close(&current_file);
}

void not_implemented(int line, const char * func, const char * file, const char * var_name, const char * macro_name)
{
    fprintf(stderr, "argument \"%s\" to macro-function \"%s\" has a type that is not yet implemented. (macro callsite is file %s line %d function %s). aborting.\n", var_name, macro_name, file, line, func);
    abort();
}



