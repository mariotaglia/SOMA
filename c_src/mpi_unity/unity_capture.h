#ifndef _UNITY_CAPTURE_H_
#define _UNITY_CAPTURE_H_
#include <hdf5.h>




enum action {FILE_READ, FILE_WRITE, IGNORE_CALLS};
enum open_mode {READ_ONLY, FAIL_IF_EXISTS, EVADE_IF_EXISTS, OVERRIDE_IF_EXISTS, NO_OP_IF_EXISTS};
enum groups {INPUTS, RESULTS, SPECIAL};

struct capture_file
{
    hid_t file_id;
    enum action act;
    MPI_Comm comm;
    int comm_rank;
    int comm_size;
};

typedef struct capture_file capture_file_t;

#define SAVE_INT_ARR(array, size, rwa_file, group) save_int(array, size, #array, rwa_file, group)
#define SAVE_INT(variable, capture_file, group) save_int(&variable, 1, #variable, capture_file, group)
int save_int (const int *val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

#define SAVE_DOUBLE_ARR(array, size, rwa_file, group) save_double(array, size, #array, rwa_file, group)
#define SAVE_DOUBLE(variable, capture_file, group) save_double(&variable, 1, #variable, capture_file, group)
int save_double (const double *val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

#define SAVE_FLOAT_ARR(array, size, rwa_file, group) save_float(array, size, #array, rwa_file, group)
#define SAVE_FLOAT(variable, capture_file, group) save_float(&variable, 1, #variable, capture_file, group)
int save_float (const float *val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

#define SAVE_UINT64_ARR(array, size, rwa_file, group) save_uint64(array, size, #array, rwa_file, group)
#define SAVE_UINT64(variable, capture_file, group) save_uint64(&variable, 1, #variable, capture_file, group)
int save_uint64 (const uint64_t *val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

#define SAVE_UINT16_ARR(array, size, rwa_file, group) save_uint16(array, size, #array, rwa_file, group)
#define SAVE_UINT16(variable, capture_file, group) save_uint16(&variable, 1, #variable, capture_file, group)
int save_uint16 (const uint16_t *val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

#define SAVE_UINT_ARR(array, size, rwa_file, group) save_unsigned_int(array, size, #array, rwa_file, group)
#define SAVE_UINT(variable, capture_file, group) save_unsigned_int(&variable, 1, #variable, capture_file, group)
int save_unsigned_int (const unsigned int * val, size_t size, const char *name, capture_file_t capture_file, enum groups group);

int read_uint16 (capture_file_t cf, enum groups group, const char *name, uint16_t *value, size_t size);
int read_uint64 (capture_file_t cf, enum groups group, const char *name, uint64_t *value, size_t size);
int read_int (capture_file_t cf, enum groups group, const char *name, int *value, size_t size);
int read_unsigned_int (capture_file_t cf, enum groups group, const char *name, unsigned int *value, size_t size);
int read_float (capture_file_t cf, enum groups group, const char *name, float *value, size_t size);
int read_double (capture_file_t cf, enum groups group, const char *name, double *value, size_t size);

capture_file_t capture_file_open(const char * name, MPI_Comm comm, enum open_mode mode);
void capture_file_close(capture_file_t * cf);


// capture-file convenience api that has a notion of a "current" file (global state).
// Has no feature that the above functions dont have, but is often shorter to use
void file_open(const char * name_suggestion, MPI_Comm comm);
void file_close();
void change_group(enum groups group);




#define SAVE(X) _Generic((X),\
double : save_double,\
float : save_float,\
int : save_int,\
uint16_t: save_uint16,\
uint64_t: save_uint64,\
unsigned int : save_unsigned_int\
)(&X,1,#X, get_current_file(), get_current_group())

#define SAVE_ARR(arr, size) _Generic((arr[0]),\
double : save_double,\
float : save_float,\
int : save_int,\
uint16_t : save_uint16,\
uint64_t : save_uint64,\
unsigned int : save_unsigned_int\
)(arr, size, #arr, get_current_file(), get_current_group())

capture_file_t get_current_file();
enum groups get_current_group();
void not_implemented(int line, const char * func, const char * file, const char * var_name, const char * macro_name);

#endif //_UNITY_CAPTURE_H_
