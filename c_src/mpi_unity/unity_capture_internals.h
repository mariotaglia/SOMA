#ifndef _UNITY_CAPTURE_INTERNALS_H_
#define _UNITY_CAPTURE_INTERNALS_H_
#include <hdf5.h>

#include "unity_capture.h"

char * append_number_to_filename(const char * oldname, int num);

int save_values (void *val, size_t size, const char *name, capture_file_t capture_file, enum groups group, hid_t h5type);
int read_values (capture_file_t cf, enum groups group, const char *name, void *value, size_t size, hid_t h5type);



#endif //_UNITY_CAPTURE_INTERNALS_H_
