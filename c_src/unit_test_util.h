#ifndef _UNIT_TEST_UTIL_H_
#define _UNIT_TEST_UTIL_H_
#include <zconf.h>

#define FOR_EACH_PARAMETER(type, function, ...) \
    type parameter_arr[] = __VA_ARGS__;\
    size_t parameter_arr_size = sizeof(parameter_arr)/sizeof(type);\
    for (size_t iter=0; iter < parameter_arr_size; iter++)\
        {\
            type parameter = parameter_arr[0];\
            function(parameter);\
        }


#endif //_UNIT_TEST_UTIL_H_
