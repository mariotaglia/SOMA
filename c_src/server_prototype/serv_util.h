#ifndef SERV_UTIL_H
#define SERV_UTIL_H


#define ABORT_ON_NULL(ptr, message) if(ptr==NULL){nullAbort(__FILE__, __LINE__, #ptr, message);}else{}


void nullAbort(const char* file, int line_no, const char* ptr_name, const char *message);


#define DPRINT(...) {int rank;\
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);\
	fprintf(stdout, "******DEBUG: rank %d line %d of file %s: ",rank, __LINE__, __FILE__);\
	fprintf(stdout, __VA_ARGS__);\
	fprintf(stdout,"\n");\
	fflush(stdout);}

#endif // SERV_UTIL_H
