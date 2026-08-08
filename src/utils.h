#include <stdio.h>

#ifndef utils_h
#define utils_h

void* xcalloc(int, size_t);
void error(const char*);
double** new_matrix(int, int);
void free_matrix(double**, int);
void print_matrix(double**, int, int);
FILE* open_file(const char*);
double* load_array(const char*, int*);
double** load_matrix(const char*, int*, int*);

#endif
