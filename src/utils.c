#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"

#define DEBUG    1
#define BUFFSIZE 1000

void* xcalloc(int n, size_t size)
{
    void* a = calloc(n, size);
    if(a == NULL)
        error("out of memory");
    return a;
}

void error(const char* msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}

double** new_matrix(int m, int n)
{
    int i;
    double** a;
    a = xcalloc(m, sizeof(double*));

    for(i = 0; i < m; i++)
    {
        a[i] = xcalloc(n, sizeof(double));
    }
    return a;
}

void free_matrix(double** a, int n)
{
    int i;
    for(i = 0; i < n; i++)
    {
        free(a[i]);
    }
    free(a);
}

void print_matrix(double** a, int n, int m)
{
    int i, j;
    for(i = 0; i < n; i++)
    {
      for(j = 0; j < m; j++)
      {
          printf("%f ", a[i][j]);
      }
      printf("\n");
    }
}

FILE* open_file(const char* path)
{
    FILE *input;

    if((input = fopen(path, "r")) == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", path);
        error("file not found");
    }
    return input;
}

double* load_array(const char* path, int* n)
{
    FILE *input;
    int i;
    double* a;

    input = open_file(path);

    if(fscanf(input, "%d", n) != 1)
        error("invalid file");

    a = xcalloc(*n, sizeof(double));
    for(i = 0; i < *n; i++)
    {
        fscanf(input, "%le", &a[i]);
    }
    fclose(input);
    return a;
}

double** load_matrix(const char* path, int* n, int *m)
{
    FILE *input;
    int i, j;
    double** a;
    double c;

    input = open_file(path);

    if(fscanf(input, "%d %d\n", n, m) != 2)
        error("invalid file");

    a = new_matrix(*n, *m);
    for(i = 0; i < *n; i++)
    {
        for(j = 0; j < *m; j++)
        {
            fscanf(input, "%le", &c);
            a[i][j] = c;
        }
    }
    fclose(input);
    return a;
}
