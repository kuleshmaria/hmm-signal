#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#include "hmm.h"

void init(hmm_t** hmm, int n_seq, double*** alphas, double*** betas, double** scaling)
{
    *hmm = new_hmm(
      "../matrices/init.txt",
      "../matrices/emission.txt",
      "../matrices/transition.txt"
    );
    *alphas = new_matrix((*hmm)->n_states, n_seq);
    *betas = new_matrix((*hmm)->n_states, n_seq);
    *scaling = xcalloc(n_seq, sizeof(double));
}

void free_all(hmm_t* hmm, double** alphas, double** betas, double* scaling)
{
    free_matrix(alphas, hmm->n_states);
    free_matrix(betas, hmm->n_states);
    free(scaling);
    free_hmm(hmm);
}

int main()
{
    hmm_t* hmm;
    double** alphas;
    double** betas;
    double* scaling;

    int seq[] = {4, 5, 5, 5, 0, 0, 1, 2, 3};
    int n_seq = 9;

    init(&hmm, n_seq, &alphas, &betas, &scaling);

    forward(alphas, scaling, hmm, seq, n_seq);
    print_matrix(alphas, hmm->n_states, n_seq);

    backward(betas, scaling, hmm, seq, n_seq);
    print_matrix(betas, hmm->n_states, n_seq);

    // free all
    free_all(hmm, alphas, betas, scaling);
    return 0;
}
