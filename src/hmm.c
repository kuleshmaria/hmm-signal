#include <stdio.h>
#include <stdlib.h>
#include "hmm.h"
#include "utils.h"

hmm_t* new_hmm(const char* init_path, const char* emission_path, const char* transition_path)
{
    hmm_t* hmm;
    int m, n;
    int a, b;

    hmm = xcalloc(1, sizeof(hmm_t));

    hmm->emission_probs = load_matrix(emission_path, &m, &n);

    hmm->init_probs = load_array(init_path, &a);
    if(a != m)
      error("init matrix has the wrong dimensions");

    hmm->transition_matrix = load_matrix(transition_path, &a, &b);
    if(a != m || b != m)
      error("transition matrix has wrong dimenions");

    hmm->n_states = m;
    hmm->n_symbols = n;
    return hmm;
}

void free_hmm(hmm_t* hmm)
{
    free(hmm->init_probs);
    free_matrix(hmm->emission_probs, hmm->n_states);
    free_matrix(hmm->transition_matrix, hmm->n_states);
    free(hmm);
}

void forward(double** alphas, double* scaling, hmm_t* hmm, int* seq, int n_seq)
{
    int i, j, k;

    scaling[0] = 0;
    for(i = 0; i < hmm->n_states; i++)
    {
        alphas[i][0] = hmm->init_probs[i] * hmm->emission_probs[i][seq[0]];
        scaling[0] += alphas[i][0];
    }
    for(i = 0; i < hmm->n_states; i++)
    {
        alphas[i][0] /= scaling[0];
    }
    for(j = 1; j < n_seq; j++)
    {
        scaling[j] = 0;
        for(i = 0; i < hmm->n_states; i++)
        {
            alphas[i][j] = 0;
            for(k = 0; k < hmm->n_states; k++)
            {
                alphas[i][j] += alphas[k][j-1]*hmm->transition_matrix[k][i];
            }
            alphas[i][j] *= hmm->emission_probs[i][seq[j]];
            scaling[j] += alphas[i][j];
        }
        for(i = 0; i < hmm->n_states; i++)
        {
            alphas[i][j] /= scaling[j];
        }
    }
}

void backward(double** betas, double* scaling, hmm_t* hmm, int* seq, int n_seq)
{
    int i, j, k;

    for(i = 0; i < hmm->n_states; i++)
    {
        betas[i][n_seq-1] = 1;
    }
    for(j = n_seq-2; j >= 0; j--)
    {
        for(i = 0; i < hmm->n_states; i++)
        {
            betas[i][j] = 0;
            for(k = 0; k < hmm->n_states; k++)
            {
                betas[i][j] += hmm->transition_matrix[i][k] * hmm->emission_probs[k][seq[j+1]] * betas[k][j+1];
            }
            betas[i][j] /= scaling[j+1];
        }
    }
}
