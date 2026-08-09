#include <stdio.h>
#include <stdlib.h>
#include "baum_welch.h"
#include "hmm.h"
#include "utils.h"

#define DEBUG 1

static void init_baum_welch(hmm_t* hmm, double*** alphas, double*** betas, double** scaling, double*** gammas, double**** xi, int n_seq)
{
    int i;
    *alphas = new_matrix(hmm->n_states, n_seq);
    *betas = new_matrix(hmm->n_states, n_seq);
    *scaling = xcalloc(n_seq, sizeof(double));
    *gammas = new_matrix(hmm->n_states, n_seq);
    *xi = xcalloc(hmm->n_states, sizeof(double**));
    for(i = 0; i < hmm->n_states; i++)
    {
        (*xi)[i] = new_matrix(hmm->n_states, n_seq-1);
    }
}

static void free_baum_welch(hmm_t* hmm, double** alphas, double** betas, double* scaling, double** gammas, double*** xi)
{
    int i;
    free_matrix(alphas, hmm->n_states);
    free_matrix(betas, hmm->n_states);
    free(scaling);
    free_matrix(gammas, hmm->n_states);
    for(i = 0; i < hmm->n_states; i++)
    {
        free_matrix(xi[i], hmm->n_states);
    }
    free(xi);
}

static void compute_gammas(double** gammas, hmm_t* hmm, double** alphas, double** betas, int n_seq)
{
    int t, i;
    double norm_sum;
    for(t = 0; t < n_seq; t++)
    {
        norm_sum = 0;
        for(i = 0; i < hmm->n_states; i++)
        {
            gammas[i][t] = alphas[i][t] * betas[i][t];
            norm_sum += gammas[i][t];
        }
        for(i = 0; i < hmm->n_states; i++)
        {
            gammas[i][t] /= norm_sum;
        }
    }
}

static void compute_xi(double*** xi, hmm_t* hmm, double** alphas, double** betas, int* seq, int n_seq)
{
    int t, i, j;
    double norm_sum;

    for(t = 0; t < n_seq - 1; t++)
    {
        norm_sum = 0;
        for(i = 0; i < hmm->n_states; i++)
        {
            for(j = 0; j < hmm->n_states; j++)
            {
                xi[i][j][t] = alphas[i][t]*hmm->transition_matrix[i][j]*betas[j][t+1]*hmm->emission_probs[j][seq[t+1]];
                norm_sum += xi[i][j][t];
            }
        }
        for(i = 0; i < hmm->n_states; i++)
        {
            for(j = 0; j < hmm->n_states; j++)
            {
                xi[i][j][t] /= norm_sum;
            }
        }
    }
}

static void update_hmm_parameters(hmm_t* hmm, double** gammas, double*** xi, int* seq, int n_seq)
{
    int i, j, t;
    double norm_sum;

    for(i = 0; i < hmm->n_states; i++)
    {
        hmm->init_probs[i] = gammas[i][0];
    }
    for(i = 0; i < hmm->n_states; i++)
    {
        norm_sum = 0;
        for(t = 0; t < n_seq-1; t++)
        {
            norm_sum += gammas[i][t];
        }
        for(j = 0; j < hmm->n_states; j++)
        {
            hmm->transition_matrix[i][j] = 0;
            for(t = 0; t < n_seq-1; t++)
            {
                hmm->transition_matrix[i][j] += xi[i][j][t];
            }
            hmm->transition_matrix[i][j] /= norm_sum;
        }
    }
    for(i = 0; i < hmm->n_states; i++)
    {
        norm_sum = 0;
        for(t = 0; t < n_seq; t++)
        {
            norm_sum += gammas[i][t];
        }
        for(j = 0; j < hmm->n_symbols; j++)
        {
            hmm->emission_probs[i][j] = 0;
            for(t = 0; t < n_seq; t++)
            {
                if(seq[t] == j)
                {
                    hmm->emission_probs[i][j] += gammas[i][t];
                }
            }
            hmm->emission_probs[i][j] /= norm_sum;
        }
    }
}

void baum_welch(hmm_t* hmm, int* seq, int n_seq, int n_iterations)
{
    double** alphas;
    double** betas;
    double* scaling;
    double** gammas;
    double*** xi;
    int k;

    init_baum_welch(hmm, &alphas, &betas, &scaling, &gammas, &xi, n_seq);

    for(k = 0; k < n_iterations; k++)
    {
        // forward-backward
        forward(alphas, scaling, hmm, seq, n_seq);
        backward(betas, scaling, hmm, seq, n_seq);

        // gamma - probability of being in state i at time t
        compute_gammas(gammas, hmm, alphas, betas, n_seq);

        // xi - probability of being in state i and j at time t and t+1
        compute_xi(xi, hmm, alphas, betas, seq, n_seq);

        // update
        update_hmm_parameters(hmm, gammas, xi, seq, n_seq);
    }
    // free
    free_baum_welch(hmm, alphas, betas, scaling, gammas, xi);
}
