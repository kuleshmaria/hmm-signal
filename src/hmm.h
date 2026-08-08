#ifndef hmm_h
#define hmm_h

typedef struct hmm_t hmm_t;
struct hmm_t
{
    int n_states;
    int n_symbols;
    double* init_probs;
    double** emission_probs;
    double** transition_matrix;
};

hmm_t* new_hmm(const char*, const char*, const char*);
void free_hmm(hmm_t*);
void forward(double**, double*, hmm_t*, int*, int);
void backward(double**, double*, hmm_t*, int*, int);

#endif
