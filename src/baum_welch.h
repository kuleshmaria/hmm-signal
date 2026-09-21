#ifndef baum_welch_h
#define baum_welch_h

#include "hmm.h"

void baum_welch(hmm_t*, int*, int, int, double);
double log_likelihood(double*, int);

#endif
