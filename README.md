## Implementation of the Baum-Welch algorithm in C

### Algorithm
Estimates HMM parameters by maximizing the likelihood given a set of observations.

init(i) - initial probability of state i

transition(i,j) - transition probability between state i and state j

emission(i,j) - probability of observation j from a state i

The log likelihood is computed as the sum of the logarithms of the alpha scaling factors.

### Computational Considerations
To avoid numerical underflow the alphas and betas are scaled.
