#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#include "hmm.h"
#include "baum_welch.h"

void init(hmm_t** hmm)
{
    *hmm = new_hmm(
      "../matrices/init.txt",
      "../matrices/emission.txt",
      "../matrices/transition.txt"
    );
}

int main()
{
    hmm_t* hmm;

    int seq[] = {4, 5, 5, 5, 0, 0, 1, 2, 3};
    int n_seq = 9;

    init(&hmm);

    baum_welch(hmm, seq, n_seq, 1);

    // free
    free_hmm(hmm);
    return 0;
}
