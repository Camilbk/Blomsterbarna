import numpy as np
import scipy


def LU(A):
    n,m = A.shape
    assert (n == m)  # square matices only

    P = np.zeros(n)
    for row in A:
        for col in A:
            for k in col:
                A[row,col] = A[row,col] - A[row,k] * A[k,col]