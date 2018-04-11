import numpy as np
import math


def LU(a):
    """
LU factorization without pivoting
Input:  matrix a
Output: L in the lower triangle of a
        U in the upper triangle of a
"""
    n, m = a.shape
    assert (n == m)  # work with square matrices only
    for j in range(n - 1):  # go column after column
        if math.fabs(a[j, j]) < 1.0E-13:
            print("Error: zero pivot in position (%d,%d)" % (j, j))
            break
        for i in range(j + 1, n):  # and then row after row
            a[i, j] = a[i, j] / a[j, j]  # compute l_ij
            for k in range(j + 1, n):
                # update a_ik - same as Gaussian elimination
                # a_ik = a_ik - l_ij * u_jk
                a[i, k] = a[i, k] - a[i, j] * a[j, k]
    return


if __name__ == "__main__":
    # Test LU factorization
    # Example 2.5 from the book
    A = np.ay([[1, 2, -1], [2, 1, -2], [-3, 1, 1]], dtype='f8')
    L = np.array([[1, 0, 0], [2, 1, 0], [-3, -7. / 3, 1]], dtype='f8')
    U = np.array([[1, 2, -1], [0, -3, 0], [0, 0, -2]], dtype='f8')
    # save A as LU overwrites it
    A1 = np.copy(A)
    LU(A1)
    # extract lower and upper triangles
    L1 = np.tril(A1, -1) + np.eye(3)
    U1 = np.triu(A1)
    print(A)
    print(L1)
    print(U1)
    # compute errors
    print(np.linalg.norm(L - L1, 'fro'), np.linalg.norm(U - U1, 'fro'), np.linalg.norm(A - np.matmul(L1, U1), 'fro'))

    # test on some random matrices
    n = 20
    L = np.tril(np.random.rand(n, n), -1) + np.eye(n)
    U = np.triu(np.random.rand(n, n))
    A = np.matmul(L, U)
    A1 = np.copy(A)
    LU(A1)
    # extract lower and upper triangles
    L1 = np.tril(A1, -1) + np.eye(n)
    U1 = np.triu(A1)
    # compute relative errors
    print(np.linalg.norm(L - L1, 'fro') / np.linalg.norm(L, 'fro'),
          np.linalg.norm(U - U1, 'fro') / np.linalg.norm(U, 'fro'),
          np.linalg.norm(A - np.matmul(L1, U1), 'fro') / np.linalg.norm(A, 'fro'))
