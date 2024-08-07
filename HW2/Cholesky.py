from math import sqrt
import numpy as np

def cholesky(A):
    n = len(A)

    # Create zero matrix for L
    L = [[0.0] * n for i in range(n)]

    # Perform the Cholesky decomposition
    for i in range(n):
        for j in range(i+1):
            tmp_sum = sum(L[i][k] * L[j][k] for k in range(j))
            
            if (i == j): # Diagonal elements
                L[i][j] = sqrt(A[i][i] - tmp_sum)
            else:
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - tmp_sum))
    return L

def solveAXb_cholesky(A,b):
    n = len(A)
    L = cholesky(A)

    Y = [[0.0] for i in range(n)]
    X = [[0.0] for i in range(n)]

    for i in range(n):
        tmp_sum = sum(L[i][k] * Y[k][0] for k in range(i))
        Y[i][0] = (b[i][0]-tmp_sum)/L[i][i]
    
    for i in range(n-1,-1,-1):
        tmp_sum = sum(L[k][i] * X[k][0] for k in range(i+1,n))
        X[i][0] = (Y[i][0] - tmp_sum)/L[i][i]
    
    return X

A = [[9, 3, 6], 
     [3, 10, 5], 
     [6, 5, 21]]

b = [[3],
     [-5],
     [16]]

print("A:")
print(np.array(A))

print("B :")
print(np.array(b))

X = solveAXb_cholesky(A,b)
print("X :")
print(np.array(X))

