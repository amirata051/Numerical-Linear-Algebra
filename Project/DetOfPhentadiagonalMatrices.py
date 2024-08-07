import numpy as np

def compute_determinant(n, beta, alpha, a, b, c):
    # Step 1: Initialize A and B sequences
    A = np.zeros(n + 2)
    B = np.zeros(n + 2)
    
    A[0] = 0
    A[1] = 1
    B[0] = 1
    B[1] = 0
    
    # Step 2: Calculate A_i(0) and B_i(0) using recurrence relations
    for i in range(2, n + 2):
        if i == 2:
            A[i] = - (b[0] * A[1] + a[0] * A[0]) / c[0]
            B[i] = - (b[0] * B[1] + a[0] * B[0]) / c[0]
        elif i == 3:
            A[i] = - (c[0] * A[i-2] + b[1] * A[i-1] + a[1] * A[i-1] + a[1] * A[i-2] + alpha[0] * A[i-3]) / beta[0]
            B[i] = - (c[0] * B[i-2] + b[1] * B[i-1] + a[1] * B[i-1] + a[1] * B[i-2] + alpha[0] * B[i-3]) / beta[0]
        else:
            A[i] = - (c[i-2] * A[i-2] + b[i-2] * A[i-1] + a[i-2] * A[i-1] + alpha[i-2] * A[i-3]) / beta[i-2]
            B[i] = - (c[i-2] * B[i-2] + b[i-2] * B[i-1] + a[i-2] * B[i-1] + alpha[i-2] * B[i-3]) / beta[i-2]

    # Step 3: Compute Q_{n+1}(0)
    Q_n1 = A[n] * B[n+1] - A[n+1] * B[n]
    
    # Step 4: Compute determinant
    determinant = - Q_n1
    
    return determinant

def create_pentadiagonal_matrix(n, beta, alpha, a, b, c):
    P = np.zeros((n, n))
    for i in range(n):
        P[i, i] = a[i]
        if i - 1 >= 0:
            P[i, i - 1] = beta[i - 1]
        if i - 2 >= 0:
            P[i, i - 2] = alpha[i - 2]
        if i + 1 < n:
            P[i, i + 1] = b[i]
        if i + 2 < n:
            P[i, i + 2] = c[i]
    return P

def test_compute_determinant():
    # Test case
    n = 6
    beta = [1, 1, 1, 1]        # beta1 to beta4
    alpha = [1, 1, 1, 1, 1]    # alpha1 to alpha5
    a = [1, 1, 1, 1, 1, 1]     # a1 to a6
    b = [1, 1, 1, 1, 1]        # b1 to b5
    c = [1, 1, 1, 1]           # c1 to c4

    # Calculate determinant using the function
    computed_determinant = compute_determinant(n, beta, alpha, a, b, c)

    # Create pentadiagonal matrix
    P = create_pentadiagonal_matrix(n, beta, alpha, a, b, c)

    # Calculate determinant using numpy
    numpy_determinant = np.linalg.det(P)

    print(f"Computed determinant: {computed_determinant}")
    print(f"Numpy determinant: {numpy_determinant}")

# Run the test
test_compute_determinant()
