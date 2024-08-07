import numpy as np
import time

def thomas_algorithm (a , b , c , f):
    n = len (f)
    alpha = [0.0 for i in range(n)]
    beta = [0.0 for i in range(n)]
    alpha[0] = b[0]

    for i in range(1,n):
        beta[i] = a[i]/alpha[i-1]
        alpha[i] = b[i] - beta[i] * c[i-1]
    
    Y = [[0.0] for i in range(n)]
    X = [[0.0] for i in range(n)]

    Y[0][0] = f[0]
    for i in range(1,n):
        Y[i][0] = f[i] - beta[i] * Y[i-1][0]

    X[-1][0] = Y[-1][0]/alpha[-1]
    for i in range(n-2,-1,-1):
        X[i][0] = (Y[i][0] - c[i] * X[i+1][0])/alpha[i]
    
    return X

a = [ -1 for i in range (1000) ]
a[0] = 0.0
b = [ 5 for i in range (1000) ]
c = [ -1 for i in range (1000) ]
c[-1] = 0.0
f = [ 1.0 for i in range (1000) ]

# a = [0 ,16, 15,-2]
# b = [4, -1, 8, 13]
# c = [-1, 2, 5, 0]
# f = [13, 65, 66, 9]

start_time = time.time()
X = thomas_algorithm (a , b , c , f)
end_time = time.time()

print (' X =' , np.array(X))

runtime = end_time - start_time
print ("Runtime : " , runtime , " seconds ")

