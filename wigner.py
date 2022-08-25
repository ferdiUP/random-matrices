import numpy as np
import matplotlib.pyplot as plt


# Defining the semi-circle law
def demicercle(t):
    return (1/(2*np.pi))*np.sqrt(4-t**2)


# Defining the Marcenko-Pastur distribution
def mp(t, c):
    if max(1-1/c, 0) == 0:
        d0 = 0
    else:
        d0 = 1-1/c
    if t.any()>=(1-np.sqrt(c))**2 and t.any()<=(1+np.sqrt(c))**2:
        d1 = np.sqrt((t-(1-np.sqrt(c))**2)*((1+np.sqrt(c))**2-t))/(2*np.pi*c*t)
    else:
        d1 = 0
    return d0+d1


# Setting the dimensions
n = 500
N = 2*n
p = 150
M = np.zeros((N, N))
W = np.zeros((N, p))
U1, U2 = np.zeros((N, N)), np.zeros((N, N))
D = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        a = np.random.normal(0, 1)
        b = np.random.normal(0, 1)
        c = np.random.normal(0, 1)
        M[i, j] = a
        M[j, i] = a
        U1[i, j] = b
        U2[i, j] = c
for i in range(N):
    for j in range(p):
        a = np.random.normal(0, 3)
        W[i, j] = a
for i in range(n):
    D[i, i] = 1
for i in range(N):
    for j in range(N):
        U2[j, i] = U2[i, j]

W = 1/p*np.dot(W, np.transpose(W))

# Getting normalized eigenvalues of M
s1 = np.linalg.eig(1/np.sqrt(N)*M)
s2 = np.linalg.eig(W)

P1 = np.dot(np.dot(U1, D), np.transpose(U1))
P2 = np.dot(np.dot(U2, D), np.transpose(U2))
s3 = np.linalg.eig(P1)
s3vsym = np.linalg.eig(P2)
s4 = np.linalg.eigvals(P2)

# Histogram plotting for convergence illutration
fig = plt.figure()
plt.hist(s1[0], 40, density=True)
x1 = np.linspace(-2, 2, 1000)
plt.plot(x1, demicercle(x1), color='red')
plt.show()
