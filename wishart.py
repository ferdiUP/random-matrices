import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# Defining the Marcenko-Pastur distribution
def mp(t, c):
    if max(1-1/c, 0) == 0:
        d0 = 0
    else:
        d0 = 1-1/c
    if (1 - np.sqrt(c))**2 <= t.any() <= (1 + np.sqrt(c))**2:
        d1 = np.sqrt((t-(1-np.sqrt(c))**2)*((1+np.sqrt(c))**2-t))/(2*np.pi*c*t)
    else:
        d1 = 0
    return d0+d1


if __name__ == "__main__":
    n = 500
    # Different rectangle sizes
    p = np.array([1500, 1000, n])
    b = np.array([2.5, 3, 4])
    c = n/p
    print(c)
    plt.figure()
    for k in range(3):
        # Matrix with independent gaussian coefficients
        B = np.random.normal(0, 1, (n, p[k]))

        # Wishart matrix
        W = 1/p[k] * B @ B.T
        s = np.linalg.eigvals(W)

        # Plots
        plot = int("13"+str(k+1))
        plt.subplot(plot)

        # Eigenvalues
        plt.hist(s, 20, density=True)

        # MP distribution
        t = np.linspace(0, b[k], 1000)
        plt.plot(t, mp(t, c[k]), 'r-')

        plt.title('$\\frac{N}{p}=$'+str(c[k]))

    plt.show()
