"""
Plot of the energy as a function of angle for a water molecule
using the standard UFF implementation of a bending potential for a
nonlinear case
"""
import numpy as np
import matplotlib.pyplot as plt


def e(theta):
    return k * (c0 + c1 * np.cos(theta) + c2 * np.cos(2*theta))


if __name__ == '__main__':

    k = 325.16182     # kcal mol rad-1
    theta0 = 1.82400  # rad

    c2 = 1 / (4 * np.sin(theta0)**2)
    c1 = -4 * c2 * np.cos(theta0)
    c0 = c2 * (2*np.cos(theta0)**2 + 1)

    print(c0, c1, c2)

    thetas = np.linspace(0.5, 3.145, num=300)

    plt.plot(thetas,
             e(thetas))
    plt.xlabel(r'$\theta$ / rad')
    plt.ylabel('E')
    plt.tight_layout()
    plt.savefig('harmonic_b_water.pdf')
