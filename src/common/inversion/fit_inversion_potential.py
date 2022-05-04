"""
Fitting for the inversion barriers using the UFF term
see eqn. 19 in ref 1. Geometries obtained from PB0-D3BJ/def2-TZVPP
optimisations of the hydrides in the gas phase using ORCA v. 5.0

References:
[1] A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A. Goddard III, and
W. M. Skiff, J. Am. Chem. Soc. 114, 25, 10024 (1992)
[2] J. D. Swalen and J. A. Ibers, J. Chern. Phys. 36, 1914 (1962)
[3] R. E. Weston, J. Am. Chern. Soc. 76, 2645 (1954)
"""
import autode as ade
import numpy as np
import matplotlib.pyplot as plt

from typing import Union


class AtomType:
    """Atom type representing a group 15 'tetrahedral' atom"""

    def __init__(self,
                 name:              str,
                 inversion_barrier: float,
                 gamma_0:           float,
                 gamma_max:         float = np.pi/2):

        self.name = name
        self.gamma_0 = gamma_0
        self.gamma_max = gamma_max
        self.barrier = inversion_barrier

        self.k = 0
        self.c0 = 0
        self.c1 = 0
        self.c2 = 0

        if inversion_barrier > 0:
            self._fit_params()

    @property
    def params(self) -> tuple:
        return self.k, self.c0, self.c1, self.c2

    def _fit_params(self) -> None:
        """Fit parameters to the inversion barrier and gamma. Requires

        E(γ_0) = 0
        E'(γ_0) = 0       (and a minimum)

        E(γ_max) = barrier
        E'(γ_max) = 0
        """
        from numpy import sin, cos

        y0 = self.gamma_0
        ymax = self.gamma_max

        c2 = (self.barrier /
              (cos(2*ymax) + 2*sin(y0)*sin(ymax)/cos(y0)
               - cos(2*y0) - 2*(sin(y0)**2)/cos(y0)))

        c1 = 2*sin(2*y0)/cos(y0)*c2

        c0 = -c1*sin(y0) - c2*cos(2*y0)
        k = (self.barrier /
             (c0 + c1*sin(ymax) + c2*cos(2*ymax)))

        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.k = k

        return None

    def plot_energy_vs_gamma(self) -> None:
        """Plot the energy as a function of the angle γ"""

        gammas = np.linspace(0, np.pi, num=300)
        energies = energy(gammas, self)

        plt.plot(57.2958 * gammas, energies, label=self.name)
        plt.xlabel('γ / degrees')
        plt.ylabel('E / kcal mol$^{-1}$')
        plt.tight_layout()
        plt.legend()
        plt.savefig(f'plot.pdf')

        return None

    def __str__(self):
        return ('InversionCentre{' +
                f'name: {self.name}, '
                f'k: {self.k:.5f}, '
                f'c0: {self.c0:.5f}, c1: {self.c1:.5f}, c2: {self.c2:.5f}'
                + '}')


def gamma(_molecule: ade.Molecule) -> float:
    """Evaluate an angle (γ) between a normal and an c-l vector for a
    central atom c, in radians"""

    v0 = _molecule.atoms.nvector(0, 1)
    v1 = _molecule.atoms.nvector(0, 2)

    v2 = np.cross(v0, v1)
    v2 /= np.linalg.norm(v2)

    return np.arccos(np.dot(v2, _molecule.atoms.nvector(0, 3)))


def energy(_gamma:     Union[float, np.ndarray],
           _atom_type: AtomType
           ) -> float:
    """Energy at a particular angle (γ) """
    k, c0, c1, c2 = _atom_type.params

    return k*(c0 + c1 * np.sin(_gamma) + c2 * np.cos(2*_gamma))


def grad(_gamma:     Union[float, np.ndarray],
         _atom_type: AtomType) -> float:
    """dE/dγ"""
    k, c0, c1, c2 = _atom_type.params

    return k*(c1 * np.cos(_gamma) - c2 * 2 * np.sin(2*_gamma))


if __name__ == "__main__":

    N_3 = AtomType('N_3',
                   inversion_barrier=0,          # UFF
                   # inversion_barrier=5.783939, # ref. 2
                   # inversion_barrier=4.66,     # PBE0-D3BJ/def2-TZVPP
                   gamma_0=gamma(ade.Molecule('H3N.xyz')))
    print(N_3)
    N_3.plot_energy_vs_gamma()

    P_3 = AtomType('P_3',
                   inversion_barrier=22,         # UFF
                   # inversion_barrier=31.54876, # ref. 3
                   # inversion_barrier=33.03,    # PBE0-D3BJ/def2-TZVPP
                   gamma_0=gamma(ade.Molecule('H3P.xyz')))
    print(P_3)
    P_3.plot_energy_vs_gamma()

    As_3 = AtomType('As_3',
                    inversion_barrier=22,      # UFF
                    # inversion_barrier=38.60, # PBE0-D3BJ/def2-TZVPP
                    gamma_0=gamma(ade.Molecule('H3As.xyz')))
    print(As_3)
    As_3.plot_energy_vs_gamma()

    Sb_3 = AtomType('Sb_3',
                    inversion_barrier=22,      # UFF
                    # inversion_barrier=44.6,  # PBE0-D3BJ/def2-TZVPP
                    gamma_0=gamma(ade.Molecule('H3Sb.xyz')))
    print(Sb_3)
    Sb_3.plot_energy_vs_gamma()

    Bi_3 = AtomType('Bi_3',
                    inversion_barrier=22,      # UFF
                    gamma_0=gamma(ade.Molecule('H3Bi.xyz')))
    print(Bi_3)
    Bi_3.plot_energy_vs_gamma()
