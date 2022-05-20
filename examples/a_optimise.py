"""Optimise a structure using the default (UFF) forcefield:"""
import mors

mol = mors.Molecule.from_xyz_file("metallocage.xyz")
mol.optimise()
mol.write_xyz_file("metallocage_opt.xyz")
