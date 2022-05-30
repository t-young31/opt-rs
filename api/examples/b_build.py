"""
Molecular structures can be built from only atomic symbols and the bond order
matrix (defined between all pairs of atoms):
"""
import optrs

mol = optrs.Molecule.from_atomic_symbols(["O", "H", "H"])

mol.set_bond_orders([0., 1., 1.,
                     1., 0., 0.,
                     1., 0., 0.])

mol.build_3d()  # Builds an rough 3D geometry
mol.optimise()  # which needs to be optimised
mol.write_xyz_file("water_opt.xyz")
