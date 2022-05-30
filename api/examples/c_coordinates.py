"""
Atomic coordinates can be set from the Python API
"""
import optrs

try:
    import autode as ade
except ModuleNotFoundError:
    exit("Failed to find autodE. Install with\n"
         "conda install autode -c conda-forge")

# Create a molecule from a list of atomic symbols (e.g. ["H", "H"])
# and set the coordinates from a flat list of coordinates [0.0, 1.0, 2.0, ...]
symbols = ade.Molecule('metallocage.xyz').atomic_symbols
coordinates = ade.Molecule('metallocage.xyz').coordinates.flatten().tolist()

mol = optrs.Molecule.from_atomic_symbols(symbols)
mol.set_coordinates(coordinates)

mol.generate_connectivty()  # Then build the connectivity (bonds, angles etc.),
mol.optimise()              # and optimise
