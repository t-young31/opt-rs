
DEG_TO_RAD = 0.017453292519943295

n_to_desc = {
    '1': 'linear',
    '2': 'trigonal',
    'R': 'resonant',
    '3': 'tetrahedral',
    '4': 'square planar',
    '5': 'trigonal bipyramidal',
    '6': 'octahedral'
}

n_to_valency = {
    '1': 2,
    '2': 3,
    'R': 3,
    '3': None,  # Depends on group
    '4': 4,
    '5': 5,
    '6': 6
}

group14_elements = ['C', 'Si', 'Ge', 'Sn', 'Pb', 'Fl']
group15_elements = ['N', 'P', 'As', 'Sb', 'Bi', 'Mc']
group16_elements = ['O', 'S', 'Se', 'Te', 'Po', 'Lv']


class AtomType:

    def __init__(self, string: str):

        items = string.split()

        if len(items) != 7:
            raise ValueError(f'Could not create an atom type from {string}. '
                             f'Must have 7 fields')

        self.name, self.r, self.theta, self.x, self.d, self.zeta, self.z = items

        self.theta = DEG_TO_RAD * float(self.theta)
        self.bridging = self.name.endswith('_b')
        self.aromatic = self.name.endswith('_R')
        self.atomic_symbol = self.name[0] if '_' in self.name else self.name[:2]

        self.valency = self._valency()

        self.oxidation_state = 0
        self._set_oxidation_state()

    def _valency(self) -> int:

        if self.name.endswith('_'):
            return 1

        if self.bridging:
            return 2

        if len(self.name) <= 2:
            return 0

        val = n_to_valency[self.name[2]]

        if isinstance(val, int):
            return val

        if self.name[2] != '3':
            raise RuntimeError

        if self.atomic_symbol in group14_elements:
            return 4

        elif self.atomic_symbol in group15_elements:
            return 3

        elif self.atomic_symbol in group16_elements:
            return 2

        return 4

    def _set_oxidation_state(self) -> None:

        try:
            val = self.name.split('+')[1]
            self.oxidation_state = int(val) if val != 'q' else 4
        except IndexError:
            pass

        return None

    @property
    def rust_struct(self) -> str:
        return ('UFFAtomType{' +
                f'name: "{self.name}", '
                f'atomic_symbol: "{self.atomic_symbol}", '
                f'bridging: {"true" if self.bridging else "false"}, '
                f'aromatic: {"true" if self.aromatic else "false"}, '
                f'valency: {self.valency}, '
                f'oxidation_state: {self.oxidation_state}, '
                f'environment: CoordinationEnvironment::None, '
                f'r: {self.r}, '
                f'theta: {self.theta:.5f}, '
                f'x: {self.x}, '
                f'd: {self.d}, '
                f'zeta: {self.zeta}, '
                f'z_eff: {self.z}'
                + '}')


if __name__ == '__main__':

    atom_types = [AtomType(line) for line in open('atom_types.txt', 'r')]

    with open('atom_types.rs', 'w') as file:

        print('use crate::ff::uff::atom_typing::{CoordinationEnvironment, '
              'UFFAtomType};\n\n\n'
              'pub(crate) const ATOM_TYPES: [UFFAtomType; 127] = [', file=file)

        for atom_type in atom_types:
            print(f'{atom_type.rust_struct}, ', file=file)

        print('];', file=file)
