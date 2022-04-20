
DEG_TO_RAD = 0.017453292519943295


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

        self.valency = 0
        self._set_valency()

        self.oxidation_state = 0
        self._set_oxidation_state()

    def _set_valency(self) -> None:

        try:
            self.valency = int(next(char.isdigit() for char
                                    in self.name.split('_')[0]))
        except StopIteration:
            pass

        if self.name.endswith('_'):
            self.valency = 1

        if self.bridging:
            self.valency = 2

        return None

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
