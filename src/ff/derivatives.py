import sympy as sym

x_i = sym.Symbol('x_i')
y_i = sym.Symbol('y_i')
z_i = sym.Symbol('z_i')

x_j = sym.Symbol('x_j')
y_j = sym.Symbol('y_j')
z_j = sym.Symbol('z_j')

x_k = sym.Symbol('x_k')
y_k = sym.Symbol('y_k')
z_k = sym.Symbol('z_k')

x_l = sym.Symbol('x_l')
y_l = sym.Symbol('y_l')
z_l = sym.Symbol('z_l')


def theta():
    """Angle between two vectors j-->i  and j-->k"""

    dxi = x_i - x_j
    dyi = y_i - y_j
    dzi = z_i - z_j
    mod_r_ij = sym.sqrt(dxi*dxi + dyi*dyi + dzi*dzi)

    dxk = x_k - x_j
    dyk = y_k - y_j
    dzk = z_k - z_j
    mod_r_ik = sym.sqrt(dxk*dxk + dyk*dyk + dzk*dzk)

    dot = dxi*dxk + dyi*dyk + dzi*dzk

    return sym.acos(dot / (mod_r_ij * mod_r_ik))


def angle_type_a(symbol):

    k = sym.Symbol('self.k_ijk')
    n = sym.Symbol('self.n')

    string = rust_expression(
        sym.diff((k/n**2) * (1 - sym.cos(n*theta())), symbol)
    )
    return f'gradient[self.{str(symbol)[-1]}].{str(symbol)[0]} += {string}'


def angle_type_b(symbol):

    k = sym.Symbol('self.k_ijk')
    c0 = sym.Symbol('self.c0__')
    c1 = sym.Symbol('self.c1__')
    c2 = sym.Symbol('self.c2__')

    theta_ = theta()
    res = sym.diff(k*(c0 + c1*sym.cos(theta_) + c2*sym.cos(2*theta_)), symbol)
    string = rust_expression(function=res)
    # string = rust_expression(function=sym.refine(res, sym.Q.real(symbol)))
    string = string.replace('__', '')

    return f'gradient[self.{str(symbol)[-1]}].{str(symbol)[0]} += {string}'


def dihedral(symbol):

    r_ij_x = x_i - x_j
    r_ij_y = y_i - y_j
    r_ij_z = z_i - z_j

    r_lk_x = x_l - x_k
    r_lk_y = y_l - y_k
    r_lk_z = z_l - z_k

    r_kj_x = x_k - x_j
    r_kj_y = y_k - y_j
    r_kj_z = z_k - z_j

    def cross_x(u_x, u_y, u_z, v_x, v_y, v_z):
        return u_y*v_z - u_z*v_y

    def cross_y(u_x, u_y, u_z, v_x, v_y, v_z):
        return u_z*v_x - u_x*v_z

    def cross_z(u_x, u_y, u_z, v_x, v_y, v_z):
        return u_x*v_y - u_y*v_x

    def mod(u_x, u_y, u_z):
        return sym.sqrt(u_x**2 + u_y**2 + u_z**2)

    def dot(u_x, u_y, u_z, v_x, v_y, v_z):
        return u_x*v_x + u_y*v_y + u_z*v_z

    v0_x = cross_x(r_ij_x, r_ij_y, r_ij_z, r_kj_x, r_kj_y, r_kj_z)
    v0_y = cross_y(r_ij_x, r_ij_y, r_ij_z, r_kj_x, r_kj_y, r_kj_z)
    v0_z = cross_z(r_ij_x, r_ij_y, r_ij_z, r_kj_x, r_kj_y, r_kj_z)
    mod_v0 = mod(v0_x, v0_y, v0_z)

    v0_x /= mod_v0
    v0_y /= mod_v0
    v0_z /= mod_v0

    v1_x = cross_x(-r_kj_x, -r_kj_y, -r_kj_z, r_lk_x, r_lk_y, r_lk_z)
    v1_y = cross_y(-r_kj_x, -r_kj_y, -r_kj_z, r_lk_x, r_lk_y, r_lk_z)
    v1_z = cross_z(-r_kj_x, -r_kj_y, -r_kj_z, r_lk_x, r_lk_y, r_lk_z)
    mod_v1 = mod(v1_x, v1_y, v1_z)
    v1_x /= mod_v1
    v1_y /= mod_v1
    v1_z /= mod_v1

    mod_r_kj = mod(r_kj_x, r_kj_y, r_kj_z)
    r_kj_x /= mod_r_kj
    r_kj_y /= mod_r_kj
    r_kj_z /= mod_r_kj

    """
    value = -np.arctan2(np.dot(np.cross(vec1, vec_xy), vec2),
                    np.dot(vec1, vec2))
    """
    v2_x = cross_x(v0_x, v0_y, v0_z, r_kj_x, r_kj_y, r_kj_z)
    v2_y = cross_y(v0_x, v0_y, v0_z, r_kj_x, r_kj_y, r_kj_z)
    v2_z = cross_z(v0_x, v0_y, v0_z, r_kj_x, r_kj_y, r_kj_z)

    phi = -sym.atan2(dot(v2_x, v2_y, v2_z, v1_x, v1_y, v1_z),
                     dot(v0_x, v0_y, v0_z, v1_x, v1_y, v1_z))

    n = sym.Symbol("self.n_phi")
    phi0 = sym.Symbol("self.phi0")

    res = (0.5 * sym.Symbol("self.v_phi")
           * (1 - sym.cos(n*phi0)*sym.cos(n*phi)))

    string = rust_expression(function=sym.diff(res, symbol))
    # string = rust_expression(function=sym.refine(res, sym.Q.real(symbol)))
    string = string.replace('__', '')

    return f'gradient[self.{str(symbol)[-1]}].{str(symbol)[0]} += {string}'


def convert_to_trait(string, name):
    """Convert a prepend function call to a suffix for a rust trait"""

    if f'{name}' not in string:
        return string

    string = string.replace(f'a{name}', 'X')
    terms = string.split(f'{name}')

    for i, term in enumerate(terms[1:]):
        count = 0

        for j, char in enumerate(term):
            if char == '(':
                count += 1

            if char == ')':
                count -= 1

            if count == 0:
                term = term[:j+1] + f'.{name}()' + term[j+1:]
                break

        terms[i+1] = term

    string = ''.join(terms)
    string = string.replace('X', f'a{name}')

    return string


def convert_to_trait_with_arg(string, name):
    """Convert a dyadic function to a trait with a single argument"""
    terms = string.split(f'{name}')

    for i, term in enumerate(terms[1:]):
        count = 0

        for j, char in enumerate(term):
            if char == '(':
                count += 1

            if char == ')':
                count -= 1

            if count == 1 and char == ',':
                term = term[:j] + f').{name}(' + term[j+1:]
                break

        terms[i+1] = term

    return ''.join(terms)


def rust_expression(function):
    """Print a valid rust expression from a sympy function"""

    string = str(function).replace('**2', '.powi(2)')
    string = string.replace('+ 1', '+ 1.')
    string = string.replace('**(3/2)', '.powf(1.5)')
    string = string.replace('/2', '/2.')
    string = string.replace('2*', '2.*')

    for fn in ('sqrt', 'acos', 'cos', 'asin', 'sin'):
        string = convert_to_trait(string, fn)

    string = convert_to_trait_with_arg(string, 'atan2')

    return string


def truncated_expression(string):
    """Truncate a bracketed expression down to something that has
    matching parentheses e.g.

    (z_i - z_j).powi(2)).p   --->   (z_i - z_j)
    """

    count = 0
    final_idx = 0

    for i, char in enumerate(string):
        count += 1 if char == '(' else (-1 if char == ')' else 0)

        if count == 0 and char == ')':
            final_idx = i+1

    if final_idx == 0:
        raise RuntimeError("Unclosed bracket")

    return string[:final_idx]


def first_avail_variable_name(_lines):
    """Generate a variable name that has not yet been used"""

    i = 0
    while any(f'let v{i}' in l for l in _lines):
        i += 1

    return f'v{i}'


def extract_variables(_lines):
    """Given a set of rust code lines extract the common variables to them"""

    variables = set()

    non_var_lines = [l for l in _lines if 'let' not in l]

    first_line = non_var_lines[0]
    start_idx = 0 if '=' not in first_line else len(first_line.split('= ')[0])
    end_idx = min([len(l) for l in non_var_lines])//2

    for i in range(start_idx, end_idx):

        for j in range(i, end_idx):
            substring = first_line[i:j]

            if not all(substring in l for l in non_var_lines):
                break

        if substring.startswith('(') or substring.startswith('v'):
            try:
                variable = truncated_expression(substring)

                if len(variable) > 5:
                    variables.add(variable)
            except RuntimeError:
                continue

    if len(variables) == 0:
        return _lines

    smallest_var = next(iter(sorted(list(variables), key=lambda x: len(x))))

    var_name = first_avail_variable_name(_lines)
    n_var_lines = sum('let' in l for l in _lines)

    for i, line in enumerate(_lines[n_var_lines:]):
        _lines[n_var_lines+i] = line.replace(smallest_var, var_name)

    if '.' not in smallest_var:
        # Remove unnecessary parentheses
        smallest_var = smallest_var[1:-1]

    _lines.insert(n_var_lines, f'let {var_name} = {smallest_var}')

    return None


def print_angle_gradient(angle_type):

    lines = []

    for x in (x_i, y_i, z_i, x_j, y_j, z_j, x_k, y_k, z_k):
        lines.append(angle_type(x))

    for _ in range(40):
        extract_variables(lines)

    return print(";\n".join(lines))


def print_dihedral_gradient():

    lines = []

    for x in (x_i, y_i, z_i, x_j, y_j, z_j, x_k, y_k, z_k, x_l, y_l, z_l):
        lines.append(dihedral(x))

    for _ in range(1):
        extract_variables(lines)

    return print(";\n".join(lines))


if __name__ == '__main__':

    # print_angle_gradient(angle_type=angle_type_a)
    # print_angle_gradient(angle_type=angle_type_b)
    print_dihedral_gradient()
