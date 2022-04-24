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


def rust_expression(function):
    """Print a valid rust expression from a sympy function"""

    string = str(function).replace('**2', '.powi(2)')
    string = string.replace('+ 1', '+ 1.')
    string = string.replace('**(3/2)', '.powf(1.5)')
    string = string.replace('/2', '/2.')
    string = string.replace('2*', '2.*')

    for fn in ('sqrt', 'acos', 'cos', 'asin', 'sin'):
        string = convert_to_trait(string, fn)

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
    end_idx = min([len(l) for l in non_var_lines])

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
        # Remove unnesacerray parentheses
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


if __name__ == '__main__':

    print_angle_gradient(angle_type=angle_type_a)
    # print_angle_gradient(angle_type=angle_type_b)
