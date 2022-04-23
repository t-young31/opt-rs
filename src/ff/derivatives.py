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

    res = sym.diff((k/n**2) * (1 - sym.cos(n*theta())), symbol)

    print(f'gradient[self.{str(symbol)[-1]}].{str(symbol)[0]} +=',
          rust_expression(function=sym.refine(res, sym.Q.real(symbol))),
          end=';\n')

    return None


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

    print(f'gradient[self.{str(symbol)[-1]}].{str(symbol)[0]} +=',
          string,
          end=';\n')

    return None


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


if __name__ == '__main__':

    for x in (x_i, y_i, z_i, x_j, y_j, z_j, x_k, y_k, z_k):
        angle_type_a(x)
        # angle_type_b(x)
