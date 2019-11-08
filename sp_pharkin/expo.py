from sympy import init_printing, symbols, exp, solve, Eq


def c_t(unknown, **kwargs):
    init_printing(use_unicode=True)
    c_t, c_0, k, t = symbols('c_t c_0 k t')
    eq1 = Eq(c_t, c_0 * exp(-1 * k * t))
    print(eq1)
    expr = solve(eq1, unknown)
    print(expr[0])
    result = expr[0].evalf(subs=kwargs)
    print(result)
