from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity


def format_output(quantity, string, output_unit, decimals):
    if output_unit:
        quantity = quantity.to(ureg(output_unit))

    quantity = round(quantity, decimals)

    return (string, quantity.magnitude, '{!s}'.format(quantity.units), '{!s}'.format(quantity), quantity)


def generic_a_eq_b_x_c(a, b, c, names):
    if a and c:
        string = names[1]
        quantity = a / c

    if a and b:
        string = names[2]
        quantity = a / b

    if b and c:
        string = names[0]
        quantity = b * c

    return (string, quantity)
