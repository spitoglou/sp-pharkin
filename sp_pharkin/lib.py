from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def format_output(quantity, string, output_unit, decimals):
    if output_unit:
        quantity = quantity.to(ureg(output_unit))

    quantity = round(quantity, decimals)

    return (
        string,
        quantity.magnitude,
        "{!s}".format(quantity.units),
        "{!s}".format(quantity),
        quantity,
    )


def generic_a_eq_b_x_c(a, b, c, names):
    """
    Solve for unknown in equation: a = b × c

    Given any two of three variables, compute the third.

    Args:
        a, b, c: Variables (can be False/None if unknown)
        names: Tuple of (name_a, name_b, name_c) for return values

    Returns:
        Tuple of (name, quantity) where quantity is the computed result

    Raises:
        ValueError: If not exactly 2 parameters are provided
    """
    provided = sum([bool(a), bool(b), bool(c)])

    if provided != 2:
        raise ValueError(
            f"generic_a_eq_b_x_c requires exactly 2 of 3 parameters. "
            f"Got {provided}: a={a is not False}, b={b is not False}, c={c is not False}"
        )

    if a and c:
        string = names[1]
        quantity = a / c

    elif a and b:
        string = names[2]
        quantity = a / b

    elif b and c:
        string = names[0]
        quantity = b * c

    return (string, quantity)
