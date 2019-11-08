from .lib import generic_a_eq_b_x_c, format_output
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity


def salt_factor(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('delivered_drug', 0)
    b = kwargs.get('dose_of_salt', 0)
    c = kwargs.get('salt_factor', 0)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Delivered Drug', 'Dose of Salt', 'Salt Factor'])

    return format_output(quantity, string, output_unit, decimals)


def bioavailability(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('delivered_drug', 0)
    b = kwargs.get('dose_administered', 0)
    c = kwargs.get('bioavailability', 0)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Delivered Drug', 'Dose Administered', 'Bioavailability'])

    return format_output(quantity, string, output_unit, decimals)
