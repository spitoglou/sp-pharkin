from .lib import format_output, generic_a_eq_b_x_c
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity


def clearance_flow_extraction_rate(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('clearance', False)
    b = kwargs.get('Q', False)
    c = kwargs.get('E', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Clearance', 'Flow(Q)', 'Extraction Rate(E)'])

    return format_output(quantity, string, output_unit, decimals)


def clearance_elimination_rate_constant_volume(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('clearance', False)
    b = kwargs.get('K', False)
    c = kwargs.get('volume', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Clearance', 'Elimination Rate Constant(K)', 'Volume'])

    return format_output(quantity, string, output_unit, decimals)


def average_clearance_weight(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('clearance', False)
    b = kwargs.get('average_clearance', False)
    c = kwargs.get('weight', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Clearance', 'Average Clearance', 'Weight'])

    return format_output(quantity, string, output_unit, decimals)
