from .lib import format_output, generic_a_eq_b_x_c
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


def volume_of_distribution_weight(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('volume_of_distribution', False)
    b = kwargs.get(
        'mean_volume_of_distribution', False)
    c = kwargs.get('weight', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Volume of Distribution', 'Mean Volume of Distribution', 'Weight'])

    return format_output(quantity, string, output_unit, decimals)


def dose_concentration_volume(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('dose', False)
    b = kwargs.get(
        'concentration', False)
    c = kwargs.get('volume', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Dose', 'Concentration', 'Volume'])

    return format_output(quantity, string, output_unit, decimals)


def target_concentration(min, max):
    result = (Q_(min) + Q_(max)) / 2
    return ('Target Concentration', result.magnitude, '{!s}'.format(result.units), '{!s}'.format(result), result)


def rate_of_elimination_mass_k(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('rate_of_elimination', False)
    b = kwargs.get(
        'mass', False)
    c = kwargs.get('K', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Rate of Elimination', 'Mass', 'Elimination Rate Constant(K)'])

    return format_output(quantity, string, output_unit, decimals)


def half_life_k(**kwargs):
    import math
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = math.log(2)
    b = kwargs.get(
        'K', False)
    c = kwargs.get('half_life', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Ln(2)', 'Elimination Rate Constant(K)', 'Half-Life'])

    return format_output(quantity, string, output_unit, decimals)


def extraction_rate(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get('c_diff', False)
    b = kwargs.get('E', False)
    c = kwargs.get('c_in', False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ['Concentration Difference(C_in - C_out)', 'Extraction Ratio(E)', 'Input Concentration(C_in)'])

    return format_output(quantity, string, output_unit, decimals)
