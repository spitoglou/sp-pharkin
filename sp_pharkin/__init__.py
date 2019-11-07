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


if __name__ == '__main__':

    v = average_clearance_weight(
        average_clearance='0.04L/(hour*kilogram)',
        weight='75 kilogram'
    )
    print(v)

    v = average_clearance_weight(
        clearance='3L/hour',
        weight='75 kilogram'
    )
    print(v)

    print('4.1')
    v = half_life_k(
        K='0.0131 /hour',
    )
    print(v)

    print('4.2')
    v = clearance_flow_extraction_rate(
        Q='1.2L/min',
        E=0.04,
        output_unit='L/hour'
    )
    print(v)

    print('4.3')
    v = half_life_k(
        half_life='12hour',
        decimals=3
    )
    print(v)

    print('4.4')
    v = clearance_elimination_rate_constant_volume(
        volume='25L',
        K='0.03/hour',
        output_unit='mL/hour'
    )
    print(v)

    print('4.5')
    k = clearance_elimination_rate_constant_volume(
        volume='175L',
        clearance='0.43L/min',
        decimals=10,
        output_unit='1/hour'
    )
    print(k)
    v = half_life_k(
        K=k[4],
        decimals=10,
        output_unit='hour'
    )
    print(v)

    print('4.6')
    k = clearance_elimination_rate_constant_volume(
        volume='12L',
        clearance='12.5mL/min',
        decimals=4,
        output_unit='1/hour'
    )
    print(k)

    print('4.6')
    k = average_clearance_weight(
        average_clearance='1.5mL/min/kilogram',
        weight='85kilogram',
        output_unit='L/hour'
    )
    print(k)
