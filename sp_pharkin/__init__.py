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


if __name__ == '__main__':
    s = salt_factor(
        output_unit='ug',
        salt_factor=.8,
        delivered_drug='400 mg'
    )
    print(s)
    s = bioavailability(
        output_unit='ng',
        bioavailability=.8,
        dose_administered='400 mg'
    )
    print(s)
    s = volume_of_distribution_weight(
        output_unit='L',
        mean_volume_of_distribution='0.5 mL/kilogram',
        weight='80 kilogram'
    )
    print(s)

    print('3.1')
    tar_con = target_concentration('150 ug/L', '250 ug/L')
    c1 = tar_con[4]
    v1 = volume_of_distribution_weight(
        mean_volume_of_distribution='0.72 L/kilogram',
        weight='65 kilogram'
    )[4]
    s = dose_concentration_volume(
        concentration=c1,
        volume=v1,
        output_unit='mg',
        decimals=1
    )
    print(s)

    print('3.2')
    v1 = volume_of_distribution_weight(
        mean_volume_of_distribution='0.91 L/kilogram',
        weight='80 kilogram'
    )[4]
    tar_con = target_concentration('20 ug/L', '40 ug/L')
    c1 = tar_con[4]
    c2 = dose_concentration_volume(
        dose='5mg',
        volume=v1,
        output_unit='ug/L',
        decimals=0
    )[4]
    print(c2, c1)

    print('3.3')
    v = dose_concentration_volume(
        dose='200ug',
        concentration='12.5ng/mL',
        output_unit='L'
    )
    print(v)

    print('3.4')
    d = salt_factor(
        dose_of_salt='5mg',
        salt_factor=0.7
    )[4]
    c = dose_concentration_volume(
        dose=d,
        volume='70L',
        output_unit='ng/mL'
    )
    print(c)

    print('3.5')
    v = dose_concentration_volume(
        dose='0.5mg',
        concentration='20ng/mL',
        output_unit='L'
    )
    print(v)

    print('3.6')
    tar_con = target_concentration('400 ng/mL', '700 ng/mL')
    c1 = tar_con[4]
    v1 = volume_of_distribution_weight(
        mean_volume_of_distribution='0.45 L/kilogram',
        weight='70 kilogram'
    )[4]
    d1 = dose_concentration_volume(
        volume=v1,
        concentration=c1
    )[4]
    d2 = salt_factor(
        delivered_drug=d1,
        salt_factor=0.75,
        output_unit='mg'
    )
    print(d2)

    r = rate_of_elimination_mass_k(
        rate_of_elimination='1mg/hr',
        mass='10mg'
    )
    print(r)
    t = half_life_k(
        K=r[4]
    )
    print(t)

    e = extraction_rate(
        c_in=10,
        c_diff=6
    )
    print(e)

    e = extraction_rate(
        c_diff='6 mg/L',
        E=0.6
    )
    print(e)

    cl = clearance_flow_extraction_rate(
        Q='2L/min',
        E=0.5
    )
    print(cl)

    Q = clearance_flow_extraction_rate(
        clearance='1L/min',
        E=0.5
    )
    print(Q)

    cl = clearance_elimination_rate_constant_volume(
        K='0.1/hour',
        volume='50L'
    )
    print(cl)

    v = clearance_elimination_rate_constant_volume(
        K='0.1/hour',
        clearance='5L/hour'
    )
    print(v)
