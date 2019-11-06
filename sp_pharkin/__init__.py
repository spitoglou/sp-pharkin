from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity


def salt_factor(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    delivered_drug = kwargs.get('delivered_drug', 0)
    dose_of_salt = kwargs.get('dose_of_salt', 0)
    salt_factor = kwargs.get('salt_factor', 0)

    if delivered_drug and salt_factor:
        string = 'Dose of Salt'
        quantity = delivered_drug / salt_factor

    if delivered_drug and dose_of_salt:
        string = 'Salt Factor'
        quantity = delivered_drug / dose_of_salt

    if dose_of_salt and salt_factor:
        string = 'Delivered Drug'
        quantity = dose_of_salt * salt_factor

    if output_unit:
        quantity = round(quantity.to(ureg(output_unit)), decimals)

    return (string, quantity.magnitude, '{!s}'.format(quantity.units), '{!s}'.format(quantity), quantity)


def bioavailability(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    delivered_drug = kwargs.get('delivered_drug', 0)
    dose_administered = kwargs.get('dose_administered', 0)
    bioavailability = kwargs.get('bioavailability', 0)

    if delivered_drug and bioavailability:
        string = 'Dose Administered'
        quantity = delivered_drug / bioavailability

    if delivered_drug and dose_administered:
        string = 'Bioavailability'
        quantity = delivered_drug / dose_administered

    if dose_administered and bioavailability:
        string = 'Delivered Drug'
        quantity = dose_administered * bioavailability

    if output_unit:
        quantity = round(quantity.to(ureg(output_unit)), decimals)

    return (string, quantity.magnitude, '{!s}'.format(quantity.units), '{!s}'.format(quantity), quantity)


def volume_of_distribution_weight(**kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    volume_of_distribution = kwargs.get('volume_of_distribution', False)
    mean_volume_of_distribution = kwargs.get(
        'mean_volume_of_distribution', False)
    weight = kwargs.get('weight', False)

    if volume_of_distribution and mean_volume_of_distribution:
        string = 'Weight'
        quantity = volume_of_distribution / mean_volume_of_distribution

    if volume_of_distribution and weight:
        string = 'Mean Volume of Distribution'
        quantity = volume_of_distribution / weight

    if mean_volume_of_distribution and weight:
        string = 'Volume of Distribution'
        quantity = mean_volume_of_distribution * weight

    if output_unit:
        quantity = round(quantity.to(ureg(output_unit)), decimals)

    return (string, quantity.magnitude, '{!s}'.format(quantity.units), '{!s}'.format(quantity), quantity)


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
