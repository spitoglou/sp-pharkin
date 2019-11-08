def test_chapter_3():
    from sp_pharkin import target_concentration, dose_concentration_volume, volume_of_distribution_weight
    from sp_pharkin.reduction_factors import salt_factor
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
    assert s[3] == '9.4 milligram'

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
    assert c2 > c1

    print('3.3')
    v = dose_concentration_volume(
        dose='200ug',
        concentration='12.5ng/mL',
        output_unit='L'
    )
    assert v[3] == '16.0 liter'

    print('3.4')
    d = salt_factor(
        dose_of_salt='5mg',
        salt_factor=0.7
    )
    print(d)
    c = dose_concentration_volume(
        dose=d[3],
        volume='70L',
        output_unit='ng/mL'
    )
    assert c[3] == '50.0 nanogram / milliliter'

    print('3.5')
    v = dose_concentration_volume(
        dose='0.5mg',
        concentration='20ng/mL',
        output_unit='L'
    )
    assert v[3] == '25.0 liter'

    print('3.6')
    tar_con = target_concentration('400 ng/mL', '700 ng/mL')
    c1 = tar_con[3]
    v1 = volume_of_distribution_weight(
        mean_volume_of_distribution='0.45 L/kilogram',
        weight='70 kilogram'
    )[3]
    d1 = dose_concentration_volume(
        volume=v1,
        concentration=c1
    )[3]
    d2 = salt_factor(
        delivered_drug=d1,
        salt_factor=0.75,
        output_unit='mg'
    )
    assert d2[3] == '23.1 milligram'
