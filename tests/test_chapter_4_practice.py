def test_chapter_4_practice_questions():
    import sp_pharkin.clearance as sppkcl
    import sp_pharkin as sppk
    print('4.1')
    v = sppk.half_life_k(
        K='0.0131 /hour',
    )
    assert v[3] == '52.91 hour'

    print('4.2')
    v = sppkcl.clearance_flow_extraction_rate(
        Q='1.2L/min',
        E=0.04,
        output_unit='L/hour'
    )
    assert v[3] == '2.88 liter / hour'

    print('4.3')
    v = sppk.half_life_k(
        half_life='12hour',
        decimals=3
    )
    assert v[3] == '0.058 / hour'

    print('4.4')
    v = sppkcl.clearance_elimination_rate_constant_volume(
        volume='25L',
        K='0.03/hour',
        output_unit='mL/hour'
    )
    assert v[3] == '750.0 milliliter / hour'

    print('4.5')
    k = sppkcl.clearance_elimination_rate_constant_volume(
        volume='175L',
        clearance='0.43L/min',
        decimals=10,
        output_unit='1/hour'
    )
    print(k)
    v = sppk.half_life_k(
        K=k[3],
        output_unit='hour'
    )
    assert v[3] == '4.7 hour'

    print('4.6')
    k = sppkcl.clearance_elimination_rate_constant_volume(
        volume='12L',
        clearance='12.5mL/min',
        decimals=4,
        output_unit='1/hour'
    )
    assert k[3] == '0.0625 / hour'

    print('4.7')
    k = sppkcl.average_clearance_weight(
        average_clearance='1.5mL/min/kilogram',
        weight='85kilogram',
        output_unit='L/hour'
    )
    print(k)
    assert k[3] == '7.65 liter / hour'
