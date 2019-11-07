

def test_elimination_rate_constant():
    from sp_pharkin import rate_of_elimination_mass_k

    result = rate_of_elimination_mass_k(
        rate_of_elimination='1mg/hr',
        mass='10mg'
    )
    assert result[0] == 'Elimination Rate Constant(K)'
    assert result[1] == 0.1
    assert result[2] == '1 / hour'


def test_half_life():
    from sp_pharkin import half_life_k, rate_of_elimination_mass_k
    r = rate_of_elimination_mass_k(
        rate_of_elimination='1mg/hr',
        mass='10mg'
    )
    result = half_life_k(
        K=r[4]
    )
    assert result[0] == 'Half-Life'
    assert result[1] == 6.93
    assert result[2] == 'hour'


def test_extraction_ratio():
    from sp_pharkin import extraction_rate
    result = extraction_rate(
        c_in=10,
        c_diff=6
    )
    assert result[0] == 'Extraction Ratio(E)'
    assert result[1] == 0.6
    assert result[2] == 'dimensionless'

    result = extraction_rate(
        c_diff='6 mg/L',
        E=0.6
    )
    assert result[0] == 'Input Concentration(C_in)'
    assert result[1] == 10.0
    assert result[2] == 'milligram / liter'


def test_clearance():
    from sp_pharkin import clearance_flow_extraction_rate, clearance_elimination_rate_constant_volume
    result = clearance_flow_extraction_rate(
        Q='2L/min',
        E=0.5
    )
    assert result[0] == 'Clearance'
    assert result[1] == 1.0
    assert result[2] == 'liter / minute'

    result = clearance_flow_extraction_rate(
        clearance='1L/min',
        E=0.5
    )
    assert result[0] == 'Flow(Q)'
    assert result[1] == 2.0
    assert result[2] == 'liter / minute'

    result = clearance_elimination_rate_constant_volume(
        K='0.1/hour',
        volume='50L'
    )
    assert result[0] == 'Clearance'
    assert result[1] == 5.0
    assert result[2] == 'liter / hour'

    result = clearance_elimination_rate_constant_volume(
        K='0.1/hour',
        clearance='5L/hour'
    )
    assert result[0] == 'Volume'
    assert result[1] == 50.0
    assert result[2] == 'liter'
