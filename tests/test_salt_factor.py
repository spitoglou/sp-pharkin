

def test_salt_factor1():
    from sp_pharkin import salt_factor
    result = salt_factor(
        delivered_drug='400 mg',
        salt_factor=0.8
    )
    assert result[1] == 500
    assert result[2] == 'milligram'


def test_salt_factor2():
    from sp_pharkin import salt_factor
    result = salt_factor(
        delivered_drug='400 mg',
        dose_of_salt='500mg'
    )
    assert result[1] == 0.8
    assert result[2] == 'dimensionless'


def test_salt_factor3():
    from sp_pharkin import salt_factor
    result = salt_factor(
        salt_factor=0.8,
        dose_of_salt='500mg'
    )
    assert result[1] == 400
    assert result[2] == 'milligram'
