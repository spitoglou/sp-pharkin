

def test_bioavailability1():
    from sp_pharkin import bioavailability
    result = bioavailability(
        delivered_drug='400 mg',
        bioavailability=0.8
    )
    assert result[1] == 500
    assert result[2] == 'milligram'


def test_bioavailability2():
    from sp_pharkin import bioavailability
    result = bioavailability(
        delivered_drug='400 mg',
        dose_administered='500mg'
    )
    assert result[1] == 0.8
    assert result[2] == 'dimensionless'


def test_bioavailability3():
    from sp_pharkin import bioavailability
    result = bioavailability(
        bioavailability=0.8,
        dose_administered='500mg'
    )
    assert result[1] == 400
    assert result[2] == 'milligram'
