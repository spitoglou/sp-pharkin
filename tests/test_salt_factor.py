

def test_salt_factor1():
    from sp_pharkin import salt_factor
    assert salt_factor(
        delivered_drug=400,
        salt_factor=0.8
    )[1] == 500


def test_salt_factor2():
    from sp_pharkin import salt_factor
    assert salt_factor(
        delivered_drug=400,
        dose_of_salt=500
    )[1] == 0.8


def test_salt_factor3():
    from sp_pharkin import salt_factor
    assert salt_factor(
        salt_factor=0.8,
        dose_of_salt=500
    )[1] == 400
