import pytest


@pytest.mark.parametrize(
    "test_case",
    [
        {
            "test_id": "salt_factor_from_dose",
            "kwargs": {"delivered_drug": "400 mg", "salt_factor": 0.8},
            "expected_label": "Dose of Salt",
            "expected_value": 500,
            "expected_unit": "milligram",
        },
        {
            "test_id": "salt_factor_from_delivered",
            "kwargs": {"delivered_drug": "400 mg", "dose_of_salt": "500mg"},
            "expected_label": "Salt Factor",
            "expected_value": 0.8,
            "expected_unit": "dimensionless",
        },
        {
            "test_id": "salt_factor_from_dose_salt",
            "kwargs": {"salt_factor": 0.8, "dose_of_salt": "500mg"},
            "expected_label": "Delivered Drug",
            "expected_value": 400,
            "expected_unit": "milligram",
        },
    ],
    ids=lambda x: x["test_id"],
)
def test_salt_factor(test_case):
    """
    Parametrized test for salt_factor function.

    Tests all three calculation scenarios:
    - Calculate Dose of Salt from Delivered Drug and Salt Factor
    - Calculate Salt Factor from Delivered Drug and Dose of Salt
    - Calculate Delivered Drug from Salt Factor and Dose of Salt
    """
    from sp_pharkin import salt_factor

    result = salt_factor(**test_case["kwargs"])

    assert result[0] == test_case["expected_label"]
    assert result[1] == test_case["expected_value"]
    assert result[2] == test_case["expected_unit"]
