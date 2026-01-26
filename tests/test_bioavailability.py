import pytest


@pytest.mark.parametrize(
    "test_case",
    [
        {
            "test_id": "bioavailability_from_dose",
            "kwargs": {"delivered_drug": "400 mg", "bioavailability": 0.8},
            "expected_label": "Dose Administered",
            "expected_value": 500,
            "expected_unit": "milligram",
        },
        {
            "test_id": "bioavailability_from_delivered",
            "kwargs": {"delivered_drug": "400 mg", "dose_administered": "500mg"},
            "expected_label": "Bioavailability",
            "expected_value": 0.8,
            "expected_unit": "dimensionless",
        },
        {
            "test_id": "bioavailability_from_dose_admin",
            "kwargs": {"bioavailability": 0.8, "dose_administered": "500mg"},
            "expected_label": "Delivered Drug",
            "expected_value": 400,
            "expected_unit": "milligram",
        },
    ],
    ids=lambda x: x["test_id"],
)
def test_bioavailability(test_case):
    """
    Parametrized test for bioavailability function.

    Tests all three calculation scenarios:
    - Calculate Dose Administered from Delivered Drug and Bioavailability
    - Calculate Bioavailability from Delivered Drug and Dose Administered
    - Calculate Delivered Drug from Bioavailability and Dose Administered
    """
    from sp_pharkin import bioavailability

    result = bioavailability(**test_case["kwargs"])

    assert result[0] == test_case["expected_label"]
    assert result[1] == test_case["expected_value"]
    assert result[2] == test_case["expected_unit"]
