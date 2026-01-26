"""
Parametrized core pharmacokinetics tests.

Consolidates repetitive test scenarios from multiple test files using pytest parametrize.
Tests salt factor, bioavailability, clearance, and volume calculations with realistic
pharmaceutical and pharmacokinetic parameters.
"""

import pytest


class TestSaltFactorParametrized:
    """Parametrized tests for salt factor calculations."""

    @pytest.mark.parametrize(
        "kwargs,expected_label,expected_value,expected_unit",
        [
            # Test: Calculate Dose of Salt from Delivered Drug and Salt Factor
            (
                {"delivered_drug": "400 mg", "salt_factor": 0.86},
                "Dose of Salt",
                465,
                "milligram",
            ),
            (
                {"delivered_drug": "250 mg", "salt_factor": 0.88},
                "Dose of Salt",
                284,
                "milligram",
            ),
            (
                {"delivered_drug": "500 mg", "salt_factor": 0.90},
                "Dose of Salt",
                556,
                "milligram",
            ),
            # Test: Calculate Salt Factor from Delivered Drug and Dose of Salt
            (
                {"delivered_drug": "400 mg", "dose_of_salt": "465mg"},
                "Salt Factor",
                0.86,
                "dimensionless",
            ),
            (
                {"delivered_drug": "250 mg", "dose_of_salt": "284mg"},
                "Salt Factor",
                0.88,
                "dimensionless",
            ),
            (
                {"delivered_drug": "500 mg", "dose_of_salt": "556mg"},
                "Salt Factor",
                0.90,
                "dimensionless",
            ),
            # Test: Calculate Delivered Drug from Salt Factor and Dose of Salt
            (
                {"salt_factor": 0.86, "dose_of_salt": "465mg"},
                "Delivered Drug",
                400,
                "milligram",
            ),
            (
                {"salt_factor": 0.88, "dose_of_salt": "284mg"},
                "Delivered Drug",
                250,
                "milligram",
            ),
            (
                {"salt_factor": 0.90, "dose_of_salt": "556mg"},
                "Delivered Drug",
                500,
                "milligram",
            ),
        ],
        ids=[
            "salt_factor_dose_0.86_400mg",
            "salt_factor_dose_0.88_250mg",
            "salt_factor_dose_0.90_500mg",
            "salt_factor_from_delivered_0.86",
            "salt_factor_from_delivered_0.88",
            "salt_factor_from_delivered_0.90",
            "delivered_from_salt_0.86",
            "delivered_from_salt_0.88",
            "delivered_from_salt_0.90",
        ],
    )
    def test_salt_factor(self, kwargs, expected_label, expected_value, expected_unit):
        """
        Parametrized test for salt_factor function with multiple scenarios.

        Tests three calculation modes:
        - Calculate Dose of Salt from Delivered Drug and Salt Factor
        - Calculate Salt Factor from Delivered Drug and Dose of Salt
        - Calculate Delivered Drug from Salt Factor and Dose of Salt
        """
        from sp_pharkin import salt_factor

        result = salt_factor(**kwargs)

        assert result[0] == expected_label
        assert result[1] == pytest.approx(
            expected_value, rel=0.01
        )  # Allow 1% rounding tolerance
        assert result[2] == expected_unit
        assert len(result) == 5


class TestBioavailabilityParametrized:
    """Parametrized tests for bioavailability calculations."""

    @pytest.mark.parametrize(
        "kwargs,expected_label,expected_value,expected_unit",
        [
            # Test: Calculate Dose Administered from Delivered Drug and Bioavailability
            (
                {"delivered_drug": "200 mg", "bioavailability": 0.5},
                "Dose Administered",
                400,
                "milligram",
            ),
            (
                {"delivered_drug": "300 mg", "bioavailability": 0.75},
                "Dose Administered",
                400,
                "milligram",
            ),
            (
                {"delivered_drug": "475 mg", "bioavailability": 0.95},
                "Dose Administered",
                500,
                "milligram",
            ),
            # Test: Calculate Bioavailability from Delivered Drug and Dose Administered
            (
                {"delivered_drug": "200 mg", "dose_administered": "400mg"},
                "Bioavailability",
                0.5,
                "dimensionless",
            ),
            (
                {"delivered_drug": "300 mg", "dose_administered": "400mg"},
                "Bioavailability",
                0.75,
                "dimensionless",
            ),
            (
                {"delivered_drug": "475 mg", "dose_administered": "500mg"},
                "Bioavailability",
                0.95,
                "dimensionless",
            ),
            # Test: Calculate Delivered Drug from Bioavailability and Dose Administered
            (
                {"bioavailability": 0.5, "dose_administered": "400mg"},
                "Delivered Drug",
                200,
                "milligram",
            ),
            (
                {"bioavailability": 0.75, "dose_administered": "400mg"},
                "Delivered Drug",
                300,
                "milligram",
            ),
            (
                {"bioavailability": 0.95, "dose_administered": "500mg"},
                "Delivered Drug",
                475,
                "milligram",
            ),
        ],
        ids=[
            "bioavail_dose_0.5_200mg",
            "bioavail_dose_0.75_300mg",
            "bioavail_dose_0.95_475mg",
            "bioavail_from_delivered_0.5",
            "bioavail_from_delivered_0.75",
            "bioavail_from_delivered_0.95",
            "delivered_from_bioavail_0.5",
            "delivered_from_bioavail_0.75",
            "delivered_from_bioavail_0.95",
        ],
    )
    def test_bioavailability(
        self, kwargs, expected_label, expected_value, expected_unit
    ):
        """
        Parametrized test for bioavailability function with multiple scenarios.

        Tests three calculation modes:
        - Calculate Dose Administered from Delivered Drug and Bioavailability
        - Calculate Bioavailability from Delivered Drug and Dose Administered
        - Calculate Delivered Drug from Bioavailability and Dose Administered
        """
        from sp_pharkin import bioavailability

        result = bioavailability(**kwargs)

        assert result[0] == expected_label
        assert result[1] == pytest.approx(
            expected_value, rel=0.01
        )  # Allow 1% rounding tolerance
        assert result[2] == expected_unit
        assert len(result) == 5


class TestClearanceParametrized:
    """Parametrized tests for clearance calculations using Q × E formula."""

    @pytest.mark.parametrize(
        "kwargs,expected_label,expected_value,expected_unit",
        [
            # Test: Calculate Clearance from Flow and Extraction Ratio (Cl = Q × E)
            ({"Q": "5L/min", "E": 0.4}, "Clearance", 2.0, "liter / minute"),
            ({"Q": "3L/min", "E": 0.6}, "Clearance", 1.8, "liter / minute"),
            ({"Q": "2L/min", "E": 0.5}, "Clearance", 1.0, "liter / minute"),
            # Test: Calculate Flow from Clearance and Extraction Ratio
            ({"clearance": "2.0L/min", "E": 0.4}, "Flow(Q)", 5.0, "liter / minute"),
            ({"clearance": "1.8L/min", "E": 0.6}, "Flow(Q)", 3.0, "liter / minute"),
            ({"clearance": "1.0L/min", "E": 0.5}, "Flow(Q)", 2.0, "liter / minute"),
            # Test: Calculate Extraction Ratio from Clearance and Flow
            (
                {"clearance": "2.0L/min", "Q": "5L/min"},
                "Extraction Rate(E)",
                0.4,
                "dimensionless",
            ),
            (
                {"clearance": "1.8L/min", "Q": "3L/min"},
                "Extraction Rate(E)",
                0.6,
                "dimensionless",
            ),
            (
                {"clearance": "1.0L/min", "Q": "2L/min"},
                "Extraction Rate(E)",
                0.5,
                "dimensionless",
            ),
        ],
        ids=[
            "clearance_Q_5_E_0.4",
            "clearance_Q_3_E_0.6",
            "clearance_Q_2_E_0.5",
            "flow_from_clearance_0.4",
            "flow_from_clearance_0.6",
            "flow_from_clearance_0.5",
            "extraction_ratio_from_Cl_0.4",
            "extraction_ratio_from_Cl_0.6",
            "extraction_ratio_from_Cl_0.5",
        ],
    )
    def test_clearance_flow_extraction_rate(
        self, kwargs, expected_label, expected_value, expected_unit
    ):
        """
        Parametrized test for clearance calculated from flow and extraction ratio (Cl = Q × E).

        Tests three calculation modes:
        - Calculate Clearance from Flow and Extraction Ratio
        - Calculate Flow from Clearance and Extraction Ratio
        - Calculate Extraction Ratio from Clearance and Flow
        """
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        result = clearance_flow_extraction_rate(**kwargs)

        assert result[0] == expected_label
        assert result[1] == pytest.approx(
            expected_value, rel=0.01
        )  # Allow 1% rounding tolerance
        assert result[2] == expected_unit
        assert len(result) == 5


class TestClearanceKVParametrized:
    """Parametrized tests for clearance calculations using K × V formula."""

    @pytest.mark.parametrize(
        "kwargs,expected_label,expected_value,expected_unit",
        [
            # Test: Calculate Clearance from Elimination Rate Constant and Volume (Cl = K × V)
            ({"K": "0.08/hour", "volume": "60L"}, "Clearance", 4.8, "liter / hour"),
            ({"K": "0.12/hour", "volume": "50L"}, "Clearance", 6.0, "liter / hour"),
            ({"K": "0.1/hour", "volume": "50L"}, "Clearance", 5.0, "liter / hour"),
            # Test: Calculate Elimination Rate Constant from Clearance and Volume
            (
                {"clearance": "4.8L/hour", "volume": "60L"},
                "Elimination Rate Constant(K)",
                0.08,
                "1 / hour",
            ),
            (
                {"clearance": "6.0L/hour", "volume": "50L"},
                "Elimination Rate Constant(K)",
                0.12,
                "1 / hour",
            ),
            (
                {"clearance": "5.0L/hour", "volume": "50L"},
                "Elimination Rate Constant(K)",
                0.1,
                "1 / hour",
            ),
            # Test: Calculate Volume from Clearance and Elimination Rate Constant
            ({"K": "0.08/hour", "clearance": "4.8L/hour"}, "Volume", 60.0, "liter"),
            ({"K": "0.12/hour", "clearance": "6.0L/hour"}, "Volume", 50.0, "liter"),
            ({"K": "0.1/hour", "clearance": "5.0L/hour"}, "Volume", 50.0, "liter"),
        ],
        ids=[
            "clearance_K_0.08_V_60",
            "clearance_K_0.12_V_50",
            "clearance_K_0.1_V_50",
            "elimination_K_from_Cl_0.08",
            "elimination_K_from_Cl_0.12",
            "elimination_K_from_Cl_0.1",
            "volume_from_Cl_0.08",
            "volume_from_Cl_0.12",
            "volume_from_Cl_0.1",
        ],
    )
    def test_clearance_elimination_rate_constant_volume(
        self, kwargs, expected_label, expected_value, expected_unit
    ):
        """
        Parametrized test for clearance calculated from elimination rate constant and volume (Cl = K × V).

        Tests three calculation modes:
        - Calculate Clearance from Elimination Rate Constant and Volume
        - Calculate Elimination Rate Constant from Clearance and Volume
        - Calculate Volume from Clearance and Elimination Rate Constant
        """
        from sp_pharkin.clearance import clearance_elimination_rate_constant_volume

        result = clearance_elimination_rate_constant_volume(**kwargs)

        assert result[0] == expected_label
        assert result[1] == pytest.approx(
            expected_value, rel=0.01
        )  # Allow 1% rounding tolerance
        assert result[2] == expected_unit
        assert len(result) == 5


class TestExtractionRatioParametrized:
    """Parametrized tests for extraction ratio calculations."""

    @pytest.mark.parametrize(
        "kwargs,expected_label,expected_value,expected_unit",
        [
            # Test: Calculate Extraction Ratio from Input and Difference Concentrations
            ({"c_in": 10, "c_diff": 6}, "Extraction Ratio(E)", 0.6, "dimensionless"),
            ({"c_in": 20, "c_diff": 12}, "Extraction Ratio(E)", 0.6, "dimensionless"),
            ({"c_in": 15, "c_diff": 7.5}, "Extraction Ratio(E)", 0.5, "dimensionless"),
            # Test: Calculate Input Concentration from Difference and Extraction Ratio
            (
                {"c_diff": "6 mg/L", "E": 0.6},
                "Input Concentration(C_in)",
                10.0,
                "milligram / liter",
            ),
            (
                {"c_diff": "12 mg/L", "E": 0.6},
                "Input Concentration(C_in)",
                20.0,
                "milligram / liter",
            ),
            (
                {"c_diff": "7.5 mg/L", "E": 0.5},
                "Input Concentration(C_in)",
                15.0,
                "milligram / liter",
            ),
        ],
        ids=[
            "extraction_ratio_c_in_10_c_diff_6",
            "extraction_ratio_c_in_20_c_diff_12",
            "extraction_ratio_c_in_15_c_diff_7.5",
            "input_concentration_from_E_0.6",
            "input_concentration_from_E_0.6_2",
            "input_concentration_from_E_0.5",
        ],
    )
    def test_extraction_ratio(
        self, kwargs, expected_label, expected_value, expected_unit
    ):
        """
        Parametrized test for extraction ratio calculations.

        Tests two calculation modes:
        - Calculate Extraction Ratio from Input and Difference Concentrations
        - Calculate Input Concentration from Difference Concentration and Extraction Ratio
        """
        from sp_pharkin import extraction_rate

        result = extraction_rate(**kwargs)

        assert result[0] == expected_label
        assert result[1] == pytest.approx(
            expected_value, rel=0.01
        )  # Allow 1% rounding tolerance
        assert result[2] == expected_unit
        assert len(result) == 5
