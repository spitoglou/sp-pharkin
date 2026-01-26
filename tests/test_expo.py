"""
Comprehensive tests for exponential decay solver module

Tests cover solving for:
- C(t): Concentration at time t
- C₀: Initial concentration
- k: Elimination rate constant
- t: Time elapsed

With various unit combinations and validation of return formats.
"""

import pytest
import math


class TestSolveForCt:
    """Test solving for concentration at time t: C(t) = C₀·e^(-kt)"""

    def test_c_t_basic_mg_per_L(self):
        """Test basic C(t) calculation with mg/L units"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour", t="2 hour")

        # C(t) = 100 * e^(-0.5 * 2) = 100 * e^(-1) ≈ 36.79
        assert result[0] == "Concentration at Time t (C(t))"
        assert abs(result[1] - 36.79) < 0.01
        assert result[2] == "milligram / liter"
        assert isinstance(result[3], str)
        assert result[4] is not None

    def test_c_t_zero_time(self):
        """Test C(t) when t=0 should equal C₀"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.3 1/hour", t="0 hour")

        # C(t) = 100 * e^(-0.3 * 0) = 100 * e^(0) = 100
        assert abs(result[1] - 100.0) < 0.01

    def test_c_t_micrograms_per_mL(self):
        """Test C(t) with ug/mL units"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="50 ug/mL", k="0.1 1/hour", t="5 hour")

        # C(t) = 50 * e^(-0.1 * 5) = 50 * e^(-0.5) ≈ 30.33
        assert round(result[1], 2) == 30.33
        assert result[2] == "microgram / milliliter"

    def test_c_t_with_output_unit_conversion(self):
        """Test C(t) with output unit conversion (mg/L to ug/mL)"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(
            c_0="100 mg/L", k="0.5 1/hour", t="2 hour", output_unit="ug/mL"
        )

        # C(t) ≈ 36.79 mg/L = 36.79 ug/mL
        assert round(result[1], 2) == 36.79
        assert result[2] == "microgram / milliliter"

    def test_c_t_with_decimals_precision(self):
        """Test C(t) with custom decimal precision"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour", t="2 hour", decimals=4)

        assert result[1] == round(36.7879, 4)

    def test_c_t_missing_parameters(self):
        """Test C(t) raises error with missing parameters"""
        from sp_pharkin import solve_for_c_t

        with pytest.raises(ValueError, match="c_0, k, and t are required"):
            solve_for_c_t(
                c_0="100 mg/L",
                k="0.5 1/hour",
                # missing t
            )

    def test_c_t_return_tuple_format(self):
        """Test C(t) returns proper 5-tuple format"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour", t="2 hour")

        assert isinstance(result, tuple)
        assert len(result) == 5
        assert isinstance(result[0], str)  # name
        assert isinstance(result[1], (int, float))  # magnitude
        assert isinstance(result[2], str)  # units
        assert isinstance(result[3], str)  # formatted string
        assert result[4] is not None  # quantity


class TestSolveForC0:
    """Test solving for initial concentration: C₀ = C(t) / e^(-kt)"""

    def test_c_0_basic_mg_per_L(self):
        """Test basic C₀ calculation with mg/L units"""
        from sp_pharkin import solve_for_c_0

        result = solve_for_c_0(c_t="36.79 mg/L", k="0.5 1/hour", t="2 hour")

        # C₀ = 36.79 / e^(-0.5 * 2) = 36.79 / e^(-1) ≈ 100
        assert result[0] == "Initial Concentration (C₀)"
        assert abs(result[1] - 100.0) < 0.2
        assert result[2] == "milligram / liter"

    def test_c_0_from_c_t_calculation(self):
        """Test C₀ calculated from previous C(t) result"""
        from sp_pharkin import solve_for_c_0

        result = solve_for_c_0(c_t="30.33 ug/mL", k="0.1 1/hour", t="5 hour")

        # C₀ = 30.33 / e^(-0.1 * 5) ≈ 50
        assert abs(result[1] - 50.0) < 0.2

    def test_c_0_with_output_unit_conversion(self):
        """Test C₀ with output unit conversion"""
        from sp_pharkin import solve_for_c_0

        result = solve_for_c_0(
            c_t="50 mg/L", k="0.3 1/hour", t="3 hour", output_unit="mg/L"
        )

        # C₀ = 50 / e^(-0.3 * 3) ≈ 122.98
        assert abs(result[1] - 122.98) < 0.2
        assert result[2] == "milligram / liter"

    def test_c_0_with_decimals_precision(self):
        """Test C₀ with custom decimal precision"""
        from sp_pharkin import solve_for_c_0

        result = solve_for_c_0(c_t="36.79 mg/L", k="0.5 1/hour", t="2 hour", decimals=1)

        assert result[1] == 100.0

    def test_c_0_missing_parameters(self):
        """Test C₀ raises error with missing parameters"""
        from sp_pharkin import solve_for_c_0

        with pytest.raises(ValueError, match="c_t, k, and t are required"):
            solve_for_c_0(
                c_t="36.79 mg/L",
                k="0.5 1/hour",
                # missing t
            )

    def test_c_0_return_tuple_format(self):
        """Test C₀ returns proper 5-tuple format"""
        from sp_pharkin import solve_for_c_0

        result = solve_for_c_0(c_t="36.79 mg/L", k="0.5 1/hour", t="2 hour")

        assert isinstance(result, tuple)
        assert len(result) == 5
        assert result[0] == "Initial Concentration (C₀)"


class TestSolveForK:
    """Test solving for elimination rate constant: k = -ln(C(t) / C₀) / t"""

    def test_k_basic_calculation(self):
        """Test basic k calculation"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(c_0="100 mg/L", c_t="36.79 mg/L", t="2 hour")

        # k = -ln(36.79 / 100) / 2 = -ln(0.3679) / 2 ≈ 0.5
        assert result[0] == "Elimination Rate Constant (k)"
        assert round(result[1], 2) == 0.5
        assert result[2] == "1 / hour"

    def test_k_from_expo_values(self):
        """Test k calculated from known exponential decay"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(c_0="50 ug/mL", c_t="30.33 ug/mL", t="5 hour")

        # k = -ln(30.33 / 50) / 5 ≈ 0.1
        assert round(result[1], 2) == 0.1

    def test_k_half_life_relationship(self):
        """Test k calculation when C(t) = C₀/2 (half-life)"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(
            c_0="100 mg/L", c_t="50 mg/L", t="1 hour"  # half-life of 1 hour
        )

        # k = -ln(0.5) / 1 = ln(2) ≈ 0.693
        assert round(result[1], 2) == 0.69

    def test_k_with_output_unit_conversion(self):
        """Test k with different time unit"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(
            c_0="100 mg/L", c_t="36.79 mg/L", t="7200 second", output_unit="1/hour"
        )

        # k ≈ 0.5 1/hour (2 hours = 7200 seconds)
        assert round(result[1], 2) == 0.5
        assert result[2] == "1 / hour"

    def test_k_missing_parameters(self):
        """Test k raises error with missing parameters"""
        from sp_pharkin import solve_for_k

        with pytest.raises(ValueError, match="c_0, c_t, and t are required"):
            solve_for_k(
                c_0="100 mg/L",
                c_t="36.79 mg/L",
                # missing t
            )

    def test_k_return_tuple_format(self):
        """Test k returns proper 5-tuple format"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(c_0="100 mg/L", c_t="36.79 mg/L", t="2 hour")

        assert isinstance(result, tuple)
        assert len(result) == 5
        assert result[0] == "Elimination Rate Constant (k)"


class TestSolveForT:
    """Test solving for time elapsed: t = -ln(C(t) / C₀) / k"""

    def test_t_basic_calculation(self):
        """Test basic t calculation"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(c_0="100 mg/L", c_t="36.79 mg/L", k="0.5 1/hour")

        # t = -ln(36.79 / 100) / 0.5 ≈ 2
        assert result[0] == "Time Elapsed (t)"
        assert round(result[1], 2) == 2.0
        assert result[2] == "hour"

    def test_t_from_expo_values(self):
        """Test t calculated from known exponential decay"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(c_0="50 ug/mL", c_t="30.33 ug/mL", k="0.1 1/hour")

        # t = -ln(30.33 / 50) / 0.1 ≈ 5
        assert round(result[1], 2) == 5.0

    def test_t_half_life_calculation(self):
        """Test t calculation for half-life (C(t) = C₀/2)"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(c_0="100 mg/L", c_t="50 mg/L", k="0.693 1/hour")  # ln(2)

        # t = -ln(0.5) / 0.693 ≈ 1
        assert round(result[1], 2) == 1.0

    def test_t_with_output_unit_conversion(self):
        """Test t with output unit conversion to minutes"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(
            c_0="100 mg/L", c_t="36.79 mg/L", k="0.5 1/hour", output_unit="minute"
        )

        # t ≈ 2 hours = 120 minutes
        assert abs(result[1] - 120.0) < 0.1
        assert result[2] == "minute"

    def test_t_missing_parameters(self):
        """Test t raises error with missing parameters"""
        from sp_pharkin import solve_for_t

        with pytest.raises(ValueError, match="c_0, c_t, and k are required"):
            solve_for_t(
                c_0="100 mg/L",
                c_t="36.79 mg/L",
                # missing k
            )

    def test_t_return_tuple_format(self):
        """Test t returns proper 5-tuple format"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(c_0="100 mg/L", c_t="36.79 mg/L", k="0.5 1/hour")

        assert isinstance(result, tuple)
        assert len(result) == 5
        assert result[0] == "Time Elapsed (t)"


class TestExpoIntegration:
    """Integration tests combining multiple expo functions"""

    def test_round_trip_c_0_to_c_t(self):
        """Test calculating C(t) then back to C₀"""
        from sp_pharkin import solve_for_c_t, solve_for_c_0

        # Calculate C(t)
        result1 = solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour", t="2 hour")
        c_t_value = f"{result1[1]} {result1[2]}"

        # Calculate C₀ from C(t)
        result2 = solve_for_c_0(c_t=c_t_value, k="0.5 1/hour", t="2 hour")

        # Should return to ~100
        assert round(result2[1], 1) == 100.0

    def test_round_trip_k_calculation(self):
        """Test calculating k from C₀, C(t), and t"""
        from sp_pharkin import solve_for_k

        result = solve_for_k(c_0="100 mg/L", c_t="36.79 mg/L", t="2 hour")

        # k should be ~0.5
        assert round(result[1], 1) == 0.5

    def test_round_trip_t_calculation(self):
        """Test calculating t from C₀, C(t), and k"""
        from sp_pharkin import solve_for_t

        result = solve_for_t(c_0="100 mg/L", c_t="36.79 mg/L", k="0.5 1/hour")

        # t should be ~2
        assert round(result[1], 1) == 2.0

    def test_clinical_scenario_drug_elimination(self):
        """Test realistic clinical scenario: drug elimination from body"""
        from sp_pharkin import solve_for_c_t, solve_for_t

        # Initial dose: 500 mg/L
        # Elimination constant: 0.2 per hour (typical half-life ~3.5 hours)
        # Find concentration after 8 hours
        result_c_t = solve_for_c_t(c_0="500 mg/L", k="0.2 1/hour", t="8 hour")

        # C(t) = 500 * e^(-0.2 * 8) ≈ 100.95
        assert abs(result_c_t[1] - 100.95) < 0.5

        # When does concentration drop to 100 mg/L?
        result_t = solve_for_t(c_0="500 mg/L", c_t="100 mg/L", k="0.2 1/hour")

        # t = -ln(100/500) / 0.2 ≈ 8.05 hours
        assert abs(result_t[1] - 8.05) < 0.1

    def test_multiple_time_points(self):
        """Test calculating concentration at multiple time points"""
        from sp_pharkin import solve_for_c_t

        c_0_value = "100 mg/L"
        k_value = "0.3 1/hour"

        results = []
        for t_hours in [1, 2, 3, 4, 5]:
            result = solve_for_c_t(c_0=c_0_value, k=k_value, t=f"{t_hours} hour")
            results.append(result[1])

        # Verify decreasing concentration over time
        assert results[0] > results[1] > results[2] > results[3] > results[4]


class TestExpoEdgeCases:
    """Test edge cases and error conditions"""

    def test_c_t_zero_concentration(self):
        """Test C(t) approaches zero over time"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="1 1/hour", t="10 hour")

        # C(t) = 100 * e^(-10) ≈ 0.0045
        assert result[1] < 1.0

    def test_very_small_k(self):
        """Test with very small elimination constant"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.01 1/hour", t="5 hour")

        # C(t) = 100 * e^(-0.01 * 5) ≈ 95.12
        assert round(result[1], 2) == 95.12

    def test_large_k(self):
        """Test with large elimination constant"""
        from sp_pharkin import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="5 1/hour", t="2 hour")

        # C(t) = 100 * e^(-5 * 2) ≈ 0.0000454
        assert result[1] < 0.001

    def test_numeric_string_input(self):
        """Test with numeric strings instead of Quantity strings"""
        from sp_pharkin import solve_for_k

        # This should work because Q_() handles numeric conversion
        result = solve_for_k(c_0="100 mg/L", c_t="50 mg/L", t="1 hour")

        # k = -ln(0.5) / 1 ≈ 0.693
        assert round(result[1], 2) == 0.69
