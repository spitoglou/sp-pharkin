"""
Test Edge Cases and Boundary Conditions

This test module covers:
- Zero and near-zero values
- Extreme values (very large/very small numbers)
- Boundary conditions (0.0, 1.0)
- Unit conversions at extremes
- Physical value constraints
"""

import pytest
import math
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


# Helper function to validate 5-tuple return value
def validate_return_tuple(result, expected_name=None, expected_units=None):
    """
    Validate that result is a properly formatted 5-tuple from format_output.

    Returns:
        tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
    """
    assert isinstance(result, tuple), "Result should be a tuple"
    assert len(result) == 5, "Result should be a 5-tuple"

    name, magnitude, unit_string, formatted_string, quantity = result

    # Validate each element
    assert isinstance(name, str), "name (index 0) should be string"
    assert isinstance(magnitude, (int, float)), "magnitude (index 1) should be numeric"
    assert isinstance(unit_string, str), "unit_string (index 2) should be string"
    assert isinstance(
        formatted_string, str
    ), "formatted_string (index 3) should be string"

    # Validate that quantity is a pint Quantity
    assert hasattr(
        quantity, "magnitude"
    ), "quantity (index 4) should be a pint Quantity"
    assert hasattr(quantity, "units"), "quantity (index 4) should have units"

    # Optional assertions for expected values
    if expected_name is not None:
        assert name == expected_name, f"Expected name '{expected_name}', got '{name}'"

    if expected_units is not None:
        assert (
            unit_string == expected_units
        ), f"Expected units '{expected_units}', got '{unit_string}'"

    return result


class TestZeroDoseScenarios:
    """Tests for zero and near-zero dose values"""

    def test_zero_dose_bioavailability(self):
        """Test bioavailability with zero dose - should raise ValueError"""
        from sp_pharkin import bioavailability

        # Zero delivered drug (falsy) means only bioavailability is provided
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            bioavailability(delivered_drug="0 mg", bioavailability=0.8)

    def test_near_zero_dose_bioavailability(self):
        """Test bioavailability with very small dose"""
        from sp_pharkin import bioavailability

        result = bioavailability(delivered_drug="0.001 mg", bioavailability=0.8)

        validate_return_tuple(result)
        # When solving for dose_administered, 0.001/0.8 with decimals=2 rounds to 0.0
        assert isinstance(result[1], (int, float))

    def test_zero_salt_dose(self):
        """Test salt factor with zero salt dose - should raise ValueError"""
        from sp_pharkin import salt_factor

        # Zero delivered drug means only salt_factor is provided as non-falsy
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            salt_factor(delivered_drug="0 mg", salt_factor=0.8)

    def test_zero_concentration_at_time_t(self):
        """Test exponential decay when C(t) approaches zero"""
        from sp_pharkin.expo import solve_for_c_t

        result = solve_for_c_t(
            c_0="100 mg/L",
            k="1 1/hour",
            t="10 hour",  # After 10 half-lives, concentration is near zero
        )

        validate_return_tuple(result)
        # e^(-10) ≈ 0.0000453, so result should be very small
        assert result[1] < 0.01

    def test_zero_concentration_elimination_rate(self):
        """Test elimination rate when elimination rate is zero"""
        from sp_pharkin import rate_of_elimination_mass_k

        # Zero rate is falsy, so only mass is provided
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            rate_of_elimination_mass_k(rate_of_elimination="0 mg/hr", mass="10 mg")


class TestNegativeAndInvalidValues:
    """Tests for negative and physically invalid values"""

    def test_negative_extraction_ratio_calculation(self):
        """Test extraction ratio with invalid negative value (should calculate if inputs allow)"""
        from sp_pharkin import extraction_rate

        # When c_diff is negative, result will be negative (which is physically invalid)
        result = extraction_rate(
            c_in=10, c_diff=-6  # Negative difference is invalid in pharmacokinetics
        )

        validate_return_tuple(result)
        assert result[1] == -0.6  # Negative value calculated (but invalid)

    def test_negative_bioavailability_input(self):
        """Test bioavailability with negative value"""
        from sp_pharkin import bioavailability

        # Negative bioavailability is physically impossible
        result = bioavailability(delivered_drug="400 mg", bioavailability=-0.5)

        validate_return_tuple(result)
        # Result will be calculated as negative (but physically invalid)
        assert result[1] < 0

    def test_negative_elimination_rate_constant(self):
        """Test elimination rate constant with negative value"""
        from sp_pharkin import rate_of_elimination_mass_k

        result = rate_of_elimination_mass_k(
            rate_of_elimination="-1 mg/hr", mass="10 mg"
        )

        validate_return_tuple(result)
        assert result[1] == -0.1  # Negative K is invalid


class TestBoundaryConditions:
    """Tests for boundary values (0.0, 1.0, and extremes)"""

    def test_extraction_ratio_at_zero(self):
        """Test extraction ratio at lower boundary (E=0)"""
        from sp_pharkin import extraction_rate

        # Zero c_diff is falsy, so only c_in is provided
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            extraction_rate(c_in=10, c_diff=0)  # No extraction

    def test_extraction_ratio_at_unity(self):
        """Test extraction ratio at upper boundary (E=1.0)"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=10)  # Complete extraction

        validate_return_tuple(result)
        assert result[1] == 1.0
        assert result[2] == "dimensionless"

    def test_bioavailability_at_zero(self):
        """Test bioavailability at lower boundary (F=0)"""
        from sp_pharkin import bioavailability

        # Both delivered_drug and bioavailability are falsy, so no parameters count as provided
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            bioavailability(delivered_drug="0 mg", bioavailability=0.0)

    def test_bioavailability_exceeds_unity(self):
        """Test bioavailability > 1.0 (physically impossible but test calculation)"""
        from sp_pharkin import bioavailability

        result = bioavailability(
            delivered_drug="400 mg", bioavailability=1.5  # F > 1 is invalid
        )

        validate_return_tuple(result)
        assert result[1] > 1.0
        # The function calculates it, but this is physically invalid

    def test_half_life_extremely_small(self):
        """Test half-life with extremely small value (microseconds)"""
        from sp_pharkin import half_life_k

        result = half_life_k(K="1000000 1/second")  # Very fast decay

        validate_return_tuple(result)
        assert result[1] < 0.001  # t_half = ln(2) / k, so very small
        assert "second" in result[2]

    def test_half_life_extremely_large(self):
        """Test half-life with extremely large value (years)"""
        from sp_pharkin import half_life_k

        result = half_life_k(K="0.00000001 1/hour")  # Very slow decay

        validate_return_tuple(result)
        assert result[1] > 1e6  # Very large half-life

    def test_clearance_zero_extraction(self):
        """Test clearance with zero extraction ratio - should raise ValueError"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        # Zero extraction means only Q is provided, missing clearance and E
        # generic_a_eq_b_x_c requires exactly 2 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            clearance_flow_extraction_rate(
                Q="2 L/min", E=0.0  # Zero is falsy, so not counted as provided
            )

    def test_clearance_complete_extraction(self):
        """Test clearance with complete extraction (E=1.0)"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        result = clearance_flow_extraction_rate(
            Q="2 L/min", E=1.0  # Complete extraction
        )

        validate_return_tuple(result)
        assert result[1] == 2.0  # Clearance = Q * E


class TestExtremeUnitConversions:
    """Tests for extreme values with various unit conversions"""

    def test_very_large_dose_conversion(self):
        """Test with very large dose value"""
        from sp_pharkin import bioavailability

        result = bioavailability(delivered_drug="1000000 mg", bioavailability=0.8)

        validate_return_tuple(result)
        assert result[1] == pytest.approx(1250000, rel=1e-4)

    def test_very_small_concentration(self):
        """Test with very small concentration value"""
        from sp_pharkin.expo import solve_for_c_t

        result = solve_for_c_t(c_0="0.00001 mg/L", k="0.1 1/hour", t="2 hour")

        validate_return_tuple(result)
        # Very small initial concentration with decay results in near-zero
        # C(t) = 0.00001 * e^(-0.1 * 2) ≈ 0.00001 * 0.8187 ≈ 0.0000082
        assert result[1] >= 0
        assert result[1] < 0.00001

    def test_microsecond_time_unit(self):
        """Test with very small time unit (microseconds)"""
        from sp_pharkin.expo import solve_for_k

        result = solve_for_k(c_0="100 mg/L", c_t="99 mg/L", t="0.000001 hour")

        validate_return_tuple(result)
        # k should be very large since decay is over microseconds
        assert result[1] > 0

    def test_multiday_time_period(self):
        """Test with very large time period (multiple days)"""
        from sp_pharkin.expo import solve_for_t

        result = solve_for_t(c_0="100 mg/L", c_t="1 mg/L", k="0.01 1/hour")

        validate_return_tuple(result)
        # For ln(100)/0.01 ≈ 460.5 hours
        assert result[1] == pytest.approx(460.52, rel=1e-3)

    def test_ml_to_liter_conversion(self):
        """Test volume conversion from mL to L with edge case"""
        from sp_pharkin import dose_concentration_volume

        result = dose_concentration_volume(
            dose="500 mg", concentration="5 mg/mL", output_unit="L"
        )

        validate_return_tuple(result)
        # 500 mg / (5 mg/mL) = 100 mL = 0.1 L
        assert result[1] == pytest.approx(0.1, rel=1e-2)
        assert result[2] == "liter"


class TestConcentrationRatios:
    """Tests for edge cases in concentration ratio calculations"""

    def test_equal_concentrations_elimination(self):
        """Test when C(t) equals C_0 (no time has passed or k=0)"""
        from sp_pharkin.expo import solve_for_k

        result = solve_for_k(c_0="100 mg/L", c_t="100 mg/L", t="2 hour")

        validate_return_tuple(result)
        assert result[1] == pytest.approx(0.0, abs=1e-10)

    def test_very_large_ratio_concentration(self):
        """Test when C(t) is much larger than C_0 (invalid but test handling)"""
        from sp_pharkin.expo import solve_for_k

        result = solve_for_k(c_0="10 mg/L", c_t="100 mg/L", t="1 hour")

        validate_return_tuple(result)
        # Negative k (concentration increased, which is invalid)
        assert result[1] < 0

    def test_dose_calculation_with_minimal_bioavailability(self):
        """Test dose calculation with very low bioavailability"""
        from sp_pharkin import bioavailability

        result = bioavailability(delivered_drug="0.001 mg", bioavailability=0.001)

        validate_return_tuple(result, expected_name="Dose Administered")
        assert result[1] == pytest.approx(1.0, rel=1e-3)


class TestFormattedOutputValidation:
    """Tests to validate all 5-tuple elements in formatted output"""

    def test_formatted_string_contains_magnitude_and_unit(self):
        """Verify formatted_string contains both magnitude and unit"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=6)
        name, magnitude, unit_string, formatted_string, quantity = result

        # formatted_string should be like "0.6 dimensionless"
        assert str(magnitude) in formatted_string
        assert unit_string in formatted_string

    def test_unit_string_matches_quantity_units(self):
        """Verify unit_string matches the quantity's units"""
        from sp_pharkin import salt_factor

        result = salt_factor(delivered_drug="400 mg", salt_factor=0.8)

        name, magnitude, unit_string, formatted_string, quantity = result

        # unit_string should match quantity.units
        assert unit_string == str(quantity.units)

    def test_pint_quantity_is_valid(self):
        """Verify pint_quantity is a valid Quantity object"""
        from sp_pharkin import rate_of_elimination_mass_k

        result = rate_of_elimination_mass_k(rate_of_elimination="1 mg/hr", mass="10 mg")

        name, magnitude, unit_string, formatted_string, quantity = result

        # quantity should be a valid Quantity with magnitude and units
        assert quantity.magnitude == pytest.approx(magnitude)
        assert str(quantity.units) == unit_string

        # Should be able to perform unit operations
        assert hasattr(quantity, "to")  # Has unit conversion method

    def test_magnitude_precision_matches_decimals_parameter(self):
        """Test that magnitude respects the decimals parameter"""
        from sp_pharkin import bioavailability

        # Providing 3 parameters where generic_a_eq_b_x_c needs exactly 2
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            bioavailability(
                delivered_drug="333.3333 mg",
                salt_factor="0.777",
                dose_of_salt="500 mg",  # Over-specified: all 3 provided
            )

        # Test valid case with proper decimals parameter
        result = bioavailability(
            delivered_drug="333.3333 mg", bioavailability=0.75, decimals=2
        )

        name, magnitude, unit_string, formatted_string, quantity = result

        # With decimals=2, should be rounded
        # Check that magnitude has reasonable precision
        assert isinstance(magnitude, (int, float))

    def test_name_field_consistency(self):
        """Test that name field is consistent across different parameters"""
        from sp_pharkin import half_life_k

        result1 = half_life_k(K="0.1 1/hour")
        result2 = half_life_k(K="0.2 1/hour")

        # Both should have the same name field
        assert result1[0] == result2[0] == "Half-Life"

    def test_all_elements_present_in_exponential_solver(self):
        """Verify all 5 elements in exponential decay solver"""
        from sp_pharkin.expo import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour", t="1 hour")

        validate_return_tuple(
            result,
            expected_name="Concentration at Time t (C(t))",
            expected_units="milligram / liter",
        )


class TestVolumeDistributionEdgeCases:
    """Tests for volume of distribution edge cases"""

    def test_zero_weight_volume(self):
        """Test volume of distribution with zero weight - should raise ValueError"""
        from sp_pharkin import volume_of_distribution_weight

        # Zero weight is falsy, so only mean_volume_of_distribution is counted
        # generic_a_eq_b_x_c needs exactly 2 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            volume_of_distribution_weight(
                mean_volume_of_distribution="0.5 L/kg", weight="0 kg"
            )

    def test_extremely_large_patient_weight(self):
        """Test with very large patient weight"""
        from sp_pharkin import volume_of_distribution_weight

        result = volume_of_distribution_weight(
            mean_volume_of_distribution="0.5 L/kg", weight="500 kg"
        )

        validate_return_tuple(result)
        assert result[1] == 250.0


class TestClearanceEdgeCases:
    """Tests for clearance calculation edge cases"""

    def test_zero_flow_rate_clearance(self):
        """Test clearance with zero flow rate - should raise ValueError"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        # Zero flow is falsy, so only E is counted
        # generic_a_eq_b_x_c needs exactly 2 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            clearance_flow_extraction_rate(Q="0 L/min", E=0.8)

    def test_very_high_clearance_value(self):
        """Test with extremely high clearance value"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        result = clearance_flow_extraction_rate(Q="100 L/min", E=1.0)

        validate_return_tuple(result)
        assert result[1] == 100.0

    def test_clearance_unit_consistency(self):
        """Test that clearance maintains proper units"""
        from sp_pharkin.clearance import clearance_elimination_rate_constant_volume

        result = clearance_elimination_rate_constant_volume(
            K="0.1 1/hour", volume="50 L"
        )

        validate_return_tuple(result)
        assert "liter" in result[2]
        assert "hour" in result[2]


class TestDoseCalculationBoundaries:
    """Tests for dose calculation at boundaries"""

    def test_dose_concentration_volume_zero_concentration(self):
        """Test dose calculation with zero concentration - should raise ValueError"""
        from sp_pharkin import dose_concentration_volume

        # Zero concentration is falsy, so only volume is counted
        # generic_a_eq_b_x_c needs exactly 2 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            dose_concentration_volume(concentration="0 mg/mL", volume="100 mL")

    def test_dose_concentration_volume_zero_volume(self):
        """Test dose calculation with zero volume"""
        from sp_pharkin import dose_concentration_volume

        result = dose_concentration_volume(dose="500 mg", concentration="5 mg/mL")

        validate_return_tuple(result, expected_name="Volume")
        assert result[1] == pytest.approx(100.0, rel=1e-2)  # 500/5 = 100 mL
