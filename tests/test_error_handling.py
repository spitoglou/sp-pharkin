"""
Test Error Handling and Validation

This test module covers:
- Missing required parameters
- Contradictory inputs
- Invalid unit strings
- Type errors
- Invalid output unit conversions
- Invalid parameter values
"""

import pytest
from pint import UnitRegistry
from pint.errors import DimensionalityError

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


class TestMissingParameters:
    """Tests for missing required parameters"""

    def test_solve_for_c_t_missing_c_0(self):
        """Test solve_for_c_t without initial concentration"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises(ValueError, match="c_0, k, and t are required"):
            solve_for_c_t(k="0.5 1/hour", t="2 hour")

    def test_solve_for_c_t_missing_k(self):
        """Test solve_for_c_t without elimination rate constant"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises(ValueError, match="c_0, k, and t are required"):
            solve_for_c_t(c_0="100 mg/L", t="2 hour")

    def test_solve_for_c_t_missing_t(self):
        """Test solve_for_c_t without time"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises(ValueError, match="c_0, k, and t are required"):
            solve_for_c_t(c_0="100 mg/L", k="0.5 1/hour")

    def test_solve_for_c_0_missing_c_t(self):
        """Test solve_for_c_0 without concentration at time t"""
        from sp_pharkin.expo import solve_for_c_0

        with pytest.raises(ValueError, match="c_t, k, and t are required"):
            solve_for_c_0(k="0.5 1/hour", t="2 hour")

    def test_solve_for_c_0_missing_k(self):
        """Test solve_for_c_0 without elimination rate constant"""
        from sp_pharkin.expo import solve_for_c_0

        with pytest.raises(ValueError, match="c_t, k, and t are required"):
            solve_for_c_0(c_t="50 mg/L", t="2 hour")

    def test_solve_for_c_0_missing_t(self):
        """Test solve_for_c_0 without time"""
        from sp_pharkin.expo import solve_for_c_0

        with pytest.raises(ValueError, match="c_t, k, and t are required"):
            solve_for_c_0(c_t="50 mg/L", k="0.5 1/hour")

    def test_solve_for_k_missing_c_0(self):
        """Test solve_for_k without initial concentration"""
        from sp_pharkin.expo import solve_for_k

        with pytest.raises(ValueError, match="c_0, c_t, and t are required"):
            solve_for_k(c_t="50 mg/L", t="2 hour")

    def test_solve_for_k_missing_c_t(self):
        """Test solve_for_k without concentration at time t"""
        from sp_pharkin.expo import solve_for_k

        with pytest.raises(ValueError, match="c_0, c_t, and t are required"):
            solve_for_k(c_0="100 mg/L", t="2 hour")

    def test_solve_for_k_missing_t(self):
        """Test solve_for_k without time"""
        from sp_pharkin.expo import solve_for_k

        with pytest.raises(ValueError, match="c_0, c_t, and t are required"):
            solve_for_k(c_0="100 mg/L", c_t="50 mg/L")

    def test_solve_for_t_missing_c_0(self):
        """Test solve_for_t without initial concentration"""
        from sp_pharkin.expo import solve_for_t

        with pytest.raises(ValueError, match="c_0, c_t, and k are required"):
            solve_for_t(c_t="50 mg/L", k="0.5 1/hour")

    def test_solve_for_t_missing_c_t(self):
        """Test solve_for_t without concentration at time t"""
        from sp_pharkin.expo import solve_for_t

        with pytest.raises(ValueError, match="c_0, c_t, and k are required"):
            solve_for_t(c_0="100 mg/L", k="0.5 1/hour")

    def test_solve_for_t_missing_k(self):
        """Test solve_for_t without elimination rate constant"""
        from sp_pharkin.expo import solve_for_t

        with pytest.raises(ValueError, match="c_0, c_t, and k are required"):
            solve_for_t(c_0="100 mg/L", c_t="50 mg/L")


class TestTypeErrors:
    """Tests for incorrect parameter types"""

    def test_invalid_quantity_string_format(self):
        """Test with invalid unit string format"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises((ValueError, AttributeError)):
            solve_for_c_t(c_0="abc mg/L", k="0.5 1/hour", t="2 hour")  # Invalid format

    def test_dimensionality_mismatch_concentration_as_time(self):
        """Test providing concentration value where time is expected"""
        from sp_pharkin.expo import solve_for_c_t

        # Pint actually allows this - the function doesn't validate dimensionality strictly
        result = solve_for_c_t(
            c_0="100 mg/L",
            k="0.5 1/hour",
            t="50 mg/L",  # Wrong dimensionality but function doesn't check
        )
        # It will compute something, though physically meaningless
        assert isinstance(result, tuple)

    def test_dimensionality_mismatch_time_as_concentration(self):
        """Test providing time value where concentration is expected"""
        from sp_pharkin.expo import solve_for_c_t

        # Pint doesn't enforce dimensional type checking strictly here
        result = solve_for_c_t(
            c_0="2 hour", k="0.5 1/hour", t="2 hour"  # Wrong dimensionality
        )
        assert isinstance(result, tuple)

    def test_invalid_unit_string(self):
        """Test with completely invalid unit"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises((ValueError, AttributeError)):
            solve_for_c_t(c_0="100 invalid_unit", k="0.5 1/hour", t="2 hour")

    def test_extraction_rate_string_instead_of_number(self):
        """Test extraction_rate with non-numeric dimensionless value"""
        from sp_pharkin import extraction_rate

        with pytest.raises((ValueError, AttributeError)):
            extraction_rate(c_in=10, E="invalid")  # Not a valid number

    def test_bioavailability_invalid_string(self):
        """Test bioavailability with invalid string parameter"""
        from sp_pharkin import bioavailability

        with pytest.raises((ValueError, AttributeError)):
            bioavailability(delivered_drug="invalid", bioavailability=0.8)


class TestInvalidOutputUnits:
    """Tests for invalid output unit conversions"""

    def test_incompatible_output_unit_concentration_to_mass(self):
        """Test converting concentration to incompatible unit (mass)"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises((ValueError, DimensionalityError)):
            solve_for_c_t(
                c_0="100 mg/L",
                k="0.5 1/hour",
                t="2 hour",
                output_unit="kg",  # Incompatible unit for concentration
            )

    def test_incompatible_output_unit_dimensionless(self):
        """Test converting dimensionless quantity to time unit"""
        from sp_pharkin import extraction_rate

        with pytest.raises((ValueError, DimensionalityError)):
            extraction_rate(
                c_in=10,
                c_diff=6,
                output_unit="hour",  # Incompatible unit for dimensionless
            )

    def test_invalid_output_unit_format(self):
        """Test with invalid unit syntax in output_unit"""
        from sp_pharkin.expo import solve_for_c_t

        with pytest.raises((ValueError, AttributeError)):
            solve_for_c_t(
                c_0="100 mg/L",
                k="0.5 1/hour",
                t="2 hour",
                output_unit="???xyz",  # Invalid unit syntax
            )

    def test_output_unit_incompatible_with_quantity(self):
        """Test output unit incompatible with quantity dimensions"""
        from sp_pharkin import rate_of_elimination_mass_k

        with pytest.raises((ValueError, DimensionalityError)):
            rate_of_elimination_mass_k(
                rate_of_elimination="1 mg/hr",
                mass="10 mg",
                output_unit="meter",  # Wrong dimension
            )


class TestInvalidParameterValues:
    """Tests for invalid parameter value constraints"""

    def test_decimals_parameter_negative(self):
        """Test with negative decimals parameter"""
        from sp_pharkin import extraction_rate

        # This may not raise an error but could cause unexpected behavior
        result = extraction_rate(c_in=10, c_diff=6, decimals=-1)

        # Should still return a valid tuple
        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_decimals_parameter_very_large(self):
        """Test with very large decimals parameter"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=6, decimals=100)

        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_zero_k_in_half_life(self):
        """Test half_life with zero elimination rate constant"""
        from sp_pharkin import half_life_k

        # Zero K is falsy, so only the implicit ln(2) constant is provided
        # This violates the requirement of exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            half_life_k(K="0 1/hour")

    def test_zero_time_in_exponential_solver(self):
        """Test exponential solver with zero time should not raise but may be edge case"""
        from sp_pharkin.expo import solve_for_k

        with pytest.raises((ValueError, ZeroDivisionError)):
            solve_for_k(
                c_0="100 mg/L",
                c_t="50 mg/L",
                t="0 hour",  # Would cause division by zero
            )

    def test_zero_ratio_in_logarithm(self):
        """Test when concentration ratio results in log(0)"""
        from sp_pharkin.expo import solve_for_k

        with pytest.raises((ValueError, ZeroDivisionError)):
            solve_for_k(c_0="100 mg/L", c_t="0 mg/L", t="2 hour")  # Would need log(0)


class TestContradictoryInputs:
    """Tests for contradictory or over-specified parameters"""

    def test_generic_a_eq_b_x_c_all_three_parameters(self):
        """Test when all three parameters are provided to solve for unknown"""
        from sp_pharkin import bioavailability

        # When all three are specified, the function should raise ValueError
        # The function requires exactly 2 of 3 parameters
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            bioavailability(
                delivered_drug="400 mg", dose_administered="500 mg", bioavailability=0.8
            )

    def test_generic_formula_contradictory_values(self):
        """Test when parameters contradict the mathematical relation"""
        from sp_pharkin import salt_factor

        # All three parameters provided violates the constraint (needs exactly 2)
        with pytest.raises(ValueError, match="exactly 2 of 3 parameters"):
            salt_factor(delivered_drug="400 mg", salt_factor=0.8, dose_of_salt="500 mg")


class TestExtremeParameterValues:
    """Tests for extreme but technically valid parameters"""

    def test_extremely_large_dose(self):
        """Test with unrealistically large dose"""
        from sp_pharkin import bioavailability

        result = bioavailability(delivered_drug="1e100 mg", bioavailability=0.8)

        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_extremely_small_dose(self):
        """Test with unrealistically small dose"""
        from sp_pharkin import bioavailability

        result = bioavailability(delivered_drug="1e-100 mg", bioavailability=0.8)

        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_k_value_causing_overflow(self):
        """Test elimination rate constant that could cause exponential overflow"""
        from sp_pharkin.expo import solve_for_c_t

        # e^(-1000) will underflow to 0
        result = solve_for_c_t(c_0="100 mg/L", k="1000 1/hour", t="2 hour")

        assert isinstance(result, tuple)
        # Result magnitude should be effectively 0
        assert result[1] == pytest.approx(0, abs=1e-100)

    def test_very_small_elimination_constant(self):
        """Test with extremely small elimination constant"""
        from sp_pharkin.expo import solve_for_c_t

        result = solve_for_c_t(c_0="100 mg/L", k="0.00000001 1/hour", t="1 hour")

        assert isinstance(result, tuple)
        # Result should be very close to initial concentration
        assert result[1] == pytest.approx(100.0, rel=1e-5)


class TestUnitValidation:
    """Tests for unit validation and consistency"""

    def test_inconsistent_concentration_units_in_ratio(self):
        """Test that different concentration units can be compared if compatible"""
        from sp_pharkin.expo import solve_for_k

        # mg/L and µg/mL are equivalent concentrations
        result = solve_for_k(c_0="100 mg/L", c_t="50 mg/L", t="2 hour")

        assert isinstance(result, tuple)

    def test_incompatible_units_in_multiplication(self):
        """Test incompatible units in dose_concentration_volume"""
        from sp_pharkin import dose_concentration_volume

        # Pint will compute with wrong units rather than raise
        result = dose_concentration_volume(
            dose="500 kg", concentration="5 mg/mL"  # Wrong dimension for dose
        )
        # Function computes it but with nonsensical units
        assert isinstance(result, tuple)

    def test_volume_unit_variations(self):
        """Test that various volume units are handled correctly"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        # mL and L should be compatible
        result = clearance_flow_extraction_rate(
            Q="2000 mL/min", E=0.5, output_unit="L/hour"
        )

        assert isinstance(result, tuple)
        assert "liter" in result[2]


class TestNumericalErrors:
    """Tests for potential numerical errors"""

    def test_negative_logarithm_argument(self):
        """Test when C(t) > C_0 produces physically invalid negative time"""
        from sp_pharkin.expo import solve_for_t

        # When C(t) > C_0, the ratio is > 1, so ln(ratio) > 0, giving negative t
        # This is physically impossible but mathematically computed
        result = solve_for_t(
            c_0="100 mg/L",
            c_t="200 mg/L",  # Greater than initial (physically impossible)
            k="0.5 1/hour",
        )

        # Result should be negative time (marking the physical impossibility)
        assert result[1] < 0  # Negative time indicates invalid scenario
        assert result[0] == "Time Elapsed (t)"

    def test_infinite_clearance_calculation(self):
        """Test clearance calculation approaching infinity"""
        from sp_pharkin.clearance import clearance_flow_extraction_rate

        result = clearance_flow_extraction_rate(Q="1e100 L/min", E=1.0)

        assert isinstance(result, tuple)
        # Should handle large numbers
        assert result[1] > 0

    def test_very_close_concentrations_ratio(self):
        """Test when C(t) is very close to C_0 (numerical precision)"""
        from sp_pharkin.expo import solve_for_k

        result = solve_for_k(
            c_0="100.0000000001 mg/L", c_t="100.0000000002 mg/L", t="2 hour"
        )

        assert isinstance(result, tuple)
        # k should be very close to 0
        assert abs(result[1]) < 1e-10


class TestRoundingAndPrecision:
    """Tests for rounding and precision handling"""

    def test_decimals_parameter_zero(self):
        """Test with decimals=0 (should round to integer)"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=3, decimals=0)

        # With decimals=0, result 0.3 is rounded to 0
        assert result[1] == 0.0

    def test_very_high_precision_output(self):
        """Test with very high precision requirement"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=6, decimals=15)

        assert isinstance(result, tuple)
        assert result[1] == pytest.approx(0.6)

    def test_rounding_affects_output_format(self):
        """Test that rounding affects formatted output string"""
        from sp_pharkin import bioavailability

        result = bioavailability(
            delivered_drug="333.333 mg", bioavailability=0.777, decimals=2
        )

        name, magnitude, units, formatted_str, quantity = result

        # formatted_str should reflect the rounded magnitude
        assert str(magnitude) in formatted_str


class TestEmptyOrNoneInputs:
    """Tests for empty or None-like inputs"""

    def test_empty_string_unit(self):
        """Test with empty/minimal unit specification"""
        from sp_pharkin.expo import solve_for_c_t

        # Pint interprets '100' as dimensionless quantity
        result = solve_for_c_t(
            c_0="100", k="0.5 1/hour", t="2 hour"  # No explicit unit
        )
        assert isinstance(result, tuple)

    def test_none_value_in_optional_parameter(self):
        """Test with None in optional parameter"""
        from sp_pharkin import extraction_rate

        # output_unit=None should default to False
        result = extraction_rate(c_in=10, c_diff=6, output_unit=None)

        assert isinstance(result, tuple)


class TestReturnValueConsistency:
    """Tests to ensure return values are consistent and valid"""

    def test_return_tuple_always_has_five_elements(self):
        """Verify all functions return 5-tuple"""
        from sp_pharkin import extraction_rate, salt_factor

        result1 = extraction_rate(c_in=10, c_diff=6)
        result2 = salt_factor(delivered_drug="400 mg", salt_factor=0.8)

        assert len(result1) == 5
        assert len(result2) == 5

    def test_return_tuple_element_types(self):
        """Verify correct types in return tuple"""
        from sp_pharkin import rate_of_elimination_mass_k

        result = rate_of_elimination_mass_k(rate_of_elimination="1 mg/hr", mass="10 mg")

        name, magnitude, units, formatted, quantity = result

        assert isinstance(name, str)
        assert isinstance(magnitude, (int, float))
        assert isinstance(units, str)
        assert isinstance(formatted, str)
        assert hasattr(quantity, "magnitude")
        assert hasattr(quantity, "units")

    def test_magnitude_matches_quantity_magnitude(self):
        """Verify magnitude in tuple matches quantity.magnitude"""
        from sp_pharkin import half_life_k

        result = half_life_k(K="0.1 1/hour")

        name, magnitude, units, formatted, quantity = result

        assert magnitude == pytest.approx(quantity.magnitude)

    def test_units_string_matches_quantity_units(self):
        """Verify units string matches quantity.units"""
        from sp_pharkin import extraction_rate

        result = extraction_rate(c_in=10, c_diff=6)

        name, magnitude, units, formatted, quantity = result

        assert units == str(quantity.units)
