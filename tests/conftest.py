"""
Shared test fixtures and utilities for pharmacokinetics tests.

This module centralizes all pytest fixtures for pharmacokinetics tests,
including basic utilities and re-exports from test_fixtures.py.
"""

import pytest
from pint import UnitRegistry

# Re-export fixtures from test_fixtures module
from .test_fixtures import (
    # Dose fixtures
    standard_dose,
    high_dose,
    low_dose,
    pediatric_dose,
    dose_range,
    # Concentration fixtures
    standard_concentration,
    high_concentration,
    low_concentration,
    concentration_range,
    # Unit fixtures
    mass_units,
    volume_units,
    time_units,
    combined_units,
    # Extraction ratio fixtures
    low_extraction_ratio,
    moderate_extraction_ratio,
    high_extraction_ratio,
    very_high_extraction_ratio,
    extraction_ratio_range,
    # Salt factor fixtures
    low_salt_factor,
    moderate_salt_factor,
    high_salt_factor,
    salt_factor_range,
    # Bioavailability fixtures
    poor_bioavailability,
    moderate_bioavailability,
    high_bioavailability,
    complete_bioavailability,
    bioavailability_scenarios,
    # Pharmacokinetic parameter fixtures
    elimination_rate_constant_values,
    volume_of_distribution_values,
    blood_flow_rates,
    renal_clearance_values,
    # Complex scenario fixtures
    drug_administration_scenarios,
    clearance_calculation_scenarios,
    realistic_drug_parameters,
)


@pytest.fixture(scope="session")
def unit_registry():
    """Provide a unit registry for consistent unit handling across tests."""
    ureg = UnitRegistry()
    return ureg


@pytest.fixture
def Q_(unit_registry):
    """Provide Quantity constructor for unit-aware calculations."""
    return unit_registry.Quantity


class TupleValidator:
    """Helper class for validating 5-tuple return values from pharmacokinetics functions."""

    @staticmethod
    def validate_result(result, expected_label, expected_value, expected_unit):
        """
        Validate a pharmacokinetics function result tuple.

        Args:
            result: Tuple returned from pharmacokinetics function (label, value, unit, quantity, magnitude)
            expected_label: Expected result label/name
            expected_value: Expected numeric value
            expected_unit: Expected unit string

        Returns:
            bool: True if all assertions pass
        """
        assert len(result) == 5, f"Expected 5-tuple, got {len(result)}-tuple"
        assert (
            result[0] == expected_label
        ), f"Label mismatch: expected {expected_label}, got {result[0]}"
        assert (
            result[1] == expected_value
        ), f"Value mismatch: expected {expected_value}, got {result[1]}"
        assert (
            result[2] == expected_unit
        ), f"Unit mismatch: expected {expected_unit}, got {result[2]}"
        return True


@pytest.fixture
def validate_tuple():
    """Provide tuple validator for test assertions."""
    return TupleValidator.validate_result


# Test data fixtures for common pharmacokinetics scenarios


@pytest.fixture
def salt_factor_test_cases():
    """Parametrized test data for salt_factor function."""
    return [
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
    ]


@pytest.fixture
def bioavailability_test_cases():
    """Parametrized test data for bioavailability function."""
    return [
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
    ]
