"""
Reusable test fixtures for pharmacokinetics tests.

Provides common test data fixtures for doses, concentrations, units, extraction ratios,
and bioavailability scenarios for use across test modules.
"""

import pytest


# ============================================================================
# Dose Fixtures
# ============================================================================


@pytest.fixture
def standard_dose():
    """Standard pharmaceutical dose of 500 mg."""
    return "500 mg"


@pytest.fixture
def high_dose():
    """High pharmaceutical dose of 1000 mg."""
    return "1000 mg"


@pytest.fixture
def low_dose():
    """Low pharmaceutical dose of 100 mg."""
    return "100 mg"


@pytest.fixture
def pediatric_dose():
    """Pediatric pharmaceutical dose of 50 mg."""
    return "50 mg"


@pytest.fixture
def dose_range():
    """Range of common pharmaceutical doses."""
    return {
        "low": "100 mg",
        "standard": "500 mg",
        "high": "1000 mg",
        "pediatric": "50 mg",
        "geriatric": "250 mg",
    }


# ============================================================================
# Concentration Fixtures
# ============================================================================


@pytest.fixture
def standard_concentration():
    """Standard plasma concentration of 10 mg/L."""
    return "10 mg/L"


@pytest.fixture
def high_concentration():
    """High plasma concentration of 20 mg/L."""
    return "20 mg/L"


@pytest.fixture
def low_concentration():
    """Low plasma concentration of 5 mg/L."""
    return "5 mg/L"


@pytest.fixture
def concentration_range():
    """Range of common plasma concentrations."""
    return {
        "low": "5 mg/L",
        "standard": "10 mg/L",
        "high": "20 mg/L",
        "therapeutic": "8-12 mg/L",
        "toxic": "25 mg/L",
    }


# ============================================================================
# Unit Fixtures
# ============================================================================


@pytest.fixture
def mass_units():
    """Common mass units in pharmacokinetics."""
    return {
        "milligram": "mg",
        "microgram": "mcg",
        "gram": "g",
        "kilogram": "kg",
    }


@pytest.fixture
def volume_units():
    """Common volume units in pharmacokinetics."""
    return {
        "liter": "L",
        "milliliter": "mL",
        "microliter": "mcL",
    }


@pytest.fixture
def time_units():
    """Common time units in pharmacokinetics."""
    return {
        "hour": "hour",
        "minute": "minute",
        "second": "second",
        "day": "day",
    }


@pytest.fixture
def combined_units():
    """Common combined units in pharmacokinetics."""
    return {
        "mass_per_volume": "mg/L",
        "volume_per_time": "L/hour",
        "mass_per_time": "mg/hour",
        "concentration_units": "mg/L",
    }


# ============================================================================
# Extraction Ratio Fixtures
# ============================================================================


@pytest.fixture
def low_extraction_ratio():
    """Low extraction ratio (0.3)."""
    return 0.3


@pytest.fixture
def moderate_extraction_ratio():
    """Moderate extraction ratio (0.5)."""
    return 0.5


@pytest.fixture
def high_extraction_ratio():
    """High extraction ratio (0.7)."""
    return 0.7


@pytest.fixture
def very_high_extraction_ratio():
    """Very high extraction ratio (0.9)."""
    return 0.9


@pytest.fixture
def extraction_ratio_range():
    """Range of extraction ratios for different scenarios."""
    return {
        "low": 0.2,
        "low_moderate": 0.3,
        "moderate": 0.5,
        "high_moderate": 0.7,
        "high": 0.85,
        "very_high": 0.95,
    }


# ============================================================================
# Salt Factor Fixtures
# ============================================================================


@pytest.fixture
def low_salt_factor():
    """Low salt factor (0.86)."""
    return 0.86


@pytest.fixture
def moderate_salt_factor():
    """Moderate salt factor (0.88)."""
    return 0.88


@pytest.fixture
def high_salt_factor():
    """High salt factor (0.90)."""
    return 0.90


@pytest.fixture
def salt_factor_range():
    """Range of typical pharmaceutical salt factors."""
    return {
        "low": 0.80,
        "moderate_low": 0.86,
        "moderate": 0.88,
        "moderate_high": 0.90,
        "high": 0.95,
    }


# ============================================================================
# Bioavailability Fixtures
# ============================================================================


@pytest.fixture
def poor_bioavailability():
    """Poor bioavailability (0.5 or 50%)."""
    return 0.5


@pytest.fixture
def moderate_bioavailability():
    """Moderate bioavailability (0.75 or 75%)."""
    return 0.75


@pytest.fixture
def high_bioavailability():
    """High bioavailability (0.95 or 95%)."""
    return 0.95


@pytest.fixture
def complete_bioavailability():
    """Complete bioavailability (1.0 or 100%)."""
    return 1.0


@pytest.fixture
def bioavailability_scenarios():
    """Range of bioavailability scenarios."""
    return {
        "poor": {
            "value": 0.5,
            "description": "Poor bioavailability, significant first-pass metabolism",
            "example": "Nitroglycerin (sublingual)",
        },
        "moderate_low": {
            "value": 0.6,
            "description": "Moderate-low bioavailability",
            "example": "Aspirin",
        },
        "moderate": {
            "value": 0.75,
            "description": "Moderate bioavailability",
            "example": "Acetaminophen",
        },
        "moderate_high": {
            "value": 0.85,
            "description": "Moderate-high bioavailability",
            "example": "Ibuprofen",
        },
        "high": {
            "value": 0.95,
            "description": "High bioavailability",
            "example": "Most IV drugs",
        },
        "complete": {
            "value": 1.0,
            "description": "Complete bioavailability (IV administration)",
            "example": "Intravenous drugs",
        },
    }


# ============================================================================
# Pharmacokinetic Parameter Fixtures
# ============================================================================


@pytest.fixture
def elimination_rate_constant_values():
    """Common elimination rate constant values."""
    return {
        "slow": 0.05,
        "slow_moderate": 0.08,
        "moderate": 0.1,
        "moderate_fast": 0.15,
        "fast": 0.2,
    }


@pytest.fixture
def volume_of_distribution_values():
    """Common volume of distribution values in liters."""
    return {
        "small": 10,
        "moderate_small": 30,
        "moderate": 50,
        "moderate_large": 70,
        "large": 100,
    }


@pytest.fixture
def blood_flow_rates():
    """Common blood flow rates in L/min."""
    return {
        "low": 1.0,
        "moderate_low": 2.0,
        "moderate": 3.0,
        "moderate_high": 5.0,
        "high": 6.0,
    }


@pytest.fixture
def renal_clearance_values():
    """Common renal clearance values in mL/min."""
    return {
        "low": 30,
        "low_moderate": 60,
        "moderate": 90,
        "moderate_high": 120,
        "high": 150,
    }


# ============================================================================
# Complex Scenario Fixtures
# ============================================================================


@pytest.fixture
def drug_administration_scenarios():
    """Complex drug administration scenarios with multiple parameters."""
    return {
        "oral_poor_bioavail": {
            "dose_administered": "500 mg",
            "bioavailability": 0.5,
            "expected_delivered": 250,
            "description": "Oral administration with poor bioavailability",
        },
        "oral_moderate_bioavail": {
            "dose_administered": "500 mg",
            "bioavailability": 0.75,
            "expected_delivered": 375,
            "description": "Oral administration with moderate bioavailability",
        },
        "oral_high_bioavail": {
            "dose_administered": "500 mg",
            "bioavailability": 0.95,
            "expected_delivered": 475,
            "description": "Oral administration with high bioavailability",
        },
        "salt_form_0.86": {
            "dose_of_salt": "500 mg",
            "salt_factor": 0.86,
            "expected_delivered": 430,
            "description": "Salt form with factor 0.86",
        },
        "salt_form_0.90": {
            "dose_of_salt": "500 mg",
            "salt_factor": 0.90,
            "expected_delivered": 450,
            "description": "Salt form with factor 0.90",
        },
    }


@pytest.fixture
def clearance_calculation_scenarios():
    """Complex clearance calculation scenarios."""
    return {
        "flow_extraction": {
            "flow": "5 L/min",
            "extraction_ratio": 0.4,
            "expected_clearance": 2.0,
            "description": "Clearance from flow and extraction ratio",
        },
        "elimination_volume": {
            "elimination_rate_constant": "0.1 /hour",
            "volume": "50 L",
            "expected_clearance": 5.0,
            "description": "Clearance from elimination rate constant and volume",
        },
        "low_clearance": {
            "flow": "3 L/min",
            "extraction_ratio": 0.2,
            "expected_clearance": 0.6,
            "description": "Low clearance drug",
        },
        "high_clearance": {
            "flow": "5 L/min",
            "extraction_ratio": 0.9,
            "expected_clearance": 4.5,
            "description": "High clearance drug",
        },
    }


# ============================================================================
# Test Data Set Fixtures
# ============================================================================


@pytest.fixture
def realistic_drug_parameters():
    """Realistic drug parameters for pharmacokinetic calculations."""
    return {
        "paracetamol": {
            "bioavailability": 0.75,
            "vd": 50,
            "cl": 5.0,
            "ke": 0.1,
            "t_half": 2.0,
        },
        "ibuprofen": {
            "bioavailability": 0.80,
            "vd": 10,
            "cl": 0.04,
            "ke": 0.004,
            "t_half": 170,
        },
        "warfarin": {
            "bioavailability": 0.95,
            "vd": 9,
            "cl": 0.002,
            "ke": 0.0001,
            "t_half": 40,
        },
    }
