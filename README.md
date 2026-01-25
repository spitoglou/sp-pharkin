# sp-pharkin

A Python library for basic pharmacokinetics calculations with unit-aware arithmetic.

Implements formulas from *"Pharmacokinetics"* by Philip Rowe, with practice questions from the book converted into tests.

## Features

- **Unit-aware calculations** using [Pint](https://pint.readthedocs.io/) - pass quantities as strings like `'10 mg/L'`
- **Flexible equation solving** - given any two variables, compute the third
- **Consistent output format** - all functions return a standardized 5-tuple
- **Symbolic math support** - exponential decay solver using SymPy

## Installation

```bash
# Clone the repository
git clone https://github.com/spitoglou/sp-pharkin.git
cd sp-pharkin

# Install with uv
uv sync
```

## Quick Start

```python
from sp_pharkin import half_life_k, dose_concentration_volume
from sp_pharkin.clearance import clearance_flow_extraction_rate

# Calculate half-life from elimination rate constant
result = half_life_k(K='0.1/hour')
print(result[3])  # '6.93 hour'

# Calculate dose from concentration and volume
result = dose_concentration_volume(
    concentration='200 ug/L',
    volume='50 L',
    output_unit='mg'
)
print(result[3])  # '10.0 milligram'

# Calculate clearance from flow and extraction rate
result = clearance_flow_extraction_rate(
    Q='1.2 L/min',
    E=0.7,
    output_unit='L/hour'
)
print(result[3])  # '50.4 liter / hour'
```

## Output Format

All functions return a 5-tuple:

```python
(name, magnitude, unit_string, formatted_string, pint_quantity)
```

| Index | Name | Description |
|-------|------|-------------|
| 0 | `name` | Human-readable name of the result |
| 1 | `magnitude` | Numeric value |
| 2 | `unit_string` | Unit as a string |
| 3 | `formatted_string` | Complete formatted result |
| 4 | `pint_quantity` | Pint Quantity object for chaining calculations |

## Available Functions

### Core Functions (`sp_pharkin.functions`)

| Function | Equation | Description |
|----------|----------|-------------|
| `volume_of_distribution_weight()` | VD = Mean_VD × Weight | Volume of distribution from mean values and body weight |
| `dose_concentration_volume()` | Dose = C × V | Solve for dose, concentration, or volume |
| `target_concentration()` | (Min + Max) / 2 | Therapeutic target as midpoint of range |
| `rate_of_elimination_mass_k()` | Rate = Mass × K | Elimination rate calculations |
| `half_life_k()` | t½ = ln(2) / K | Half-life and elimination constant |
| `extraction_rate()` | ΔC = E × C_in | Extraction ratio calculations |

### Clearance Functions (`sp_pharkin.clearance`)

| Function | Equation | Description |
|----------|----------|-------------|
| `clearance_flow_extraction_rate()` | Cl = Q × E | Clearance from flow and extraction |
| `clearance_elimination_rate_constant_volume()` | Cl = K × V | Clearance from K and volume |
| `average_clearance_weight()` | Cl = Avg_Cl × Weight | Patient-specific clearance |

### Reduction Factors (`sp_pharkin.reduction_factors`)

| Function | Equation | Description |
|----------|----------|-------------|
| `salt_factor()` | Drug = Dose × SF | Adjust for salt formulations |
| `bioavailability()` | Drug = Dose × F | Adjust for absorption efficiency |

### Exponential Decay (`sp_pharkin.expo`)

| Function | Equation | Description |
|----------|----------|-------------|
| `c_t()` | C(t) = C₀ × e^(-kt) | Symbolic solver for decay equations |

## Optional Parameters

All functions accept these optional keyword arguments:

- `output_unit` - Convert result to specified unit (e.g., `'mg'`, `'L/hour'`)
- `decimals` - Round result to N decimal places (default: 2)

```python
result = half_life_k(
    K='0.1/hour',
    output_unit='minute',
    decimals=1
)
print(result[3])  # '415.9 minute'
```

## Chaining Calculations

Use the Pint Quantity (index 4) to chain calculations:

```python
from sp_pharkin import volume_of_distribution_weight, dose_concentration_volume

# Step 1: Calculate volume of distribution
vd = volume_of_distribution_weight(
    mean_volume_of_distribution='0.72 L/kg',
    weight='65 kg'
)[4]  # Extract pint_quantity

# Step 2: Use it to calculate required dose
dose = dose_concentration_volume(
    concentration='200 ug/L',
    volume=vd,
    output_unit='mg'
)
print(dose[3])  # '9.36 milligram'
```

## Development

```bash
# Install dependencies
uv sync

# Run all tests
uv run pytest

# Run specific test file
uv run pytest tests/test_chapter4.py

# Run specific test
uv run pytest tests/test_chapter4.py::test_clearance

# Lint code
uv run ruff check .

# Format code
uv run ruff format .
```

## Dependencies

- **pint** - Unit handling and dimensional analysis
- **sympy** - Symbolic mathematics for exponential equations

## License

MIT
