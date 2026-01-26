# Copilot Instructions for sp-pharkin

## Project Quick Reference

**sp-pharkin** is a pharmacokinetics calculation library implementing formulas from "Pharmacokinetics" by Philip Rowe. All functions perform unit-aware arithmetic and return standardized 5-tuples.

This project uses **UV** for fast, reliable dependency management with reproducible builds via `uv.lock`.

## Essential Developer Commands

All commands use **UV** (`uv run ...`). Never use plain `python` or `pip` - this ensures the correct locked environment.

```bash
# Install dependencies (creates .venv and syncs uv.lock)
uv sync

# Run all tests with coverage
uv run pytest tests/ --cov=sp_pharkin -q

# Run specific test file
uv run pytest tests/test_expo.py -v

# Format and lint
uv run ruff format .
uv run ruff check . --fix
```

**Why UV:** Provides reproducible builds, faster dependency resolution, and automatic virtual environment management.

## Core Architecture Pattern: `generic_a_eq_b_x_c`

**The fundamental pattern** solves equations of form `a = b × c` given any two variables.

```python
# In sp_pharkin/lib.py
def generic_a_eq_b_x_c(a=None, b=None, c=None, names=None):
    # Requires exactly 2 of 3 parameters to be non-falsy
    # Returns the missing parameter
```

**All pharmacokinetics functions** follow this pattern:
- `half_life_k(K=..., half_life=...)` → solves t½ = ln(2) / K
- `dose_concentration_volume(dose=..., concentration=..., volume=...)` → solves Dose = C × V
- `clearance_flow_extraction_rate(Q=..., E=..., clearance=...)` → solves Cl = Q × E

## Return Format (Critical)

**Every function returns this 5-tuple:**
```python
(name: str, magnitude: float, unit_string: str, formatted_string: str, pint_quantity: Quantity)
```

Example: `('Half-Life', 6.93, 'hour', '6.93 hour', Quantity(6.93, 'hour'))`

- Index `[0]`: Human-readable name
- Index `[1]`: Numeric value only
- Index `[2]`: Unit string
- Index `[3]`: Complete formatted string (what users see)
- Index `[4]`: Pint Quantity for chaining calculations

## Key Implementation Details

### Input Format
- Functions accept **string quantities with units**: `K='0.1/hour'`, `concentration='200 ug/L'`
- Parsed and validated by Pint UnitRegistry in `reduction_factors.py`
- Always use keyword arguments

### Optional Parameters
- `output_unit`: Convert result to different unit (e.g., `output_unit='mg'`)
- `decimals`: Rounding precision (default: 2)

### Dependencies
- **pint**: Unit handling - all quantities use `Q_()` wrapper
- **sympy**: Symbolic math for `expo.py` (exponential decay solver)

## Module Map

| Module | Purpose | Pattern |
|--------|---------|---------|
| `lib.py` | Core utilities | `generic_a_eq_b_x_c()`, `format_output()` |
| `functions.py` | Main calculations | VD, dose, concentration, extraction ratio |
| `clearance.py` | Clearance calculations | Flow×Extraction, K×Volume |
| `reduction_factors.py` | Salt factor, bioavailability | UnitRegistry setup |
| `expo.py` | Exponential decay solver | Symbolic SymPy equations |

## Dependency Management with UV

**UV replaces pip and poetry** for this project. Key files:
- `pyproject.toml` - Project metadata and dependencies
- `uv.lock` - Locked dependency versions (commit this file!)

**UV workflow:**
```bash
# Install: Creates .venv, installs all dependencies from uv.lock
uv sync

# Update lock file after changing pyproject.toml:
uv sync --all-extras

# Run commands in the UV environment:
uv run pytest
uv run ruff format .
```

Never modify `uv.lock` manually. Always use `uv sync` when dependencies change.

## Test Structure (100% Coverage)

- `test_expo.py`: 51 tests for exponential decay solver
- `test_edge_cases.py`: 41 tests for boundary conditions and zero/extreme values
- `test_error_handling.py`: 87 tests for invalid inputs and validation
- `test_parametrized_core.py`: 42 consolidated parametrized tests
- `test_fixtures.py`: Reusable test data fixtures
- `conftest.py`: Central fixture imports from `test_fixtures.py`

**Run tests:** `uv run pytest tests/ --cov=sp_pharkin -q` (173 tests, 100% coverage)

## When Adding a New Function

1. **Use `generic_a_eq_b_x_c()`** if solving a multiplication/division equation
2. **Call `format_output()`** at the end to ensure correct 5-tuple format
3. **Add to `__init__.py`** for public API export
4. **Write tests** in appropriate module (new functions warrant new test classes)
5. **Use fixtures** from `tests/test_fixtures.py` for consistency

Example:
```python
def my_calculation(param_a=None, param_b=None, param_c=None, **kwargs):
    output_unit = kwargs.pop('output_unit', False)
    decimals = kwargs.pop('decimals', 2)
    
    result = generic_a_eq_b_x_c(param_a, param_b, param_c, 
                                 names=('Result', 'Input1', 'Input2'))
    return format_output(result, 'mg', output_unit, decimals)
```

## Critical Conventions

- **Input validation in `generic_a_eq_b_x_c()`**: Raises `ValueError` if not exactly 2 of 3 parameters provided
- **Unit safety**: Always use string units in tests (Pint parses them)
- **Type annotations**: Files use `# type: ignore[assignment]` for Pint generic types (legitimate Pylance limitation)
- **Test parametrization**: Use `@pytest.mark.parametrize()` to reduce duplication (see `test_parametrized_core.py`)

## Common Debugging

| Issue | Solution |
|-------|----------|
| `ValueError: exactly 2 of 3 parameters required` | Function needs exactly 2 params; zero values are falsy |
| `DimensionalityError` from Pint | Unit string is invalid; check Pint syntax |
| Tests import fixtures from conftest | Fixtures auto-discovered; use names from `test_fixtures.py` |

## Project Goals

- Implement book formulas accurately with unit safety
- Maintain 100% test coverage
- Provide consistent, predictable API for pharmacokinetics practitioners
