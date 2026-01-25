# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

sp-pharkin is a Python library for basic pharmacokinetics calculations. It implements formulas from "Pharmacokinetics" by Philip Rowe, with practice questions from the book converted into tests.

## Commands

```bash
# Install dependencies
uv sync

# Run all tests
uv run pytest

# Run a single test file
uv run pytest tests/test_chapter4.py

# Run a specific test
uv run pytest tests/test_chapter4.py::test_clearance

# Lint code
uv run ruff check .

# Format code
uv run ruff format .
```

## Dependencies

- **pint**: Unit handling library - all quantities use `Q_()` wrapper for unit-aware calculations
- **sympy**: Symbolic math for exponential decay equations (in `expo.py`)

## Architecture

### Core Pattern: `generic_a_eq_b_x_c`

Most calculation functions follow a pattern for solving equations of form `a = b × c`. Given any two variables, the third is computed. Functions use keyword arguments and return a 5-tuple:
```python
(name, magnitude, unit_string, formatted_string, pint_quantity)
```

### Module Structure

- **lib.py**: Core utilities - `format_output()` for result formatting, `generic_a_eq_b_x_c()` for the multiplication/division pattern
- **functions.py**: Main pharmacokinetics functions (volume of distribution, dose/concentration/volume, half-life, extraction rate)
- **reduction_factors.py**: Salt factor and bioavailability calculations
- **clearance.py**: Drug clearance calculations (flow × extraction rate, K × volume)
- **expo.py**: Symbolic exponential decay solver using sympy

### Usage Pattern

Functions accept string quantities with units (parsed by pint):
```python
from sp_pharkin import half_life_k
result = half_life_k(K='0.1/hour')  # Returns half-life in hours
```

Optional kwargs: `output_unit` (convert result), `decimals` (rounding precision, default 2)
