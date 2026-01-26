"""
Exponential Decay Solver Module

Solves for variables in the exponential decay equation:
C(t) = C₀·e^(-kt)

Where:
- C(t): Concentration at time t
- C₀: Initial concentration
- k: Elimination rate constant
- t: Time elapsed
"""

from .lib import format_output
from pint import UnitRegistry
import math

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def solve_for_c_t(**kwargs):
    """
    Solve for C(t): concentration at time t

    C(t) = C₀·e^(-kt)

    Args:
        c_0: Initial concentration (e.g., '100 mg/L')
        k: Elimination rate constant (e.g., '0.5 1/hour')
        t: Time elapsed (e.g., '2 hour')
        output_unit: Optional output unit (e.g., 'mg/L')
        decimals: Number of decimal places (default: 2)

    Returns:
        tuple: (name, magnitude, units, formatted_str, quantity)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    c_0 = kwargs.get("c_0", None)
    k = kwargs.get("k", None)
    t = kwargs.get("t", None)

    if c_0 is None or k is None or t is None:
        raise ValueError("c_0, k, and t are required to solve for c_t")

    # C(t) = C₀·e^(-kt)
    quantity = c_0 * math.exp(-1 * (k * t).magnitude)

    return format_output(
        quantity, "Concentration at Time t (C(t))", output_unit, decimals
    )


def solve_for_c_0(**kwargs):
    """
    Solve for C₀: initial concentration

    C₀ = C(t) / e^(-kt)

    Args:
        c_t: Concentration at time t (e.g., '50 mg/L')
        k: Elimination rate constant (e.g., '0.5 1/hour')
        t: Time elapsed (e.g., '2 hour')
        output_unit: Optional output unit (e.g., 'mg/L')
        decimals: Number of decimal places (default: 2)

    Returns:
        tuple: (name, magnitude, units, formatted_str, quantity)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    c_t = kwargs.get("c_t", None)
    k = kwargs.get("k", None)
    t = kwargs.get("t", None)

    if c_t is None or k is None or t is None:
        raise ValueError("c_t, k, and t are required to solve for c_0")

    # C₀ = C(t) / e^(-kt) = C(t) · e^(kt)
    # Convert c_t to proper concentration units before calculation
    quantity = c_t / math.exp(-1 * (k * t).magnitude)

    return format_output(quantity, "Initial Concentration (C₀)", output_unit, decimals)


def solve_for_k(**kwargs):
    """
    Solve for k: elimination rate constant

    k = -ln(C(t) / C₀) / t

    Args:
        c_0: Initial concentration (e.g., '100 mg/L')
        c_t: Concentration at time t (e.g., '50 mg/L')
        t: Time elapsed (e.g., '2 hour')
        output_unit: Optional output unit (e.g., '1/hour')
        decimals: Number of decimal places (default: 2)

    Returns:
        tuple: (name, magnitude, units, formatted_str, quantity)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    c_0 = kwargs.get("c_0", None)
    c_t = kwargs.get("c_t", None)
    t = kwargs.get("t", None)

    if c_0 is None or c_t is None or t is None:
        raise ValueError("c_0, c_t, and t are required to solve for k")

    # k = -ln(C(t) / C₀) / t
    ratio = (c_t / c_0).magnitude
    k_value = -math.log(ratio) / t.magnitude

    # Create the quantity with proper units (1/time)
    quantity = k_value / t.units

    return format_output(
        quantity, "Elimination Rate Constant (k)", output_unit, decimals
    )


def solve_for_t(**kwargs):
    """
    Solve for t: time elapsed

    t = -ln(C(t) / C₀) / k

    Args:
        c_0: Initial concentration (e.g., '100 mg/L')
        c_t: Concentration at time t (e.g., '50 mg/L')
        k: Elimination rate constant (e.g., '0.5 1/hour')
        output_unit: Optional output unit (e.g., 'hour')
        decimals: Number of decimal places (default: 2)

    Returns:
        tuple: (name, magnitude, units, formatted_str, quantity)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    c_0 = kwargs.get("c_0", None)
    c_t = kwargs.get("c_t", None)
    k = kwargs.get("k", None)

    if c_0 is None or c_t is None or k is None:
        raise ValueError("c_0, c_t, and k are required to solve for t")

    # t = -ln(C(t) / C₀) / k
    ratio = (c_t / c_0).magnitude
    t_value = -math.log(ratio) / k.magnitude

    # Create the quantity with proper units (1/k_units)
    quantity = t_value * (1 / k.units)

    return format_output(quantity, "Time Elapsed (t)", output_unit, decimals)
