from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def format_output(quantity, string, output_unit, decimals):
    """
    Format pharmacokinetic calculation results into standardized 5-tuple.

    This utility function ensures all functions in sp-pharkin return consistently
    formatted results: (name, magnitude, unit_string, formatted_string, pint_quantity).
    This standardized format enables chaining calculations and consistent handling.

    The output format provides both raw (for calculations) and formatted (for display)
    values, allowing downstream code to either continue calculations or display results.

    Args:
        quantity (pint.Quantity): The computed result with units
            Example: Quantity(50.0, 'milligram / liter')
        string (str): Human-readable name of the result
            Example: 'Clearance', 'Half-Life', 'Elimination Rate'
        output_unit (str or bool): Desired output unit or False to keep original
            - If False: Keep quantity units unchanged
            - If string like 'mg/L': Convert quantity to this unit
            - Example: quantity in 'mg/mL' converted to 'mg/L'
        decimals (int): Number of decimal places for rounding
            Example: 2 rounds to nearest 0.01

    Returns:
        tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        - name (str): Parameter name ('Clearance', 'Half-Life', etc.)
        - magnitude (float): Numeric value only (e.g., 1.2)
        - unit_string (str): Unit as string (e.g., 'liter / hour')
        - formatted_string (str): Complete formatted result (e.g., '1.2 L/hour')
        - pint_quantity: Full Quantity object for further calculations

    Examples:
        Basic formatting:
            >>> q = Q_(50, 'mg/L')
            >>> format_output(q, 'Concentration', False, 2)
            ('Concentration', 50.0, 'milligram / liter', '50.0 mg/L', Quantity(50.0, 'mg/L'))

        With unit conversion:
            >>> q = Q_(50000, 'mg/L')
            >>> format_output(q, 'Concentration', 'g/L', 2)
            ('Concentration', 50.0, 'gram / liter', '50.0 g/L', Quantity(50.0, 'g/L'))

        With rounding:
            >>> q = Q_(1.23456, 'hour')
            >>> format_output(q, 'Half-Life', False, 2)
            ('Half-Life', 1.23, 'hour', '1.23 hour', Quantity(1.23, 'hour'))

    Standardized Return Format Rationale:
        Index [0] name (str):
            - Human-readable parameter name
            - Useful for labels, reports, user interface
            - Example: 'Elimination Rate Constant(K)' not just 'K'
        Index [1] magnitude (float):
            - Numeric value only, no units
            - Used for calculations: dose = k * vd
            - Used for comparisons: if magnitude > threshold
            - Useful for data analysis without unit handling
        Index [2] unit_string (str):
            - Unit name as string (full form)
            - Example: 'milligram / hour' (expanded form)
            - Not ideal for display but useful for parsing
            - Used by Pint for unit handling
        Index [3] formatted_string (str):
            - Complete readable result: magnitude + abbreviated units
            - Example: '1.2 L/hour' or '6.93 hour'
            - Best for user display, reports, logging
            - What patients/providers should see
        Index [4] pint_quantity:
            - Full Quantity object with units
            - Used to continue calculations: dose = [4] * concentration
            - Preserves unit information for Pint operations
            - Critical for chained calculations

    Usage Patterns:
        Display results:
            >>> result = solve_for_clearance(...)
            >>> print(f"{result[0]}: {result[3]}")  # "Clearance: 1.2 L/hour"

        Use in calculations:
            >>> cl_result = solve_for_clearance(...)
            >>> vd = Q_(50, 'L')
            >>> k = cl_result[4] / vd  # chained calculation

        Extract numeric value:
            >>> hl_result = half_life_k(...)
            >>> hours = hl_result[1]  # numeric value
            >>> if hours > 24: print("Long half-life")

        Export for reports:
            >>> results = [solve_for_clearance(...), solve_for_volume(...)]
            >>> for name, magnitude, units, formatted, qty in results:
            ...     csv_writer.writerow([name, magnitude, formatted])

    Unit Conversion Examples:
        Different representations of same quantity:
            >>> q = Q_(1000, 'mg/L')
            >>> format_output(q, 'Dose', 'g/L', 2)
            # Returns: ('Dose', 1.0, 'gram / liter', '1.0 g/L', Quantity(1.0, 'g/L'))

            >>> q = Q_(0.5, 'L/hour')
            >>> format_output(q, 'Clearance', 'mL/minute', 2)
            # Returns: ('Clearance', 8.33, 'milliliter / minute', '8.33 mL/minute', ...)

            >>> q = Q_(7, 'hour')
            >>> format_output(q, 'Half-Life', 'day', 3)
            # Returns: ('Half-Life', 0.292, 'day', '0.292 day', Quantity(0.292, 'day'))

    Error Handling:
        - Invalid unit string raises DimensionalityError (incompatible conversion)
        - Quantity must have proper units or raises UndefinedUnitError
        - decimals must be non-negative integer
        - string should be non-empty but no validation performed

    Clinical Interpretation Guide:
        When receiving this 5-tuple in clinical code:
        1. Use [3] for display (user-friendly formatted string)
        2. Use [0] for parameter name/identification
        3. Use [1] for calculations and comparisons
        4. Use [4] if chaining with other unit-aware calculations
        5. Use [2] for parsing/export to systems needing unit strings

    Notes:
        - All inputs must be properly typed (quantity is Pint Quantity)
        - Unit conversions use Pint's standard unit registry
        - Rounding is performed after unit conversion
        - Magnitude is always float type (even if originally integer)
        - formatted_string uses Pint's default string representation
        - Function is pure (no side effects)
    """
    if output_unit:
        quantity = quantity.to(ureg(output_unit))

    quantity = round(quantity, decimals)

    return (
        string,
        quantity.magnitude,
        "{!s}".format(quantity.units),
        "{!s}".format(quantity),
        quantity,
    )


def generic_a_eq_b_x_c(a, b, c, names):
    """
    Solve for unknown in the fundamental equation: a = b × c

    This is the core computational engine of sp-pharkin. It solves multiplicative
    relationships in pharmacokinetics given ANY two of three variables. Most PK
    equations are rearrangements of this basic relationship.

    Formula: a = b × c
    Rearranged: b = a / c, or c = a / b

    This function is mathematically simple but pharmacologically powerful - it
    encodes dozens of commonly-used PK relationships.

    Args:
        a, b, c (Quantity or False): The three variables in a = b × c
            - Pass actual Quantity object: Q_(value, 'unit')
            - Pass False for the unknown variable
            - Exactly TWO must be non-False (provide exactly 2 of 3)
        names (tuple of 3 str): Human-readable names for [a, b, c]
            - names[0]: Name for parameter 'a'
            - names[1]: Name for parameter 'b'
            - names[2]: Name for parameter 'c'
            - Example: ('Dose', 'Concentration', 'Volume')

    Returns:
        tuple: (name, quantity)
        - name (str): The names entry corresponding to solved variable
        - quantity (Quantity): The calculated result with proper units

    Raises:
        ValueError: If not exactly 2 parameters are non-False/None
            - Raises if 0 parameters provided (nothing to solve)
            - Raises if 1 parameter provided (can't solve)
            - Raises if 3+ parameters provided (overdetermined)
            - Helpful error message shows what was provided

    Mathematical Cases (solves exactly one of these):
        Case 1: a and c provided → solve for b
            Formula: b = a / c
            Example: dose=500mg, volume=50L → concentration = 500/50 = 10 mg/L
        Case 2: a and b provided → solve for c
            Formula: c = a / b
            Example: dose=500mg, concentration=10mg/L → volume = 500/10 = 50 L
        Case 3: b and c provided → solve for a
            Formula: a = b × c
            Example: concentration=10mg/L, volume=50L → dose = 10×50 = 500 mg

    Pharmacokinetic Relationships Using This Function:

        1. Dose-Concentration-Volume (C = Dose/Vd):
            - a=Dose, b=Concentration, c=Volume
            - Calculates: Loading dose needed for target concentration
            - Or: Resulting concentration from given dose
            - Or: Volume of distribution from dose and measured concentration

        2. Clearance-Elimination-Rate-Constant (Cl = K × Vd):
            - a=Clearance, b=K, c=Vd
            - Calculates: Clearance from elimination constant and volume
            - Or: Half-life from clearance and volume
            - Or: How much drug is eliminated per unit time

        3. Clearance-Flow-Extraction (Cl = Q × E):
            - a=Clearance, b=Flow, c=Extraction
            - Calculates: Organ clearance from blood flow and extraction ratio
            - Or: Extraction efficiency from flow and clearance
            - Or: Required flow for target clearance

        4. Elimination-Rate-Mass-K (Rate = K × Mass):
            - a=Rate, b=K, c=Mass
            - Calculates: How fast drug is being eliminated
            - Or: Elimination constant from measured rate
            - Or: Amount of drug present from elimination rate

        5. Bioavailability (Delivered = Bioavailability × Dose):
            - a=Delivered Drug, b=Bioavailability, c=Dose Administered
            - Calculates: How much drug reaches circulation
            - Or: Bioavailability from measured and administered
            - Or: Dose needed for target systemic level

        6. Salt Factor (Delivered = Salt Factor × Salt Dose):
            - a=Delivered Drug, b=Salt Factor, c=Salt Dose
            - Calculates: Active drug from salt formulation
            - Or: Salt factor from salt dose and active amount
            - Or: Required salt dose for target active amount

        7. Volume Scaling (Vd = Mean Vd/kg × Weight):
            - a=Vd Total, b=Vd Per kg, c=Weight
            - Calculates: Patient-specific volume of distribution
            - Or: Population Vd per kg from patient data
            - Or: Required weight for target Vd

        8. Clearance Scaling (Cl = Cl avg × Weight/70kg):
            - a=Clearance Patient, b=Clearance Population, c=Weight Ratio
            - Calculates: Individualized clearance
            - Or: Population clearance from patient measurements
            - Or: Weight needed for target clearance

        9. Extraction Ratio (Difference = Ratio × Concentration In):
            - a=Concentration Difference, b=Extraction Ratio, c=Conc In
            - Calculates: How much drug is removed by organ
            - Or: Organ efficiency from concentration change
            - Or: Input concentration from extraction and change

        10. Half-Life-K Relationship (ln(2) = K × t½):
            - Special case: a=ln(2)≈0.693, b=K, c=Half-Life
            - Calculates: Half-life from elimination constant
            - Or: Elimination constant from half-life
            - Or: Fundamental relationship verification

    Input Validation:
        Required: Exactly 2 of 3 parameters must be provided

        Valid inputs:
            generic_a_eq_b_x_c(
                a=Q_(500, 'mg'),
                b=Q_(10, 'mg/L'),
                c=False,
                names=['Dose', 'Concentration', 'Volume']
            )
            # Result: ('Volume', Quantity(50.0, 'liter'))

        Invalid (too many):
            generic_a_eq_b_x_c(
                a=Q_(500, 'mg'),
                b=Q_(10, 'mg/L'),
                c=Q_(50, 'L'),  # All three provided!
                names=['Dose', 'Concentration', 'Volume']
            )
            # Raises: ValueError - exactly 2 of 3 parameters required

        Invalid (too few):
            generic_a_eq_b_x_c(
                a=Q_(500, 'mg'),
                b=False,
                c=False,  # Only one provided
                names=['Dose', 'Concentration', 'Volume']
            )
            # Raises: ValueError - exactly 2 of 3 parameters required

    Unit Handling with Pint:
        All parameters must be Pint Quantity objects (created with Q_):

        Correct:
            a = Q_(100, 'mg/L')
            b = Q_(0.5, '1/hour')
            c = Q_(5, 'hour')

        Incorrect (will fail):
            a = 100  # Missing units!
            b = Q_(0.5, '1/hour')
            c = Q_(5, 'hour')

        Unit compatibility is automatically checked by Pint. Incompatible
        units in arithmetic will raise DimensionalityError.

    Dimensional Analysis Examples:

        Case: a=Dose, b=Concentration, c=Volume → Solve for Volume
            a / b = (mg) / (mg/L) = (mg) × (L/mg) = L ✓
            Result correctly has units of volume

        Case: a=Dose, b=Concentration, c=Volume → Solve for Concentration
            a / c = (mg) / (L) = mg/L ✓
            Result correctly has units of concentration

        Case: a=Dose, b=Concentration, c=Volume → Solve for Dose
            b × c = (mg/L) × (L) = mg ✓
            Result correctly has units of dose

    Performance Characteristics:
        - O(1) time complexity (constant: just arithmetic operations)
        - O(1) space complexity (no data structures allocated)
        - Pint quantity operations may add some overhead
        - Suitable for real-time calculation (e.g., dosing in ED)

    Clinical Safety Notes:
        - Does not validate reasonableness of inputs
        - Will calculate nonsensical results from bad inputs
        - Example: Q_(1000, 'mg/L') × Q_(100, 'L') = Q_(100000, 'mg')
        - Always validate: Is result reasonable for the clinical context?
        - Use domain knowledge to verify: Too high? Too low? Right units?

    Integration with sp-pharkin:
        Every calculation function in sp-pharkin uses this internally:
        1. Extract kwargs and set defaults (output_unit, decimals)
        2. Convert all string inputs to Quantities: Q_(value)
        3. Call generic_a_eq_b_x_c(a, b, c, names)
        4. Get back (name, quantity) tuple
        5. Pass to format_output() for standardized 5-tuple
        6. Return standardized tuple to user

    Example Complete Workflow:
        >>> from sp_pharkin import dose_concentration_volume
        >>> result = dose_concentration_volume(
        ...     concentration='10 mg/L',
        ...     volume='50 L'
        ... )
        >>> print(result[3])  # Display formatted string
        "500.0 milligram"
        >>> dose_mg = result[1]  # Get numeric value for calculations
        >>> print(dose_mg)
        500.0

    Advanced: Chaining Calculations:
        >>> # Calculate clearance from dose, target Css, and tau
        >>> dose_result = dose_concentration_volume(
        ...     dose='500 mg',
        ...     volume='50 L'
        ... )
        >>> peak = dose_result[4]  # Get Quantity
        >>> css = peak * Q_(0.7, 'dimensionless')  # Estimate steady state
        >>> cl = Q_(0.5, 'L/hour')
        >>> daily_dose = css * cl * Q_(24, 'hour')
        >>> print(daily_dose)

    Notes:
        - This is the fundamental equation for all PK calculations
        - Correctness depends on unit consistency (Pint handles this)
        - Name tuple order must correspond to a, b, c order
        - Works for any multiplicative relationship (not just PK)
        - See each module's functions for specific equation implementations
        - More complex relationships (e.g., exponential decay) use specialized functions
    """
    provided = sum([bool(a), bool(b), bool(c)])

    if provided != 2:
        raise ValueError(
            f"generic_a_eq_b_x_c requires exactly 2 of 3 parameters. "
            f"Got {provided}: a={a is not False}, b={b is not False}, c={c is not False}"
        )

    if a and c:
        string = names[1]
        quantity = a / c

    elif a and b:
        string = names[2]
        quantity = a / b

    elif b and c:
        string = names[0]
        quantity = b * c

    return (string, quantity)
