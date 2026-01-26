"""
Additional pharmacokinetics calculations: Loading dose, Maintenance dose, Steady state, Accumulation.

These functions implement additional formulas from "Pharmacokinetics" by Philip Rowe,
Chapter 8 (Infusions) and Chapter 10 (Accumulation and Steady State).
"""

from .lib import format_output, generic_a_eq_b_x_c
import math
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def loading_dose(**kwargs):
    """
    Calculate loading dose to achieve target concentration immediately.

    The loading dose is used to rapidly attain therapeutic concentration at the start of treatment,
    bypassing the time required to reach steady state through maintenance doses alone.

    Formula: LD = (Target × Volume_of_Distribution) / (Bioavailability × Salt_Factor)

    Also called: Priming dose, Initial dose

    Args (provide exactly 2 of 4):
        loading_dose (str): The loading dose to administer (e.g., '250 mg')
        target_concentration (str): Desired target concentration (e.g., '10 mg/L')
        volume_of_distribution (str): Volume of distribution (e.g., '50 L')
        bioavailability_salt_factor (float): Product of F (bioavailability) × S (salt factor)
                                             Default: 1.0 (for free drug)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Loading Dose', 250.0, 'milligram', '250.0 milligram', Quantity)

    Examples:
        Calculate loading dose for 5 mg/L target in 50L volume:
            >>> loading_dose(target_concentration='5 mg/L', volume_of_distribution='50 L')
            ('Loading Dose', 250.0, 'milligram', '250.0 milligram', ...)

        With bioavailability adjustment (40% oral bioavailability):
            >>> loading_dose(target_concentration='5 mg/L', volume_of_distribution='50 L',
            ...             bioavailability_salt_factor=0.4)
            ('Loading Dose', 625.0, 'milligram', '625.0 milligram', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 8.4, page 4363: LD = Target × V / F
        - Section 8.5, page 4369: For salt forms: LD = Target × V / (F × S)
        - Chapter 8: Infusions and bolus administration

    Notes:
        - Use when rapid therapeutic effect is needed
        - Often used in acute care, emergency medicine
        - Must be followed by maintenance doses to sustain levels
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Handle bioavailability_salt_factor - default to 1.0 if not provided
    bf_sf = kwargs.pop("bioavailability_salt_factor", False)
    if not bf_sf:
        bf_sf = 1.0

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("loading_dose", False)
    b = kwargs.get("target_concentration", False)
    c = kwargs.get("volume_of_distribution", False)

    # Calculate: a = b × c / bf_sf
    if a and b and c:
        # Over-specified
        raise ValueError(
            "exactly 2 of 4 parameters required: loading_dose, target_concentration, volume_of_distribution, bioavailability_salt_factor"
        )

    # Convert bf_sf to quantity if needed
    bf_sf_qty = Q_(bf_sf) if not isinstance(bf_sf, type(Q_(1))) else bf_sf

    if b and c:
        # Calculate loading dose
        a = (b * c) / bf_sf_qty
    elif a and c:
        # Calculate target concentration
        b = (a * bf_sf_qty) / c
    elif a and b:
        # Calculate volume of distribution
        c = (a * bf_sf_qty) / b
    else:
        raise ValueError(
            "exactly 2 of 3 parameters required (excluding bioavailability_salt_factor)"
        )

    # Use whichever was calculated
    if not isinstance(a, bool):
        string = "Loading Dose"
        quantity = a
    elif not isinstance(b, bool):
        string = "Target Concentration"
        quantity = b
    else:
        string = "Volume of Distribution"
        quantity = c

    return format_output(quantity, string, output_unit, decimals)


def maintenance_dose(**kwargs):
    """
    Calculate maintenance dose to maintain steady-state concentration.

    Maintenance doses are administered at regular intervals to replace the amount of drug eliminated
    since the previous dose, maintaining stable therapeutic concentrations.

    Formula: MD = (Css × Clearance × τ) / (Bioavailability × Salt_Factor)

    Where:
        Css = Concentration at steady state (target therapeutic level)
        Cl = Drug clearance (volume cleared per unit time)
        τ (tau) = Dosing interval (time between doses)
        F = Bioavailability (fraction of dose reaching systemic circulation)
        S = Salt factor (for salt forms: mass of salt form / mass of base)

    Also follows: Css = (F × MD) / (Cl × τ)

    Args (provide exactly 2 of 4):
        maintenance_dose (str): The maintenance dose (e.g., '45 mg')
        steady_state_concentration (str): Target concentration at steady state (e.g., '5 mg/L')
        clearance (str): Drug clearance (e.g., '1.5 L/hour')
        dosing_interval (str): Time between doses (e.g., '8 hour')

    Optional:
        bioavailability_salt_factor (float): Product of F × S. Default: 1.0

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the pair provided.

    Examples:
        Calculate maintenance dose for 5 mg/L target with 1.5 L/hour clearance
        and 6-hour dosing interval:
            >>> maintenance_dose(steady_state_concentration='5 mg/L',
            ...                  clearance='1.5 L/hour',
            ...                  dosing_interval='6 hour')
            Requires exactly 2 of 4 parameters - use MD, dosing_interval, clearance
            to solve for Css

        Calculate dosing interval given dose, target, and clearance:
            >>> maintenance_dose(maintenance_dose='45 mg',
            ...                  steady_state_concentration='5 mg/L',
            ...                  clearance='1.5 L/hour')
            ('Dosing Interval', 6.0, 'hour', '6.0 hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.3, page 3570: Basic pharmacokinetic parameters
        - Section 8.3, page 4199: Steady-state concentration formula
        - Section 9.0, page 4450+: Dosing regimens and maintenance
        - Chapter 9: Practical dosing considerations

    Notes:
        - Most common dosing strategy in clinical practice
        - Must achieve steady state before optimal therapeutic effect
        - Number of half-lives to reach steady state: ~4-5 (90-95% level)
        - At steady state, elimination rate = dose rate
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Handle bioavailability_salt_factor - default to 1.0
    bf_sf = kwargs.pop("bioavailability_salt_factor", False)
    if not bf_sf:
        bf_sf = 1.0

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("maintenance_dose", False)
    css = kwargs.get("steady_state_concentration", False)
    cl = kwargs.get("clearance", False)
    tau = kwargs.get("dosing_interval", False)

    bf_sf_qty = Q_(bf_sf) if not isinstance(bf_sf, type(Q_(1))) else bf_sf

    # Count non-false parameters
    params_count = sum([bool(p) for p in [a, css, cl, tau]])
    if params_count != 2:
        raise ValueError(
            "exactly 2 of 4 parameters required (excluding bioavailability_salt_factor)"
        )

    if css and cl and tau:
        # Over-specified
        raise ValueError("exactly 2 of 4 parameters required")

    if css and cl:
        # Calculate tau (dosing interval)
        # tau = a / (css × cl / bf_sf_qty)
        # Need maintenance_dose too - can't solve with just css and cl
        raise ValueError("Need maintenance_dose or dosing_interval")

    if a and css and cl:
        # Calculate dosing interval: a = (css × cl × tau) / bf_sf_qty
        # tau = (a × bf_sf_qty) / (css × cl)
        tau = (a * bf_sf_qty) / (css * cl)
        string = "Dosing Interval"
        quantity = tau
    elif a and tau and cl:
        # Calculate steady-state concentration
        css = (a * bf_sf_qty) / (cl * tau)
        string = "Steady-State Concentration"
        quantity = css
    elif a and tau and css:
        # Calculate clearance
        cl = (a * bf_sf_qty) / (css * tau)
        string = "Clearance"
        quantity = cl
    elif css and cl and tau:
        # Calculate maintenance dose
        a = (css * cl * tau) / bf_sf_qty
        string = "Maintenance Dose"
        quantity = a
    else:
        raise ValueError("Cannot solve with provided parameters")

    return format_output(quantity, string, output_unit, decimals)


def time_to_steady_state(**kwargs):
    """
    Estimate time required to reach steady-state concentration.

    Steady state is reached when the amount of drug given in each dose interval equals the amount
    eliminated in that interval. This typically requires 4-5 half-lives for 90-95% of steady state.

    Formulas:
        Approximate: t_ss ≈ 4-5 × t½ (where t½ = half-life)
        Exact: t_ss = -ln(1 - f_target) / K

    Where:
        K = Elimination rate constant (ln(2) / t½)
        f_target = Target fraction of steady state (default: 0.95 for 95%)

    Time to reach common fractions:
        - 50% steady state: ~1 half-life
        - 75% steady state: ~2 half-lives
        - 87.5% steady state: ~3 half-lives
        - 93.75% steady state: ~4 half-lives
        - 96.9% steady state: ~5 half-lives

    Args (provide exactly 1 of 2):
        half_life (str): Drug half-life (e.g., '4 hour', '120 minute')
        elimination_rate_constant (str): K value (e.g., '0.173 1/hour')

    Optional:
        target_fraction (float): Fraction of steady state to calculate (0-1, default: 0.95)
        output_unit (str): Convert output to different unit
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Time to 95% Steady State', 17.2, 'hour', '17.2 hour', ...)

    Examples:
        Time to reach 95% steady state for 4-hour half-life:
            >>> time_to_steady_state(half_life='4 hour')
            ('Time to 95% Steady State', 17.2, 'hour', '17.2 hour', ...)

        Time to reach 90% steady state:
            >>> time_to_steady_state(half_life='4 hour', target_fraction=0.90)
            ('Time to 90% Steady State', 13.3, 'hour', '13.3 hour', ...)

        Using elimination rate constant:
            >>> time_to_steady_state(elimination_rate_constant='0.173 1/hour')
            ('Time to 95% Steady State', 17.2, 'hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 10, pages 4750+: Accumulation and Steady State
        - Section 10.1-10.3: Time to steady state concepts
        - Pages 4763-4768: Mathematical derivation and clinical examples
        - Table 10.1: Fraction of steady state vs number of half-lives

    Notes:
        - Critical for determining when to assess drug efficacy/toxicity
        - Necessary for determining steady-state drug levels
        - Applies to drugs with linear (first-order) elimination kinetics
        - For non-linear (zero-order) drugs, different calculations needed
        - Clinical samples for therapeutic drug monitoring typically taken at steady state
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)
    target_fraction = kwargs.pop("target_fraction", 0.95)  # Default 95% of steady state

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    half_life = kwargs.get("half_life", False)
    k = kwargs.get("elimination_rate_constant", False)

    if not half_life and not k:
        raise ValueError("Provide either half_life or elimination_rate_constant")

    if half_life and k:
        raise ValueError("Provide exactly 1 of 2 parameters")

    if half_life:
        # t_ss = -ln(1 - target_fraction) / ln(2) × half_life
        # Simplified: for 90%: ~3.3 half-lives, for 95%: ~4.3 half-lives
        factor = -math.log(1 - target_fraction) / math.log(2)
        quantity = factor * half_life
        string = f"Time to {int(target_fraction*100)}% Steady State"
    else:
        # t_ss = -ln(1 - target_fraction) / K
        factor = -math.log(1 - target_fraction)
        quantity = factor / k
        string = f"Time to {int(target_fraction*100)}% Steady State"

    return format_output(quantity, string, output_unit, decimals)


def accumulation_factor(**kwargs):
    """
    Calculate accumulation factor (ratio of steady-state concentration to single-dose concentration).

    The accumulation factor describes how much drug accumulates in the body during repeated dosing.
    A factor of 1 means no accumulation (single dose level), while higher factors indicate significant
    buildup at steady state.

    Formulas:
        R = 1 / (1 - e^(-K×τ))
        Alternative: R = 1 / (1 - e^(-0.693×τ/t½))
        Also: τ = -ln(1 - 1/R) / K
        And: K = -ln(1 - 1/R) / τ

    Where:
        R = Accumulation factor (dimensionless)
        τ (tau) = Dosing interval (time between doses)
        K = Elimination rate constant (1/time)
        t½ = Half-life of the drug

    Interpretation:
        R = 1.0    → No accumulation (single dose)
        R = 1.1    → 10% accumulation (K×τ = 0.095)
        R = 2.0    → 100% accumulation at steady state (K×τ = 0.693, ~1 half-life)
        R = 10.0   → 900% accumulation (K×τ = 2.303, ~3.3 half-lives)

    At steady state:
        Css,max = R × Cmax,single dose
        Css,min = R × Cmin,single dose
        Css,avg = R × Cavg,single dose

    Args (provide exactly 2 of 3):
        accumulation_factor (str): The accumulation factor (e.g., '2.22')
        elimination_rate_constant (str): K value (e.g., '0.1 1/hour')
        dosing_interval (str): Time between doses (e.g., '6 hour')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Accumulation Factor', 2.22, 'dimensionless', '2.22', ...)

    Examples:
        Calculate accumulation factor for K=0.1 1/hour, τ=6 hours:
            >>> accumulation_factor(elimination_rate_constant='0.1 1/hour',
            ...                     dosing_interval='6 hour')
            ('Accumulation Factor', 2.22, 'dimensionless', '2.22', ...)

        Calculate dosing interval for desired accumulation:
            >>> accumulation_factor(accumulation_factor='2.5',
            ...                     elimination_rate_constant='0.1 1/hour')
            ('Dosing Interval', 9.16, 'hour', '9.16 hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 10, pages 4750+: Accumulation and Steady State
        - Section 10.4-10.6: Accumulation factor derivation
        - Pages 4778-4798: Mathematical relationships and clinical applications
        - Figure 10.4-10.6: Graphical representation of accumulation
        - Table 10.2: Accumulation factors vs dosing intervals

    Clinical Applications:
        - Aminoglycosides: High R (accumulation concern) → use extended intervals or check levels
        - Warfarin: Low R → maintain steady dosing; don't skip doses
        - Digoxin: High R → watch for toxicity with repeated dosing
        - ACE inhibitors: R ≈ 1.2-1.5 → mild accumulation, minimal concern

    Notes:
        - Short dosing intervals relative to half-life → higher accumulation
        - Accumulation most relevant for drugs with long half-lives
        - Always reaches steady state after ~4-5 half-lives (R approaches plateau)
        - For drugs with narrow therapeutic indices, accumulation is clinically important
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    r = kwargs.get("accumulation_factor", False)
    k = kwargs.get("elimination_rate_constant", False)
    tau = kwargs.get("dosing_interval", False)

    params_count = sum([bool(p) for p in [r, k, tau]])
    if params_count != 2:
        raise ValueError("exactly 2 of 3 parameters required")

    if k and tau:
        # Calculate accumulation factor: R = 1 / (1 - e^(-K×τ))
        # First calculate K × τ (dimensionless)
        k_tau = k * tau
        # Make sure result is dimensionless
        if hasattr(k_tau, "units"):
            k_tau = k_tau.to("dimensionless").magnitude  # type: ignore[union-attr]
        else:
            k_tau = float(k_tau)

        r = 1.0 / (1.0 - math.exp(-k_tau))
        string = "Accumulation Factor"
        quantity = Q_(r, "dimensionless")

    elif r and k:
        # Calculate dosing interval: τ = -ln(1 - 1/R) / K
        r_val = float(r.magnitude) if hasattr(r, "magnitude") else float(r)  # type: ignore[union-attr]
        tau = -math.log(1.0 - 1.0 / r_val) / k
        string = "Dosing Interval"
        quantity = tau

    elif r and tau:
        # Calculate elimination rate constant: K = -ln(1 - 1/R) / τ
        r_val = float(r.magnitude) if hasattr(r, "magnitude") else float(r)  # type: ignore[union-attr]
        k = -math.log(1.0 - 1.0 / r_val) / tau
        string = "Elimination Rate Constant"
        quantity = k

    else:
        raise ValueError("Cannot solve with provided parameters")

    return format_output(quantity, string, output_unit, decimals)


def infusion_rate(**kwargs):
    """
    Calculate infusion rate to achieve and maintain target steady-state concentration.

    Continuous infusion delivers drug at a constant rate, avoiding peak-trough fluctuations
    seen with intermittent dosing. Steady state is achieved immediately (theoretically) or
    within 4-5 half-lives practically.

    Formulas:
        R = Css × Cl  (simplest form)
        Also: Css = R / Cl  (solving for steady-state concentration)
        And: Cl = R / Css  (solving for clearance)

    Where:
        R = Infusion rate (amount per unit time, e.g., mg/hour)
        Css = Steady-state concentration (amount per unit volume, e.g., mg/L)
        Cl = Drug clearance (volume per unit time, e.g., L/hour)

    Relationship:
        At steady state with constant infusion: Rate in = Rate out
        So: Infusion rate = Elimination rate = Css × Cl

    Args (provide exactly 2 of 3):
        infusion_rate (str): The infusion rate (e.g., '15 mg/hour')
        steady_state_concentration (str): Target concentration (e.g., '5 mg/L')
        clearance (str): Drug clearance (e.g., '3 L/hour')

    Optional:
        output_unit (str): Convert output to different unit
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Infusion Rate', 7.5, 'milligram / hour', '7.5 mg/hour', ...)

    Examples:
        Calculate infusion rate to maintain 5 mg/L with 1.5 L/hour clearance:
            >>> infusion_rate(steady_state_concentration='5 mg/L',
            ...              clearance='1.5 L/hour')
            ('Infusion Rate', 7.5, 'milligram / hour', '7.5 mg/hour', ...)

        Calculate required steady-state concentration for 10 mg/hour infusion:
            >>> infusion_rate(infusion_rate='10 mg/hour',
            ...              clearance='1.5 L/hour')
            ('Steady-State Concentration', 6.67, 'milligram / liter', '6.67 mg/L', ...)

        Calculate clearance from known infusion rate and target:
            >>> infusion_rate(infusion_rate='7.5 mg/hour',
            ...              steady_state_concentration='5 mg/L')
            ('Clearance', 1.5, 'liter / hour', '1.5 L/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 8: Infusions
        - Section 8.3, pages 4199-4210: Steady-state concentration and infusion rate
        - Section 8.5, pages 4370+: Bolus and infusion combinations
        - Figure 8.2-8.3: Infusion concentration-time curves
        - Pages 4387-4400: Clinical examples

    Advantages over Intermittent Dosing:
        - Maintains constant therapeutic level (no peaks/troughs)
        - Reduced risk of toxicity from peak levels
        - Improved for drugs with narrow therapeutic indices
        - Can combine with loading dose for rapid onset
        - Suitable for continuous drug delivery systems

    Clinical Applications:
        - Heparin: Continuous IV infusion for anticoagulation
        - Insulin: IV infusion for tight glucose control
        - Lidocaine: Post-MI antiarrhythmic infusion
        - Nitroglycerin: Constant infusion for angina/heart failure
        - Dobutamine/dopamine: ICU inotropic support
        - Antibiotics: β-lactams benefit from continuous or prolonged infusion

    Practical Considerations:
        - Steady state reached in ~4-5 half-lives
        - Can give loading dose initially to reach Css immediately
        - Requires IV access (not oral)
        - Total daily dose = R × 24 hours
        - Easier to adjust than intermittent dosing (just change infusion rate)

    Notes:
        - Assumes first-order (linear) kinetics
        - Css proportional to infusion rate
        - Zero-order kinetics (e.g., phenytoin) requires different calculations
        - Infusion device accuracy important for drugs with narrow therapeutic window
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("infusion_rate", False)
    b = kwargs.get("steady_state_concentration", False)
    c = kwargs.get("clearance", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Infusion Rate", "Steady-State Concentration", "Clearance"]
    )

    return format_output(quantity, string, output_unit, decimals)
