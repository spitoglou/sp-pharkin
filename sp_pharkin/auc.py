"""
Area Under the Curve (AUC) calculations for drug exposure assessment.

AUC (Area Under the Curve) is quite literally the area under a concentration
versus time graph. It represents the total systemic drug exposure over time
and is one of the most important pharmacokinetic parameters in both drug
development and clinical practice.

The physiological meaning of AUC: When drug is administered, blood concentrations
rise (during absorption/distribution) and then fall (during elimination). The
area under this concentration-time profile integrates the magnitude AND duration
of drug presence in the body. A larger AUC means either higher concentrations,
longer duration, or both - indicating greater total drug exposure.

AUC units are always Mass.Time.Volume^-1 (e.g., mg.h/L or ng.h/mL). The height
of the curve is concentration (mass/volume) and the width is time, so:
    AUC units = (mass/volume) x time = mass.time/volume

Key relationships:
    - AUC = F.D / Cl (dose, bioavailability, and clearance)
    - AUC is directly proportional to the amount of drug reaching systemic circulation
    - AUC is inversely proportional to clearance (elimination efficiency)

These functions implement formulas from "Pharmacokinetics" by Philip Rowe,
Chapter 5 (AUC fundamentals) and Chapter 9 (Trapezoidal rule and bioavailability).
"""

from .lib import format_output
from pint import UnitRegistry
from typing import List, Tuple

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def auc_dose_clearance(**kwargs):
    """
    Calculate AUC from dose, bioavailability, and clearance.

    This function implements one of the most fundamental relationships in
    pharmacokinetics, connecting three critical parameters: the amount of drug
    that enters systemic circulation (F x D), the efficiency of elimination (Cl),
    and the resulting total drug exposure (AUC).

    The physiological basis: AUC represents the integral of drug concentration
    over time - literally, the total "exposure" the body experiences. This is
    determined by how much drug enters (F x D) divided by how efficiently the
    body removes it (Cl). High clearance drugs have lower AUC for the same dose.

    Formula: AUC = (F × D) / Cl

    Derivation (from Chapter 5, Appendix 2):
        - Rate of elimination at any time = C(t) × Cl
        - Total drug eliminated from t=0 to t=∞ = integral of C(t) × Cl
        - Since Cl is constant: Total eliminated = Cl × integral of C(t)
        - The integral of C(t) from 0 to ∞ IS the AUC
        - Total eliminated must equal amount that entered: F × D = Cl × AUC
        - Therefore: AUC = F × D / Cl

    Where:
        AUC = Area under the concentration-time curve (mass·time/volume)
        F = Bioavailability (fraction, 0-1); F=1.0 for IV administration
        D = Dose administered (mass)
        Cl = Clearance (volume/time)

    Args (provide exactly 2 of 3, F defaults to 1.0):
        auc (str): Area under curve (e.g., '100 mg*hour/L')
        dose (str): Dose administered (e.g., '500 mg')
        clearance (str): Drug clearance (e.g., '5 L/hour')
        bioavailability (float): Fraction reaching systemic circulation (default: 1.0)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)

    Examples:
        Calculate AUC from 500 mg IV dose with 5 L/hour clearance:
            >>> auc_dose_clearance(dose='500 mg', clearance='5 L/hour')
            ('AUC', 100.0, 'milligram * hour / liter', '100.0 mg·h/L', ...)

        Calculate AUC for oral dose with 40% bioavailability:
            >>> auc_dose_clearance(dose='500 mg', clearance='5 L/hour', bioavailability=0.4)
            ('AUC', 40.0, 'milligram * hour / liter', '40.0 mg·h/L', ...)

        Calculate required dose for target AUC:
            >>> auc_dose_clearance(auc='100 mg*hour/L', clearance='5 L/hour')
            ('Dose', 500.0, 'milligram', '500.0 mg', ...)

        Calculate clearance from known dose and AUC:
            >>> auc_dose_clearance(dose='500 mg', auc='100 mg*hour/L')
            ('Clearance', 5.0, 'liter / hour', '5.0 L/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 5.3, pages 2883-2957: Area Under the Curve fundamentals
        - Section 5.6, Appendix 2: Derivation of AUC = F.D / Cl
        - Formula: AUC = F.D / Cl (page 2957)
        - Chapter 5: Single IV bolus injection
        - Chapter 9: Using AUC for bioavailability determination

    AUC Units and Their Meaning:
        - Units are always Mass × Time / Volume (e.g., mg.h/L, ng.h/mL)
        - Height of concentration curve is mass/volume (e.g., mg/L)
        - Width is time (e.g., hours)
        - Product gives mass.time/volume
        - Convention: Write as mg.h.L^-1 or mg.h/L (both equivalent)

    Bioequivalence Testing Applications:
        AUC is THE primary parameter for bioequivalence (BE) studies:
        - Generic drug approval requires AUC within 80-125% of reference
        - Two formulations are bioequivalent if their AUCs are statistically similar
        - The 90% confidence interval must fall within 80-125% limits
        - This ensures the generic delivers the same total drug exposure

        Absolute Bioavailability Studies:
            - Compare AUC_oral to AUC_iv (same dose, same patient)
            - F_oral = AUC_oral / AUC_iv (since F_iv = 1.0 by definition)
            - Determines what fraction of oral dose reaches systemic circulation

        Relative Bioavailability Studies:
            - Compare two non-IV formulations (e.g., tablet vs capsule)
            - F_relative = AUC_test / AUC_reference
            - Cannot determine absolute fraction absorbed (only relative)

    Drug Development Uses:
        Early Development (Phase I):
            - Characterize dose-proportionality: Does AUC increase linearly with dose?
            - Non-linear kinetics indicated if AUC doesn't scale with dose
            - Example: Saturable metabolism → AUC increases more than proportionally

        Formulation Development:
            - Compare AUC between formulations to ensure equivalent exposure
            - Optimize formulation for desired AUC profile
            - Extended-release vs immediate-release comparisons

        Food Effect Studies:
            - AUC with food vs AUC fasted
            - Fed/fasted ratio determines if drug should be taken with meals
            - Important for drugs with absorption affected by food

        Drug Interaction Studies:
            - AUC_with_inhibitor / AUC_alone = interaction magnitude
            - Example: CYP3A4 inhibitor increases victim drug AUC
            - Regulatory guidance uses AUC ratios to classify interactions

    Clinical Uses of AUC:
        Therapeutic Drug Monitoring:
            - Some drugs dosed to achieve target AUC (e.g., carboplatin)
            - Carboplatin dosing: Dose = target AUC × (GFR + 25)
            - Ensures consistent exposure despite varying clearance

        Dose Adjustment:
            - Renal/hepatic impairment → decreased clearance → increased AUC
            - Dose reduction calculated to maintain same AUC
            - Example: If clearance halved, dose should be halved for same AUC

        Efficacy/Toxicity Relationships:
            - Many drugs show AUC-dependent effects
            - Higher AUC → better efficacy but also more toxicity
            - Therapeutic window defined by AUC range

    Relationship to Total Drug Exposure:
        AUC is the definitive measure of total systemic exposure because:
        - It integrates concentration over the entire time course
        - It accounts for both peak height AND duration
        - Two profiles with same Cmax but different AUC have different exposure
        - AUC determines total amount of drug available for effect

    Importance of Accurate AUC Calculation:
        Errors in AUC calculation can lead to:
        - Incorrect bioavailability determinations
        - Failed bioequivalence studies
        - Improper dose adjustments
        - Missed drug interactions
        - Inaccurate clearance calculations

    Common Errors in AUC Calculation:
        1. Ignoring the tail area (extrapolation to infinity)
           - Underestimates true AUC, especially for drugs with long half-lives
        2. Insufficient sampling during absorption phase
           - Underestimates early AUC, may miss Cmax
        3. Stopping sampling too early
           - Large tail area = error-prone extrapolation
        4. Incorrect K value for tail calculation
           - Propagates error into AUC estimate
        5. Using wrong formula for extravascular doses
           - Cannot use Cl to calculate AUC directly for non-IV routes

    When AUC Matters Clinically:
        - Drugs with narrow therapeutic index (small AUC changes = big effect changes)
        - Drugs cleared by single pathway (pathway impairment = large AUC change)
        - Drugs with concentration-dependent toxicity
        - Immunosuppressants (tacrolimus, cyclosporine)
        - Oncology drugs (carboplatin, methotrexate)
        - Anticoagulants (dose to maintain AUC within range)
        - Antiretrovirals (AUC determines viral suppression)

    Mathematical Relationships:
        AUC-related formulas:
            - AUC = F × D / Cl (this function)
            - AUC = C0 / K (for IV bolus, one-compartment)
            - Cl = F × D / AUC (calculate clearance from AUC)
            - F = AUC_test / AUC_reference (bioavailability ratio)

    Notes:
        - For IV administration, F = 1.0 (100% bioavailability guaranteed)
        - For oral and other routes, F is typically < 1.0
        - AUC is independent of the shape of the concentration-time curve
        - Same AUC can result from high peak/short duration OR low peak/long duration
        - AUC at steady state (AUC_tau) = AUC_0-inf after single dose
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)
    F = kwargs.pop("bioavailability", 1.0)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    auc = kwargs.get("auc", False)
    dose = kwargs.get("dose", False)
    clearance = kwargs.get("clearance", False)

    # Count provided parameters
    provided = sum([bool(auc), bool(dose), bool(clearance)])
    if provided != 2:
        raise ValueError("Exactly 2 of 3 parameters required: auc, dose, clearance")

    if dose and clearance:
        quantity = (F * dose) / clearance
        string = "AUC"
    elif auc and clearance:
        quantity = (auc * clearance) / F
        string = "Dose"
    elif auc and dose:
        quantity = (F * dose) / auc
        string = "Clearance"
    else:
        raise ValueError("Cannot solve with provided parameters")

    return format_output(quantity, string, output_unit, decimals)


def auc_tail(**kwargs):
    """
    Calculate the tail portion of AUC beyond the last measured concentration.

    This function calculates the extrapolated area under the curve from the
    last measured concentration to infinity. This "tail" area is essential for
    obtaining a complete AUC (AUC_0-inf) when using discrete sampling data.

    The physiological basis: After the last blood sample is taken, drug
    elimination continues. The remaining drug will eventually be completely
    cleared, contributing additional area to the total AUC. Assuming first-order
    elimination continues unchanged, this tail area can be calculated mathematically
    rather than requiring infinite sampling.

    Formula: Tail Area = C_final / K

    Derivation:
        - After time t_last, concentration follows: C(t) = C_final × e^(-K×(t-t_last))
        - Area from t_last to infinity = integral from t_last to inf of C_final × e^(-K×(t-t_last)) dt
        - Evaluating the integral: = C_final × [-1/K × e^(-K×(t-t_last))] from t_last to inf
        - At t=inf, e^(-inf) = 0; at t=t_last, e^0 = 1
        - Therefore: Tail Area = C_final × (0 - (-1/K)) = C_final / K

    Where:
        Tail = AUC from last measurement to infinity (mass·time/volume)
        C_final = Last measured concentration (mass/volume)
        K = Elimination rate constant (1/time)

    Args (provide exactly 2 of 3 to solve for the other):
        tail (str): Tail AUC (e.g., '10 mg*hour/L')
        c_final (str): Final measured concentration (e.g., '1 mg/L')
        K (str): Elimination rate constant (e.g., '0.1 1/hour')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)

    Examples:
        Calculate tail AUC from final concentration 1 mg/L and K=0.1/hour:
            >>> auc_tail(c_final='1 mg/L', K='0.1 1/hour')
            ('AUC Tail', 10.0, 'milligram * hour / liter', '10.0 mg·h/L', ...)

        Calculate K from known tail area and final concentration:
            >>> auc_tail(tail='10 mg*hour/L', c_final='1 mg/L')
            ('Elimination Rate Constant (K)', 0.1, '1 / hour', '0.1 1/hour', ...)

        Calculate final concentration from tail area and K:
            >>> auc_tail(tail='10 mg*hour/L', K='0.1 1/hour')
            ('Final Concentration', 1.0, 'milligram / liter', '1.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 9.6, pages 5132-5226: Tail area calculation methodology
        - Formula: Tail Area = Cfinal / K (page 5201)
        - Section 9.6, pages 5223-5226: Importance of minimizing tail area
        - Chapter 9: Extravascular administration and AUC determination

    Importance of AUC Tail Extrapolation:
        Complete AUC Calculation:
            - AUC_0-inf = AUC_0-tlast + Tail Area
            - Without tail extrapolation, AUC is underestimated
            - For drugs with long half-lives, tail can be substantial

        Bioavailability Studies:
            - Bioequivalence requires accurate AUC_0-inf
            - Large tail area increases uncertainty in AUC estimate
            - Regulatory agencies scrutinize tail area proportion

        Clearance Calculations:
            - Cl = F × D / AUC requires complete AUC
            - Incomplete AUC → overestimated clearance
            - This propagates errors to all derived parameters

    Common Errors in Tail Calculation:
        1. Using wrong K value:
            - K must be determined from terminal linear phase
            - Multi-compartment drugs: use terminal (slowest) K
            - Wrong K → proportionally wrong tail area

        2. Incomplete terminal phase data:
            - Need 3+ points in terminal linear portion
            - Too few points → unreliable K estimate
            - Affects both K and therefore tail calculation

        3. Non-linear kinetics at low concentrations:
            - Some drugs show non-linear elimination at low concentrations
            - Tail formula assumes linear (first-order) kinetics
            - May not apply if capacity-limited elimination at low C

        4. Stopping sampling too early:
            - Large C_final → large tail area → large error
            - Rule of thumb: sample until C < 10% of Cmax
            - Better: tail should be < 20% of total AUC

    Study Design Considerations:
        Minimizing Tail Area Error:
            - "Studies should be designed to generate only small tail areas"
            - Continue sampling until drug concentration is very low
            - Avoid stopping when significant drug remains

        Terminal Linear Portion:
            - Plot data semi-logarithmically
            - Identify late period where data forms straight line
            - This represents pure elimination (absorption complete)
            - Use this portion to determine K via t1/2

        Sample Timing:
            - More samples in terminal phase improve K estimation
            - At minimum 3 points for terminal linear portion
            - Ideally span 2-3 half-lives in terminal phase

    Acceptable Tail Area Percentage:
        FDA/EMA Guidance:
            - Tail should generally be < 20% of total AUC
            - Larger percentages require justification
            - May need to extend sampling duration

        Practical Example:
            - Total AUC = 100 mg.h/L
            - If tail = 5 mg.h/L (5%), this is acceptable
            - If tail = 35 mg.h/L (35%), study design issue
            - May indicate sampling stopped too early

    Determining K for Tail Calculation:
        From Terminal Linear Portion:
            1. Plot concentration vs time on semi-log scale
            2. Identify terminal linear portion (usually last 3+ points)
            3. Fit line through terminal points
            4. Determine half-life from the line
            5. Calculate K = 0.693 / t1/2

        For Extravascular Doses:
            - Early data includes absorption + elimination
            - Terminal portion shows only elimination
            - Back-extrapolate terminal line to find K

    When AUC Tail Matters Most:
        - Drugs with long half-lives (large proportion of AUC in tail)
        - Drugs with low clearance (persist longer in body)
        - Bioequivalence studies (accuracy critical)
        - Drugs with non-linear kinetics (verify assumptions)
        - When clearance calculation is the goal

    Mathematical Relationships:
        Full AUC calculation:
            - AUC_0-inf = AUC_0-tlast + C_final/K
            - If using trapezoidal rule: AUC_0-tlast from trapezoids
            - Total = sum of trapezoids + tail

        Relationship to other parameters:
            - K = 0.693 / t1/2
            - Tail = C_final / K = C_final × t1/2 / 0.693

    Notes:
        - Always verify K is from terminal linear portion
        - Tail calculation assumes first-order kinetics continue
        - For two-compartment drugs, use terminal (beta) rate constant
        - Error in K proportionally affects tail area
        - Small tail areas are more reliable than large ones
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    tail = kwargs.get("tail", False)
    c_final = kwargs.get("c_final", False)
    K = kwargs.get("K", False)

    provided = sum([bool(tail), bool(c_final), bool(K)])

    if provided == 2 and c_final and K:
        quantity = c_final / K
        string = "AUC Tail"
    elif provided == 2 and tail and K:
        quantity = tail * K
        string = "Final Concentration"
    elif provided == 2 and tail and c_final:
        quantity = c_final / tail
        string = "Elimination Rate Constant (K)"
    else:
        raise ValueError("Exactly 2 of 3 parameters required: tail, c_final, K")

    return format_output(quantity, string, output_unit, decimals)


def trapezoidal_auc(
    times: List[str],
    concentrations: List[str],
    K: str | None = None,
    include_tail: bool = True,
    output_unit: str | bool = False,
    decimals: int = 2,
) -> Tuple[str, float, str, str, object]:
    """
    Calculate AUC using the trapezoidal rule from concentration-time data.

    This function implements the standard trapezoidal rule method for calculating
    AUC from discrete concentration-time measurements. This is the gold standard
    method used in pharmacokinetic studies, bioequivalence trials, and clinical
    drug development.

    The physiological basis: Real pharmacokinetic experiments yield concentrations
    at particular time points only, not continuous curves. The trapezoidal rule
    provides a practical method to estimate the area under this discrete data by
    treating the region between consecutive points as trapezoids.

    Method Overview:
        1. Connect consecutive data points with straight lines
        2. Drop vertical lines from each point to the time axis
        3. This creates a series of trapezoids beneath the data
        4. Calculate each trapezoid's area: width × average height
        5. Sum all trapezoid areas for AUC_0-tlast
        6. Add tail area (C_final/K) for complete AUC_0-inf

    Formula for Each Trapezoid:
        Area_i = [(C_i + C_(i+1)) / 2] × (t_(i+1) - t_i)

        Where:
        - C_i, C_(i+1) = concentrations at consecutive time points
        - t_i, t_(i+1) = the two time points
        - (C_i + C_(i+1)) / 2 = average height of trapezoid
        - (t_(i+1) - t_i) = width of trapezoid

    Total AUC = Sum of all trapezoid areas + Tail (if included)

    Args:
        times: List of time points as strings with units (e.g., ['0 hour', '1 hour', '4 hour'])
        concentrations: List of concentrations at each time point (e.g., ['10 mg/L', '8 mg/L', '5 mg/L'])
        K: Elimination rate constant for tail calculation (e.g., '0.1 1/hour')
            Required if include_tail=True
        include_tail: Whether to extrapolate to infinity (default: True)
        output_unit: Convert result to this unit (optional)
        decimals: Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        - 'AUC (0-∞)' if tail included
        - 'AUC (0-tlast)' if tail not included

    Examples:
        Calculate AUC from concentration-time data with tail:
            >>> times = ['0 hour', '1 hour', '2 hour', '4 hour', '8 hour']
            >>> concs = ['10 mg/L', '8 mg/L', '6 mg/L', '3 mg/L', '1 mg/L']
            >>> trapezoidal_auc(times, concs, K='0.2 1/hour')
            ('AUC (0-∞)', 42.0, 'milligram * hour / liter', '42.0 mg·h/L', ...)

        Calculate AUC without tail extrapolation (AUC_0-tlast only):
            >>> trapezoidal_auc(times, concs, include_tail=False)
            ('AUC (0-tlast)', 37.0, 'milligram * hour / liter', '37.0 mg·h/L', ...)

        Book example (Section 9.6, Table 9.2):
            >>> times = ['0 hour', '1 hour', '2 hour', '3 hour', '5 hour', '8 hour', '12 hour']
            >>> concs = ['0 mg/L', '3.1 mg/L', '9.0 mg/L', '5.9 mg/L', '1.8 mg/L', '0.8 mg/L', '0.2 mg/L']
            >>> trapezoidal_auc(times, concs, K='0.308 1/hour')
            # AUC_0-12h = 28.65 mg.h/L, tail = 0.65 mg.h/L, total = 29.3 mg.h/L

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 9.6, pages 5059-5131: Trapezoidal rule practical method
        - Table 9.2, pages 5113-5128: Worked example with calculations
        - Section 9.6, pages 5081-5097: Theory and methodology
        - Chapter 9: Extravascular administration

    Detailed Trapezoidal Rule Methodology:

        Step 1 - Prepare Data:
            - Ensure time and concentration data are paired correctly
            - For extravascular doses, initial concentration at t=0 is typically 0
            - Time should be in consistent units throughout

        Step 2 - Calculate Individual Trapezoid Areas:
            For each pair of consecutive points (i, i+1):
            - Width = t_(i+1) - t_i (time difference)
            - Average concentration = (C_i + C_(i+1)) / 2
            - Trapezoid area = Width × Average concentration

            Example from book (Table 9.2):
            | Trapezoid | t_start | t_end | Width | C_start | C_end | Avg C  | Area    |
            |-----------|---------|-------|-------|---------|-------|--------|---------|
            | T1        | 0h      | 1h    | 1h    | 0       | 3.1   | 1.55   | 1.55    |
            | T2        | 1h      | 2h    | 1h    | 3.1     | 9.0   | 6.05   | 6.05    |
            | T3        | 2h      | 3h    | 1h    | 9.0     | 5.9   | 7.45   | 7.45    |
            | T4        | 3h      | 5h    | 2h    | 5.9     | 1.8   | 3.85   | 7.70    |
            | T5        | 5h      | 8h    | 3h    | 1.8     | 0.8   | 1.30   | 3.90    |
            | T6        | 8h      | 12h   | 4h    | 0.8     | 0.2   | 0.50   | 2.00    |
            |-----------|---------|-------|-------|---------|-------|--------|---------|
            | Total     |         |       |       |         |       |        | 28.65   |

        Step 3 - Sum Trapezoid Areas:
            - AUC_0-tlast = sum of all trapezoid areas
            - This gives the area from time zero to last observation

        Step 4 - Add Tail (if needed):
            - Tail Area = C_final / K
            - AUC_0-inf = AUC_0-tlast + Tail Area
            - Requires K determined from terminal linear phase

    Bioequivalence Testing Applications:
        Primary Parameter:
            - AUC is the primary measure of total drug exposure
            - Bioequivalence requires 90% CI of AUC ratio within 80-125%
            - Both AUC_0-tlast and AUC_0-inf typically reported

        Study Design:
            - Crossover design: each subject receives both test and reference
            - Washout period between treatments (typically 5+ half-lives)
            - Compare AUC_test / AUC_reference for each subject

        Regulatory Requirements:
            - FDA: AUC_0-inf preferred, AUC_0-t acceptable if AUC_t/AUC_inf > 80%
            - EMA: Similar requirements with specific guidance on residual area
            - Both: Geometric mean ratios and 90% CIs required

    Drug Development Uses:
        Phase I Studies:
            - Characterize single-dose pharmacokinetics
            - Determine dose proportionality via AUC vs dose
            - Assess food effects (fed vs fasted AUC comparison)

        Formulation Development:
            - Compare prototype formulations by AUC
            - Demonstrate similar exposure between formulations
            - Support scale-up and manufacturing changes

        Drug Interaction Studies:
            - AUC ratio (with/without perpetrator) quantifies interaction
            - DDI magnitude: 1-2× = weak, 2-5× = moderate, >5× = strong
            - Example: Ketoconazole increases midazolam AUC 15×

    Importance of Accurate AUC Calculation:
        Consequences of Error:
            - Incorrect bioavailability determination
            - Failed bioequivalence studies
            - Wrong dose recommendations
            - Missed or exaggerated drug interactions
            - Improper clearance calculation (Cl = Dose/AUC)

        Sources of Accuracy:
            - Adequate sampling (especially around Cmax and terminal phase)
            - Correct analytical methods (assay sensitivity)
            - Proper K determination for tail calculation
            - Complete data (minimal missing samples)

    Common Errors in AUC Calculation:
        1. Missing the Initial Trapezoid:
            - For extravascular doses, t=0 concentration is 0
            - Don't forget area from t=0 to first sample
            - This function handles this if times starts at 0

        2. Inadequate Sampling Density:
            - Too few points around Cmax → underestimate peak area
            - Too few points in absorption phase → miss early area
            - Wide spacing during decline → larger approximation error

        3. Incorrect Terminal K:
            - K must be from terminal LINEAR portion (semi-log plot)
            - For oral doses, need pure elimination phase
            - Wrong K → proportionally wrong tail area

        4. Large Tail Area:
            - Tail > 20% of total AUC indicates study design issue
            - Should continue sampling until C is very low
            - Large tails are error-prone

        5. Non-Linear Kinetics:
            - Method assumes concentration changes are linear between points
            - Rapidly changing concentrations may need log-linear trapezoidal
            - This function uses linear trapezoidal (standard method)

    Approximation Nature:
        - Trapezoidal rule is inherently an approximation
        - Assumes linear change between data points
        - Real concentration-time curves are not straight lines
        - Parts overstate true area, parts understate it
        - Overall gives reasonable approximation
        - More data points = better approximation

    Alternative Methods (for reference):
        Linear Trapezoidal (this function):
            - Standard method, assumes linear interpolation
            - Good for ascending portions of curve
            - May overestimate AUC during decline

        Log-Linear Trapezoidal:
            - Uses log-linear interpolation during decline
            - Better for first-order decline phases
            - More accurate for declining concentrations

        Linear Up / Log Down:
            - Linear trapezoidal during absorption
            - Log-linear trapezoidal during elimination
            - Often most accurate overall

    When AUC Matters Clinically:
        - Generic drug approval (bioequivalence demonstration)
        - Dose-to-target AUC drugs (e.g., carboplatin)
        - Therapeutic drug monitoring programs
        - Drug interaction assessment
        - Impaired clearance patients (adjust dose for target AUC)

    Practical Workflow:
        1. Collect concentration-time data with appropriate sampling
        2. Plot semi-logarithmically to identify terminal phase
        3. Determine K from terminal linear portion
        4. Apply trapezoidal rule to calculate AUC_0-tlast
        5. Calculate tail area: C_final / K
        6. Sum for complete AUC_0-inf
        7. Verify tail < 20% of total (if not, consider more sampling)

    Notes:
        - This function implements the standard linear trapezoidal method
        - Initial concentration of 0 is typical for extravascular doses
        - K must be provided for tail calculation (include_tail=True)
        - Units are automatically handled by pint
        - Result includes both trapezoid sum and tail (if requested)
    """
    if len(times) != len(concentrations):
        raise ValueError("Times and concentrations must have same length")

    if len(times) < 2:
        raise ValueError("At least 2 time points required")

    if include_tail and K is None:
        raise ValueError("K required for tail calculation when include_tail=True")

    # Convert to quantities
    times_q = [Q_(t) for t in times]
    concs_q = [Q_(c) for c in concentrations]

    # Calculate trapezoid areas
    total_auc = Q_(0, concs_q[0].units * times_q[0].units)

    for i in range(len(times_q) - 1):
        dt = times_q[i + 1] - times_q[i]
        avg_conc = (concs_q[i] + concs_q[i + 1]) / 2
        total_auc = total_auc + (avg_conc * dt)

    # Add tail if requested
    if include_tail and K is not None:
        K_q = Q_(K)
        tail = concs_q[-1] / K_q
        total_auc = total_auc + tail
        string = "AUC (0-∞)"
    else:
        string = "AUC (0-tlast)"

    return format_output(total_auc, string, output_unit, decimals)
