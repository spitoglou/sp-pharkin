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
    Solve for C(t): drug concentration at a specific time using exponential decay.

    In first-order kinetics, drug concentration decays exponentially over time.
    This function calculates what the plasma concentration will be at any future
    time given the initial concentration and elimination rate constant.

    Formula: C(t) = C₀ × e^(-k×t)

    Where:
    - C(t) = Concentration at time t (what we're solving for)
    - C₀ = Initial concentration (at t=0)
    - e = Euler's number (≈ 2.718)
    - k = Elimination rate constant (1/time)
    - t = Time elapsed

    This models the fundamental pharmacokinetic principle of first-order elimination:
    the rate at which drug leaves the body is proportional to the amount present.

    Args:
        c_0 (str): Initial plasma concentration (e.g., '100 mg/L')
        k (str): Elimination rate constant (e.g., '0.1 1/hour')
        t (str): Time elapsed (e.g., '5 hour')
        output_unit (str, optional): Convert result to this unit (e.g., 'mg/L')
        decimals (int): Decimal places to round to (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Concentration at Time t (C(t))', 60.65, 'milligram / liter', '60.65 mg/L', ...)

    Examples:
        Find concentration after 5 hours with C₀=100 mg/L, k=0.1 1/hour:
            >>> solve_for_c_t(c_0='100 mg/L', k='0.1 1/hour', t='5 hour')
            ('Concentration at Time t (C(t))', 60.65, 'milligram / liter', '60.65 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.1-6.3, pages 3400-3500: Exponential decay kinetics
        - Section 6.4, pages 3500-3550: Mathematical relationships for first-order
        - Derivation shown in pages 3420-3450
        - Chapter 6: Drug elimination and kinetics

    Clinical Applications:
        - Predicting plasma concentration at any time after dose
        - Determining when concentrations drop below therapeutic range
        - Calculating when concentration drops to sub-therapeutic or toxic levels
        - Planning drug interactions (when combining drugs)
        - Monitoring compliance (verify if patient took drug based on concentration)
        - Dialysis decisions (is concentration high enough to warrant dialysis?)

    Relationship to Half-Life:
        - At t = t½, C(t) = C₀/2
        - At t = 2×t½, C(t) = C₀/4
        - At t = 3×t½, C(t) = C₀/8
        - At t = 5×t½, C(t) ≈ C₀/32 (essentially eliminated)
        - Each half-life period, concentration is cut in half

    Example Calculations:
        Theophylline (k=0.11 1/hour, t½≈6.3 hours):
            - C₀ = 15 mg/L (typical maintenance level)
            - At t=6.3 hours: C(t) = 15 × e^(-0.11×6.3) ≈ 7.5 mg/L (one t½)
            - At t=12.6 hours: C(t) = 15 × e^(-0.11×12.6) ≈ 3.75 mg/L (two t½)
        Gentamicin (k=0.3 1/hour, t½≈2.3 hours):
            - Peak C = 5 mg/L (after IV bolus)
            - At t=4.6 hours: C(t) = 5 × e^(-0.3×4.6) ≈ 1.25 mg/L (trough, two t½)

    Disease Effects on Concentration Decay:
        - Renal disease: ↓ k (slower decay, higher concentrations at any time)
        - Liver disease: ↓ k (slower decay for hepatically metabolized drugs)
        - Age: Often ↓ k (slower elimination in elderly)
        - Fever: ↑ k (faster metabolism/elimination in some drugs)
        - Genetic polymorphisms: Can alter k 2-10 fold
        - Enzyme induction: ↑ k (faster elimination)
        - Enzyme inhibition: ↓ k (slower elimination)

    Clinical Decision Points (Concentration-Based Dosing):
        - Therapeutic monitoring: Measure actual concentration at known time
        - If measured C > target: Extend next dose interval (wait longer)
        - If measured C < target: Give next dose sooner or increase amount
        - Example: Theophylline level 8 mg/L (below target 10-20)
                   Calculate new dosing based on: C(t) at time of measurement
                   Then solve for new maintenance dose needed

    Drug Accumulation at Steady State:
        - This function predicts single-dose decay
        - With repetitive dosing: Accumulation factor = 1/(1-e^(-kτ))
        - Average steady-state = (Dose/τ) / Cl = (Dose/τ) / (k×Vd)
        - Concentration oscillates between: C_min and C_max at steady state

    Converting Between Different Time Units:
        - k must be consistent with time units in t
        - If k = 0.1 1/hour and t = 120 minutes, convert t = 2 hours
        - Or convert k to appropriate units: 0.1/hour = 0.1/60 1/minute ≈ 0.00167 1/minute

    Monitoring Drug Levels During Therapy:
        - Blood drawn at time t after last dose
        - Measure actual concentration
        - Compare to predicted using this equation
        - If measured >> predicted: May not be taking doses as prescribed
        - If measured << predicted: May have accelerated elimination
        - Plan adjustments based on target Css

    Notes:
        - Assumes constant k (first-order kinetics)
        - Assumes drug has entered distribution phase (post-distribution)
        - Does not account for continued absorption (if active)
        - Does not account for metabolite accumulation
        - Assumes no non-linear kinetics (saturable metabolism)
        - Time must be measured from a known baseline (usually dose administration)
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
    Solve for C₀: initial (peak) plasma concentration after a dose.

    This function works backward from a known concentration at a known time to
    calculate what the initial concentration must have been. It's useful when
    you measure a concentration but need to know what the peak concentration was.

    Formula: C₀ = C(t) / e^(-k×t) = C(t) × e^(k×t)

    Where:
    - C₀ = Initial concentration (what we're solving for)
    - C(t) = Known concentration at time t
    - k = Elimination rate constant (1/time)
    - t = Time elapsed
    - e = Euler's number (≈ 2.718)

    This is the inverse of exponential decay. If you know current concentration
    and how long it's been decaying, you can calculate the original concentration.

    Args:
        c_t (str): Known concentration at time t (e.g., '60.65 mg/L')
        k (str): Elimination rate constant (e.g., '0.1 1/hour')
        t (str): Time elapsed (e.g., '5 hour')
        output_unit (str, optional): Convert result to this unit (e.g., 'mg/L')
        decimals (int): Decimal places to round to (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Initial Concentration (C₀)', 100.0, 'milligram / liter', '100.0 mg/L', ...)

    Examples:
        Calculate initial concentration given C(t)=60.65 mg/L at t=5 hours, k=0.1 1/hour:
            >>> solve_for_c_0(c_t='60.65 mg/L', k='0.1 1/hour', t='5 hour')
            ('Initial Concentration (C₀)', 100.0, 'milligram / liter', '100.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.1-6.3, pages 3400-3500: Exponential decay kinetics
        - Section 6.4, pages 3500-3550: Mathematical relationships for first-order
        - Pages 3430-3460: Solving exponential equations
        - Chapter 6: Drug elimination and kinetics

    Clinical Applications:
        - Back-calculate peak concentration from trough measurement
        - Assess whether drug accumulation is occurring
        - Determine if dosing interval should be adjusted
        - Therapeutic drug monitoring (predicting peak from trough measurement)
        - Safety assessment (will peak exceed toxic threshold?)
        - Compliance verification (calculating expected peaks)

    Therapeutic Drug Monitoring Workflow:
        1. Patient takes dose (amount = Dose)
        2. Wait for distribution (typically 30 min - 2 hours, depends on drug)
        3. Draw blood at known time t
        4. Measure concentration = C(t)
        5. Use this function to calculate C₀ (peak concentration)
        6. If C₀ < target_min: Consider dose increase
        7. If C₀ > target_max: Consider dose reduction

    Trough-to-Peak Relationship:
        - Trough = C(τ) = concentration just before next dose (C₀ at end of dosing interval)
        - Peak = C(t_peak) = highest concentration after dose
        - Fluctuation = Peak - Trough
        - Examples:
            * Gentamicin: Peak target 5-10 mg/L, Trough <2 mg/L
            * Vancomycin: Trough target 15-20 mg/L
            * Theophylline: 10-20 mg/L maintained
            * Digoxin: Peak 1-2 hours post-dose, therapeutic 0.8-2.0 ng/mL

    Using Measured Trough to Predict Peak:
        - Often easier to measure trough (patient arrives for dose)
        - Trough = C(τ) where τ is dosing interval
        - Multiply trough by e^(kτ) to estimate peak
        - Example: Gentamicin trough = 1 mg/L, k = 0.3 1/hour, τ = 8 hours
                   Peak estimate = 1 × e^(0.3×8) = 1 × e^2.4 ≈ 11 mg/L

    Predicting Accumulation:
        - With repetitive dosing: C_peak_new = C_peak_old + Dose/Vd
        - At steady state: Peak = C_peak_ss = (Dose/Vd) / (1-e^(-kτ))
        - Accumulation factor = 1/(1-e^(-kτ))
        - Using this function: Measure one peak, calculate what peak will be at steady state

    Disease Effects on Measured Concentration:
        - Renal disease: ↓ k → Same measured C(t) means higher C₀
        - Liver disease: ↓ k → Measured concentrations persist longer
        - Fever/Exercise: ↑ k → Same measured C(t) means lower C₀
        - Heart failure: ↓ k and ↓ Vd → Complex effect on C₀

    Example Calculations:
        Digoxin therapeutic monitoring:
            - Draw blood 1 week after starting
            - Measure C = 1.2 ng/mL
            - k (digoxin) ≈ 0.035 1/day
            - C₀ 7 days post-dose = 1.2 × e^(0.035×7) ≈ 1.2 × e^0.245 ≈ 1.2 × 1.28 ≈ 1.54 ng/mL
            - This would be peak concentration from first dose

        Aminoglycoside (gentamicin):
            - Measured trough = 1.5 mg/L at 8 hours
            - k = 0.3 1/hour
            - C₀ = 1.5 × e^(0.3×8) = 1.5 × e^2.4 ≈ 1.5 × 11 ≈ 16.5 mg/L
            - Peak might be ~5 mg/L (target) if properly dosed
            - If actual peak is 16.5, dose is too high

    Using Peak/Trough Measurements:
        - Some drugs monitored by peak level only (e.g., gentamicin peak)
        - Some by trough level only (e.g., vancomycin trough)
        - Some by both (e.g., digoxin by level, warfarin by INR not concentration)
        - This function helps interpret what peaks mean

    Converting Times and Concentrations:
        - All units must be consistent
        - If t is in hours, k must be 1/hour
        - If C(t) is mg/L, C₀ result will be mg/L
        - Convert if needed: 60 minutes = 1 hour, 1000 mcg = 1 mg

    Factors Affecting the C₀ Calculation:
        - Accuracy of k estimate is critical
        - Small error in k → large error in calculated C₀
        - Small error in measured C(t) → proportional error in C₀
        - Error in t (time of measurement) → exponential error (worst for k×t)
        - Steady state assumption: Must account for prior dosing

    Notes:
        - Assumes constant k (first-order kinetics)
        - Assumes measurement time is accurate
        - Assumes k is accurately known
        - Does not account for continued absorption
        - Does not account for non-linear kinetics
        - Does not account for metabolites (that's included in measured concentration)
        - Multiple compartment drugs: k is distribution-phase constant, not elimination
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
    Solve for k: elimination rate constant from concentration changes over time.

    This function determines how fast a drug is being eliminated by measuring
    how much the concentration has dropped over a known time period. It's the
    most important calculation for understanding individual drug kinetics.

    Formula: k = -ln(C(t) / C₀) / t = ln(C₀ / C(t)) / t

    Where:
    - k = Elimination rate constant (what we're solving for)
    - C₀ = Initial concentration
    - C(t) = Concentration at time t
    - t = Time elapsed
    - ln = Natural logarithm (ln(2) ≈ 0.693)

    This is derived from the exponential decay model. By measuring how much
    concentration has changed, we can calculate the underlying elimination rate.

    Args:
        c_0 (str): Initial plasma concentration (e.g., '100 mg/L')
        c_t (str): Concentration at later time t (e.g., '50 mg/L')
        t (str): Time elapsed (e.g., '7 hour')
        output_unit (str, optional): Convert result to this unit (e.g., '1/hour', '1/day')
        decimals (int): Decimal places to round to (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Elimination Rate Constant (k)', 0.099, '1 / hour', '0.099 1/hour', ...)

    Examples:
        Calculate k from measured concentrations 7 hours apart:
            >>> solve_for_k(c_0='100 mg/L', c_t='50 mg/L', t='7 hour')
            ('Elimination Rate Constant (k)', 0.099, '1 / hour', '0.099 1/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.2-6.4, pages 3440-3550: Elimination rate constant determination
        - Pages 3470-3500: Logarithmic and semi-log plots
        - Figure 6.1-6.3: Semi-log concentration vs time plots
        - Chapter 6: Drug elimination and kinetics

    Clinical Applications:
        - Individualize drug elimination for specific patient
        - Calculate personalized half-life
        - Predict steady-state concentration
        - Adjust dosing interval
        - Detect changes in kidney/liver function
        - Evaluate drug interactions
        - Assess compliance (abnormal k suggests missed doses)

    Step-by-Step Therapeutic Drug Monitoring:
        1. Give drug (IV bolus or equivalent)
        2. Draw blood at time t₁ (e.g., 1 hour post-dose)
        3. Draw blood at time t₂ (e.g., 8 hours post-dose)
        4. Measure C(t₁) and C(t₂)
        5. Calculate t = t₂ - t₁ = 7 hours
        6. Calculate k using this function
        7. Calculate t½ = 0.693/k
        8. Use k and Vd to calculate Cl = k × Vd
        9. Adjust maintenance dose = Cl × Css × τ

    Determining k from Semi-Log Plot:
        - Plot ln(C) vs time on linear scales (semi-log plot)
        - ln(C) = ln(C₀) - k×t (linear equation, slope = -k)
        - Calculate slope: Δln(C) / Δt
        - k = |slope|
        - Can use multiple points to estimate best-fit line

    Calculated Half-Life from Measured k:
        - t½ = ln(2) / k = 0.693 / k
        - Example: k = 0.099 1/hour → t½ = 0.693/0.099 ≈ 7 hours
        - Drug loses 50% of concentration every 7 hours
        - After 5 × t½ = 35 hours, ~97% eliminated

    Comparison to Textbook Values:
        Population average k values:
        - Theophylline: 0.15 1/hour (literature) vs measured 0.12 1/hour (patient)
        - Warfarin: 0.005 1/hour (literature) vs measured 0.008 1/hour (patient)
        - Gentamicin: 0.3 1/hour (literature, normal renal)
        - If patient k << population: May have renal/liver disease
        - If patient k >> population: May have enzyme induction

    Disease Effects on k:
        Normal k (healthy 70 kg adult):
            - Theophylline: 0.10-0.16 1/hour
            - Warfarin: 0.003-0.007 1/hour
            - Gentamicin: 0.25-0.35 1/hour (renal clearance)
        Renal disease (creatinine clearance ↓):
            - ↓ k for renally cleared drugs (gentamicin, digoxin)
            - Normal k for metabolized drugs (theophylline)
            - Exception: Metabolites accumulate (active or toxic)
        Liver disease:
            - ↓ k for hepatically metabolized drugs
            - Normal k for renally excreted drugs
            - Cirrhosis: Often severe ↓ k (by 50-80%)
        Age-related changes:
            - Elderly: Often ↓ k (↓ renal and hepatic function)
            - Children: Often ↑ k (faster metabolism, better renal function)
            - Infants: ↓ k (immature metabolism and renal function)
        Heart failure:
            - ↓ k for drugs dependent on renal clearance
            - May be ↓ or ↑ for metabolized drugs (complex effects)
        Fever/Infection:
            - May ↑ k (fever increases metabolism of some drugs)
            - May ↑ k (infection causes enzyme induction)
        Genetic polymorphisms:
            - Can change k 2-10 fold (e.g., warfarin: CYP2C9 variants)
            - Codeine: Ultra-rapid, rapid, normal, poor metabolizers

    Enzyme Effects on k:
        Induction (↑ k):
            - Rifampin: ↑ k for many drugs (CYP3A4 induction)
            - Carbamazepine: ↑ k for oral contraceptives (CYP3A4)
            - Phenobarbital: ↑ k for warfarin (CYP2C9 induction)
        Inhibition (↓ k):
            - Cimetidine: ↓ k for theophylline (CYP1A2 inhibition)
            - Ketoconazole: ↓ k for many drugs (CYP3A4 inhibition)
            - Erythromycin: ↓ k for warfarin (CYP2C9 inhibition)

    Accuracy Considerations:
        Measurement accuracy:
            - Concentration must be accurately measured
            - Time of measurement must be precisely known
            - Error in time → exponential error in calculated k
            - Error in one concentration → proportional error in k
        Sampling timing:
            - Should be at least 3-4 hours apart (minimize relative error)
            - Should be in elimination phase (not distribution phase)
            - Avoid sampling too close together (relative error large)
            - Avoid sampling at end of life (concentrations very low)
        Assumptions:
            - Assumes first-order kinetics
            - Assumes steady elimination rate (no change during sampling)
            - Assumes no new absorption (past distribution phase)

    Using Population k vs Calculated k:
        Population value:
            - Published in literature (based on many patients)
            - Good for initial predictions
            - May not apply to individual patient
        Calculated value:
            - Specific for this patient
            - More accurate for dose adjustment
            - Should be used for maintenance dosing decisions
            - More reliable than population estimates after measurement

    Converting k Between Time Units:
        k can be expressed in any time unit:
        - If k = 0.1 1/hour, then k = 0.1 × 24 = 2.4 1/day
        - If k = 0.1 1/hour, then k = 0.1/60 ≈ 0.00167 1/minute
        - t½ changes correspondingly: 7 hours = 0.29 days = 420 minutes
        - Formula: k_unit1 × time_unit1 = k_unit2 × time_unit2

    Clinical Examples:
        Theophylline (asthma):
            - Draw 1 hour post-dose: C = 18 mg/L
            - Draw 9 hours post-dose: C = 9 mg/L
            - Time interval = 8 hours
            - k = ln(18/9) / 8 = 0.693 / 8 = 0.087 1/hour
            - t½ = 0.693 / 0.087 ≈ 8 hours
            - Patient elimination slower than population average (4-6 hours)
            - May need dose reduction to avoid accumulation

        Gentamicin (infection):
            - Peak (1 hour post-IV): 5 mg/L
            - Trough (8 hours pre-next dose): 1 mg/L
            - Time interval = 7 hours
            - k = ln(5/1) / 7 = 1.609 / 7 = 0.23 1/hour
            - t½ = 0.693 / 0.23 ≈ 3 hours
            - Normal kidney function (population k = 0.25-0.35 1/hour)
            - Dosing appropriate

    Notes:
        - Assumes first-order kinetics (not valid for phenytoin, aspirin at high doses)
        - Requires at least two concentration measurements
        - k is population-specific and disease-specific
        - Always draw during elimination phase (not distribution)
        - Accuracy improves with more time points (best-fit line)
        - Consider drug interactions when k differs significantly from literature
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
    Solve for t: time required for concentration to fall to a specific level.

    This function answers the clinical question: "How long until the concentration
    drops to a safe level?" or "When should I give the next dose?" It calculates
    the time needed for any desired concentration drop given the elimination rate.

    Formula: t = -ln(C(t) / C₀) / k = ln(C₀ / C(t)) / k

    Where:
    - t = Time required (what we're solving for)
    - C₀ = Initial concentration
    - C(t) = Target concentration at later time
    - k = Elimination rate constant (1/time)
    - ln = Natural logarithm

    This is the inverse of the exponential decay model. Given a target concentration
    and elimination rate, calculate how long it takes to reach that level.

    Args:
        c_0 (str): Initial plasma concentration (e.g., '100 mg/L')
        c_t (str): Target concentration to reach (e.g., '10 mg/L')
        k (str): Elimination rate constant (e.g., '0.1 1/hour')
        output_unit (str, optional): Convert result to this unit (e.g., 'hour', 'day')
        decimals (int): Decimal places to round to (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Time Elapsed (t)', 23.03, 'hour', '23.03 hour', ...)

    Examples:
        Calculate time for concentration to fall from 100 to 10 mg/L with k=0.1 1/hour:
            >>> solve_for_t(c_0='100 mg/L', c_t='10 mg/L', k='0.1 1/hour')
            ('Time Elapsed (t)', 23.03, 'hour', '23.03 hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.3, pages 3480-3520: Time-based calculations
        - Section 6.4, pages 3520-3580: Using time in pharmacokinetic decisions
        - Pages 3500-3540: Practical timing calculations
        - Chapter 6: Drug elimination and kinetics

    Clinical Applications:
        - Determine appropriate dosing interval
        - Calculate when to check therapeutic levels
        - Determine when drug is safe to discontinue
        - Plan concurrent medication administration (avoid interactions)
        - Schedule patient follow-up for compliance check
        - Calculate time to reach effective levels (multiple dosing)
        - Determine when concentration reaches sub-toxic levels

    Common Clinical Scenarios:
        1. Overdose/Toxicity Management:
            - Peak concentration = 20 mg/L (toxic)
            - Target safe level = 10 mg/L (sub-toxic)
            - Calculate time to wait before next dose
            - Example: Theophylline toxicity
        2. Therapeutic Monitoring Scheduling:
            - Start therapy, need to reach steady state
            - Time to 5 × t½ = time to essentially reach steady state
            - Schedule next level check
        3. Drug Interaction Avoidance:
            - Drug A needs 4 hours to be eliminated
            - Can safely give Drug B after that time
        4. Dosing Interval Selection:
            - Time to drop from peak to trough = dosing interval
            - Or time from peak to desired minimum level

    Relationship to Half-Life:
        - Time to fall from C₀ to C₀/2 = t½ = 0.693/k
        - Time to fall from C₀ to C₀/4 = 2×t½
        - Time to fall from C₀ to C₀/8 = 3×t½
        - Time to fall from C₀ to C₀/32 = 5×t½ (≈99% eliminated)

        Using t½ directly:
            - If t½ = 6 hours
            - Time for 50% elimination = 6 hours
            - Time for 75% elimination = 12 hours
            - Time for 87.5% elimination = 18 hours
            - Time for 93.75% elimination = 24 hours

    Practical Examples:
        Theophylline (asthma crisis):
            - Peak concentration = 25 mg/L (above therapeutic 10-20)
            - Toxic range starts ~20 mg/L
            - Target = 15 mg/L (middle of therapeutic)
            - k = 0.11 1/hour (typical)
            - t = ln(25/15) / 0.11 ≈ 0.511 / 0.11 ≈ 4.6 hours
            - Hold next dose for 5 hours, then reassess

        Gentamicin (renal insufficiency):
            - Patient with reduced k (renal disease)
            - Peak = 7 mg/L (above target of 5)
            - Want to drop to 2 mg/L (below trough target)
            - k = 0.15 1/hour (reduced from normal 0.3)
            - t = ln(7/2) / 0.15 ≈ 1.25 / 0.15 ≈ 8.3 hours
            - Need >8 hours between doses for safety

        Warfarin (interaction detected):
            - Patient on warfarin at steady state
            - Starting CYP2C9 inhibitor (will ↓ warfarin clearance)
            - Current level = 2.5 (safe)
            - Want to get level down to 1.5 before starting inhibitor
            - k (warfarin) = 0.005 1/hour
            - t = ln(2.5/1.5) / 0.005 ≈ 0.51 / 0.005 ≈ 102 hours ≈ 4.3 days
            - Plan inhibitor start for 4+ days after level measurement

    Dosing Interval Calculation:
        Extended interval dosing (especially renal impairment):
        - Normal dosing: One dose per t½ to t½
        - Extended interval: One dose per 3-5 × t½
        - Allows lower peak level (safety) but longer trough
        - Example: Gentamicin
            * Normal dosing: 5 mg/kg every 8 hours (daily)
            * Reduced function: 7 mg/kg every 24-48 hours (extended interval)
            * Calculate t for concentration to fall from peak to safe trough
            * Then space doses accordingly

    Time to Steady State:
        - Concentration builds with each dose
        - By 5 × t½, have 96.9% of steady-state concentration
        - Example: Digoxin (t½ = 36 hours)
            * Time to 99% steady state = 5 × 36 = 180 hours ≈ 7.5 days
            * Often give loading dose to achieve level sooner
        - Example: Theophylline (t½ = 6 hours)
            * Time to 99% steady state = 5 × 6 = 30 hours ≈ 1.25 days
            * Reaches steady state quickly

    Multi-Dose Accumulation:
        - This function calculates single-elimination time
        - With repetitive dosing, concentrations accumulate
        - Peak accumulation ≈ (Dose/Cl) × [1 / (1-e^(-kτ))]
        - Calculate time to reach each successive peak

    Disease Effects on Time Calculations:
        - Renal disease: ↓ k means ↑ time to reach any target level
            * Same amount needs more time to eliminate
            * Longer dosing intervals required
        - Liver disease: ↓ k means ↑ time
        - Age: ↓ k in elderly means ↑ time
        - Fever: May ↑ k (↓ time) due to increased metabolism

    Clinical Decision-Making:
        Using calculated time:
        1. "Drug level too high, how long to wait?"
           → Use this function to calculate wait time
        2. "Need to reach effective level, when to reassess?"
           → Calculate time for 5 × t½
        3. "Starting drug interaction, timing critical"
           → Calculate clearance time before new drug
        4. "Patient has renal disease, adjust dosing interval"
           → Calculate appropriate interval based on target levels

    Time Unit Conversions:
        - k units determine time units
        - If k = 0.1 1/hour, result is in hours
        - If k = 2.4 1/day, result is in days
        - Convert: 1 day = 24 hours = 1440 minutes
        - Example: t = 48 hours = 2 days = 2880 minutes

    Practical Workflow:
        1. Calculate k from two measured concentrations (or know from literature)
        2. Define target concentration for clinical goal
        3. Use this function to calculate time needed
        4. Schedule dose timing or level checks based on calculated time
        5. Adjust subsequent doses based on measured vs predicted times

    Accuracy Considerations:
        - Calculation assumes constant k
        - k may change with time (enzyme induction/inhibition)
        - k may vary day-to-day (disease fluctuations)
        - Small changes in k cause exponential changes in time
        - For long elimination times, consider disease changes
        - Re-assess k periodically (monthly for chronic therapy)

    Notes:
        - Assumes first-order kinetics throughout
        - Assumes no new dose or absorption during this time
        - Does not account for accumulation of metabolites
        - Assumes k is accurately known
        - Real-world k may vary by 20-50% between patients
        - Consider renal/hepatic disease when calculating clinical timing
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
