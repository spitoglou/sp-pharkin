"""
Core pharmacokinetics calculations for dose, distribution, and elimination.

These functions implement fundamental pharmacokinetic relationships from
"Pharmacokinetics" by Philip Rowe, particularly Chapter 5 (Distribution) and
Chapter 6 (Elimination).
"""

from .lib import format_output, generic_a_eq_b_x_c
from pint import UnitRegistry
import math

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def volume_of_distribution_weight(**kwargs):
    """
    Calculate volume of distribution adjusted for body weight.

    Volume of distribution (Vd) is the theoretical volume in which a drug is distributed
    to produce the observed plasma concentration after absorption. It reflects the extent
    of drug distribution to tissues relative to plasma.

    Formula: Vd = Mean_Vd × Weight

    Also used to scale population parameters to individual patients:
        Vd_patient = Vd_population × (Patient_weight / Reference_weight)

    Args (provide exactly 2 of 3):
        volume_of_distribution (str): Calculated volume for this patient (e.g., '50 L')
        mean_volume_of_distribution (str): Population average (e.g., '0.7 L/kg')
        weight (str): Patient body weight (e.g., '70 kg')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Volume of Distribution', 50.0, 'liter', '50.0 liter', ...)

    Examples:
        Calculate Vd for 75 kg patient with mean Vd of 0.7 L/kg:
            >>> volume_of_distribution_weight(mean_volume_of_distribution='0.7 L/kg',
            ...                               weight='75 kg')
            ('Volume of Distribution', 52.5, 'liter', '52.5 liter', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 5.1-5.3, pages 3100-3200: Volume of distribution definition and calculation
        - Table 5.1: Population Vd values for common drugs
        - Chapter 5: Distribution pharmacokinetics

    Clinical Interpretation:
        - Vd < 0.1 L/kg: Confined to intravascular space (e.g., heparin)
        - Vd 0.1-0.4 L/kg: Primarily plasma and interstitial fluid
        - Vd 0.4-1.0 L/kg: Distributed throughout body water
        - Vd > 1.0 L/kg: Highly tissue-bound (e.g., digoxin, chloroquine)

    Notes:
        - Vd is NOT a real physical volume, but a mathematical concept
        - Reflects both plasma protein binding and tissue distribution
        - Important for calculating loading doses
        - Varies with disease state (edema, dehydration, ascites)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("volume_of_distribution", False)
    b = kwargs.get("mean_volume_of_distribution", False)
    c = kwargs.get("weight", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Volume of Distribution", "Mean Volume of Distribution", "Weight"]
    )

    return format_output(quantity, string, output_unit, decimals)


def dose_concentration_volume(**kwargs):
    """
    Calculate dose, concentration, or volume relationship in body.

    This fundamental relationship connects three core pharmacokinetic variables:
    how much drug (dose), how concentrated it is (plasma concentration), and
    over what volume it's distributed (volume of distribution).

    Formula: Dose = Concentration × Volume (after complete absorption/distribution)

    More generally: Amount in body = Concentration × Vd

    This is the basis for calculating loading doses and initial concentrations.

    Args (provide exactly 2 of 3):
        dose (str): Amount of drug administered (e.g., '500 mg')
        concentration (str): Plasma concentration achieved (e.g., '10 mg/L')
        volume (str): Volume of distribution (e.g., '50 L')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the pair provided.
        Example: ('Dose', 500.0, 'milligram', '500.0 milligram', ...)

    Examples:
        Calculate dose needed for 10 mg/L concentration in 50 L volume:
            >>> dose_concentration_volume(concentration='10 mg/L', volume='50 L')
            ('Dose', 500.0, 'milligram', '500.0 milligram', ...)

        Calculate resulting concentration from 500 mg in 50 L:
            >>> dose_concentration_volume(dose='500 mg', volume='50 L')
            ('Concentration', 10.0, 'milligram / liter', '10.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 5.1, page 3100+: Fundamental concentration relationship
        - Section 6.1, page 3400+: Dose-concentration relationships
        - Chapters 5-6: Distribution principles

    Clinical Significance:
        - Used to calculate initial IV bolus doses
        - Basis for therapeutic drug monitoring
        - Essential for understanding peak concentrations
        - Used in drug interactions analysis

    Notes:
        - Assumes drug distributes instantaneously to Vd (after IV bolus)
        - For oral drugs, account for bioavailability: Absorbed dose = IV dose × F
        - Accounts for all drug in body, both bound and free
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("dose", False)
    b = kwargs.get("concentration", False)
    c = kwargs.get("volume", False)

    string, quantity = generic_a_eq_b_x_c(a, b, c, ["Dose", "Concentration", "Volume"])

    return format_output(quantity, string, output_unit, decimals)


def target_concentration(min, max):
    """
    Calculate target (therapeutic) concentration from therapeutic window.

    The therapeutic window (or therapeutic range) defines the plasma concentration
    range where a drug is expected to be clinically effective with acceptable
    safety. Target concentration is typically the midpoint of this range for
    most dosing calculations.

    Formula: Target = (Minimum + Maximum) / 2

    Args:
        min (str): Lower bound of therapeutic range (e.g., '50 mg/L')
        max (str): Upper bound of therapeutic range (e.g., '150 mg/L')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Target Concentration', 100.0, 'milligram / liter', '100.0 mg/L', ...)

    Examples:
        Midpoint between 50 and 150 mg/L:
            >>> target_concentration('50 mg/L', '150 mg/L')
            ('Target Concentration', 100.0, 'milligram / liter', '100.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 11.1-11.2, pages 5000+: Therapeutic drug monitoring
        - Table 11.1: Therapeutic ranges for common drugs
        - Chapter 11: Clinical applications and dosing

    Therapeutic Windows (Examples):
        - Theophylline: 10-20 mg/L (wide window, variable metabolism)
        - Warfarin: Monitored by INR, not concentration
        - Digoxin: 0.8-2.0 ng/mL (narrow window, toxicity risk)
        - Phenytoin: 10-20 mg/L (non-linear kinetics)
        - Gentamicin: Peak 5-10 mg/L, trough <2 mg/L (nephrotoxicity)
        - Lithium: 0.6-1.2 mEq/L (very narrow, toxicity risk)
        - Aminophylline: 10-20 mg/L (similar to theophylline)

    Notes:
        - Target varies by indication and patient factors
        - Some drugs need different targets for different conditions
        - Genetic factors affect ideal concentration (pharmacogenomics)
        - Disease states (renal failure, liver disease) may shift therapeutic range
        - Drug interactions can shift effective therapeutic level
    """
    min_q = Q_(min)
    max_q = Q_(max)
    result = (min_q + max_q) / 2  # type: ignore[operator]

    return (
        "Target Concentration",
        result.magnitude,
        "{!s}".format(result.units),
        "{!s}".format(result),
        result,
    )


def rate_of_elimination_mass_k(**kwargs):
    """
    Calculate drug elimination rate, mass, or elimination rate constant.

    The elimination rate describes how quickly a drug is removed from the body.
    In first-order kinetics (most drugs), the rate is proportional to the amount
    of drug present - the higher the concentration, the faster it's eliminated.

    Formula: Elimination Rate = Mass × K

    Where:
        Elimination Rate = Amount of drug eliminated per unit time (e.g., mg/hour)
        Mass = Amount of drug in body (e.g., mg)
        K = First-order elimination rate constant (e.g., 1/hour)

    Also: dC/dt = -K × C (rate of concentration change over time)

    Args (provide exactly 2 of 3):
        rate_of_elimination (str): Amount eliminated per unit time (e.g., '10 mg/hour')
        mass (str): Amount of drug in body (e.g., '500 mg')
        K (str): Elimination rate constant (e.g., '0.02 1/hour')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the pair provided.
        Example: ('Elimination Rate', 10.0, 'milligram / hour', '10.0 mg/hour', ...)

    Examples:
        Calculate elimination rate for 500 mg with K=0.02 1/hour:
            >>> rate_of_elimination_mass_k(mass='500 mg', K='0.02 1/hour')
            ('Elimination Rate', 10.0, 'milligram / hour', '10.0 mg/hour', ...)

        Calculate K from rate and mass:
            >>> rate_of_elimination_mass_k(rate_of_elimination='10 mg/hour', mass='500 mg')
            ('Elimination Rate Constant(K)', 0.02, '1 / hour', '0.02 1/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.2-6.3, pages 3400-3500: Elimination kinetics
        - Section 6.4, page 3550+: First-order elimination mathematics
        - Chapter 6: Drug elimination processes

    Relationship to Other Parameters:
        K = 0.693 / t½ (half-life relationship)
        Clearance = K × Vd
        Half-life = 0.693 / K
        t = ln(C₀/C) / K (time to reach concentration C)

    Clinical Interpretation:
        - Higher K → Faster elimination → Shorter half-life → More frequent dosing
        - Lower K → Slower elimination → Longer half-life → Less frequent dosing
        - K unchanged by dose amount (first-order kinetics property)
        - K may change in disease (renal failure increases K for renally cleared drugs)

    Notes:
        - Assumes first-order kinetics (linear: rate ∝ concentration)
        - Some drugs show zero-order kinetics at therapeutic doses (e.g., phenytoin, ethanol)
        - K varies between individuals due to genetics, disease, age
        - Measured from plasma concentration vs time curve (slope = -K)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("rate_of_elimination", False)
    b = kwargs.get("mass", False)
    c = kwargs.get("K", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Rate of Elimination", "Mass", "Elimination Rate Constant(K)"]
    )

    return format_output(quantity, string, output_unit, decimals)


def half_life_k(**kwargs):
    """
    Calculate half-life or elimination rate constant relationship.

    Half-life (t½) is the time required for plasma concentration (or total body
    amount) of a drug to decrease by 50%. It's determined exclusively by the
    elimination rate constant (K) and is independent of dose.

    Formula: ln(2) = K × t½
    Rearranged: t½ = ln(2) / K = 0.693 / K
    Also: K = ln(2) / t½ = 0.693 / t½

    Clinical Significance of Half-life:
        - Determines dosing frequency for maintenance doses
        - Time to reach 50% steady state: 1 t½
        - Time to reach 87.5% steady state: 3 t½
        - Time to reach 93.75% steady state: 4 t½
        - Time to reach 96.9% steady state: 5 t½
        - Time to eliminate 87.5% of dose: 3 t½

    Args (provide exactly 2 of 3):
        K (str): Elimination rate constant (e.g., '0.1 1/hour')
        half_life (str): Half-life of drug (e.g., '6.9 hour')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Half-Life', 6.93, 'hour', '6.93 hour', ...)

    Examples:
        Calculate half-life from K:
            >>> half_life_k(K='0.1 1/hour')
            ('Half-Life', 6.93, 'hour', '6.93 hour', ...)

        Calculate K from half-life:
            >>> half_life_k(half_life='6.93 hour')
            ('Elimination Rate Constant(K)', 0.1, '1 / hour', '0.1 1/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.3, pages 3470-3490: Half-life and elimination kinetics
        - Section 6.4, pages 3500-3550: Mathematical relationships
        - Table 6.1: Half-lives of common drugs
        - Chapter 6: Elimination processes

    Half-life Examples (from Table 6.1):
        - Very short: Penicillin (0.5 hr), Acetaminophen (2-3 hr)
        - Short: Theophylline (3-4 hr), Gentamicin (2-3 hr)
        - Moderate: Warfarin (35-45 hr), Phenytoin (20-30 hr)
        - Long: Digoxin (36-40 hr), Chloroquine (weeks)
        - Very long: Vitamin A (several months), DDT (years)

    Dosing Interval Guidelines:
        - Standard: Dosing interval ≈ t½ (maintains 50-100% of Css)
        - More frequent: q(t½/2) (maintains 75-100%)
        - Less frequent: q(2×t½) (maintenance but lower levels)
        - Extended interval: q(5×t½) (used to accumulate to Css faster)

    Clinical Implications:
        - Short t½ drugs need frequent dosing (e.g., penicillin q4-6h)
        - Long t½ drugs allow once-daily dosing (e.g., warfarin once daily)
        - Long t½ drugs show significant accumulation (e.g., digoxin)
        - Starting therapy: Load then maintain, or just maintain
        - Stopping therapy: Takes 5 t½ to eliminate from system

    Notes:
        - t½ is constant for a given drug in a given patient (independent of dose)
        - t½ changes in disease affecting clearance (renal/hepatic disease)
        - Genetic polymorphisms can change t½ by 2-10 fold (e.g., warfarin, codeine)
        - Age affects t½ for some drugs (clearance decreases with age)
        - First-order kinetics: t½ independent of concentration
        - Zero-order kinetics: t½ depends on concentration (non-linear behavior)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = math.log(2)
    b = kwargs.get("K", False)
    c = kwargs.get("half_life", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Ln(2)", "Elimination Rate Constant(K)", "Half-Life"]
    )

    return format_output(quantity, string, output_unit, decimals)


def extraction_rate(**kwargs):
    """
    Calculate extraction ratio, concentration difference, or input concentration.

    The extraction ratio (E) represents the fraction of drug extracted/eliminated
    by an organ during a single pass through its blood supply. It's a measure of
    organ efficiency in removing drug from blood.

    Formula: Extraction Ratio = (Concentration_in - Concentration_out) / Concentration_in

    Or: C_diff = E × C_in

    Where:
        E = Extraction Ratio (0-1, or 0-100%)
        C_in = Concentration entering organ (arterial)
        C_out = Concentration leaving organ (venous)
        C_diff = C_in - C_out (amount extracted)

    Interpretation:
        - E = 0.2 (20%) → Low extraction → Drug poorly removed per pass
        - E = 0.5 (50%) → Moderate extraction → Significant removal
        - E = 0.8 (80%) → High extraction → Very efficient removal
        - E → 1.0 (>90%) → Very high extraction → Limited by blood flow

    Relationship to Clearance:
        Clearance = Q × E  (where Q = organ blood flow)
        - High E, low Q → Metabolism-limited clearance
        - Low E, high Q → Flow-limited clearance
        - E approaching 1.0 → Clearance limited by blood flow

    Args (provide exactly 2 of 3):
        c_diff (str): Difference in concentration (e.g., '5 mg/L')
        E (str): Extraction ratio, 0-1 (e.g., '0.2')
        c_in (str): Input (arterial) concentration (e.g., '10 mg/L')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Extraction Ratio(E)', 0.5, 'dimensionless', '0.5', ...)

    Examples:
        Calculate extraction ratio from concentrations:
            >>> extraction_rate(c_diff='5 mg/L', c_in='10 mg/L')
            ('Extraction Ratio(E)', 0.5, 'dimensionless', '0.5', ...)

        Calculate concentration drop across organ:
            >>> extraction_rate(E='0.5', c_in='10 mg/L')
            ('Concentration Difference(C_in - C_out)', 5.0, 'mg/L', '5.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 7.1-7.2, pages 3700-3800: Hepatic and renal clearance
        - Section 7.3, pages 3800-3900: Extraction ratio and clearance
        - Table 7.1: Extraction ratios for common drugs and organs
        - Chapter 7: Organ-specific elimination

    Extraction Ratios for Common Drugs (Liver):
        - High E (>0.5): Propranolol, morphine, lidocaine, nitroglycerin
        - Moderate E (0.2-0.5): Phenytoin, warfarin, diazepam
        - Low E (<0.2): Theophylline, acetaminophen, aspirin

    Organs and Their Extraction Ratios (Range):
        - Liver: 0.1-0.9 (varies by drug)
        - Kidney: 0.5-1.0 (for actively secreted drugs)
        - Lung: 0.5-0.9 (for some drugs on first pass)
        - Brain: Very low E (blood-brain barrier)

    First-Pass Metabolism:
        - Drugs with high hepatic E show first-pass effect (reduced oral bioavailability)
        - Example: Propranolol (E=0.7) has F=0.2-0.3 (only 20-30% reaches systemic circulation)
        - Solution: Give IV or use alternative routes (sublingual, transdermal)

    Clinical Applications:
        - Predicts hepatic disease impact on clearance
        - Used in calculating drug interactions
        - Important for assessing first-pass metabolism
        - Guides dosing for patients with liver/kidney disease

    Notes:
        - E cannot exceed 1.0 (can't extract more than present)
        - E depends on both intrinsic clearance and organ blood flow
        - Disease affecting liver/kidney function changes E
        - Enzyme induction increases E (increases clearance)
        - Enzyme inhibition decreases E (decreases clearance, increases levels)
        - Age-related changes in extractive organs affect E
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("c_diff", False)
    b = kwargs.get("E", False)
    c = kwargs.get("c_in", False)

    string, quantity = generic_a_eq_b_x_c(
        a,
        b,
        c,
        [
            "Concentration Difference(C_in - C_out)",
            "Extraction Ratio(E)",
            "Input Concentration(C_in)",
        ],
    )

    return format_output(quantity, string, output_unit, decimals)
