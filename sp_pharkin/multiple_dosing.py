"""
Multiple IV dosing pharmacokinetics calculations for steady-state concentrations.

These functions implement formulas from "Pharmacokinetics" by Philip Rowe,
Chapter 10 (Multiple Dosing and Steady State), focusing on intermittent IV
bolus dosing at regular intervals.

Key concepts:
- Css,max: Peak concentration at steady state (immediately after dose)
- Css,min: Trough concentration at steady state (just before next dose)
- Fluctuation: Degree of concentration variation between doses
"""

from .lib import format_output
import math
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def css_max(**kwargs):
    """
    Calculate peak (maximum) concentration at steady state for multiple IV doses.

    At steady state with repeated IV bolus dosing, Css,max represents the peak
    plasma concentration achieved immediately after each dose. This is when
    the drug level is highest during each dosing interval.

    The physiological basis: After each IV bolus dose, the drug concentration
    spikes to its maximum. At steady state, drug remaining from previous doses
    adds to each new dose, creating a predictable peak that balances elimination
    over the dosing interval. The peak concentration determines the risk of
    concentration-dependent toxicity while ensuring adequate therapeutic effect.

    Formula: Css,max = (D/V) * 1/(1 - e^(-K*tau))

    Or equivalently: Css,max = C0 * R

    Where:
        D = Dose administered (e.g., mg)
        V = Volume of distribution (e.g., L)
        K = Elimination rate constant (e.g., 1/hour)
        tau = Dosing interval (e.g., hours)
        C0 = Initial concentration from single dose = D/V
        R = Accumulation factor = 1/(1 - e^(-K*tau))

    The term (D/V) represents the concentration immediately after a single dose,
    and 1/(1 - e^(-K*tau)) is the accumulation factor that accounts for drug
    remaining from previous doses.

    Args (provide exactly 4 of 5):
        css_max (str): Peak steady-state concentration (e.g., '20 mg/L')
        dose (str): Dose administered per dosing interval (e.g., '500 mg')
        volume (str): Volume of distribution (e.g., '50 L')
        K (str): Elimination rate constant (e.g., '0.1 1/hour')
        dosing_interval (str): Time between doses, tau (e.g., '8 hour')

    Optional:
        output_unit (str): Convert output to different unit
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Peak Steady-State Concentration (Css,max)', 20.0, 'mg/L', '20.0 mg/L', ...)

    Examples:
        Calculate Css,max from dose, volume, K, and dosing interval:
            >>> css_max(dose='500 mg', volume='50 L', K='0.1 1/hour',
            ...         dosing_interval='8 hour')
            ('Peak Steady-State Concentration (Css,max)', 18.18, ...)

        Calculate required dose for target Css,max:
            >>> css_max(css_max='20 mg/L', volume='50 L', K='0.1 1/hour',
            ...         dosing_interval='8 hour')
            ('Dose', 550.0, 'milligram', '550.0 milligram', ...)

        Calculate volume from known Css,max, dose, K, and interval:
            >>> css_max(css_max='20 mg/L', dose='550 mg', K='0.1 1/hour',
            ...         dosing_interval='8 hour')
            ('Volume of Distribution', 50.0, 'liter', '50.0 liter', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 10: Multiple Dosing and Steady State
        - Section 10.3.2: Peak and trough concentrations (pages 99-100)
        - Pages 5651-5714: Mathematical derivation and clinical examples
        - Figure 10.4: Average, Maximum and Minimum concentrations at steady state
        - Appendix: Derivation of Css,max = D/V x 1/(1 - e^(-K*tau))

    Drug-Specific Examples and Therapeutic Ranges:

        Aminoglycosides (Gentamicin, Tobramycin, Amikacin):
            The classic case where Css,max monitoring is essential. These drugs
            exhibit concentration-dependent killing - higher peaks = better
            bactericidal effect. However, prolonged high levels cause toxicity.

            Gentamicin/Tobramycin (conventional dosing, q8h):
                - Target Css,max: 5-10 mg/L (varies by infection severity)
                - Mild-moderate infection: 5-6 mg/L peak acceptable
                - Severe/life-threatening: 8-10 mg/L peak required
                - Typical dose: 1.5-2 mg/kg q8h
                - Vd: ~0.25 L/kg (distributes to extracellular fluid)
                - K: varies with renal function (normal: ~0.3 h^-1)

            Gentamicin/Tobramycin (extended-interval/once-daily dosing):
                - Target Css,max: 15-25 mg/L (intentionally high)
                - Allows drug-free interval for renal recovery
                - Dose: 5-7 mg/kg once daily
                - Requires normal renal function

            Amikacin:
                - Target Css,max: 20-30 mg/L (conventional)
                - Target Css,max: 56-64 mg/L (extended-interval)
                - Dose: 15 mg/kg/day divided q8-12h

            Example calculation (from book, page 5699):
                Patient with K = 0.03 h^-1 (poor renal function)
                Dose = 80 mg, Vd = 20 L, tau = 8 hours
                Css,max = 80/20 x 1/(1 - e^(-0.03 x 8))
                       = 4 x 1/(1 - 0.787)
                       = 4 x 4.69 = 18.76 mg/L
                This is toxic! Requires extended interval (see book section 10.3.3)

        Vancomycin:
            While trough monitoring is standard, peaks may be monitored in
            some settings (though less commonly now with AUC-based dosing).
            - Historical target Css,max: 25-40 mg/L (1-2 hours post-infusion)
            - Modern practice: AUC/MIC targeting (400-600 for MRSA)
            - Dose: 15-20 mg/kg q8-12h
            - Vd: ~0.4-1.0 L/kg
            - K: depends on renal function

        Theophylline:
            Narrow therapeutic index drug requiring careful peak monitoring.
            - Therapeutic range: 10-20 mg/L (some sources: 5-15 mg/L)
            - Toxic effects begin >20 mg/L (seizures, arrhythmias)
            - Often given as aminophylline (salt factor S = 0.8)
            - Vd: ~0.48 L/kg (35 L for 70 kg patient)
            - K: ~0.04-0.08 h^-1 (highly variable)
            - Half-life: ~9 hours (range 4-16, affected by smoking, age, disease)

            Css,max considerations:
                - Smoking: Increases K by 50-100% (lower peaks for same dose)
                - Heart failure: Decreases K by 50% (higher peaks, toxicity risk)
                - Liver disease: Decreases K significantly
                - Slow-release formulations preferred to minimize fluctuation

        Digoxin:
            For digoxin, average concentration (Css,av) is typically sufficient
            rather than peak/trough monitoring. However, understanding Css,max
            is relevant for timing of blood samples.
            - Therapeutic range: 0.8-2.0 microgram/L (Css,av)
            - Target mid-range: 1.4 microgram/L
            - Oral bioavailability: ~70% (tablets)
            - Very long half-life: ~42 hours
            - Blood samples taken at least 6 hours post-dose (distribution phase)
            - Toxicity: >2.0 microgram/L associated with arrhythmias, GI symptoms

        Phenytoin:
            Exhibits saturable (non-linear) kinetics - peak monitoring complex.
            - Therapeutic range: 10-20 mg/L (total), 1-2 mg/L (free)
            - Toxic effects >20 mg/L (nystagmus, ataxia, confusion)
            - Non-linear kinetics mean small dose changes cause large peak changes
            - Vd: ~0.65 L/kg
            - This equation applies only at low concentrations (first-order region)

        Lithium:
            Used in bipolar disorder with narrow therapeutic index.
            - Therapeutic range: 0.6-1.2 mEq/L (maintenance)
            - Acute mania: 0.8-1.2 mEq/L
            - Toxic: >1.5 mEq/L
            - Samples taken 12 hours post-dose (standardized trough)
            - Half-life: ~18-24 hours

    Clinical Interpretation of Css,max:

        Peak Too High (Toxicity Risk):
            - Aminoglycosides: Ototoxicity (vestibular/cochlear damage)
            - Theophylline: Seizures, tachyarrhythmias, nausea/vomiting
            - Vancomycin: Infusion reactions ("red man syndrome")
            - Phenytoin: Nystagmus, ataxia, cognitive impairment

            Actions when peak too high:
            1. Reduce dose size (maintains interval)
            2. Extend dosing interval (maintains dose)
            3. Switch to sustained-release formulation
            4. Consider alternative agent

        Peak Too Low (Therapeutic Failure Risk):
            - Aminoglycosides: Inadequate bactericidal effect
            - Theophylline: Uncontrolled bronchospasm
            - Vancomycin: Treatment failure, resistance development

            Actions when peak too low:
            1. Increase dose size
            2. Verify drug is being administered correctly
            3. Check for drug interactions affecting Vd or K
            4. Assess volume status (edema increases Vd)

    Disease State Effects on Css,max:

        Renal Impairment:
            - Aminoglycosides: K decreases dramatically
            - With normal dose/interval: Css,max increases (toxicity)
            - With reduced GFR: Extend interval, not just reduce dose
            - Book example: K = 0.03 h^-1 patient required 48-hour intervals

        Hepatic Impairment:
            - Theophylline: K decreases (reduced metabolism)
            - Phenytoin: Free fraction increases (reduced protein binding)
            - Expect higher Css,max for same dose
            - Reduce dose by 50% or more in severe cirrhosis

        Heart Failure:
            - Reduces organ blood flow and clearance
            - Theophylline: K decreases by ~50%
            - Digoxin: Vd may decrease (less tissue perfusion)
            - Loading doses may need reduction

        Critical Illness/Sepsis:
            - Increased Vd (capillary leak, fluid resuscitation)
            - May need higher initial doses for adequate peaks
            - Variable K (organ dysfunction)
            - Frequent TDM essential

        Burns:
            - Dramatically increased Vd for aminoglycosides
            - May need 2-3x normal doses for adequate peaks
            - Altered protein binding

        Obesity:
            - Aminoglycosides: Use ideal or adjusted body weight
            - Theophylline: Use ideal body weight
            - Actual body weight overestimates Vd for hydrophilic drugs

        Pregnancy:
            - Increased Vd (expanded plasma volume)
            - Increased K (enhanced renal clearance)
            - May need higher doses for equivalent peaks

        Pediatrics:
            - Neonates: Reduced K (immature renal/hepatic function)
            - Children: Often higher K than adults (faster metabolism)
            - Age-specific dosing essential

        Elderly:
            - Reduced K (declining renal function)
            - May have reduced Vd (less lean body mass)
            - Start with lower doses, monitor closely

    Relationship to Other Parameters:

        Connection to Css,min (trough):
            - Css,max = Css,min + D/V (fundamental relationship)
            - The dose increment (D/V) is the difference between peak and trough
            - At steady state: amount eliminated = dose administered

        Connection to Css,avg (average):
            - Css,max = Css,avg * K * tau / (1 - e^(-K*tau))
            - Css,max > Css,avg (always, due to post-dose spike)
            - The larger the fluctuation, the greater the difference

        Peak-to-Trough Ratio:
            - Css,max / Css,min = e^(K*tau)
            - Ratio depends only on K and tau, not on dose
            - Ratio = 2 when tau = t_half (one half-life between doses)
            - Ratio = 4 when tau = 2 * t_half

        Accumulation Factor (R):
            - R = 1/(1 - e^(-K*tau)) = Css,max / C0
            - R approaches 1 when tau >> t_half (no accumulation)
            - R approaches infinity when tau << t_half (extensive accumulation)
            - R = 2 when tau = t_half

        Connection to Half-Life:
            - t_half = 0.693 / K
            - Steady state reached after 4-5 half-lives
            - Longer half-life = more accumulation for same interval

    When to Use This Calculation Clinically:

        Primary Indications:
            1. Aminoglycoside dosing optimization
            2. Theophylline toxicity assessment
            3. Designing dosing regimens for new patients
            4. Adjusting doses after measured levels

        This Function Is Appropriate When:
            - Drug is given as IV bolus (instantaneous administration)
            - Drug follows one-compartment model
            - First-order elimination kinetics apply
            - Steady state has been reached (4-5 half-lives on therapy)

        Limitations (When NOT to Use):
            - Drugs given as IV infusion (different kinetics)
            - Two-compartment drugs during distribution phase
            - Non-linear (saturable) kinetics (e.g., high-dose phenytoin)
            - Before steady state is reached
            - Extravascular administration (oral, IM) - use different equations

    Common Clinical Scenarios:

        Scenario 1: New Aminoglycoside Patient
            Given: Weight 70 kg, normal renal function
            Goal: Calculate initial gentamicin regimen for Css,max 8 mg/L
            Approach:
                - Estimate Vd = 0.25 L/kg x 70 kg = 17.5 L
                - Estimate K = 0.3 h^-1 (normal renal function)
                - Choose tau = 8 hours (standard q8h dosing)
                - Solve for dose to achieve target Css,max

        Scenario 2: Renally Impaired Patient
            Given: CrCl 30 mL/min, current regimen causing high peaks
            Goal: Adjust dosing to achieve safe Css,max
            Approach:
                - Recalculate K based on reduced renal function
                - Consider extending interval (book recommends this approach)
                - May need 24-48 hour intervals for severe impairment
                - Book example: 110 mg q48h achieved targets with K = 0.03 h^-1

        Scenario 3: Theophylline Toxicity
            Given: Patient with Css,max 25 mg/L (toxic), on q12h dosing
            Goal: Reduce Css,max to 15 mg/L
            Approach:
                - Calculate dose reduction needed
                - Consider slow-release formulation to reduce fluctuation
                - Evaluate contributing factors (smoking cessation, new drugs)

        Scenario 4: Loading Dose Calculation
            Given: Need to achieve therapeutic Css,max immediately
            Goal: Calculate appropriate loading dose
            Approach:
                - LD = Target Css,max x Vd / F
                - Loading dose gets patient to steady-state levels immediately
                - Follow with maintenance doses to sustain levels

    Mathematical Notes:

        Derivation (from book Appendix, page 6155):
            At steady state, the decline during the dosing interval exactly equals
            the rise caused by the next dose.
            - Rise from dose = D/V
            - Decline during interval = Css,max x (1 - e^(-K*tau))
            - At steady state: Css,max x (1 - e^(-K*tau)) = D/V
            - Rearranging: Css,max = D/V x 1/(1 - e^(-K*tau))

        Units Consideration:
            - K (1/time) x tau (time) must be dimensionless
            - Always ensure units cancel in the exponential term
            - Example: K = 0.03 h^-1 x tau = 8 h = 0.24 (dimensionless)

    Notes:
        - Assumes instantaneous IV bolus administration
        - Assumes complete distribution before next dose
        - First-order elimination kinetics throughout
        - Steady state reached after approximately 4-5 half-lives
        - For drugs with significant distribution phase, true peak may be higher
          immediately after dose, then decline during distribution
        - Peak sampling for aminoglycosides: 30-60 min post-dose (after distribution)
        - This equation has restricted applicability but is essential for
          aminoglycoside dosing in clinical practice
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    css_max_val = kwargs.get("css_max", False)
    dose = kwargs.get("dose", False)
    volume = kwargs.get("volume", False)
    k = kwargs.get("K", False)
    tau = kwargs.get("dosing_interval", False)

    # Count provided parameters
    params = [css_max_val, dose, volume, k, tau]
    params_count = sum([bool(p) for p in params])

    if params_count != 4:
        raise ValueError(
            f"exactly 4 of 5 parameters required (css_max, dose, volume, K, dosing_interval). "
            f"Got {params_count}."
        )

    # Calculate K * tau (dimensionless)
    if k and tau:
        k_tau = k * tau
        if hasattr(k_tau, "units"):
            k_tau_val = k_tau.to("dimensionless").magnitude
        else:
            k_tau_val = float(k_tau)
    else:
        k_tau_val = None

    # Accumulation factor: R = 1 / (1 - e^(-K*tau))
    if k_tau_val is not None:
        accumulation_factor = 1.0 / (1.0 - math.exp(-k_tau_val))

    if not css_max_val:
        # Calculate Css,max = (D/V) * R
        c0 = dose / volume  # Initial concentration from single dose
        quantity = c0 * accumulation_factor
        string = "Peak Steady-State Concentration (Css,max)"

    elif not dose:
        # Calculate dose: D = Css,max * V / R
        quantity = css_max_val * volume / accumulation_factor
        string = "Dose"

    elif not volume:
        # Calculate volume: V = D / (Css,max / R) = D * R / Css,max
        quantity = dose * accumulation_factor / css_max_val
        string = "Volume of Distribution"

    elif not k:
        # Calculate K from Css,max, dose, volume, tau
        # Css,max = (D/V) / (1 - e^(-K*tau))
        # 1 - e^(-K*tau) = (D/V) / Css,max
        # e^(-K*tau) = 1 - (D/V) / Css,max
        # -K*tau = ln(1 - (D/V) / Css,max)
        # K = -ln(1 - (D/V) / Css_max) / tau
        c0 = dose / volume
        ratio = c0 / css_max_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)

        if ratio_val >= 1.0:
            raise ValueError(
                "Invalid parameters: D/V must be less than Css,max for valid K calculation"
            )

        k_val = -math.log(1.0 - ratio_val) / tau
        quantity = k_val
        string = "Elimination Rate Constant (K)"

    elif not tau:
        # Calculate tau from Css,max, dose, volume, K
        # Css,max = (D/V) / (1 - e^(-K*tau))
        # 1 - e^(-K*tau) = (D/V) / Css,max
        # e^(-K*tau) = 1 - (D/V) / Css_max
        # -K*tau = ln(1 - (D/V) / Css_max)
        # tau = -ln(1 - (D/V) / Css_max) / K
        c0 = dose / volume
        ratio = c0 / css_max_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)

        if ratio_val >= 1.0:
            raise ValueError(
                "Invalid parameters: D/V must be less than Css,max for valid tau calculation"
            )

        tau_val = -math.log(1.0 - ratio_val) / k
        quantity = tau_val
        string = "Dosing Interval (tau)"

    else:
        raise ValueError("Cannot determine which parameter to solve for")

    return format_output(quantity, string, output_unit, decimals)


def css_min(**kwargs):
    """
    Calculate trough (minimum) concentration at steady state for multiple IV doses.

    At steady state with repeated IV bolus dosing, Css,min represents the trough
    plasma concentration just before the next dose is administered. This is when
    the drug level is lowest during each dosing interval.

    The physiological basis: After each dose, drug is continuously eliminated
    from the body. At steady state, the trough represents the lowest point in
    the concentration cycle - the moment when the next dose is due. The trough
    concentration is critical for assessing both therapeutic adequacy (is there
    enough drug present?) and toxicity risk (has the drug accumulated?).

    Formulas:
        Css,min = Css,max * e^(-K*tau)
        Or: Css,min = Css,max - D/V

    Where:
        Css,max = Peak steady-state concentration
        K = Elimination rate constant
        tau = Dosing interval
        D = Dose administered
        V = Volume of distribution

    The relationship Css,min = Css,max - D/V reflects that the difference between
    peak and trough is exactly the concentration increment from one dose (D/V).
    This is a fundamental relationship that always holds at steady state.

    Args (provide exactly 2 of 3 for primary formula, or 3 of 4 for alternative):
        css_min (str): Trough steady-state concentration (e.g., '5 mg/L')
        css_max (str): Peak steady-state concentration (e.g., '20 mg/L')
        K (str): Elimination rate constant (e.g., '0.1 1/hour')
        dosing_interval (str): Time between doses, tau (e.g., '8 hour')

        Alternative parameters (if css_max not provided):
        dose (str): Dose administered (e.g., '500 mg')
        volume (str): Volume of distribution (e.g., '50 L')

    Optional:
        output_unit (str): Convert output to different unit
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Trough Steady-State Concentration (Css,min)', 8.18, ...)

    Examples:
        Calculate Css,min from Css,max, K, and dosing interval:
            >>> css_min(css_max='20 mg/L', K='0.1 1/hour', dosing_interval='8 hour')
            ('Trough Steady-State Concentration (Css,min)', 8.99, ...)

        Calculate Css,min using dose and volume (Css,min = Css,max - D/V):
            >>> css_min(css_max='20 mg/L', dose='500 mg', volume='50 L')
            ('Trough Steady-State Concentration (Css,min)', 10.0, ...)

        Calculate Css,max from known Css,min:
            >>> css_min(css_min='9 mg/L', K='0.1 1/hour', dosing_interval='8 hour')
            ('Peak Steady-State Concentration (Css,max)', 20.0, ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 10: Multiple Dosing and Steady State
        - Section 10.3.2: Peak and trough concentrations (pages 99-100)
        - Section 10.3.3: Aminoglycoside use with poor renal function (page 100)
        - Pages 5706-5714: Css,min = Css,max - D/V derivation and examples
        - Figure 10.4: Average, Maximum and Minimum concentrations at steady state
        - Appendix (page 6185): Derivation of Css,min = Css,max - D/V

    Drug-Specific Examples and Therapeutic Ranges:

        Aminoglycosides (Gentamicin, Tobramycin, Amikacin):
            The trough is CRITICAL for aminoglycoside monitoring because:
            - Nephrotoxicity correlates with sustained trough elevation
            - Ototoxicity (irreversible) correlates with prolonged exposure
            - Goal: Allow drug-free interval for renal tubular cell recovery

            Gentamicin/Tobramycin (conventional dosing, q8h):
                - Target Css,min: <2 mg/L (strict interpretation)
                - Many clinicians accept <1 mg/L as safer target
                - Troughs >2 mg/L: Significantly increased nephrotoxicity risk
                - Typical patient: Css,min should be 0.5-1.5 mg/L

            Gentamicin/Tobramycin (extended-interval/once-daily dosing):
                - Target Css,min: <1 mg/L or undetectable
                - Goal is near-complete elimination before next dose
                - Allows 8-12 hour drug-free interval
                - Preferred in patients with good renal function

            Amikacin:
                - Target Css,min: <5 mg/L (conventional)
                - Target Css,min: <2.5 mg/L (extended-interval)
                - Less nephrotoxic than gentamicin at equivalent doses

            Book example (page 5710):
                Css,max = 18.76 mg/L, D/V = 4 mg/L
                Css,min = 18.76 - 4 = 14.76 mg/L
                This is dangerously high! Patient needs extended interval.

        Vancomycin:
            Vancomycin is the classic example where TROUGH monitoring is standard.
            - Therapeutic Css,min: 10-20 mg/L (traditional targets)
            - Severe MRSA infections: 15-20 mg/L trough recommended
            - Less serious infections: 10-15 mg/L may suffice
            - Css,min <10 mg/L: Risk of treatment failure, resistance
            - Css,min >20 mg/L: Increased nephrotoxicity risk

            Modern AUC-based monitoring:
                - AUC/MIC ratio 400-600 preferred for MRSA
                - Troughs of 15-20 mg/L approximate this target
                - Two-level sampling increasingly used

            Timing of trough sample:
                - Within 30 minutes before next dose
                - Must be at true steady state (after 4th dose typically)

            Loading dose consideration:
                - Vancomycin has long half-life (~6-12 hours)
                - Loading dose (25-30 mg/kg) recommended for serious infections
                - First trough meaningful only after steady state reached

        Theophylline:
            For theophylline, both peak and trough matter due to narrow index.
            - Therapeutic range: 10-20 mg/L (entire range)
            - Typical Css,min target: >10 mg/L (avoid breakthrough symptoms)
            - Trough <5 mg/L: Poor bronchodilation, symptoms likely
            - Slow-release formulations minimize peak-trough difference

            Nighttime symptoms consideration:
                - If Css,min too low overnight, patient may wake with bronchospasm
                - May need bedtime dosing or sustained-release at night
                - Extended-release theophylline helps maintain adequate troughs

        Digoxin:
            Unique considerations for digoxin sampling:
            - Therapeutic range: 0.8-2.0 microgram/L (Css,av typically used)
            - Trough samples: Take at least 6 hours post-dose
            - The 6-hour delay is due to prolonged distribution phase
            - Earlier samples show falsely elevated (distribution-phase) levels
            - Very long half-life (~42 hours) means minimal fluctuation

            Book note (page 7079):
                "With digoxin, we do not generally need to concern ourselves with
                the peak and trough levels. It is sufficient to control Css,av
                within a suitable range."

        Phenytoin:
            Trough monitoring preferred due to complex kinetics.
            - Therapeutic range: 10-20 mg/L (total), 1-2 mg/L (free)
            - Trough samples standard practice
            - Saturable kinetics mean trough less predictable
            - Small dose changes can cause large trough changes
            - Albumin levels affect free fraction (adjust for hypoalbuminemia)

        Lithium:
            Standardized trough monitoring at 12 hours post-dose.
            - Therapeutic range: 0.6-1.2 mEq/L (maintenance)
            - Acute mania: 0.8-1.2 mEq/L
            - Toxic: >1.5 mEq/L (narrow margin!)
            - 12-hour post-dose sample is standard worldwide
            - Consistency in timing essential for interpretation

        Carbamazepine:
            - Therapeutic range: 4-12 mg/L
            - Trough monitoring preferred
            - Auto-induction: K increases over 2-4 weeks of therapy
            - Initial troughs may differ from steady-state troughs

        Valproic Acid:
            - Therapeutic range: 50-100 mg/L
            - Trough monitoring standard
            - Highly protein-bound (90%); free levels in renal/hepatic disease

    Clinical Interpretation of Css,min:

        Trough Too High (Toxicity/Accumulation Risk):
            - Aminoglycosides: Nephrotoxicity, ototoxicity (irreversible)
            - Vancomycin: Nephrotoxicity (combined with aminoglycosides = worse)
            - Theophylline: Persistent nausea, tachycardia, seizure risk
            - Phenytoin: Ataxia, nystagmus, cognitive impairment
            - Lithium: Tremor, polyuria, encephalopathy

            Clinical actions when trough too high:
            1. EXTEND dosing interval (most effective for aminoglycosides)
            2. Reduce dose size (maintains interval)
            3. Hold dose(s) until trough in range
            4. Assess renal function (often declining)
            5. Check for new interacting medications

        Trough Too Low (Sub-therapeutic/Failure Risk):
            - Aminoglycosides: Infection not controlled (but toxicity still possible)
            - Vancomycin: Treatment failure, emergence of resistance (VRSA risk)
            - Theophylline: Breakthrough bronchospasm, especially nocturnal
            - Antiepileptics: Breakthrough seizures
            - Lithium: Mood destabilization

            Clinical actions when trough too low:
            1. Increase dose size
            2. Shorten dosing interval (within limits)
            3. Verify adherence
            4. Check for drug interactions (enzyme induction)
            5. Assess volume status (increased Vd dilutes drug)

    Disease State Effects on Css,min:

        Renal Impairment (CRITICAL for many drugs):
            - Aminoglycosides: Dramatically prolonged half-life
                K normal: ~0.3 h^-1 (t_half ~2-3 hours)
                K severe CKD: ~0.03 h^-1 (t_half ~23 hours)
            - Vancomycin: Half-life increases from 6 to 24+ hours
            - Digoxin: Reduced renal clearance (70% of elimination)
            - Result: Css,min rises progressively if interval not extended

            Book approach for aminoglycosides in renal failure (Section 10.3.3):
                - Calculate patient's half-life based on K
                - Determine how many half-lives needed for peak-to-trough decline
                - Set interval accordingly (may need 24-48 hours or longer)
                - Example: K = 0.03 h^-1, need Css,max/Css,min = 4
                - This requires 2 half-lives = 2 x 23 = 46 hours between doses

        Hepatic Impairment:
            - Theophylline: Reduced metabolism, higher troughs
            - Phenytoin: Reduced protein binding, more free drug
            - Lidocaine: Significantly reduced clearance
            - Adjust for Child-Pugh class (reduce doses 30-50%)

        Heart Failure:
            - Reduced organ perfusion affects clearance
            - Theophylline: Clearance reduced ~50%
            - Lidocaine: Reduced volume of distribution and clearance
            - Higher troughs for same dose

        Critical Illness:
            - Highly variable kinetics
            - Fluid shifts alter Vd
            - Organ dysfunction alters clearance
            - Frequent monitoring essential

        Hypoalbuminemia:
            - Phenytoin, valproic acid: Increased free fraction
            - Total trough may appear therapeutic but free is toxic
            - Use free drug levels or apply correction formulas

        Obesity:
            - Use ideal/adjusted body weight for Vd estimates
            - Actual weight overestimates Vd for hydrophilic drugs
            - May result in higher-than-expected troughs

        Pediatrics:
            - Neonates: Extended half-lives (immature elimination)
            - Children: Often shorter half-lives than adults
            - Age-appropriate dosing intervals essential

        Elderly:
            - Reduced renal function (even with "normal" creatinine)
            - Estimate CrCl using Cockcroft-Gault
            - Generally expect higher troughs for same regimen

    Relationship to Other Parameters:

        Connection to Css,max (peak):
            - Css,min = Css,max - D/V (fundamental at steady state)
            - Css,min = Css,max * e^(-K*tau) (decay relationship)
            - The dose increment (D/V) connects peak and trough

        Peak-to-Trough Ratio:
            - Css,max / Css,min = e^(K*tau)
            - If tau = t_half: ratio = 2 (50% decline between doses)
            - If tau = 2 x t_half: ratio = 4 (75% decline)
            - Extended-interval aminoglycosides: ratio may be 10-20+

        Connection to Fluctuation:
            - Fluctuation = (Css,max - Css,min) / Css,avg x 100%
            - Low Css,min relative to Css,max means high fluctuation
            - Controlled-release formulations raise Css,min (reduce fluctuation)

        Relationship to Half-Life:
            - Longer half-life = less decline between doses = higher Css,min
            - t_half = 0.693 / K
            - For same dose and interval, longer t_half = higher Css,min

    When to Use This Calculation Clinically:

        Primary Indications for Trough Monitoring:
            1. Vancomycin - Standard of care for all patients
            2. Aminoglycosides - Essential for toxicity prevention
            3. Lithium - Standardized 12-hour trough
            4. Antiepileptics - Phenytoin, carbamazepine, valproic acid
            5. Theophylline - If on chronic therapy

        This Function Is Appropriate When:
            - Drug is given as IV bolus (instantaneous administration)
            - Drug follows one-compartment model
            - First-order elimination kinetics apply
            - Steady state has been reached

        When to Sample for Trough:
            - Aminoglycosides: Within 30 minutes before next dose
            - Vancomycin: Within 30 minutes before 4th or 5th dose
            - Digoxin: At least 6 hours post-dose (distribution)
            - Phenytoin: Anytime at steady state (minimal fluctuation)
            - Lithium: Exactly 12 hours post-dose (standardized)

    Common Clinical Scenarios:

        Scenario 1: Aminoglycoside Trough Elevated
            Measured: Css,min = 3.5 mg/L (target <2 mg/L)
            Current: Gentamicin 120 mg q8h
            Approach:
                - Calculate how much trough needs to fall (3.5 -> 1.5)
                - Extend interval to q12h or q24h (preferred over dose reduction)
                - Recheck trough before 3rd dose on new regimen

        Scenario 2: Vancomycin Trough Subtherapeutic
            Measured: Css,min = 8 mg/L (target 15-20 mg/L for MRSA)
            Current: Vancomycin 1g q12h
            Approach:
                - Increase dose to 1.25-1.5g q12h
                - OR shorten interval to q8h (if renal function adequate)
                - Recheck trough at new steady state

        Scenario 3: Renal Function Declining
            Baseline: CrCl 80 mL/min, now CrCl 40 mL/min
            Current: Aminoglycoside regimen producing Css,min 1.8 mg/L
            Expected: Css,min will rise as K decreases
            Approach:
                - Recalculate K based on new renal function
                - Adjust interval before troughs become dangerous
                - Consider switch to extended-interval if not already

        Scenario 4: Calculating Interval for Target Trough
            Given: Want Css,max = 8 mg/L, Css,min = 1 mg/L
            Ratio: Css,max / Css,min = 8
            Approach:
                - e^(K*tau) = 8
                - K*tau = ln(8) = 2.08
                - For K = 0.1 h^-1: tau = 20.8 hours (~24h in practice)

    Mathematical Notes:

        Two Equivalent Formulas:
            1. Css,min = Css,max * e^(-K*tau)
               - Shows exponential decay from peak to trough
               - Useful when you know Css,max and kinetics

            2. Css,min = Css,max - D/V
               - Shows that dose increment (D/V) is peak-trough difference
               - Useful when you know both concentrations and dose

        Derivation (from book Appendix, page 6193):
            At steady state, the rise caused by each dose exactly equals the
            decline during the dosing interval:
            - Css,max = Css,min + D/V (rise equals increment)
            - Rearranging: Css,min = Css,max - D/V

    Notes:
        - Trough level is the most commonly monitored parameter for many drugs
        - Timing of sample is critical - must be true pre-dose
        - Longer dosing interval relative to half-life produces lower trough
        - Extended-interval dosing deliberately creates very low troughs
        - Trough helps assess accumulation, especially in renal impairment
        - Always interpret trough in context of peak and clinical response
        - Some drugs (digoxin) require delayed sampling due to distribution
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    css_min_val = kwargs.get("css_min", False)
    css_max_val = kwargs.get("css_max", False)
    k = kwargs.get("K", False)
    tau = kwargs.get("dosing_interval", False)
    dose = kwargs.get("dose", False)
    volume = kwargs.get("volume", False)

    # Calculate K * tau if available
    if k and tau:
        k_tau = k * tau
        if hasattr(k_tau, "units"):
            k_tau_val = k_tau.to("dimensionless").magnitude
        else:
            k_tau_val = float(k_tau)
        decay_factor = math.exp(-k_tau_val)
    else:
        decay_factor = None

    # Determine which formula to use and what to solve for

    # Formula 1: Css,min = Css,max * e^(-K*tau)
    if css_max_val and k and tau and not css_min_val:
        quantity = css_max_val * decay_factor
        string = "Trough Steady-State Concentration (Css,min)"

    elif css_min_val and k and tau and not css_max_val:
        # Css,max = Css,min / e^(-K*tau) = Css,min * e^(K*tau)
        quantity = css_min_val / decay_factor
        string = "Peak Steady-State Concentration (Css,max)"

    elif css_min_val and css_max_val and tau and not k:
        # e^(-K*tau) = Css,min / Css,max
        # -K*tau = ln(Css,min / Css,max)
        # K = -ln(Css,min / Css_max) / tau
        ratio = css_min_val / css_max_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)

        if ratio_val <= 0 or ratio_val >= 1:
            raise ValueError(
                "Invalid parameters: Css,min must be positive and less than Css,max"
            )

        k_val = -math.log(ratio_val) / tau
        quantity = k_val
        string = "Elimination Rate Constant (K)"

    elif css_min_val and css_max_val and k and not tau:
        # e^(-K*tau) = Css,min / Css,max
        # -K*tau = ln(Css,min / Css_max)
        # tau = -ln(Css,min / Css_max) / K
        ratio = css_min_val / css_max_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)

        if ratio_val <= 0 or ratio_val >= 1:
            raise ValueError(
                "Invalid parameters: Css,min must be positive and less than Css,max"
            )

        tau_val = -math.log(ratio_val) / k
        quantity = tau_val
        string = "Dosing Interval (tau)"

    # Formula 2: Css,min = Css,max - D/V
    elif css_max_val and dose and volume and not css_min_val:
        c0 = dose / volume
        quantity = css_max_val - c0
        string = "Trough Steady-State Concentration (Css,min)"

    elif css_min_val and dose and volume and not css_max_val:
        c0 = dose / volume
        quantity = css_min_val + c0
        string = "Peak Steady-State Concentration (Css,max)"

    elif css_min_val and css_max_val and volume and not dose:
        # D/V = Css,max - Css,min
        # D = (Css,max - Css,min) * V
        quantity = (css_max_val - css_min_val) * volume
        string = "Dose"

    elif css_min_val and css_max_val and dose and not volume:
        # D/V = Css,max - Css,min
        # V = D / (Css,max - Css,min)
        quantity = dose / (css_max_val - css_min_val)
        string = "Volume of Distribution"

    else:
        raise ValueError(
            "Invalid parameter combination. Provide: "
            "(css_max, K, dosing_interval) to calculate css_min, or "
            "(css_min, K, dosing_interval) to calculate css_max, or "
            "(css_max, dose, volume) to calculate css_min, or similar valid combinations."
        )

    return format_output(quantity, string, output_unit, decimals)


def fluctuation(**kwargs):
    """
    Calculate the extent of concentration fluctuation at steady state.

    Fluctuation quantifies the degree of variation between peak and trough
    concentrations during a dosing interval at steady state. It is expressed
    as a percentage relative to the average steady-state concentration.

    The physiological basis: With intermittent dosing, drug concentrations
    inevitably rise and fall. The extent of this "swing" between peak and
    trough matters clinically - too much fluctuation means the patient
    experiences both potentially toxic peaks and sub-therapeutic troughs.
    Fluctuation helps quantify this problem and guides decisions about
    dose frequency and formulation choice.

    Formula: Fluctuation = ((Css,max - Css,min) / Css,avg) * 100%

    Simplified (since Css,max - Css,min = D/V at steady state):
        Fluctuation = (D/V) / Css,avg * 100%

    Where:
        Css,max = Peak steady-state concentration
        Css,min = Trough steady-state concentration
        Css,avg = Average steady-state concentration
        D = Dose
        V = Volume of distribution

    The term (D/V) represents the concentration increment from each dose,
    which equals the peak-to-trough difference at steady state.

    Interpretation:
        - Very low fluctuation (<25%): Near-constant levels (infusion-like)
        - Low fluctuation (25-50%): Well-controlled variation
        - Moderate fluctuation (50-100%): Typical for most drugs
        - High fluctuation (100-200%): Significant swings
        - Very high fluctuation (>200%): Extreme variation (extended-interval dosing)

    Args (provide appropriate parameters for your calculation):
        fluctuation (str): Fluctuation percentage (e.g., '80' for 80%)
        css_max (str): Peak steady-state concentration (e.g., '20 mg/L')
        css_min (str): Trough steady-state concentration (e.g., '8 mg/L')
        css_avg (str): Average steady-state concentration (e.g., '12.5 mg/L')
        dose (str): Dose administered (e.g., '500 mg')
        volume (str): Volume of distribution (e.g., '50 L')

    Optional:
        output_unit (str): Convert output to different unit (for non-% results)
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Fluctuation', 96.0, 'percent', '96.0 percent', ...)

    Examples:
        Calculate fluctuation from Css,max, Css,min, and Css,avg:
            >>> fluctuation(css_max='20 mg/L', css_min='8 mg/L', css_avg='12.5 mg/L')
            ('Fluctuation', 96.0, 'percent', '96.0 percent', ...)

        Calculate fluctuation using D/V and Css,avg:
            >>> fluctuation(dose='500 mg', volume='50 L', css_avg='12.5 mg/L')
            ('Fluctuation', 80.0, 'percent', '80.0 percent', ...)

        Calculate Css,avg from known fluctuation and peak/trough:
            >>> fluctuation(fluctuation='96', css_max='20 mg/L', css_min='8 mg/L')
            ('Average Steady-State Concentration (Css,avg)', 12.5, ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 10: Multiple Dosing and Steady State
        - Section 10.6: Extent of fluctuation in drug concentrations (pages 103-107)
        - Section 10.6.1: Half-life effects on fluctuation
        - Section 10.6.2: Dose division to reduce fluctuation
        - Section 10.6.3: Slow release dosage forms
        - Section 10.6.4: Practical steps to restrict fluctuation
        - Figure 10.7: Effect of dose division on fluctuation
        - Figure 10.8: Fast and slow release formulations
        - Figure 10.9: Control of fluctuation by dose division vs. slow release

    Drug-Specific Examples and Therapeutic Considerations:

        Theophylline:
            Classic example where fluctuation control is essential.
            - Therapeutic range: 10-20 mg/L (narrow!)
            - High fluctuation causes:
                * Peak toxicity: nausea, tachycardia, seizures (>20 mg/L)
                * Trough failure: breakthrough bronchospasm (<10 mg/L)

            Without fluctuation control:
                - Half-life ~9 hours (highly variable)
                - Once-daily dosing: Extreme fluctuation (peaks toxic, troughs ineffective)
                - Book quote (page 6061): "This is why drugs such as theophylline
                  are made available in slow release formulations."

            Strategies for theophylline (from book Section 10.6.4):
                1. Increase dose division:
                   - 200 mg once daily -> 100 mg q12h -> 50 mg q6h
                   - Each division reduces fluctuation
                   - q6h dosing achieves smooth envelope but poor compliance

                2. Use slow-release formulations:
                   - 100 mg SR twice daily matches 50 mg IR four times daily
                   - Better patient compliance with equivalent control
                   - Book Figure 10.9 demonstrates this equivalence

            Target fluctuation for theophylline:
                - Ideally <50% to stay within 10-20 mg/L range
                - Slow-release formulations essential for most patients

        Aminoglycosides (Gentamicin, Tobramycin, Amikacin):
            Unique case where HIGH fluctuation is deliberately sought.

            Conventional q8h dosing:
                - Moderate fluctuation (~100-150%)
                - Peak 5-10 mg/L, Trough <2 mg/L
                - Fluctuation = (10-2)/6 x 100 = 133%

            Extended-interval (once-daily) dosing:
                - VERY HIGH fluctuation (>500%)
                - Peak 15-25 mg/L, Trough <1 mg/L (often undetectable)
                - Fluctuation = (25-0.5)/8 x 100 = 306%

            Why high fluctuation is beneficial for aminoglycosides:
                1. Concentration-dependent killing: Higher peaks = better efficacy
                2. Post-antibiotic effect: Killing continues after drug gone
                3. Toxicity time-dependent: Low troughs allow renal recovery
                4. Net result: Better efficacy + less toxicity with high fluctuation

            Book example (pages 5715-5720):
                Patient with K = 0.03 h^-1 needed 48-hour dosing interval
                This creates extreme fluctuation but achieves both:
                - Adequate peak for efficacy (5-10 mg/L)
                - Low trough for safety (<2 mg/L)

        Digoxin:
            Very low fluctuation due to long half-life.
            - Half-life: ~42 hours
            - Once-daily dosing: Minimal fluctuation
            - Fluctuation typically <25%
            - Peak and trough nearly equal
            - Css,avg (0.8-2.0 microgram/L) is sufficient target

            Book note (page 7079): "With digoxin, we do not generally need to
            concern ourselves with the peak and trough levels."

        Phenytoin:
            Complicated by non-linear kinetics.
            - At therapeutic levels: Low fluctuation (long apparent half-life)
            - Near saturation: Fluctuation increases unpredictably
            - Controlled-release formulations available (Dilantin Kapseals)
            - Fluctuation less relevant than maintaining trough in range

        Lithium:
            Moderate fluctuation acceptable.
            - Half-life: 18-24 hours
            - Once or twice daily dosing standard
            - Fluctuation ~50-75% with once daily
            - Slow-release formulations reduce GI side effects and fluctuation

        Antihypertensives (Beta-blockers, CCBs):
            Fluctuation affects 24-hour BP control.
            - High fluctuation: BP controlled only part of day
            - Extended-release formulations preferred for:
                * Once-daily convenience
                * Smooth 24-hour BP control
                * Reduced morning surge risk

        Opioid Analgesics:
            Fluctuation affects pain control and side effects.
            - Immediate-release: High fluctuation
                * Peak: Sedation, respiratory depression, euphoria
                * Trough: Pain breakthrough
            - Extended-release: Low fluctuation
                * Smoother pain control
                * More stable respiratory function
                * Less abuse potential (debated)

    Clinical Interpretation and Decision-Making:

        When High Fluctuation Is Problematic:

            Narrow Therapeutic Index Drugs:
                - Theophylline: 10-20 mg/L range too narrow for high fluctuation
                - Phenytoin: 10-20 mg/L with toxicity just above
                - Lithium: 0.6-1.2 mEq/L with toxicity >1.5
                - Warfarin: Small changes in level = major INR changes
                - Digoxin: 0.8-2.0 microgram/L (though long t_half helps)

            Peak-Related Toxicity:
                - Aminoglycosides: Ototoxicity (though high peaks desired)
                - Vancomycin: Infusion reactions at high concentrations
                - Theophylline: Seizures, arrhythmias at high peaks
                - Some antipsychotics: Extrapyramidal effects at peaks

            Trough-Related Failure:
                - Antibiotics: Sub-MIC levels allow bacterial regrowth
                - Antiepileptics: Breakthrough seizures in troughs
                - Bronchodilators: Nocturnal symptoms (theophylline)
                - Antihypertensives: Early morning BP surge

        When High Fluctuation Is Acceptable or Desirable:

            Concentration-Dependent Effects:
                - Aminoglycosides: Higher peak = better killing
                - Fluoroquinolones: Peak/MIC ratio matters more than time
                - Daptomycin: Concentration-dependent bactericidal activity

            Post-Antibiotic Effect:
                - Drug effect persists after concentration drops
                - Aminoglycosides, fluoroquinolones exhibit this
                - Allows drug-free interval without efficacy loss

            Toxicity Time-Dependent:
                - Low trough allows tissue recovery
                - Aminoglycoside nephrotoxicity correlates with trough
                - Extended intervals = better safety profile

            Wide Therapeutic Index:
                - Penicillins: Large safety margin
                - Cephalosporins: Generally safe even with high peaks
                - Higher fluctuation tolerable

    Disease State Effects on Optimal Fluctuation:

        Renal Impairment:
            - Extended intervals increase fluctuation
            - For aminoglycosides: This is desirable (extended-interval approach)
            - For most drugs: May need to accept higher fluctuation
            - Alternative: Use continuous infusion if fluctuation intolerable

        Hepatic Impairment:
            - Reduced K means longer half-life
            - Longer half-life naturally reduces fluctuation
            - Same dose/interval = lower fluctuation
            - May not need formulation change

        Critical Illness:
            - Highly variable pharmacokinetics
            - Continuous infusion often preferred (zero fluctuation)
            - Example: Continuous vancomycin infusion in ICU

        Malabsorption:
            - Oral slow-release may not work as designed
            - Drug released but not absorbed as intended
            - May need IV or immediate-release with frequent dosing

    Strategies to Control Fluctuation:

        1. Increase Dose Division (Book Section 10.6.2):
            From book Figure 10.7:
            - 200 mg once daily: Extreme fluctuation (toxic peaks, ineffective troughs)
            - 100 mg twice daily: Moderate fluctuation (acceptable)
            - 50 mg four times daily: Low fluctuation (best control)

            Trade-off: Compliance decreases as frequency increases
            - Once daily: >90% adherence typical
            - Twice daily: ~80% adherence
            - Four times daily: <50% adherence often

        2. Use Slow-Release Formulations (Book Section 10.6.3):
            From book Figure 10.8:
            - Same dose, different formulation
            - Slow release: Lower, later peaks; higher troughs
            - Effect: Smoother concentration profile

            Advantages:
            - Reduced fluctuation with same dosing frequency
            - Better patient compliance (fewer doses)
            - Reduced peak-related side effects

            Book quote (page 6062): "We can suppress fluctuations in concentration
            either by increasing dose division or by using a slow release dosage
            formulation. The latter is more convenient for patients."

        3. Continuous Infusion (Fluctuation = 0%):
            - Eliminates all fluctuation
            - Reserved for:
                * ICU settings (continuous monitoring possible)
                * Drugs where fluctuation is particularly harmful
                * Short-term high-intensity therapy
            - Examples: Continuous vancomycin, aminophylline, heparin

        4. Shorten Half-Life (Usually Not Practical):
            - Enzyme induction (not controllable)
            - Dialysis (only for select situations)
            - Not a routine clinical strategy

    Relationship to Dosing Interval and Half-Life:

        Mathematical Relationship:
            - Fluctuation increases as tau/t_half increases
            - tau = t_half: Peak/trough = 2:1, fluctuation ~67%
            - tau = 2 x t_half: Peak/trough = 4:1, fluctuation ~150%
            - tau = 3 x t_half: Peak/trough = 8:1, fluctuation ~233%

        The "Goldilocks" Half-Life (Book Section 10.5-10.6):
            From book (page 5902-5905):
            "Drug companies would prefer to avoid drugs with either a long half-life
            (concentrations take too long to accumulate) or a short half-life
            (concentrations fluctuate wildly). They like it to be just right!"

            - Very long half-life (>24h): Slow accumulation, needs loading dose
            - Very short half-life (<4h): High fluctuation, frequent dosing needed
            - Ideal half-life (~8-12h): Reasonable accumulation, manageable fluctuation

        Choosing Dosing Interval:
            - tau < t_half: Low fluctuation but frequent dosing
            - tau = t_half: Moderate fluctuation, reasonable dosing
            - tau > t_half: High fluctuation, infrequent dosing

    When to Calculate Fluctuation Clinically:

        Primary Indications:
            1. Designing dosing regimens for narrow therapeutic index drugs
            2. Deciding between immediate-release vs. sustained-release
            3. Troubleshooting patients with peak toxicity or trough failure
            4. Comparing different dosing strategies
            5. Educating patients about timing of side effects

        This Function Helps Answer:
            - "Why does this patient have side effects after each dose?"
            - "Why does the patient have symptoms before the next dose?"
            - "Would a sustained-release formulation help this patient?"
            - "How would changing dosing frequency affect drug levels?"

    Common Clinical Scenarios:

        Scenario 1: Theophylline Patient with Peak Toxicity
            Current: Aminophylline 300 mg q8h (S=0.8, so theophylline 240 mg q8h)
            Symptoms: Nausea and tremor 2-3 hours after each dose (peak)
            Measured: Css,max 22 mg/L, Css,min 11 mg/L, Css,avg 15 mg/L
            Fluctuation = (22-11)/15 x 100 = 73%

            Options:
            1. Switch to slow-release 240 mg q8h (reduces fluctuation)
            2. Change to 120 mg q4h (better fluctuation but poor compliance)
            3. Reduce dose and accept lower average (may lose efficacy)

        Scenario 2: Aminoglycoside Regimen Comparison
            Patient: 70 kg, normal renal function
            Option A: Gentamicin 120 mg q8h (conventional)
            Option B: Gentamicin 420 mg q24h (extended-interval)

            Option A fluctuation:
                Css,max ~8 mg/L, Css,min ~1.5 mg/L, Css,avg ~4 mg/L
                Fluctuation = (8-1.5)/4 x 100 = 163%

            Option B fluctuation:
                Css,max ~20 mg/L, Css,min ~0.5 mg/L, Css,avg ~4 mg/L
                Fluctuation = (20-0.5)/4 x 100 = 488%

            Despite higher fluctuation, Option B is preferred:
            - Higher peak for better killing
            - Longer drug-free interval for renal recovery

        Scenario 3: Converting to Slow-Release
            Current: Drug X 50 mg q6h (4 times daily)
            Goal: Reduce to twice daily dosing
            Required: Slow-release 100 mg q12h formulation

            Both achieve same Css,avg but SR has:
            - Lower Css,max (less peak toxicity)
            - Higher Css,min (less trough failure)
            - Lower fluctuation
            - Better compliance

    Mathematical Notes:

        Alternative Fluctuation Measures:
            This function uses: (Css,max - Css,min) / Css,avg x 100%

            Other definitions exist:
            - Peak-to-trough ratio: Css,max / Css,min
            - Swing: Css,max - Css,min (absolute difference)
            - % Decline: (Css,max - Css,min) / Css,max x 100%

        Relationship to Other Parameters:
            - Fluctuation = D/V / Css,avg x 100% (since Css,max - Css,min = D/V)
            - Fluctuation depends on dose, volume, and clearance
            - Larger D/V (bigger dose, smaller volume) = higher fluctuation
            - Higher Css,avg (more accumulation) = lower fluctuation

    Notes:
        - Fluctuation is dimensionless (expressed as percentage)
        - Lower fluctuation is generally preferred for patient comfort/safety
        - Exception: Aminoglycosides benefit from high fluctuation
        - Fluctuation calculation assumes steady state has been reached
        - Slow-release formulations are the most practical way to reduce fluctuation
        - Continuous infusion eliminates fluctuation entirely but requires IV access
        - Patient compliance often determines acceptable fluctuation level
        - Some fluctuation is inevitable with intermittent dosing
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    fluct = kwargs.get("fluctuation", False)
    css_max_val = kwargs.get("css_max", False)
    css_min_val = kwargs.get("css_min", False)
    css_avg_val = kwargs.get("css_avg", False)
    dose = kwargs.get("dose", False)
    volume = kwargs.get("volume", False)

    # Formula 1: Fluctuation = ((Css,max - Css,min) / Css,avg) * 100
    if css_max_val and css_min_val and css_avg_val and not fluct:
        diff = css_max_val - css_min_val
        ratio = diff / css_avg_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)
        quantity = Q_(ratio_val * 100, "percent")
        string = "Fluctuation"

    # Formula 2: Fluctuation = (D/V) / Css,avg * 100
    elif dose and volume and css_avg_val and not fluct:
        c0 = dose / volume
        ratio = c0 / css_avg_val
        if hasattr(ratio, "units"):
            ratio_val = ratio.to("dimensionless").magnitude
        else:
            ratio_val = float(ratio)
        quantity = Q_(ratio_val * 100, "percent")
        string = "Fluctuation"

    # Solve for Css,avg from fluctuation and peak/trough
    elif fluct and css_max_val and css_min_val and not css_avg_val:
        # Css,avg = (Css,max - Css,min) / (fluctuation / 100)
        diff = css_max_val - css_min_val
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        quantity = diff / (fluct_val / 100.0)
        string = "Average Steady-State Concentration (Css,avg)"

    # Solve for Css,avg from fluctuation and D/V
    elif fluct and dose and volume and not css_avg_val:
        # Css,avg = (D/V) / (fluctuation / 100)
        c0 = dose / volume
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        quantity = c0 / (fluct_val / 100.0)
        string = "Average Steady-State Concentration (Css,avg)"

    # Solve for Css,max from fluctuation, Css,min, and Css,avg
    elif fluct and css_min_val and css_avg_val and not css_max_val:
        # (Css,max - Css,min) / Css,avg = fluctuation / 100
        # Css,max - Css,min = Css,avg * fluctuation / 100
        # Css,max = Css,min + Css,avg * fluctuation / 100
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        diff = css_avg_val * (fluct_val / 100.0)
        quantity = css_min_val + diff
        string = "Peak Steady-State Concentration (Css,max)"

    # Solve for Css,min from fluctuation, Css,max, and Css,avg
    elif fluct and css_max_val and css_avg_val and not css_min_val:
        # (Css,max - Css,min) / Css,avg = fluctuation / 100
        # Css,max - Css,min = Css,avg * fluctuation / 100
        # Css,min = Css,max - Css,avg * fluctuation / 100
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        diff = css_avg_val * (fluct_val / 100.0)
        quantity = css_max_val - diff
        string = "Trough Steady-State Concentration (Css,min)"

    # Solve for dose from fluctuation, volume, and Css,avg
    elif fluct and volume and css_avg_val and not dose:
        # (D/V) / Css,avg = fluctuation / 100
        # D/V = Css_avg * fluctuation / 100
        # D = V * Css_avg * fluctuation / 100
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        quantity = volume * css_avg_val * (fluct_val / 100.0)
        string = "Dose"

    # Solve for volume from fluctuation, dose, and Css,avg
    elif fluct and dose and css_avg_val and not volume:
        # (D/V) / Css_avg = fluctuation / 100
        # D/V = Css_avg * fluctuation / 100
        # V = D / (Css_avg * fluctuation / 100)
        fluct_val = fluct.magnitude if hasattr(fluct, "magnitude") else float(fluct)
        quantity = dose / (css_avg_val * (fluct_val / 100.0))
        string = "Volume of Distribution"

    else:
        raise ValueError(
            "Invalid parameter combination. Provide: "
            "(css_max, css_min, css_avg) to calculate fluctuation, or "
            "(dose, volume, css_avg) to calculate fluctuation, or "
            "(fluctuation, css_max, css_min) to calculate css_avg, or similar valid combinations."
        )

    return format_output(quantity, string, output_unit, decimals)
