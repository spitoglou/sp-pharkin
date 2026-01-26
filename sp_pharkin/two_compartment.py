"""
Two-compartment model pharmacokinetics calculations.

The two-compartment model describes drugs that distribute between a central
compartment (blood/plasma and highly perfused organs) and a peripheral
compartment (less perfused tissues). This creates biphasic concentration-time
curves with an initial rapid distribution phase followed by a slower
elimination phase.

These functions implement formulas from "Pharmacokinetics" by Philip Rowe,
Chapter 7 (Two-Compartment Model).

How to Identify Two-Compartment Behavior:
    The practical method for distinguishing a one-compartment from a two-compartment
    drug is to inspect a semi-logarithmic graph of concentration versus time:
    - One compartment: Semi-log plot forms a simple STRAIGHT LINE
    - Two compartment: Semi-log plot shows a distinctive DOG-LEG BEND

    The bend occurs because concentrations fall faster in the initial distribution
    period (both redistribution AND elimination removing drug from blood) and slower
    in the later elimination period (elimination partially offset by drug returning
    from tissues).

Understanding the Two Phases:
    Distribution Phase (Early, governed by alpha):
        - Drug rapidly moves from blood into tissues
        - Both redistribution AND elimination deplete central compartment
        - Results in RAPID decline in blood concentration
        - Duration: typically complete within 4-5 distribution half-lives

    Elimination Phase (Late, governed by beta):
        - Equilibrium established between compartments
        - Drug returns from tissues as elimination continues
        - Results in SLOW decline in blood concentration
        - This phase determines dosing interval for most drugs

Digoxin as the Classic Two-Compartment Drug:
    Digoxin is the textbook example of a two-compartment drug with critical
    clinical implications:

    - Digoxin is polar (limited lipid solubility) and enters tissues SLOWLY
      despite cardiac muscle being well-perfused
    - The cardiac muscle (site of action) is in the SECOND compartment
    - Blood samples are from the FIRST compartment
    - Early blood levels do NOT correlate with therapeutic effect!

    Clinical Consequence:
        - A blood sample taken at 20 minutes post-dose may show very high
          (apparently toxic) levels, yet the therapeutic response is minimal
          because tissue concentrations are still low
        - It is almost impossible to interpret samples taken in the first few
          hours when blood levels are FALLING but tissue concentrations are RISING
        - From about 6 hours onwards, a steady ratio exists between blood and
          tissue concentrations, making TDM meaningful

    RULE: Blood samples for therapeutic drug monitoring of digoxin should
    NOT be taken less than 6 hours post-dose.

When One-Compartment Model is Adequate:
    More lipid-soluble drugs distribute quickly into well-perfused tissues,
    showing minimal distribution phase. For these drugs:
    - Blood levels correlate with therapeutic effect at ALL times
    - Semi-log concentration-time plot is approximately linear
    - One-compartment model provides adequate description
    - Blood sampling timing is less critical

    Examples where one-compartment may suffice:
    - Highly lipophilic drugs with rapid tissue equilibration
    - When only terminal-phase sampling is planned
    - When distribution phase is very short relative to dosing interval

Common Drugs Following Two-Compartment Kinetics:
    Cardiac Glycosides:
        - Digoxin: t1/2,alpha ~35 min, t1/2,beta ~36-48 hours
          (classic example, must wait 6+ hours for TDM)

    Local Anesthetics:
        - Lidocaine: t1/2,alpha ~8 min, t1/2,beta ~1.5-2 hours
          (rapid CNS distribution, watch for toxicity)

    Aminoglycosides:
        - Gentamicin: t1/2,alpha ~15-30 min, t1/2,beta ~2-3 hours
        - Tobramycin, Amikacin: similar profiles
        (peak/trough monitoring requires timing awareness)

    Glycopeptides:
        - Vancomycin: significant distribution phase
        (trough samples preferred for TDM)

    Benzodiazepines:
        - Diazepam: rapid CNS distribution, long terminal half-life

    Opioids:
        - Fentanyl: rapid CNS distribution, redistribution-based duration

Relationship Between Micro-constants and Hybrid Constants:
    Micro-constants (model-based, not directly measurable):
        - k12: Rate constant from central to peripheral compartment
        - k21: Rate constant from peripheral to central compartment
        - K (k10): Elimination rate constant from central compartment

    Hybrid/Macro-constants (observable from concentration-time data):
        - alpha (α): Distribution rate constant (faster, larger value)
        - beta (β): Elimination rate constant (slower, smaller value)
        - A: Intercept coefficient for distribution phase
        - B: Intercept coefficient for elimination phase

    Mathematical Relationships:
        α + β = k12 + k21 + K    (sum of roots)
        α × β = k21 × K          (product of roots)

    The hybrid constants are what we measure from the biphasic curve;
    the micro-constants describe the underlying physiological processes.

Key Parameters:
    - alpha (α): Distribution rate constant (faster, initial phase)
    - beta (β): Elimination rate constant (slower, terminal phase)
    - A, B: Intercept coefficients for each exponential phase
    - k12: Rate constant from central to peripheral compartment
    - k21: Rate constant from peripheral to central compartment
    - K (or k10): Elimination rate constant from central compartment
    - V1: Volume of central compartment
    - Vss: Volume of distribution at steady state

Reference:
    Rowe, P. Pharmacokinetics
    - Section 7.1: The two-compartment model structure
    - Section 7.2: Drug concentrations in blood (first compartment)
    - Section 7.3: Determining how many compartments a drug occupies
    - Section 7.4: Drug concentration in the second compartment
    - Section 7.5: Two-compartment systems and therapeutic drug monitoring
                   for digoxin (pages 71-72)
"""

from .lib import format_output, Q_
import math


def concentration_two_compartment(**kwargs):
    """
    Calculate plasma concentration at time t for a two-compartment model.

    After IV bolus administration, drugs following two-compartment kinetics
    show a biphasic decline in plasma concentration. The initial rapid decline
    (distribution phase, governed by α) reflects drug moving from blood to
    tissues. The slower terminal decline (elimination phase, governed by β)
    reflects drug being eliminated from the body.

    Formula: C(t) = A × e^(-α×t) + B × e^(-β×t)

    Where:
        - C(t): Concentration at time t
        - A: Intercept coefficient for distribution phase (extrapolated to t=0)
        - B: Intercept coefficient for elimination phase (extrapolated to t=0)
        - α (alpha): Distribution rate constant (faster, larger value)
        - β (beta): Elimination rate constant (slower, smaller value)
        - t: Time after administration

    At t=0: C(0) = A + B (initial concentration)

    Args:
        A (str): Distribution phase intercept coefficient (e.g., '50 mg/L')
        B (str): Elimination phase intercept coefficient (e.g., '30 mg/L')
        alpha (str): Distribution rate constant (e.g., '2.0 1/hour')
        beta (str): Elimination rate constant (e.g., '0.3 1/hour')
        t (str): Time after administration (e.g., '2 hour')
        output_unit (str, optional): Desired output unit (e.g., 'mg/L')
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Concentration at time t', 15.5, 'milligram / liter', '15.5 mg/L', ...)

    Examples:
        Calculate concentration at 2 hours:
            >>> concentration_two_compartment(
            ...     A='50 mg/L', B='30 mg/L',
            ...     alpha='2.0 1/hour', beta='0.3 1/hour',
            ...     t='2 hour'
            ... )
            ('Concentration at time t', ..., 'milligram / liter', ...)

        Digoxin concentration example (typical parameters):
            >>> # Digoxin: alpha ~1.2/hour, beta ~0.02/hour (t1/2,beta ~36 hours)
            >>> concentration_two_compartment(
            ...     A='1.5 microgram/L', B='1.0 microgram/L',
            ...     alpha='1.2 1/hour', beta='0.02 1/hour',
            ...     t='6 hour'
            ... )
            # Returns concentration after distribution phase is complete

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.1-7.3: Biphasic concentration-time curves
        - Section 7.5: Therapeutic drug monitoring considerations for digoxin

    How to Identify Two-Compartment Behavior:
        Plot ln(concentration) vs time. If the curve shows a distinct "dog-leg"
        bend rather than a straight line, the drug follows two-compartment
        kinetics. The initial steep slope reflects alpha (distribution), and
        the later gentle slope reflects beta (elimination).

    Clinical Interpretation of Parameters:
        - α >> β: Clear separation of distribution and elimination phases
        - A > B: Significant distribution into tissues (drug moves extensively
          from blood to peripheral compartment)
        - A < B: Drug remains predominantly in central compartment
        - A + B = C0: Initial concentration at time zero

    Timing and Phase Considerations:
        Distribution Phase (early times):
            - When t < 2-3 × (1/α), concentration is dominated by A × e^(-αt)
            - Blood levels fall rapidly as drug distributes to tissues
            - For digoxin: first ~2-3 hours post-dose

        Elimination Phase (later times):
            - When t > 5 × (1/α), concentration approaches B × e^(-βt)
            - Distribution is complete; decline reflects true elimination
            - For digoxin: after ~6 hours post-dose

    Digoxin-Specific Clinical Implications:
        Digoxin is the classic problem drug for two-compartment kinetics:
        - Cardiac muscle (site of action) is in the SECOND compartment
        - Blood (sampling site) is in the FIRST compartment
        - During distribution phase, blood levels are HIGH but tissue
          levels (and therapeutic effect) are still LOW
        - A sample at 20 minutes may show apparently toxic levels, yet
          the patient has minimal therapeutic response
        - CRITICAL: For digoxin TDM, wait at least 6 hours post-dose
          before sampling to ensure meaningful correlation with effect

    When to Use This Function vs One-Compartment Model:
        Use two-compartment model when:
        - Semi-log concentration-time plot shows curved (not straight) line
        - Drug is polar/water-soluble (slow tissue distribution)
        - Early samples show dramatically different levels than predicted
          by terminal-phase extrapolation
        - The drug is digoxin, aminoglycosides, vancomycin, or similar

        One-compartment may suffice when:
        - Drug is highly lipophilic (rapid tissue equilibration)
        - Only terminal-phase sampling is planned
        - Distribution phase is very short relative to sampling times

    Common Two-Compartment Drugs:
        Drug            t1/2,alpha      t1/2,beta       TDM Timing
        ---------------------------------------------------------------
        Digoxin         ~35 min         36-48 hours     Wait 6+ hours
        Lidocaine       ~8 min          1.5-2 hours     After loading
        Gentamicin      ~15-30 min      2-3 hours       Peak: 30min post
                                                        Trough: pre-dose
        Vancomycin      ~30 min         4-6 hours       Trough preferred

    Notes:
        - By convention, α > β (alpha is always the faster rate constant)
        - A and B must have same concentration units
        - α and β must have same reciprocal time units
        - Requires all 5 parameters (A, B, alpha, beta, t) to calculate
        - For practical TDM, avoid sampling during the distribution phase
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    A = Q_(kwargs["A"])
    B = Q_(kwargs["B"])
    alpha = Q_(kwargs["alpha"])
    beta = Q_(kwargs["beta"])
    t = Q_(kwargs["t"])

    # Ensure dimensionless exponents by converting rate*time
    alpha_t = (alpha * t).to("dimensionless").magnitude
    beta_t = (beta * t).to("dimensionless").magnitude

    # C(t) = A × e^(-α×t) + B × e^(-β×t)
    concentration = A * math.exp(-alpha_t) + B * math.exp(-beta_t)

    return format_output(
        concentration, "Concentration at time t", output_unit, decimals
    )


def alpha_beta_from_micro(**kwargs):
    """
    Calculate macroscopic rate constants (α, β) from microscopic rate constants.

    The two-compartment model has three micro-constants (k12, k21, K) that
    describe drug movement between compartments. These combine to give two
    macro-constants (α, β) that describe the observable biphasic decline.

    This function bridges the gap between physiological understanding (micro-
    constants that describe actual drug movement) and observable data (hybrid
    constants that can be measured from concentration-time curves).

    Formulas:
        α + β = k12 + k21 + K       (sum of hybrid constants)
        α × β = k21 × K             (product of hybrid constants)

    Solving the quadratic equation:
        α, β = [(k12 + k21 + K) ± sqrt((k12 + k21 + K)² - 4×k21×K)] / 2

    Where:
        - k12: Rate constant from central (1) to peripheral (2) compartment
        - k21: Rate constant from peripheral (2) to central (1) compartment
        - K (k10): Elimination rate constant from central compartment
        - α: Larger (distribution) rate constant - HYBRID constant
        - β: Smaller (elimination) rate constant - HYBRID constant

    Args:
        k12 (str): Transfer rate constant, compartment 1 to 2 (e.g., '1.5 1/hour')
        k21 (str): Transfer rate constant, compartment 2 to 1 (e.g., '0.8 1/hour')
        K (str): Elimination rate constant from central compartment (e.g., '0.5 1/hour')
        output_unit (str, optional): Desired output unit for rate constants
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        tuple: Two 5-tuples containing (α, β) results
        (
            ('Alpha (distribution rate constant)', ..., '1 / hour', ...),
            ('Beta (elimination rate constant)', ..., '1 / hour', ...)
        )

    Examples:
        Calculate α and β from micro-constants:
            >>> alpha_result, beta_result = alpha_beta_from_micro(
            ...     k12='1.5 1/hour', k21='0.8 1/hour', K='0.5 1/hour'
            ... )
            >>> print(alpha_result[3])  # Alpha formatted
            >>> print(beta_result[3])   # Beta formatted

        Digoxin-like parameters (slow tissue distribution):
            >>> # k12 < k21 suggests drug favors central compartment
            >>> # but eventually equilibrates with tissues
            >>> alpha_result, beta_result = alpha_beta_from_micro(
            ...     k12='0.5 1/hour', k21='0.3 1/hour', K='0.02 1/hour'
            ... )

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.2: Relationship between micro and macro constants
        - Figure 3.4: Basis of the two-compartmental model

    Understanding Micro-constants vs Hybrid Constants:
        Micro-constants (model-based, not directly measurable):
            - k12: Describes rate of drug LEAVING central compartment for tissues
                   Higher k12 = faster distribution to peripheral compartment
            - k21: Describes rate of drug RETURNING to central compartment
                   Higher k21 = faster return from tissues
            - K (k10): Describes rate of drug ELIMINATION from central compartment
                   Determines overall drug removal from body

        Hybrid/Macro-constants (observable from concentration-time data):
            - α (alpha): Reflects the DISTRIBUTION phase
                   Dominated by k12 (distribution to tissues)
                   Always the LARGER of the two hybrid constants
            - β (beta): Reflects the ELIMINATION phase
                   Reflects true elimination rate modified by redistribution
                   Always the SMALLER of the two hybrid constants

        Key Insight: While α is called the "distribution" rate constant and
        β the "elimination" rate constant, both hybrid constants are actually
        complex combinations of ALL three micro-constants. They are not
        simple surrogates for k12 and K.

    Mathematical Derivation:
        The characteristic equation for the two-compartment model is:
        λ² - (k12 + k21 + K)λ + k21×K = 0

        Using the quadratic formula:
        λ = [sum ± sqrt(sum² - 4×product)] / 2

        The two roots are α (larger) and β (smaller).

        Note that k12 appears ONLY in the sum (α + β), not in the product (α × β).
        This means k12 affects the relative magnitudes of α and β but not their
        product. Conversely, k21 and K appear in BOTH sum and product.

    Clinical Interpretation of Micro-constants:
        k12 > k21:
            - Drug favors peripheral (tissue) compartment
            - Extensive distribution, high Vss relative to V1
            - Example: Lipophilic drugs binding to tissue proteins

        k12 < k21:
            - Drug favors central compartment
            - Limited distribution, Vss closer to V1
            - Example: Polar drugs like digoxin (slow entry into tissues)

        Large K relative to k12, k21:
            - Rapid elimination dominates
            - Distribution phase may be obscured
            - May approximate one-compartment behavior

        Small K relative to k12, k21:
            - Distribution much faster than elimination
            - Clear biphasic curve with distinct phases
            - Classic two-compartment appearance on semi-log plot

    Digoxin Example - Why Distribution is Slow:
        Digoxin is a polar molecule with limited lipid solubility. Despite
        cardiac muscle being well-perfused, digoxin enters slowly because:
        - Poor passive diffusion across cell membranes
        - Depends on carrier-mediated transport
        - k12 is relatively small compared to highly lipophilic drugs

        This slow distribution (small k12) means the distribution phase
        extends for several hours, creating the TDM timing requirement
        of waiting 6+ hours post-dose before sampling.

    When One-Compartment Approximation is Valid:
        If distribution is very rapid (large k12 and k21) relative to
        elimination (small K), the drug quickly equilibrates between
        compartments and behaves approximately as one compartment.

        Mathematically: When α >> β, the A×e^(-αt) term decays to zero
        quickly, leaving only B×e^(-βt), which is a one-compartment equation.

    Notes:
        - α is always the larger value (distribution phase)
        - β is always the smaller value (elimination phase)
        - All rate constants must have the same time units
        - Requires all three micro-constants (k12, k21, K)
        - A negative discriminant indicates physically impossible parameters
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    k12 = Q_(kwargs["k12"])
    k21 = Q_(kwargs["k21"])
    K = Q_(kwargs["K"])

    # Sum and product of roots
    sum_roots = k12 + k21 + K  # α + β
    product_roots = k21 * K  # α × β

    # Quadratic formula: x = [-b ± sqrt(b² - 4ac)] / 2a
    # Here: λ² - (sum)λ + (product) = 0
    # So: λ = [sum ± sqrt(sum² - 4×product)] / 2

    sum_magnitude = sum_roots.magnitude
    product_magnitude = product_roots.magnitude
    unit = sum_roots.units

    discriminant = sum_magnitude**2 - 4 * product_magnitude

    if discriminant < 0:
        raise ValueError(
            "Invalid micro-constants: discriminant is negative. "
            "Check that k12, k21, and K values are physically reasonable."
        )

    sqrt_discriminant = math.sqrt(discriminant)

    alpha_magnitude = (sum_magnitude + sqrt_discriminant) / 2
    beta_magnitude = (sum_magnitude - sqrt_discriminant) / 2

    alpha = Q_(alpha_magnitude, unit)
    beta = Q_(beta_magnitude, unit)

    alpha_result = format_output(
        alpha, "Alpha (distribution rate constant)", output_unit, decimals
    )
    beta_result = format_output(
        beta, "Beta (elimination rate constant)", output_unit, decimals
    )

    return (alpha_result, beta_result)


def volume_central(**kwargs):
    """
    Calculate volume of the central compartment (V1).

    The central compartment represents the blood/plasma and highly perfused
    organs where drug distributes rapidly. V1 is calculated from the dose
    and the initial concentration (C0), which equals A + B at time zero.

    Understanding the central compartment is crucial for two-compartment
    drugs like digoxin, where the site of action (cardiac muscle) is in
    the SECOND compartment, not the central compartment we sample from.

    Formula: V1 = Dose / (A + B) = Dose / C0

    Where:
        - V1: Volume of central compartment
        - Dose: Amount of drug administered (IV bolus)
        - A: Distribution phase intercept coefficient
        - B: Elimination phase intercept coefficient
        - C0: Initial concentration = A + B

    Args (provide dose AND either C0 OR both A and B):
        dose (str): Administered dose (e.g., '500 mg')
        A (str, optional): Distribution phase intercept (e.g., '50 mg/L')
        B (str, optional): Elimination phase intercept (e.g., '30 mg/L')
        C0 (str, optional): Initial concentration at t=0 (e.g., '80 mg/L')
        output_unit (str, optional): Desired output unit (e.g., 'L')
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Volume of Central Compartment (V1)', 6.25, 'liter', '6.25 L', ...)

    Examples:
        Calculate V1 from dose and intercepts:
            >>> volume_central(dose='500 mg', A='50 mg/L', B='30 mg/L')
            ('Volume of Central Compartment (V1)', 6.25, 'liter', '6.25 L', ...)

        Calculate V1 from dose and initial concentration:
            >>> volume_central(dose='500 mg', C0='80 mg/L')
            ('Volume of Central Compartment (V1)', 6.25, 'liter', '6.25 L', ...)

        Digoxin example (relatively small V1 due to polar nature):
            >>> # Digoxin V1 is typically 20-30 L (larger than plasma volume
            >>> # but smaller than Vss due to extensive tissue binding)
            >>> volume_central(dose='500 microgram', C0='20 microgram/L')
            ('Volume of Central Compartment (V1)', 25.0, 'liter', '25.0 L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.3: Volume parameters
        - Section 3.6: Volume of distribution principles

    What the Central Compartment Represents:
        The first (central) compartment includes:
        - Blood/plasma (the sampling site)
        - Highly perfused organs that equilibrate rapidly with blood
        - Tissues where drug enters quickly via passive diffusion

        For polar drugs like digoxin:
        - The central compartment is relatively small
        - Cardiac muscle (site of action) is NOT in the central compartment
        - Cardiac muscle equilibrates slowly and is part of compartment 2
        - This disconnect explains why early blood levels do not predict effect

    Clinical Interpretation:
        - V1 is typically smaller than total Vd or Vss
        - V1 reflects initial distribution space (blood + well-perfused tissues)
        - For most drugs: V1 approximates plasma volume (3-5 L) or somewhat larger
        - V1 > 10 L suggests rapid extravascular distribution even before
          the slower second compartment is reached

    Volume Hierarchy in Two-Compartment Models:
        V1 < Vss < Vd(area) or Vd(beta)

        - V1: Central compartment only (smallest)
        - Vss: Both compartments at equilibrium (intermediate)
        - Vd(area): Back-extrapolated from terminal phase (largest)

    Why V1 Matters for Digoxin:
        For digoxin TDM, understanding V1 helps explain the timing problem:

        1. Immediately post-dose: All drug is in V1 (central compartment)
           - Blood concentration is VERY HIGH (Dose/V1)
           - Tissue concentration is essentially ZERO
           - Therapeutic effect is MINIMAL

        2. During distribution (0-6 hours): Drug redistributes
           - Blood concentration FALLS rapidly
           - Tissue concentration RISES
           - Therapeutic effect INCREASES
           - Blood level does NOT correlate with effect!

        3. After distribution (6+ hours): Pseudo-equilibrium
           - Blood and tissue concentrations fall in parallel
           - Ratio between compartments is constant
           - Blood level now CORRELATES with therapeutic effect
           - This is when TDM sampling is meaningful

    Loading Dose Considerations:
        For two-compartment drugs, loading dose calculations are complex:

        - If based on V1: Achieves target C0 but may be toxic if distribution
          to site of action is slow (digoxin example)
        - If based on Vss: May under-dose initially but safer for slow
          distribution drugs
        - For digoxin: Loading doses are typically given in divided doses
          over 24 hours to allow distribution and avoid toxicity

    Digoxin-Specific Clinical Guidance:
        - V1 for digoxin: approximately 20-30 L (0.3-0.4 L/kg)
        - Vss for digoxin: approximately 500-700 L (7-10 L/kg)
        - The large difference (V1 << Vss) reflects extensive tissue binding
        - A water-soluble drug like digoxin distributes into tissues slowly
          but eventually binds extensively (large Vss)

    Other Two-Compartment Drug Examples:
        Drug            V1 (typical)    Vss (typical)
        ----------------------------------------------
        Digoxin         ~25 L           ~500 L
        Lidocaine       ~50 L           ~100 L
        Gentamicin      ~15-20 L        ~20-25 L
        Vancomycin      ~15-20 L        ~30-40 L

    Notes:
        - Must provide dose AND (C0 OR both A and B)
        - V1 is NOT the same as Vss or Vd(area)
        - V1 is used for loading dose calculations, but with caution for
          drugs with slow distribution (like digoxin)
        - For obese patients, V1 should be estimated from ideal body weight
          for water-soluble drugs like digoxin
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    dose = Q_(kwargs["dose"])

    # Get C0 either directly or from A + B
    if "C0" in kwargs:
        C0 = Q_(kwargs["C0"])
    elif "A" in kwargs and "B" in kwargs:
        A = Q_(kwargs["A"])
        B = Q_(kwargs["B"])
        C0 = A + B
    else:
        raise ValueError(
            "Must provide either 'C0' or both 'A' and 'B' parameters "
            "along with 'dose' to calculate V1."
        )

    V1 = dose / C0

    return format_output(
        V1, "Volume of Central Compartment (V1)", output_unit, decimals
    )


def volume_steady_state(**kwargs):
    """
    Calculate volume of distribution at steady state (Vss).

    Vss accounts for drug distributed in both the central and peripheral
    compartments at equilibrium. It is larger than V1 because it includes
    drug that has moved to the peripheral compartment.

    For two-compartment drugs like digoxin, Vss can be dramatically larger
    than V1 due to extensive tissue binding, even though tissue distribution
    is slow. Understanding this distinction is crucial for dosing.

    Formula: Vss = V1 × (1 + k12/k21)

    Where:
        - Vss: Volume of distribution at steady state
        - V1: Volume of central compartment
        - k12: Rate constant from central to peripheral compartment
        - k21: Rate constant from peripheral to central compartment
        - k12/k21: Ratio reflecting extent of peripheral distribution

    Alternative formula: Vss = V1 + V2
    Where V2 = V1 × (k12/k21) is the peripheral compartment volume.

    Args:
        V1 (str): Volume of central compartment (e.g., '6 L')
        k12 (str): Transfer rate constant, compartment 1 to 2 (e.g., '1.5 1/hour')
        k21 (str): Transfer rate constant, compartment 2 to 1 (e.g., '0.8 1/hour')
        output_unit (str, optional): Desired output unit (e.g., 'L')
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Volume of Distribution at Steady State (Vss)', 17.25, 'liter', ...)

    Examples:
        Calculate Vss from V1 and transfer constants:
            >>> volume_steady_state(V1='6 L', k12='1.5 1/hour', k21='0.8 1/hour')
            ('Volume of Distribution at Steady State (Vss)', 17.25, 'liter', ...)

        Digoxin example (extensive tissue binding, large Vss/V1 ratio):
            >>> # Digoxin: V1 ~25 L, but Vss ~500 L due to tissue binding
            >>> # This reflects k12/k21 ratio creating large V2
            >>> volume_steady_state(V1='25 L', k12='0.8 1/hour', k21='0.04 1/hour')
            ('Volume of Distribution at Steady State (Vss)', 525.0, 'liter', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.3: Volume parameters and steady state
        - Section 3.6: Volume of distribution and tissue binding

    Understanding Vss in Two-Compartment Models:
        Vss represents the TOTAL volume of distribution once equilibrium
        between compartments is achieved. It is the sum of:
        - V1: Central compartment (blood + rapidly equilibrating tissues)
        - V2: Peripheral compartment (slowly equilibrating tissues)

        The ratio k12/k21 determines how much larger Vss is than V1:
        - k12/k21 = 1: Vss = 2 × V1 (equal distribution)
        - k12/k21 > 1: Drug favors tissues (Vss >> V1)
        - k12/k21 < 1: Drug favors blood (Vss closer to V1)

    Clinical Interpretation:
        - Vss > V1 always (drug distributes to peripheral tissues)
        - k12 > k21: Drug favors peripheral compartment (higher Vss)
        - k12 < k21: Drug favors central compartment (lower Vss)
        - Vss is the appropriate volume for steady-state calculations
        - Used for maintenance dose calculations during continuous infusion

    Relationship to Other Volumes:
        - V1: Central compartment only (smallest)
        - Vss: Both compartments at equilibrium (intermediate)
        - Vd(area) or Vd(beta): Back-extrapolated from terminal phase (largest)
        - Relationship: V1 < Vss < Vd(area)

    Digoxin - The Classic Example of Large Vss:
        Digoxin demonstrates extreme Vss/V1 disparity:
        - V1: ~25 L (0.3-0.4 L/kg) - central compartment
        - Vss: ~500-700 L (7-10 L/kg) - total at equilibrium

        Why such a large Vss for a polar drug?
        1. Digoxin binds extensively to skeletal muscle Na+/K+-ATPase
        2. This tissue binding creates a large "reservoir" (V2)
        3. Distribution INTO tissues is slow (k12 relatively small)
        4. But once there, drug remains bound (k21 even smaller)
        5. The ratio k12/k21 is large, making Vss >> V1

        Clinical implication: A relatively small dose fills the large Vss,
        but it takes multiple half-lives (days for digoxin) to reach
        steady state.

    Vss and Therapeutic Drug Monitoring:
        Understanding Vss helps interpret TDM results for digoxin:

        During Loading:
        - Initially, drug is concentrated in V1 (high blood levels)
        - Over 6+ hours, drug redistributes to fill Vss
        - Blood levels fall as drug moves to tissues
        - Sample too early = misleadingly high levels

        At Steady State:
        - Drug is distributed throughout Vss
        - Blood concentration reflects total body content
        - Trough levels (pre-dose) give reliable TDM
        - Amount in body = Css × Vss

    Using Vss for Dosing Calculations:
        Maintenance Dose = Cl × Css × tau
            - Cl: Clearance (relates to elimination)
            - Css: Target steady-state concentration
            - tau: Dosing interval

        Time to Steady State ≈ 4-5 × t1/2,beta
            - For digoxin: 4-5 × 36 hours = 6-7.5 days
            - Explains why loading doses are used when urgent

        Loading Dose Considerations:
            - For rapid effect: base on V1 or smaller volume
            - For safety with slow distribution: give over time
            - For digoxin: divided loading over 24 hours

    Common Two-Compartment Drug Volumes:
        Drug            V1 (L)      Vss (L)     Vss/V1 Ratio
        --------------------------------------------------------
        Digoxin         ~25         ~500-700    ~20-25
        Lidocaine       ~50         ~100        ~2
        Gentamicin      ~15-20      ~20-25      ~1.3
        Vancomycin      ~15-20      ~30-40      ~2
        Diazepam        ~30         ~80-100     ~3

    Why Vss/V1 Ratio Matters:
        Small ratio (< 2): Minimal redistribution phase
            - Distribution completes quickly
            - One-compartment model may suffice
            - Example: aminoglycosides

        Large ratio (> 5): Extensive tissue distribution
            - Prolonged redistribution phase
            - Clear biphasic curve on semi-log plot
            - Sampling timing is critical
            - Example: digoxin

    Effect of Body Composition on Vss:
        For water-soluble drugs like digoxin:
        - Vss correlates with lean body mass, not fat
        - Obese patients: use ideal body weight
        - Elderly: may have reduced Vss due to muscle loss
        - Renal disease: may alter tissue binding and Vss

    Notes:
        - Requires all three parameters (V1, k12, k21)
        - k12 and k21 must have same units (they divide out)
        - Vss is model-independent (same regardless of compartment number)
        - Vss is the most physiologically meaningful volume parameter
          for steady-state calculations
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    V1 = Q_(kwargs["V1"])
    k12 = Q_(kwargs["k12"])
    k21 = Q_(kwargs["k21"])

    # k12/k21 should be dimensionless
    ratio = (k12 / k21).to("dimensionless").magnitude

    # Vss = V1 × (1 + k12/k21)
    Vss = V1 * (1 + ratio)

    return format_output(
        Vss, "Volume of Distribution at Steady State (Vss)", output_unit, decimals
    )


def distribution_half_life(**kwargs):
    """
    Calculate half-life of the distribution phase.

    The distribution half-life (t1/2,alpha) represents the time for the distribution
    process to reach 50% completion. This is the rapid initial phase where
    drug moves from blood into tissues. Understanding this half-life is
    essential for determining when to sample blood for therapeutic drug
    monitoring (TDM).

    Formula: t1/2,alpha = ln(2) / alpha = 0.693 / alpha

    Where:
        - t1/2,alpha: Distribution phase half-life
        - alpha: Distribution rate constant
        - ln(2) = 0.693

    Args:
        alpha (str): Distribution rate constant (e.g., '2.0 1/hour')
        output_unit (str, optional): Desired output unit (e.g., 'minute')
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Distribution Half-Life (t1/2,alpha)', 0.35, 'hour', '0.35 hour', ...)

    Examples:
        Calculate distribution half-life:
            >>> distribution_half_life(alpha='2.0 1/hour')
            ('Distribution Half-Life (t1/2,alpha)', 0.35, 'hour', '0.35 hour', ...)

        Digoxin example (slow distribution):
            >>> # Digoxin alpha ~1.2/hour, giving t1/2,alpha ~35 minutes
            >>> distribution_half_life(alpha='1.2 1/hour')
            ('Distribution Half-Life (t1/2,alpha)', 0.58, 'hour', ...)
            # About 35 minutes - explains the 6-hour wait for TDM

        Lidocaine example (rapid distribution):
            >>> # Lidocaine alpha ~5.0/hour, giving t1/2,alpha ~8 minutes
            >>> distribution_half_life(alpha='5.0 1/hour')
            ('Distribution Half-Life (t1/2,alpha)', 0.14, 'hour', ...)
            # About 8 minutes - rapid CNS distribution

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.4: Half-lives in two-compartment models
        - Section 7.5: TDM timing for digoxin

    Why Distribution Half-Life Matters:
        The distribution phase creates a window during which blood
        concentrations do NOT correlate with tissue concentrations
        or therapeutic effect. This is because:

        1. Drug is leaving blood rapidly (falling blood levels)
        2. Drug is entering tissues rapidly (rising tissue levels)
        3. The site of action may be in the tissue compartment
        4. Blood sampling during this phase is misleading

    The 4-5 Half-Life Rule:
        Distribution is essentially complete after 4-5 distribution half-lives:
        - 1 half-life: 50% complete
        - 2 half-lives: 75% complete
        - 3 half-lives: 87.5% complete
        - 4 half-lives: 93.75% complete
        - 5 half-lives: 96.875% complete (~97%)

        After 4-5 × t1/2,alpha, the distribution term (A × e^(-alpha×t))
        has decayed to negligible levels, and concentration follows
        the elimination phase (B × e^(-beta×t)).

    Critical TDM Timing Guidance:
        RULE: Wait at least 4-5 × t1/2,alpha before sampling for TDM

        Drug            t1/2,alpha      Wait at least
        ------------------------------------------------
        Digoxin         ~35 min         ~3-6 hours*
        Lidocaine       ~8 min          ~30-40 minutes
        Gentamicin      ~15-30 min      ~1-2 hours
        Vancomycin      ~30 min         ~2 hours

        *For digoxin, the standard recommendation is 6 hours post-dose
         to ensure complete equilibration and account for variability.

    Digoxin TDM - The Classic Example:
        For digoxin, the distribution half-life of ~35 minutes means:
        - 4 × 35 min = 2.3 hours (mathematically minimum)
        - Clinical recommendation: 6 hours minimum

        Why 6 hours instead of ~2.5 hours?
        - Provides margin for patient variability
        - Ensures tissue-blood ratio is truly constant
        - Accounts for the critical nature of digoxin toxicity
        - Standard practice in clinical laboratories

        The problem with early sampling:
        "A blood sample taken at 20 minutes post-dose may show very high
        (apparently toxic) levels, yet the therapeutic response is minimal
        because tissue concentrations are still low." - Rowe

    Clinical Interpretation:
        - Short t1/2,alpha (minutes): Rapid tissue distribution
          Example: Lidocaine rapidly enters CNS - toxicity can occur quickly
        - Long t1/2,alpha (hours): Slow equilibration with tissues
          Example: Digoxin slowly enters cardiac muscle - delayed effect
        - Distribution usually complete within 4-5 × t1/2,alpha
        - Blood sampling during distribution phase gives unreliable levels
        - Wait until post-distribution for therapeutic drug monitoring

    Clinical Examples with TDM Implications:
        Digoxin: t1/2,alpha ~35 minutes
            - Cardiac muscle is in peripheral compartment (slow entry)
            - Blood levels HIGH, tissue levels LOW early after dose
            - TDM: Wait 6+ hours post-dose
            - Early samples may show "toxic" levels with minimal effect

        Lidocaine: t1/2,alpha ~8 minutes
            - CNS is rapidly perfused (fast distribution)
            - CNS toxicity can occur within minutes of IV bolus
            - TDM: Usually after loading dose is complete
            - Distribution phase is brief but clinically important

        Aminoglycosides: t1/2,alpha ~15-30 minutes
            - Distribution into extracellular fluid
            - Peak samples: 30 minutes post-infusion
            - Trough samples: immediately pre-dose
            - Distribution phase must be considered for peak timing

        Vancomycin: t1/2,alpha ~30 minutes
            - Distribution into tissues
            - Trough samples preferred (avoid distribution phase)
            - If peak needed: at least 1-2 hours post-infusion

    When One-Compartment Model Suffices:
        If t1/2,alpha is very short relative to the dosing interval or
        sampling time, the drug effectively behaves as one-compartment:
        - The distribution phase is essentially instantaneous
        - Semi-log plot appears linear
        - Blood levels always correlate with tissue levels
        - Example: Highly lipophilic drugs with rapid equilibration

    Notes:
        - alpha must have units of 1/time (e.g., 1/hour, 1/minute)
        - Result will have time units inverse to alpha's units
        - t1/2,alpha is always shorter than t1/2,beta
        - This is the half-life that determines TDM sampling timing
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    alpha = Q_(kwargs["alpha"])

    # t½,α = 0.693 / α
    half_life = math.log(2) / alpha

    return format_output(
        half_life, "Distribution Half-Life (t1/2,alpha)", output_unit, decimals
    )


def terminal_half_life(**kwargs):
    """
    Calculate half-life of the terminal (elimination) phase.

    The terminal half-life (t1/2,beta) represents the time for plasma concentration
    to decrease by 50% during the elimination phase. This is the slower
    terminal phase after distribution is complete. This half-life determines
    dosing intervals, time to steady state, and drug washout periods.

    Formula: t1/2,beta = ln(2) / beta = 0.693 / beta

    Where:
        - t1/2,beta: Terminal/elimination phase half-life
        - beta: Elimination rate constant
        - ln(2) = 0.693

    Args:
        beta (str): Elimination rate constant (e.g., '0.3 1/hour')
        output_unit (str, optional): Desired output unit (e.g., 'hour')
        decimals (int, optional): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Terminal Half-Life (t1/2,beta)', 2.31, 'hour', '2.31 hour', ...)

    Examples:
        Calculate terminal half-life:
            >>> terminal_half_life(beta='0.3 1/hour')
            ('Terminal Half-Life (t1/2,beta)', 2.31, 'hour', '2.31 hour', ...)

        Digoxin example (long terminal half-life):
            >>> # Digoxin beta ~0.02/hour, giving t1/2,beta ~36 hours
            >>> terminal_half_life(beta='0.02 1/hour')
            ('Terminal Half-Life (t1/2,beta)', 34.66, 'hour', ...)
            # About 36 hours - explains once-daily dosing and long washout

        Lidocaine example (short terminal half-life):
            >>> # Lidocaine beta ~0.4/hour, giving t1/2,beta ~1.7 hours
            >>> terminal_half_life(beta='0.4 1/hour')
            ('Terminal Half-Life (t1/2,beta)', 1.73, 'hour', ...)
            # Requires continuous infusion for maintenance

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 7: Two-Compartment Model
        - Section 7.4: Half-lives in two-compartment models
        - Chapter 14: Creatinine clearance and digoxin dosing

    Clinical Significance of t1/2,beta:
        The terminal half-life is the PRIMARY determinant of:
        1. Dosing interval selection
        2. Time to reach steady state
        3. Time to eliminate drug after discontinuation
        4. Accumulation ratio at steady state

        This is the "half-life" typically reported in drug references and
        literature. When a drug's half-life is stated without qualification,
        it almost always refers to t1/2,beta.

    The 4-5 Half-Life Rules:
        Time to Steady State:
        - 1 half-life: 50% of steady state
        - 2 half-lives: 75% of steady state
        - 3 half-lives: 87.5% of steady state
        - 4 half-lives: 93.75% of steady state
        - 5 half-lives: 96.875% of steady state (~97%)

        Drug Washout (after discontinuation):
        - 1 half-life: 50% remaining
        - 2 half-lives: 25% remaining
        - 3 half-lives: 12.5% remaining
        - 4 half-lives: 6.25% remaining
        - 5 half-lives: 3.125% remaining (~eliminated)

    Digoxin - Clinical Decision Making:
        With t1/2,beta ~36-48 hours, digoxin has important implications:

        Dosing Interval:
        - Once daily dosing is appropriate (tau < t1/2,beta)
        - Minimal peak-trough fluctuation
        - Missed doses have prolonged effect

        Time to Steady State:
        - 4-5 × 36 hours = 6-7.5 days
        - Explains why loading doses are used for urgent digitalization
        - Full effect of dose changes takes about a week

        Drug Interactions and Dose Changes:
        - Effects of dose change take 1 week to fully manifest
        - When adding/removing interacting drugs, monitor over days
        - Renal impairment prolongs t1/2,beta further

        Toxicity Management:
        - Drug persists for days after discontinuation
        - Digoxin-specific antibody (Fab fragments) may be needed
        - Supportive care must be prolonged

    Clinical Examples with Dosing Implications:
        Digoxin: t1/2,beta ~36-48 hours
            - Dosing: Once daily
            - Steady state: 6-7 days
            - Washout: 6-7 days
            - Loading: Often used (divided over 24 hours)
            - TDM: Trough levels, at least 6 hours post-dose

        Lidocaine: t1/2,beta ~1.5-2 hours
            - Dosing: Continuous IV infusion after loading
            - Steady state: 6-10 hours
            - Washout: 6-10 hours
            - Loading: IV bolus required for immediate effect
            - TDM: During infusion, any time (no distribution phase issue)

        Gentamicin: t1/2,beta ~2-3 hours
            - Dosing: Every 8 hours (traditional) or once daily (extended)
            - Steady state: 8-15 hours
            - Washout: 8-15 hours
            - Peak/trough monitoring important
            - Renal adjustment critical (half-life prolonged in renal failure)

        Vancomycin: t1/2,beta ~4-6 hours
            - Dosing: Every 8-12 hours typically
            - Steady state: 16-30 hours
            - Trough monitoring preferred
            - Renal adjustment critical

    Comparison with One-Compartment Model:
        - In one-compartment: single half-life = 0.693/K
        - In two-compartment: t1/2,beta is analogous but reflects terminal phase
        - t1/2,beta > half-life that would be calculated from K alone

        The difference arises because in two-compartment models, drug
        returning from tissues slows the apparent elimination from blood.
        Beta reflects this combined process, not just elimination.

    When One-Compartment Approximation Works:
        For drugs with rapid distribution (short t1/2,alpha), the
        terminal phase dominates, and a one-compartment model using
        t1/2,beta may be adequate for most clinical purposes:
        - Dosing calculations
        - Steady-state predictions
        - Washout estimates

        The two-compartment model is essential when:
        - Detailed concentration predictions are needed
        - Loading dose optimization is critical
        - The site of action is in the peripheral compartment
        - TDM interpretation during distribution phase is needed

    Factors Affecting t1/2,beta:
        Increased (longer) t1/2,beta:
        - Renal impairment (for renally cleared drugs like digoxin)
        - Hepatic impairment (for hepatically cleared drugs)
        - Heart failure (reduced clearance of many drugs)
        - Drug interactions (CYP450 inhibitors, etc.)
        - Advanced age (reduced clearance)

        Decreased (shorter) t1/2,beta:
        - Enzyme induction (rifampin, phenytoin, etc.)
        - Increased renal/hepatic blood flow
        - Pediatric patients (some drugs)

    Common Two-Compartment Drug Half-Lives:
        Drug            t1/2,alpha      t1/2,beta       Dosing
        ---------------------------------------------------------------
        Digoxin         ~35 min         36-48 hours     Once daily
        Lidocaine       ~8 min          1.5-2 hours     Continuous infusion
        Gentamicin      ~15-30 min      2-3 hours       Q8h or once daily
        Vancomycin      ~30 min         4-6 hours       Q8-12h
        Diazepam        ~15-30 min      20-100 hours    Variable
        Fentanyl        ~1-2 min        2-4 hours       Infusion or PRN

    Notes:
        - beta must have units of 1/time (e.g., 1/hour, 1/minute)
        - Result will have time units inverse to beta's units
        - t1/2,beta is always longer than t1/2,alpha
        - This is the clinically relevant half-life for dosing interval selection
        - Always consider patient-specific factors that may alter t1/2,beta
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    beta = Q_(kwargs["beta"])

    # t½,β = 0.693 / β
    half_life = math.log(2) / beta

    return format_output(
        half_life, "Terminal Half-Life (t1/2,beta)", output_unit, decimals
    )
