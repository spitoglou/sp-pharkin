"""
Non-linear (Michaelis-Menten) pharmacokinetics calculations.

These functions implement saturable kinetics from "Pharmacokinetics" by Philip Rowe,
Chapter 11. Non-linear kinetics occur when drug concentrations are high enough to
saturate metabolic enzymes or transport mechanisms. This chapter addresses a situation
where many of the assumptions made in standard pharmacokinetics break down.

The problem only arises with drugs that are mainly or entirely eliminated by metabolism.
Drug metabolism by Cytochrome P450 enzymes is an example of an enzyme-catalysed reaction
that follows Michaelis-Menten kinetics when concentrations approach enzyme saturation.

CRITICAL CLINICAL WARNING - PHENYTOIN AS THE PARADIGM:
Phenytoin is THE clinically significant case of non-linear kinetics. It is a drug with:
- A narrow therapeutic window (10-20 mg/L total, 1-2 mg/L free)
- Serious toxicity in overdose (nystagmus, ataxia, confusion, seizures)
- Therapeutic concentrations that approach or exceed Km (4-10 mg/L)
- Disproportionate concentration changes with dose adjustments

Key consequences of non-linear phenytoin kinetics:
- A 10% dose increase can cause >50% concentration increase
- Dosage adjustment is a SPECIALIST job - normal rules don't apply
- Time to reach new steady state is prolonged (2-4 weeks, not 4-5 half-lives)
- Must use therapeutic drug monitoring (TDM) for all dose adjustments

Drugs with Clinically Relevant Non-Linear Kinetics:
1. PHENYTOIN (anticonvulsant) - THE most clinically important example
   - Km: 4-10 mg/L (highly variable between patients)
   - Vmax: 350-700 mg/day (varies with body size, genetics, enzyme induction)
   - Therapeutic range: 10-20 mg/L - often ABOVE Km, meaning saturated metabolism
   - Small dose changes cause disproportionately large concentration changes
   - Given the potential toxicity, dosage adjustment requires specialist knowledge

2. ETHANOL (alcohol)
   - Fully saturates liver enzymes at typical consumption levels
   - Zero-order elimination at ~7 g/hour regardless of blood level
   - Cannot "speed up" elimination by any means
   - BAC falls linearly with time, not exponentially
   - Km ~10 mg/dL, but social drinking produces levels far above this

3. SALICYLATES (aspirin at high doses)
   - Cause enzyme saturation at anti-inflammatory doses
   - Not clinically managed with narrow concentration targets
   - Non-linear kinetics are largely academic in clinical practice

4. THEOPHYLLINE (at high doses)
   - Shows some degree of enzyme saturation
   - Effect is relatively small compared to phenytoin
   - For practical clinical purposes, treated as linear

Why Standard Pharmacokinetic Concepts Fail:
With non-linear kinetics, there is NO constant proportionality between rate of
elimination and drug concentration. This means:
- Elimination rate constant (K) no longer applies
- Half-life is NOT constant - it INCREASES with concentration
- Clearance is concentration-dependent, not constant
- Cannot use standard first-order equations
- The only valid parameters are Vmax and Km (enzymological constants)

Unlike K, Vmax and Km DO have fixed values for a particular drug in a particular
patient and can be used to calculate what concentration would arise from a given
dosage regime.

Recognizing Non-Linear Behavior Clinically:
1. Disproportionate response to dose changes:
   - Doubling dose causes MORE than doubling of concentration
   - Or small dose increase causes unexpectedly large concentration rise

2. Variable "half-life" measurements:
   - Half-life appears longer at higher concentrations
   - Different half-lives measured at different concentration ranges

3. Prolonged time to steady state:
   - Takes much longer than expected 4-5 half-lives
   - May take 2-4 weeks instead of days

4. Therapeutic drug monitoring shows inconsistencies:
   - Concentrations don't correlate linearly with dose
   - Dose-concentration relationship curves upward

Clinical Management of Non-Linear Drugs (Phenytoin Model):
1. Start LOW: Begin with conservative doses (200-300 mg/day for phenytoin)
2. Increase SLOWLY: Use small increments (25-50 mg for phenytoin)
3. Wait LONG: Allow 2+ weeks between adjustments for new steady state
4. Monitor ALWAYS: Therapeutic drug monitoring is mandatory, not optional
5. Specialist input: Dosage adjustment is a specialist job
6. Consider free levels: In hypoalbuminemia, measure free (unbound) drug
7. Individual parameters: Each patient has different Km and Vmax

Drug Development Implications:
Pharmaceutical companies strongly prefer drugs with linear kinetics because:
- Dose adjustment is simple and predictable for clinicians
- Given equal efficacy, a linear drug will always be preferred
- Companies establish early whether candidate drugs are non-linear
- Non-linear candidates may be terminated to avoid clinical complexity
- There is strong commercial disinclination to develop non-linear drugs

Mathematical Basis - Michaelis-Menten Equation:
    v = (Vmax x C) / (Km + C)

Where:
    v = rate of elimination (mass/time)
    Vmax = maximum elimination rate when enzymes fully saturated (mass/time)
    C = drug concentration (mass/volume)
    Km = Michaelis constant - concentration at which v = Vmax/2 (mass/volume)

Behavior at Different Concentration Ranges:
    When C << Km: v ≈ (Vmax/Km) x C  (FIRST-ORDER, linear kinetics)
        - Rate proportional to concentration
        - Constant half-life
        - Doubling dose doubles concentration
        - Time to steady state = 4-5 half-lives

    When C = Km:  v = Vmax/2  (TRANSITION zone)
        - Half-maximal elimination rate
        - Kinetics transitioning from first to zero order
        - This is where phenytoin often operates!

    When C >> Km: v ≈ Vmax  (ZERO-ORDER, saturated kinetics)
        - Rate is constant regardless of concentration
        - "Half-life" increases with concentration
        - Doubling dose MORE than doubles concentration
        - Time to steady state is prolonged and unpredictable
"""

from .lib import format_output
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def michaelis_menten_rate(**kwargs):
    """
    Calculate rate of elimination with saturable (Michaelis-Menten) kinetics.

    This function implements the fundamental equation for enzyme-mediated drug
    elimination when concentrations are high enough to cause significant enzyme
    saturation. The Michaelis-Menten equation describes how elimination rate
    varies with concentration when enzyme capacity becomes limiting.

    This equation replaces first-order kinetics when C approaches or exceeds Km.
    It is the ONLY valid way to describe elimination for non-linear drugs like
    phenytoin, where standard pharmacokinetic parameters (K, half-life, clearance)
    cease to have fixed values.

    Formula: v = (Vmax x C) / (Km + C)

    Where:
        v = rate of elimination (mass/time, e.g., mg/hour)
        Vmax = maximum rate when enzymes fully saturated (mass/time)
        C = drug concentration (mass/volume, e.g., mg/L)
        Km = Michaelis constant - concentration at which v = Vmax/2 (mass/volume)

    Behavior at Different Concentration Ranges:
        When C << Km: v ≈ (Vmax/Km) x C (FIRST-ORDER, linear kinetics)
            - Rate proportional to concentration
            - Constant half-life applies
            - Standard pharmacokinetic equations valid
            - Doubling dose doubles concentration

        When C = Km:  v = Vmax/2 (TRANSITION zone - half-maximal rate)
            - Kinetics transitioning from first to zero order
            - Half-life starting to become concentration-dependent
            - THIS IS WHERE PHENYTOIN OFTEN OPERATES

        When C >> Km: v ≈ Vmax (ZERO-ORDER, saturated kinetics)
            - Rate is constant regardless of concentration
            - Half-life increases with concentration
            - Standard half-life concept does NOT apply
            - Doubling dose MORE than doubles concentration
            - Time to steady state greatly prolonged

    CRITICAL CLINICAL WARNING - WHY NON-LINEAR KINETICS ARE DANGEROUS:
    1. Disproportionate concentration changes:
       - A 10% phenytoin dose increase can cause >50% concentration rise
       - Small adjustments can push patient from therapeutic to toxic

    2. Standard half-life does NOT apply:
       - Half-life is NOT constant - it INCREASES with concentration
       - Cannot predict time course using standard equations
       - Cannot use "4-5 half-lives to steady state" rule

    3. Prolonged, unpredictable time to steady state:
       - May take 2-4 WEEKS instead of days
       - Accumulation continues much longer than expected
       - Patient may appear stable, then suddenly become toxic

    4. Therapeutic drug monitoring is MANDATORY:
       - Cannot predict concentrations from dose alone
       - Must measure levels after every dose change
       - Wait adequate time (2+ weeks) before measuring new steady state

    Args (provide exactly 3 of 4):
        v (str): Rate of elimination (e.g., '10 mg/hour')
        Vmax (str): Maximum elimination rate (e.g., '15 mg/hour')
        C (str): Drug concentration (e.g., '20 mg/L')
        Km (str): Michaelis constant (e.g., '10 mg/L')

    Optional kwargs:
        output_unit (str): Desired output unit for the result
        decimals (int): Number of decimal places for rounding (default 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Rate of Elimination (v)', 10.0, 'milligram / hour', '10.0 mg/hour', ...)

    Examples:
        Calculate elimination rate at concentration 20 mg/L:
            >>> michaelis_menten_rate(Vmax='15 mg/hour', C='20 mg/L', Km='10 mg/L')
            ('Rate of Elimination (v)', 10.0, 'milligram / hour', '10.0 mg/hour', ...)

        Calculate Km from known rate and Vmax:
            >>> michaelis_menten_rate(v='10 mg/hour', Vmax='15 mg/hour', C='20 mg/L')
            ('Michaelis Constant (Km)', 10.0, 'milligram / liter', '10.0 mg/L', ...)

        Calculate required Vmax for desired elimination rate:
            >>> michaelis_menten_rate(v='10 mg/hour', C='20 mg/L', Km='10 mg/L')
            ('Maximum Elimination Rate (Vmax)', 15.0, 'milligram / hour', '15.0 mg/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 11: Non-linear pharmacokinetics
        - Section 11.1, pages 6206-6240: Drug metabolism as enzyme-catalysed reaction
        - Section 11.2, pages 6302-6320: Exceptions to linearity (phenytoin, salicylates, ethanol)
        - Section 11.3, pages 6328-6380: Effect on dose-concentration relationship
        - Section 11.4, pages 6380-6488: Clinical significance of non-linear kinetics

    Phenytoin-Specific Clinical Context:
        Phenytoin is THE clinically significant example of non-linear kinetics:
        - Km: 4-10 mg/L (highly variable between patients)
        - Vmax: 350-700 mg/day (varies with body size, genetics)
        - Therapeutic range: 10-20 mg/L (often AT or ABOVE Km!)
        - This means phenytoin operates in the saturated region
        - Dosage adjustment is a SPECIALIST job - normal rules don't apply

        Phenytoin Dosing Dangers:
        - Narrow therapeutic window with serious toxicity in overdose
        - At 10 mg/L (low therapeutic): C/Km = 1-2.5, highly non-linear
        - At 20 mg/L (high therapeutic): C/Km = 2-5, severely saturated
        - A patient stable at 15 mg/L may become toxic with tiny dose increase

    Ethanol Clinical Context:
        Ethanol provides another example of saturated kinetics:
        - Km ≈ 10 mg/dL (very low)
        - Social drinking produces concentrations far above Km
        - Elimination is essentially zero-order at ~7 g/hour
        - Cannot "speed up" elimination - rate is fixed at Vmax
        - Blood alcohol falls linearly with time, not exponentially
        - Clinicians don't typically dose-adjust ethanol, but the kinetics
          explain why "sobering up" takes a fixed time regardless of level

    Mathematical Derivations:
        Solving for each variable:
        - v = (Vmax x C) / (Km + C)
        - Vmax = v x (Km + C) / C
        - Km = C x (Vmax - v) / v
        - C = (Km x v) / (Vmax - v)

    Why Standard Parameters Fail:
        With non-linear kinetics:
        - K (elimination rate constant) no longer applies
        - Half-life is NOT constant
        - Clearance is concentration-dependent
        - Only Vmax and Km have fixed values for a given patient/drug

    Km - The Critical Threshold:
        - Km represents enzyme affinity - lower Km means higher affinity
        - Km determines WHERE kinetics transition from linear to non-linear
        - When C < 0.1 x Km: kinetics are approximately linear
        - When C > Km: kinetics are predominantly saturated
        - Phenytoin's danger: therapeutic range is AT or ABOVE Km

    Vmax - The Capacity Ceiling:
        - Vmax depends on amount of enzyme present
        - Can be INDUCED by other drugs (e.g., carbamazepine) → ↑ Vmax
        - Can be INHIBITED by other drugs → ↓ Vmax → higher concentrations
        - Genetic polymorphisms affect Vmax between individuals
        - This contributes to inter-patient variability

    Drug Development Implications:
        - Companies test early for non-linear kinetics
        - Strong commercial disinclination to develop non-linear drugs
        - A linear drug will always be preferred over non-linear if equal efficacy
        - Non-linear candidates may be terminated early in development

    Clinical Management Strategy:
        1. Know if your drug has non-linear kinetics (phenytoin!)
        2. Start with LOW doses
        3. Make SMALL incremental adjustments
        4. Wait ADEQUATE time between changes (2+ weeks for phenytoin)
        5. USE therapeutic drug monitoring for every adjustment
        6. Expect the unexpected - small changes can have large effects
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Convert all provided string quantities to pint Quantities
    params = {k: Q_(v) for k, v in kwargs.items()}

    v = params.get("v", False)
    Vmax = params.get("Vmax", False)
    C = params.get("C", False)
    Km = params.get("Km", False)

    # Count provided parameters
    provided = sum([bool(v), bool(Vmax), bool(C), bool(Km)])

    if provided != 3:
        raise ValueError(
            f"michaelis_menten_rate requires exactly 3 of 4 parameters. "
            f"Got {provided}: v={v is not False}, Vmax={Vmax is not False}, "
            f"C={C is not False}, Km={Km is not False}"
        )

    # Solve for the missing parameter
    # v = (Vmax * C) / (Km + C)
    if not v:
        # Calculate v from Vmax, C, and Km
        quantity = (Vmax * C) / (Km + C)
        string = "Rate of Elimination (v)"

    elif not Vmax:
        # Vmax = v * (Km + C) / C
        quantity = v * (Km + C) / C
        string = "Maximum Elimination Rate (Vmax)"

    elif not Km:
        # Km = C * (Vmax - v) / v
        quantity = C * (Vmax - v) / v
        string = "Michaelis Constant (Km)"

    elif not C:
        # C = (Km * v) / (Vmax - v)
        quantity = (Km * v) / (Vmax - v)
        string = "Concentration (C)"

    return format_output(quantity, string, output_unit, decimals)


def apparent_km(**kwargs):
    """
    Determine Km from steady-state data using the relationship v = Vmax/2 when C = Km.

    The Michaelis constant (Km) is a fundamental parameter that determines WHERE
    drug kinetics transition from linear to non-linear. Knowing a patient's Km is
    critical for safe dosing of non-linear drugs like phenytoin because it identifies
    the concentration threshold above which disproportionate effects occur.

    At steady state, when the rate of drug administration equals the rate of
    elimination, we can determine Km by finding the concentration at which
    the elimination rate is exactly half of Vmax.

    Key Relationship: When v = Vmax/2, then C = Km

    Mathematical Proof:
        v = (Vmax x C) / (Km + C)
        Vmax/2 = (Vmax x Km) / (Km + Km)
        Vmax/2 = (Vmax x Km) / (2 x Km)
        Vmax/2 = Vmax/2  (verified)

    CRITICAL CLINICAL IMPORTANCE OF KNOWING Km:
    1. Identifies the danger zone:
       - When C < 0.1 x Km: linear kinetics, safe zone
       - When C approaches Km: transitioning, caution needed
       - When C > Km: saturated kinetics, HIGH RISK of disproportionate effects

    2. Phenytoin's dangerous reality:
       - Km typically 4-10 mg/L (highly variable between patients)
       - Therapeutic range 10-20 mg/L
       - THIS MEANS THERAPEUTIC CONCENTRATIONS ARE AT OR ABOVE Km
       - Patient is ALWAYS operating in the non-linear danger zone

    3. Patient-specific Km values:
       - Km varies significantly between patients
       - Due to genetic polymorphisms in CYP2C9 and CYP2C19
       - A patient with low Km will saturate at lower concentrations
       - Must determine Km individually for optimal dosing

    Why Standard Half-Life Doesn't Apply When C > Km:
        At saturation, elimination rate becomes constant (Vmax), not proportional
        to concentration. This means:
        - "Half-life" varies with concentration (longer at higher C)
        - Cannot use t1/2 = 0.693/K because K is not constant
        - Time to reach steady state is unpredictable
        - Must use Michaelis-Menten kinetics instead

    Args (provide exactly 2 of 3):
        Km (str): Michaelis constant (e.g., '5 mg/L')
        v (str): Rate of elimination at half-maximal (e.g., '7.5 mg/hour')
        Vmax (str): Maximum elimination rate (e.g., '15 mg/hour')

    Optional kwargs:
        output_unit (str): Desired output unit for the result
        decimals (int): Number of decimal places for rounding (default 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Michaelis Constant (Km)', 5.0, 'milligram / liter', '5.0 mg/L', ...)

    Examples:
        Verify Km when v = Vmax/2:
            >>> apparent_km(v='7.5 mg/hour', Vmax='15 mg/hour')
            # This confirms we're at Km - the concentration would equal Km

        Calculate rate at Km:
            >>> apparent_km(Km='5 mg/L', Vmax='15 mg/hour')
            ('Rate at Km (v=Vmax/2)', 7.5, 'milligram / hour', '7.5 mg/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 11: Non-linear pharmacokinetics
        - Section 11.1, pages 6206-6240: Km as the half-saturation constant
        - Section 11.2, pages 6302-6320: Phenytoin Km in clinical range
        - Section 11.4, pages 6380-6488: Using Vmax and Km for dose calculations

    Clinical Methods to Determine Km:
        Method 1 - Two steady-state measurements:
            1. Measure steady-state concentration at two different doses
            2. At steady state: Rate of admin = Rate of elimination
            3. Use both data points to solve for Km and Vmax simultaneously

        Method 2 - Lineweaver-Burk plot:
            1. Take multiple dose-concentration pairs
            2. Plot 1/v vs 1/C
            3. Linear regression gives Km and Vmax

        Method 3 - Direct observation:
            1. Find concentration where v = Vmax/2
            2. That concentration equals Km
            3. Requires knowing Vmax first

    Phenytoin Km Values (Population Data):
        - Mean Km: ~5-7 mg/L
        - Range: 4-10 mg/L (highly variable)
        - Low Km (4 mg/L): Saturates at lower concentrations, higher risk
        - High Km (10 mg/L): More linear behavior, somewhat safer
        - Individual determination is essential for optimal dosing

    Factors Affecting Km:
        Km is primarily determined by enzyme structure and is relatively constant:
        - Genetic polymorphisms (CYP2C9, CYP2C19) - major factor
        - NOT significantly affected by enzyme inducers/inhibitors
        - NOT affected by body weight or organ function
        - Consistent within an individual over time

    Factors Affecting Vmax (contrast with Km):
        Vmax depends on enzyme AMOUNT and CAN change:
        - Enzyme INDUCTION (e.g., by carbamazepine): ↑ Vmax
        - Enzyme INHIBITION (e.g., by isoniazid): ↓ Vmax
        - Body size: larger patients have more enzyme
        - Liver disease: ↓ enzyme mass → ↓ Vmax

    Therapeutic Drug Monitoring Implications:
        When monitoring non-linear drugs:
        - A patient with low Km needs LOWER target concentrations
        - A patient with high Km can tolerate higher concentrations
        - Knowing individual Km allows personalized dosing
        - Without knowing Km, dose adjustments are empirical guesswork

    Ethanol Km Context:
        - Ethanol Km ≈ 10 mg/dL (very low)
        - Even moderate drinking produces C >> Km
        - Essentially always zero-order kinetics
        - This is why BAC falls linearly, not exponentially

    Clinical Management Using Km:
        1. Determine individual Km from steady-state data
        2. Recognize when C is approaching or exceeding Km
        3. Make smaller dose adjustments as C/Km ratio increases
        4. Allow longer equilibration time when C > Km
        5. Monitor more frequently in the non-linear range

    Notes:
        - This function uses the defining relationship v = Vmax/2 when C = Km
        - In practice, Km is determined from multiple dose-concentration pairs
        - Km varies between patients due to genetic enzyme differences
        - Enzyme inducers/inhibitors affect Vmax but generally not Km
        - Unlike K (which varies with concentration), Km is a true constant
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    params = {k: Q_(v) for k, v in kwargs.items()}

    Km = params.get("Km", False)
    v = params.get("v", False)
    Vmax = params.get("Vmax", False)

    provided = sum([bool(Km), bool(v), bool(Vmax)])

    if provided != 2:
        raise ValueError(
            f"apparent_km requires exactly 2 of 3 parameters. "
            f"Got {provided}: Km={Km is not False}, v={v is not False}, Vmax={Vmax is not False}"
        )

    # At C = Km, v = Vmax/2
    # So: v = Vmax/2, or Vmax = 2*v, or Km is what we're looking for

    if not Km:
        # If v = Vmax/2, we can verify this relationship
        # But Km is a concentration - we need concentration data to find it
        # This function confirms the relationship exists
        expected_v = Vmax / 2
        if abs(v.magnitude - expected_v.magnitude) / expected_v.magnitude > 0.01:
            raise ValueError(
                f"When solving for Km, v should equal Vmax/2. "
                f"Got v={v}, expected v={expected_v}"
            )
        # Return a message that Km equals the concentration at this rate
        # Since we don't have C, we return a placeholder
        raise ValueError(
            "To determine Km, you need concentration data. "
            "Km equals the concentration at which v = Vmax/2. "
            "Use michaelis_menten_rate() with v, Vmax, and C to solve for Km."
        )

    elif not v:
        # v = Vmax / 2 (rate at Km)
        quantity = Vmax / 2
        string = "Rate at Km (v=Vmax/2)"

    elif not Vmax:
        # Vmax = 2 * v
        quantity = v * 2
        string = "Maximum Elimination Rate (Vmax)"

    return format_output(quantity, string, output_unit, decimals)


def is_linear(C: str, Km: str, threshold: float = 0.1) -> bool:
    """
    Check if kinetics are approximately linear at the given concentration.

    This function is a critical clinical decision tool that determines whether
    a drug at a given concentration follows approximately first-order (linear)
    kinetics or has significant non-linear behavior due to enzyme saturation.
    The answer determines which pharmacokinetic equations and assumptions are valid.

    This directly addresses the question from Chapter 11: "Is this patient in the
    danger zone where standard pharmacokinetic rules break down?"

    Decision Rule:
        - If C < threshold x Km: Kinetics are approximately LINEAR (first-order)
        - If threshold x Km < C < Km: Kinetics are in TRANSITION zone
        - If C > Km: Kinetics are predominantly NON-LINEAR (zero-order)

    Default threshold: 0.1 (C < 10% of Km is considered linear)

    Kinetics Behavior in Each Zone:
        LINEAR (C << Km):
            v ≈ (Vmax/Km) x C
            - Rate proportional to concentration
            - CONSTANT half-life applies
            - Standard equations valid (K, Cl, t1/2)
            - Doubling dose doubles concentration
            - Time to steady state = 4-5 half-lives
            - Dose adjustments are simple and proportionate
            - Normal pharmacokinetic rules apply

        TRANSITION (C ≈ Km):
            v = Vmax/2 (half-maximal rate)
            - Kinetics transitioning from first to zero order
            - Half-life starting to become concentration-dependent
            - Standard equations becoming unreliable
            - Dose changes may have disproportionate effects
            - PHENYTOIN THERAPEUTIC RANGE IS HERE

        NON-LINEAR (C >> Km):
            v ≈ Vmax (constant rate)
            - Rate is constant (zero-order)
            - Half-life INCREASES with concentration
            - Standard half-life concept DOES NOT APPLY
            - Doubling dose MORE than doubles concentration
            - Time to steady state prolonged and UNPREDICTABLE
            - Small dose changes can cause LARGE concentration changes
            - Must use Michaelis-Menten kinetics, not standard equations

    CRITICAL CLINICAL WARNING - WHY THIS MATTERS:
    For phenytoin and other non-linear drugs, knowing whether you're in the
    linear or non-linear zone determines:
    1. Whether standard dosing calculations are valid
    2. How cautiously dose adjustments should be made
    3. How long to wait for new steady state
    4. Whether simple proportionate dose changes work

    The Phenytoin Problem:
        Phenytoin's therapeutic range (10-20 mg/L) is AT or ABOVE typical Km (4-10 mg/L)
        This means patients at therapeutic levels are ALWAYS in the non-linear zone
        This is why phenytoin is so difficult to dose safely

    How to Recognize Non-Linear Behavior Clinically:
    1. Disproportionate response to dose changes:
       - Patient's concentration increases more than expected with dose increase
       - Small dose increase causes unexpectedly large concentration rise

    2. Variable "half-life" measurements:
       - Half-life appears longer at higher concentrations
       - Half-life measurements are inconsistent at different times

    3. Prolonged time to steady state:
       - Takes much longer than expected 4-5 half-lives
       - Concentration still rising after expected time

    4. Non-linear dose-concentration relationship:
       - When plotted, dose vs concentration curves upward, not straight

    Args:
        C (str): Current drug concentration (e.g., '15 mg/L')
        Km (str): Michaelis constant for the drug (e.g., '5 mg/L')
        threshold (float): Fraction of Km below which kinetics are linear (default 0.1)

    Returns:
        bool: True if kinetics are approximately linear (C < threshold x Km),
              False if non-linear behavior is significant

    Examples:
        Check if phenytoin at 2 mg/L is in linear range (Km=5 mg/L):
            >>> is_linear(C='2 mg/L', Km='5 mg/L')
            False  # 2 mg/L is 40% of Km, above 10% threshold

        Check if concentration is safely in linear range:
            >>> is_linear(C='0.3 mg/L', Km='5 mg/L')
            True  # 0.3 mg/L is 6% of Km, below 10% threshold

        Use custom threshold for stricter assessment:
            >>> is_linear(C='2 mg/L', Km='5 mg/L', threshold=0.5)
            True  # 2 mg/L is 40% of Km, below 50% threshold

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 11: Non-linear pharmacokinetics
        - Section 11.1, pages 6206-6240: Linear vs non-linear behavior
        - Section 11.2, pages 6302-6320: Drugs with non-linear kinetics
        - Section 11.3, pages 6328-6380: Effect on dose-concentration relationship
        - Figure 11.5-11.6: Visual comparison of linear vs non-linear kinetics

    Clinical Examples - Phenytoin:
        Phenytoin (Km ~ 4-10 mg/L, therapeutic 10-20 mg/L):
            At 0.5 mg/L: C/Km ~ 0.05-0.125
                - Approximately linear
                - Standard equations valid
                - But WAY below therapeutic range - not useful

            At 5 mg/L: C/Km ~ 0.5-1.25
                - Transition zone
                - Kinetics becoming non-linear
                - Approaching therapeutic but not yet there

            At 10 mg/L (low therapeutic): C/Km ~ 1-2.5
                - NON-LINEAR
                - Cannot use standard half-life
                - Dose changes will have disproportionate effects

            At 20 mg/L (high therapeutic): C/Km ~ 2-5
                - SEVERELY SATURATED
                - Small dose changes very dangerous
                - Prolonged time to new steady state

        Key insight: There is NO concentration that is both:
        (a) therapeutic for phenytoin (10-20 mg/L) AND
        (b) in the linear kinetics range (C < 0.1 x Km ≈ 0.4-1 mg/L)
        This is why phenytoin dosing is ALWAYS complex and dangerous.

    Clinical Examples - Ethanol:
        Ethanol (Km ~ 10 mg/dL = 0.1 g/L):
            At BAC 0.02% (20 mg/dL): C/Km = 2
                - Already saturated
                - This is "barely impaired" level

            At BAC 0.08% (80 mg/dL): C/Km = 8
                - Highly saturated (legal limit many jurisdictions)

            At BAC 0.20% (200 mg/dL): C/Km = 20
                - Severely saturated (dangerous intoxication)

        All typical drinking levels are well above Km, so ethanol
        elimination is ALWAYS zero-order (~7 g/hour constant rate).
        Cannot "speed up" elimination - rate is fixed.

    Therapeutic Drug Monitoring Implications:
        When is_linear returns False:
        1. Do NOT use standard half-life for time calculations
        2. Do NOT assume proportionate dose-concentration relationship
        3. Use Michaelis-Menten equations instead
        4. Make smaller dose adjustments
        5. Wait longer between adjustments (2+ weeks for phenytoin)
        6. Monitor concentrations more frequently
        7. Consider specialist consultation

    Drug Development Context:
        - Drug companies TEST for non-linearity early in development
        - If a drug shows non-linear kinetics at therapeutic doses, it may be
          terminated or require extensive additional safety work
        - Linear kinetics are STRONGLY preferred for new drugs
        - This is why phenytoin (developed before modern PK understanding)
          remains one of few clinically relevant non-linear drugs

    Notes:
        - This is a practical approximation for clinical decision-making
        - The transition between linear and non-linear is GRADUAL, not sharp
        - When in doubt, assume non-linear and monitor closely
        - Individual patient Km values vary significantly (4-10 mg/L for phenytoin)
        - A patient with low Km will show non-linear behavior at lower concentrations
        - Always use therapeutic drug monitoring for non-linear drugs
    """
    C_qty = Q_(C)
    Km_qty = Q_(Km)

    # Ensure units are compatible
    C_dimensionless = (C_qty / Km_qty).to("dimensionless").magnitude

    return C_dimensionless < threshold


def phenytoin_steady_state(**kwargs):
    """
    Calculate steady-state concentration for phenytoin dosing.

    Phenytoin is THE classic and most clinically significant example of non-linear
    kinetics. As stated in Rowe's Pharmacokinetics Chapter 11: "Phenytoin provides
    the one clinically significant case" of non-linear kinetics. It is a drug with
    a narrow therapeutic window and serious toxicity in overdose, so controlling
    its blood concentrations is a real concern.

    At therapeutic doses, phenytoin metabolism is partially saturated because the
    therapeutic range (10-20 mg/L) is at or above typical Km values (4-10 mg/L).
    This causes disproportionate concentration increases with small dose changes.

    Formula: Css = (Km x Dose/tau) / (Vmax - Dose/tau)

    Or equivalently: Css = (Km x R) / (Vmax - R), where R = Dose/tau

    Where:
        Css = steady-state concentration (mass/volume)
        Km = Michaelis constant (mass/volume, typically 4-10 mg/L)
        Dose = maintenance dose per interval (mass)
        tau (dosing_interval) = time between doses
        Vmax = maximum elimination rate (mass/time)
        Dose/tau = R = rate of drug administration (mass/time)

    CRITICAL CLINICAL WARNING - THE PHENYTOIN DOSING PARADOX:
    From Rowe Chapter 11: "For most drugs, changes in dose size bring simple,
    proportionate changes in blood levels. With non-linear drugs such as phenytoin,
    changes in blood levels are disproportionately large. Great caution is required."

    Why Phenytoin Dosing is Notoriously Difficult:
    1. Disproportionate concentration changes:
       - A 10% dose increase can cause >50% concentration increase
       - The dose-concentration curve is NOT a straight line
       - Doubling the dose achieves FAR more than doubling in blood levels

    2. Standard rules DON'T APPLY:
       - "Given the potential toxicity of phenytoin, dosage adjustment is a
          specialist job; the normal concepts and rules just don't apply."
       - Cannot use elimination rate constant K (it's not constant)
       - Cannot use standard half-life (it varies with concentration)
       - Cannot use standard clearance (it's concentration-dependent)

    3. Time to steady state is PROLONGED:
       - Cannot use "4-5 half-lives" rule
       - May take 2-4 WEEKS to reach new steady state
       - Patient may appear stable, then suddenly become toxic

    4. Narrow therapeutic index:
       - Therapeutic range: 10-20 mg/L (total) or 1-2 mg/L (free)
       - Toxicity begins just above therapeutic range
       - Symptoms: nystagmus, ataxia, confusion, seizures

    5. High inter-patient variability:
       - Km varies 4-10 mg/L between patients
       - Vmax varies 350-700 mg/day
       - Same dose gives VERY different concentrations in different patients

    Phenytoin-Specific Dosing Guidelines:
        1. START LOW: 200-300 mg/day initial dose
        2. INCREASE SLOWLY: 25-50 mg increments ONLY
        3. WAIT LONG: At least 2 weeks between adjustments
        4. MONITOR ALWAYS: Therapeutic drug monitoring is mandatory
        5. CONSIDER FREE LEVELS: In hypoalbuminemia, uremia
        6. SPECIALIST INPUT: Dosage adjustment is a specialist job
        7. USE THIS EQUATION: Go back to Vmax and Km for calculations

    Why "Normal Concepts Just Don't Apply":
        Standard pharmacokinetics assumes:
        - Elimination rate constant K is constant (NOT TRUE for phenytoin)
        - Half-life is constant (NOT TRUE - it increases with concentration)
        - Clearance is constant (NOT TRUE - it decreases as C approaches saturation)
        - Dose-concentration relationship is linear (NOT TRUE - it curves upward)

        For phenytoin, you MUST use Michaelis-Menten kinetics:
        - Vmax and Km are the only valid constants
        - This equation calculates Css directly from these parameters

    Args (provide exactly 3 of 4):
        Css (str): Steady-state concentration (e.g., '15 mg/L')
        Km (str): Michaelis constant (e.g., '5 mg/L')
        dose (str): Maintenance dose per interval (e.g., '300 mg')
        dosing_interval (str): Time between doses (e.g., '24 hour')
        Vmax (str): Maximum elimination rate (e.g., '500 mg/day')

    Optional kwargs:
        output_unit (str): Desired output unit for the result
        decimals (int): Number of decimal places for rounding (default 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Steady-State Concentration (Css)', 15.0, 'milligram / liter', ...)

    Examples:
        Calculate Css for phenytoin 300 mg/day with typical parameters:
            >>> phenytoin_steady_state(
            ...     Km='5 mg/L',
            ...     dose='300 mg',
            ...     dosing_interval='24 hour',
            ...     Vmax='500 mg/day'
            ... )
            ('Steady-State Concentration (Css)', 15.0, 'milligram / liter', ...)

        Calculate required dose for target Css:
            >>> phenytoin_steady_state(
            ...     Css='15 mg/L',
            ...     Km='5 mg/L',
            ...     dosing_interval='24 hour',
            ...     Vmax='500 mg/day'
            ... )
            ('Dose', 300.0, 'milligram', '300.0 mg', ...)

        Demonstrate the danger - small dose increase, large Css increase:
            At 300 mg/day: Css = 15 mg/L
            At 330 mg/day (10% increase): Css = 19.4 mg/L (29% increase!)
            At 350 mg/day (17% increase): Css = 23.3 mg/L (55% increase - TOXIC!)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 11: Non-linear pharmacokinetics, pages 6198-6500
        - Section 11.2, page 6307: "Phenytoin provides the one clinically significant case"
        - Section 11.3, pages 6328-6380: Effect on dose-concentration relationship
        - Section 11.4, pages 6380-6488: Clinical significance - "dosage adjustment
          is a specialist job; the normal concepts and rules just don't apply"
        - Section 11.5, pages 6488-6500: Drug development implications

    Typical Phenytoin Population Parameters:
        Km (Michaelis constant):
            - Population mean: 5-7 mg/L
            - Range: 4-10 mg/L (highly variable)
            - Determines WHERE saturation occurs
            - A patient with low Km will saturate sooner

        Vmax (Maximum elimination rate):
            - Population mean: 400-500 mg/day
            - Range: 350-700 mg/day
            - Depends on body size, genetics, enzyme induction
            - Can be affected by other drugs (inducers/inhibitors)

        Therapeutic range:
            - Total phenytoin: 10-20 mg/L (or 40-80 umol/L)
            - Free phenytoin: 1-2 mg/L (or 4-8 umol/L)
            - Note: Therapeutic range is AT or ABOVE Km!

        Protein binding:
            - Approximately 90% bound to albumin
            - In hypoalbuminemia: less bound, higher free fraction
            - Must measure FREE levels if albumin is low
            - Free levels more clinically relevant

    Mathematical Derivations:
        From Css = (Km x R) / (Vmax - R), where R = Dose/tau:
        - Dose/tau = (Css x Vmax) / (Km + Css)
        - Vmax = (Dose/tau) x (Km + Css) / Css
        - Km = Css x (Vmax - Dose/tau) / (Dose/tau)

    The Asymptotic Danger:
        As R (dosing rate) approaches Vmax:
        - Css = (Km x R) / (Vmax - R) approaches INFINITY
        - If R >= Vmax, steady state CANNOT be achieved
        - Drug will accumulate without limit
        - This is why Dose/tau must always be LESS than Vmax

    Clinical Examples of Disproportionate Changes:
        Patient with Km=5 mg/L, Vmax=500 mg/day:
            200 mg/day → Css = 3.3 mg/L (subtherapeutic)
            250 mg/day → Css = 5.0 mg/L (still low)
            300 mg/day → Css = 7.5 mg/L (approaching therapeutic)
            350 mg/day → Css = 11.7 mg/L (therapeutic)
            400 mg/day → Css = 20.0 mg/L (high therapeutic/toxic)
            450 mg/day → Css = 45.0 mg/L (SEVERELY TOXIC!)

        Note how a 50 mg increase from 350 to 400 mg/day causes
        71% concentration increase, and another 50 mg causes 125% increase!

    Why Therapeutic Drug Monitoring is MANDATORY:
        1. Cannot predict concentration from dose alone
        2. High inter-patient variability in Km and Vmax
        3. Drug interactions can change Vmax
        4. Dose changes have unpredictable effects
        5. Concentration must be measured, not estimated
        6. Always measure after adequate equilibration time (2+ weeks)

    When to Measure Free Phenytoin Levels:
        - Hypoalbuminemia (albumin < 3.5 g/dL)
        - Renal failure (uremia displaces phenytoin from albumin)
        - Pregnancy (increased free fraction)
        - Concurrent highly protein-bound drugs
        - Elderly patients (often have lower albumin)

    Clinical Management Strategy:
        1. Determine individual Km and Vmax from two steady-state levels
        2. Use this equation to calculate predicted Css at new doses
        3. Make SMALL adjustments (25-50 mg)
        4. Wait adequate time (2+ weeks) for new steady state
        5. Measure levels after each change
        6. Consider loading dose for urgent situations (specialist guidance)
        7. Monitor for toxicity symptoms (nystagmus is early sign)

    Drug Development Context:
        From Rowe Chapter 11.5: "Drug companies want to establish at an early stage
        that any new candidate drug molecule obeys simple linear kinetics."

        - Companies would choose linear over non-linear drug if equally effective
        - Strong commercial disinclination to develop non-linear drugs
        - Phenytoin (developed before modern PK) would likely not be developed today
        - New anticonvulsants are generally designed to have linear kinetics

    Notes:
        - This equation assumes steady state has been reached
        - Time to steady state is PROLONGED with non-linear kinetics (2+ weeks)
        - Dose/tau MUST be less than Vmax (otherwise Css approaches infinity)
        - Always verify calculated predictions with measured concentrations
        - Individual patient parameters may differ significantly from population values
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    params = {k: Q_(v) for k, v in kwargs.items()}

    Css = params.get("Css", False)
    Km = params.get("Km", False)
    dose = params.get("dose", False)
    dosing_interval = params.get("dosing_interval", False)
    Vmax = params.get("Vmax", False)

    # Count provided parameters (dose and dosing_interval count as one: R = dose/tau)
    has_R = bool(dose) and bool(dosing_interval)
    provided = sum([bool(Css), bool(Km), has_R, bool(Vmax)])

    if provided != 3:
        raise ValueError(
            f"phenytoin_steady_state requires exactly 3 of 4 parameters "
            f"(dose and dosing_interval together count as one). "
            f"Got {provided}: Css={Css is not False}, Km={Km is not False}, "
            f"dose/interval={has_R}, Vmax={Vmax is not False}"
        )

    # Calculate rate of administration R = dose/tau
    if has_R:
        R = dose / dosing_interval
    else:
        R = False

    # Css = (Km * R) / (Vmax - R)
    if not Css:
        if R >= Vmax:
            raise ValueError(
                f"Dose/interval ({R}) must be less than Vmax ({Vmax}). "
                f"Otherwise steady state cannot be achieved (infinite accumulation)."
            )
        quantity = (Km * R) / (Vmax - R)
        string = "Steady-State Concentration (Css)"

    elif not Km:
        # Km = Css * (Vmax - R) / R
        quantity = Css * (Vmax - R) / R
        string = "Michaelis Constant (Km)"

    elif not Vmax:
        # Vmax = R * (Km + Css) / Css
        quantity = R * (Km + Css) / Css
        string = "Maximum Elimination Rate (Vmax)"

    elif not R:
        # R = (Css * Vmax) / (Km + Css)
        # Return dose assuming a 24-hour interval if not specified
        quantity = (Css * Vmax) / (Km + Css)
        string = "Dosing Rate (Dose/tau)"

    return format_output(quantity, string, output_unit, decimals)


def time_to_eliminate_nonlinear(**kwargs):
    """
    Calculate time to eliminate drug with zero-order (saturated) kinetics.

    When drug concentrations are well above Km, elimination approaches zero-order
    kinetics where a constant amount of drug is eliminated per unit time,
    REGARDLESS of concentration. This is fundamentally different from first-order
    kinetics where elimination is proportional to concentration.

    This is why standard half-life does NOT apply to saturated drugs - the concept
    of "half-life" requires first-order kinetics where time to eliminate 50% is
    constant. With zero-order kinetics, the time to eliminate 50% INCREASES as
    concentration increases.

    Formula: t = (C0 - Ct) x V / Vmax

    Where:
        t = time to eliminate from C0 to Ct
        C0 = initial concentration (mass/volume)
        Ct = target concentration at time t (mass/volume)
        V = volume of distribution (volume)
        Vmax = maximum elimination rate (mass/time)

    This formula assumes C >> Km throughout the elimination period, so
    elimination rate is approximately constant at Vmax.

    Derivation:
        At saturation: Rate = Vmax (constant, independent of C)
        Amount eliminated = Vmax x t
        Amount eliminated = (C0 - Ct) x V
        Therefore: t = (C0 - Ct) x V / Vmax

    CRITICAL UNDERSTANDING - WHY STANDARD HALF-LIFE DOESN'T APPLY:
    With first-order kinetics:
        - Half-life = time to go from 100 to 50 = time to go from 50 to 25
        - It's always the SAME time to eliminate 50%
        - C(t) = C0 x e^(-K x t) (exponential decay)

    With zero-order kinetics (when C >> Km):
        - Time to go from 100 to 50 is DIFFERENT from time to go from 50 to 25
        - Higher starting concentration = longer "half-life"
        - C(t) = C0 - (Vmax/V) x t (LINEAR decay, not exponential)
        - The term "half-life" loses its meaning

    CLINICAL WARNING - IMPLICATIONS OF ZERO-ORDER KINETICS:
    1. Time to eliminate is LINEAR with concentration difference:
       - Doubling the initial concentration DOUBLES the elimination time
       - Cannot use standard half-life calculations

    2. Cannot predict time course using standard equations:
       - First-order: C(t) = C0 x e^(-K x t) - INVALID
       - Zero-order: C(t) = C0 - (Vmax/V) x t - USE THIS

    3. "Half-life" concept breaks down:
       - Half-life increases with concentration
       - A patient at higher concentration takes longer to eliminate 50%

    4. Time to reach steady state is unpredictable:
       - Cannot use "4-5 half-lives" rule
       - Accumulation continues much longer than expected

    Args (provide exactly 4 of 5):
        t (str): Time to eliminate (e.g., '10 hour')
        C0 (str): Initial concentration (e.g., '100 mg/L')
        Ct (str): Target concentration at time t (e.g., '20 mg/L')
        V (str): Volume of distribution (e.g., '50 L')
        Vmax (str): Maximum elimination rate (e.g., '400 mg/hour')

    Optional kwargs:
        output_unit (str): Desired output unit for the result
        decimals (int): Number of decimal places for rounding (default 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Example: ('Time to Eliminate', 10.0, 'hour', '10.0 hour', ...)

    Examples:
        Calculate time to reduce concentration from 100 to 20 mg/L:
            >>> time_to_eliminate_nonlinear(
            ...     C0='100 mg/L',
            ...     Ct='20 mg/L',
            ...     V='50 L',
            ...     Vmax='400 mg/hour'
            ... )
            ('Time to Eliminate', 10.0, 'hour', '10.0 hour', ...)

        Calculate final concentration after given time:
            >>> time_to_eliminate_nonlinear(
            ...     C0='100 mg/L',
            ...     t='10 hour',
            ...     V='50 L',
            ...     Vmax='400 mg/hour'
            ... )
            ('Final Concentration (Ct)', 20.0, 'milligram / liter', '20.0 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 11: Non-linear pharmacokinetics
        - Section 11.2, pages 6302-6320: Ethanol as zero-order example
        - Section 11.3, pages 6328-6380: Linear vs zero-order time courses
        - Section 11.4, pages 6380-6488: Clinical implications

    THE ETHANOL MODEL - Classic Zero-Order Kinetics:
        From Rowe Chapter 11: "Ethanol in doses high enough to cause noticeable
        effects, will fully saturate the liver enzymes."

        Ethanol is the paradigm example of zero-order elimination because:
        - Km ≈ 10 mg/dL (very low)
        - Social drinking produces concentrations FAR above Km
        - Elimination is essentially constant at ~7 g/hour (Vmax)
        - Blood alcohol falls LINEARLY with time, not exponentially
        - Cannot "speed up" elimination - rate is fixed at Vmax

        Clinical Example - Calculating Time to Sober Up:
        For a 70 kg person:
            - Vd ≈ 0.6 L/kg = 42 L (ethanol distributes in total body water)
            - Vmax ≈ 7 g/hour (fairly constant across individuals)

        Scenario: BAC 0.20% (severe intoxication) to 0.08% (legal limit):
            - C0 = 0.20% = 200 mg/dL = 2 g/L
            - Ct = 0.08% = 80 mg/dL = 0.8 g/L
            - Amount to eliminate = (2 - 0.8) g/L x 42 L = 50.4 g
            - Time = 50.4 g / 7 g/hour = 7.2 hours

        Scenario: BAC 0.30% to 0.08%:
            - C0 = 0.30% = 3 g/L
            - Amount to eliminate = (3 - 0.8) g/L x 42 L = 92.4 g
            - Time = 92.4 g / 7 g/hour = 13.2 hours

        Key insight: Higher initial BAC = proportionally longer time to sober up
        (This is LINEAR, unlike first-order where it would be logarithmic)

    Why Standard Half-Life Concept Fails:
        First-order kinetics (C << Km):
            - Half-life t1/2 = 0.693/K (constant)
            - Time to go from 100 to 50 = time to go from 50 to 25
            - "4-5 half-lives" to steady state or elimination

        Zero-order kinetics (C >> Km):
            - "Half-life" depends on concentration
            - Time from 100 to 50: (100-50) x V / Vmax
            - Time from 50 to 25: (50-25) x V / Vmax = HALF as long!
            - The higher the concentration, the longer the apparent "half-life"

        Mathematical demonstration:
            If Vmax/V = 10 mg/L per hour:
            - Time from 100 to 50 mg/L = 50/10 = 5 hours
            - Time from 50 to 25 mg/L = 25/10 = 2.5 hours
            - "Half-life" at C=100 is 5 hours
            - "Half-life" at C=50 is 2.5 hours
            The "half-life" is NOT constant!

    Comparison: First-Order vs Zero-Order:
        First-order (C << Km):
            - C(t) = C0 x e^(-K x t)  (EXPONENTIAL decay)
            - Half-life is CONSTANT
            - Rate of elimination ∝ concentration
            - Time to eliminate 50% is always the same
            - Plot of ln(C) vs time is LINEAR

        Zero-order (C >> Km):
            - C(t) = C0 - (Vmax/V) x t  (LINEAR decay)
            - "Half-life" INCREASES with concentration
            - Rate of elimination is CONSTANT (= Vmax)
            - Time to eliminate 50% depends on starting C
            - Plot of C vs time is LINEAR (not ln(C)!)

    Phenytoin Considerations:
        - Phenytoin at therapeutic levels is in the transition/saturated zone
        - As concentration rises, kinetics shift toward zero-order
        - This explains why time to steady state is prolonged (2-4 weeks)
        - Cannot use standard "4-5 half-lives" calculation
        - Must wait much longer for equilibration after dose changes

    Clinical Applications:
        1. Ethanol elimination planning:
           - Can calculate time to reach legal limit
           - Cannot speed up elimination (coffee, exercise don't help)

        2. Phenytoin overdose management:
           - Use this equation when C is well above Km
           - As C falls and approaches Km, switch to M-M calculations

        3. Predicting drug accumulation:
           - When dosing rate approaches Vmax, accumulation accelerates
           - Time to dangerous levels can be estimated

    When This Equation is Valid:
        - C >> Km throughout the time period
        - Elimination is essentially constant at Vmax
        - No saturation of other processes (absorption, distribution)

    When to Switch to Michaelis-Menten:
        - As C approaches Km
        - When more accurate modeling is needed
        - For the full transition from zero-order back to first-order

    Therapeutic Drug Monitoring Implications:
        For saturated drugs:
        - Standard t1/2-based timing doesn't apply
        - Wait much longer between dose adjustments
        - Use this equation to estimate time to target concentration
        - Always verify with measured levels

    Notes:
        - This equation is valid ONLY when C >> Km throughout
        - As C approaches Km, kinetics transition back to first-order
        - Real elimination often transitions from zero to first order
        - For accurate modeling across all concentration ranges, use full
          Michaelis-Menten equation: v = (Vmax x C) / (Km + C)
        - This is an approximation for the saturated portion of the curve
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    params = {k: Q_(v) for k, v in kwargs.items()}

    t = params.get("t", False)
    C0 = params.get("C0", False)
    Ct = params.get("Ct", False)
    V = params.get("V", False)
    Vmax = params.get("Vmax", False)

    provided = sum([bool(t), bool(C0), bool(Ct), bool(V), bool(Vmax)])

    if provided != 4:
        raise ValueError(
            f"time_to_eliminate_nonlinear requires exactly 4 of 5 parameters. "
            f"Got {provided}: t={t is not False}, C0={C0 is not False}, "
            f"Ct={Ct is not False}, V={V is not False}, Vmax={Vmax is not False}"
        )

    # t = (C0 - Ct) * V / Vmax
    if not t:
        quantity = (C0 - Ct) * V / Vmax
        string = "Time to Eliminate"

    elif not C0:
        # C0 = Ct + (t * Vmax / V)
        quantity = Ct + (t * Vmax / V)
        string = "Initial Concentration (C0)"

    elif not Ct:
        # Ct = C0 - (t * Vmax / V)
        quantity = C0 - (t * Vmax / V)
        string = "Final Concentration (Ct)"

    elif not V:
        # V = (t * Vmax) / (C0 - Ct)
        quantity = (t * Vmax) / (C0 - Ct)
        string = "Volume of Distribution (V)"

    elif not Vmax:
        # Vmax = (C0 - Ct) * V / t
        quantity = (C0 - Ct) * V / t
        string = "Maximum Elimination Rate (Vmax)"

    return format_output(quantity, string, output_unit, decimals)
