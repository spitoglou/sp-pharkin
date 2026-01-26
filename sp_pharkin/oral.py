"""
Oral (extravascular) administration pharmacokinetics calculations.

These functions implement calculations for oral/extravascular drug administration
from "Pharmacokinetics" by Philip Rowe, Chapter 9 (Extravascular Administration).

Oral administration involves both absorption and elimination processes occurring
simultaneously, creating more complex pharmacokinetic profiles than IV administration.
The key additional parameter is the absorption rate constant (Ka).
"""

import math

from pint import UnitRegistry

from .lib import format_output, generic_a_eq_b_x_c

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]

# Tolerance for comparing Ka and K (to handle floating point comparisons)
_RATE_CONSTANT_TOLERANCE = 1e-6


def absorption_rate(**kwargs):
    """
    Calculate rate of drug absorption, amount awaiting absorption, or absorption rate constant.

    During oral administration, drug in the gastrointestinal tract (the "absorption
    compartment") is absorbed into systemic circulation at a rate proportional to
    the amount remaining to be absorbed. This is first-order absorption kinetics.

    The physiological basis: As drug dissolves and becomes available in the GI tract,
    it crosses the intestinal membrane into the portal circulation. The rate of this
    transfer depends on both the amount of dissolved drug available AND the efficiency
    of absorption (determined by drug properties, formulation, and patient factors).

    Formula: Rate of Absorption = Amount Awaiting Absorption (Aa) x Ka

    Where:
        - Rate of Absorption: Amount of drug absorbed per unit time (e.g., mg/hour)
        - Aa (Amount Awaiting): Drug remaining in GI tract to be absorbed (e.g., mg)
        - Ka: Absorption rate constant (e.g., 1/hour)

    The absorption rate constant (Ka) reflects the combined effects of:
        - Drug dissolution rate from the dosage form
        - GI membrane permeability (lipophilicity, molecular size, ionization)
        - GI motility and transit time through absorptive regions
        - Splanchnic blood flow (carries absorbed drug away from absorption site)
        - Formulation factors (tablet vs solution vs capsule vs sustained-release)
        - P-glycoprotein efflux (can pump drug back into GI lumen)

    Args (provide exactly 2 of 3):
        rate_of_absorption (str): Rate at which drug enters systemic circulation
            (e.g., '50 mg/hour')
        amount_awaiting (str): Amount of drug still in GI tract awaiting absorption
            (e.g., '200 mg')
        Ka (str): Absorption rate constant (e.g., '0.25 1/hour')

    Optional kwargs:
        output_unit (str or False): Desired output unit, or False to keep original
        decimals (int): Number of decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the pair provided.
        - name (str): Parameter name ('Rate of Absorption', 'Amount Awaiting Absorption', etc.)
        - magnitude (float): Numeric value only
        - unit_string (str): Unit as string
        - formatted_string (str): Complete formatted result
        - pint_quantity: Full Quantity object for further calculations

    Examples:
        Calculate absorption rate from 200 mg awaiting with Ka=0.25/hour:
            >>> absorption_rate(amount_awaiting='200 mg', Ka='0.25 1/hour')
            ('Rate of Absorption', 50.0, 'milligram / hour', '50.0 mg/hour', ...)

        Calculate amount awaiting from rate and Ka:
            >>> absorption_rate(rate_of_absorption='50 mg/hour', Ka='0.25 1/hour')
            ('Amount Awaiting Absorption', 200.0, 'milligram', '200.0 mg', ...)

        Calculate Ka from rate and amount:
            >>> absorption_rate(rate_of_absorption='50 mg/hour', amount_awaiting='200 mg')
            ('Absorption Rate Constant (Ka)', 0.25, '1 / hour', '0.25 1/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 9: Extravascular Administration
        - Section 9.1-9.2: Absorption kinetics and the absorption rate constant
        - Section 9.3: Changing the rate of absorption via formulation
        - Pages 4690-4710: Derivation of first-order absorption model
        - Figure 9.1: Schematic of drug handling with oral administration

    Typical Ka Values by Drug and Formulation:

        Very Rapid Absorption (Ka > 3/hour, t1/2,abs < 15 min):
            - Oral solutions/elixirs: Ka = 3-6 /hour
            - Sublingual nitroglycerin: Ka = 6-12 /hour (bypasses first-pass)
            - Alcohol (ethanol): Ka = 6-10 /hour (small, lipophilic molecule)
            - Acetaminophen solution: Ka = 4-6 /hour

        Rapid Absorption (Ka = 1-3/hour, t1/2,abs = 15-45 min):
            - Aspirin tablets: Ka = 1.5-3 /hour
            - Acetaminophen tablets: Ka = 1-2 /hour
            - Ibuprofen: Ka = 1-2.5 /hour
            - Most immediate-release tablets: Ka = 1-2 /hour
            - Capsules (immediate-release): Ka = 1-2 /hour

        Moderate Absorption (Ka = 0.3-1/hour, t1/2,abs = 45 min - 2 hours):
            - Propranolol: Ka = 0.5-1 /hour
            - Digoxin: Ka = 0.5-1 /hour (limited by low lipid solubility)
            - Theophylline tablets: Ka = 0.5-1 /hour
            - Most standard tablets: Ka = 0.5-1.5 /hour
            - Enteric-coated (after coating dissolves): Ka = 0.3-0.7 /hour

        Slow Absorption (Ka < 0.3/hour, t1/2,abs > 2 hours):
            - Theophylline sustained-release: Ka = 0.1-0.3 /hour
            - Morphine extended-release: Ka = 0.1-0.2 /hour
            - Griseofulvin: Ka = 0.1-0.2 /hour (very poorly soluble)
            - Sustained-release formulations: Ka = 0.05-0.3 /hour
            - Transdermal patches: Ka = 0.02-0.1 /hour
            - Depot IM injections: Ka = 0.01-0.1 /hour

    Formulation Effects on Ka (Same Drug, Different Forms):

        Theophylline Example (from Rowe, Chapter 9):
            - Solution: Ka approximately 2-4 /hour
            - Immediate-release tablet: Ka approximately 0.5-1 /hour
            - Sustained-release: Ka approximately 0.1-0.3 /hour
            The same drug can have 10-40 fold differences in Ka based on formulation.

        General Formulation Hierarchy (fastest to slowest absorption):
            1. IV injection (instantaneous, Ka = infinity theoretically)
            2. Sublingual/buccal tablets (rapid, bypasses first-pass)
            3. Oral solutions/elixirs (no dissolution step needed)
            4. Suspensions (small particles, dissolve quickly)
            5. Capsules (gelatin shell dissolves, then drug)
            6. Immediate-release tablets (disintegration + dissolution)
            7. Enteric-coated tablets (delayed release in intestine)
            8. Sustained/extended-release (controlled release over time)
            9. Transdermal patches (slow, continuous input)
            10. Depot IM/SC injections (very slow release from depot)

    Food-Drug Interactions Affecting Absorption Rate:

        Food DELAYS Absorption (decreased Ka, increased Tmax):
            - Most drugs: Food slows gastric emptying, delaying drug delivery
              to absorptive small intestine
            - Acetaminophen: Ka reduced 50% with high-fat meal
            - Aspirin: Peak delayed 1-2 hours with food
            - Penicillins: Absorption reduced and delayed with food
            - Fluoroquinolones: Delayed by dairy (Ca++ chelation)

        Food ENHANCES Absorption (increased Ka or F):
            - Lipophilic drugs: Fat in meal enhances solubilization
            - Griseofulvin: AUC increased 2-3 fold with fatty meal
            - Propranolol: Bioavailability increased with food (reduced first-pass)
            - Carbamazepine: Better absorbed with food

        Food Has MINIMAL Effect:
            - Theophylline: Generally unaffected (take with or without food)
            - Some sustained-release: Designed for food-independent release
            - Digoxin: Minimal food effect on extent (slight delay in rate)

        Clinical Recommendation Patterns:
            - "Take on empty stomach": Drug absorption impaired by food
            - "Take with food": Drug causes GI upset or needs fat for absorption
            - "Take with or without food": Minimal clinically significant effect

    Factors Affecting Ka - Detailed Mechanisms:

        Drug Properties:
            - Lipophilicity (LogP):
              * LogP > 2: Rapid passive diffusion across membranes
              * LogP < 0: Slow absorption, may need transporters
              * Very high LogP (>5): May be too lipophilic, poor dissolution
            - Molecular weight:
              * < 500 Da: Generally good absorption potential
              * > 500 Da: Reduced passive permeability
            - pKa and ionization:
              * Weak acids (aspirin, NSAIDs): Better absorbed in stomach
              * Weak bases (propranolol): Better absorbed in intestine
              * Ionized forms: Poor membrane permeability
            - Solubility:
              * BCS Class I (high solubility, high permeability): Ka limited by permeability
              * BCS Class II (low solubility, high permeability): Ka limited by dissolution
              * BCS Class III (high solubility, low permeability): Ka limited by permeability
              * BCS Class IV (low solubility, low permeability): Variable, often low Ka

        Formulation Factors:
            - Particle size: Micronization increases surface area, faster dissolution
            - Salt form: Different salts have different dissolution rates
            - Excipients:
              * Surfactants can enhance wetting and dissolution
              * Some binders can slow disintegration
              * Enteric coatings deliberately delay release
            - Manufacturing process: Compression force, granulation method affect disintegration

        Patient/Physiological Factors:
            - Gastric emptying rate:
              * Faster emptying = faster delivery to absorptive intestine
              * Delayed by food, anticholinergics, opioids
              * Accelerated by metoclopramide, erythromycin
            - GI motility:
              * Hypermotility can reduce contact time (diarrhea)
              * Hypomotility can increase absorption time
            - GI pH:
              * Antacids, PPIs raise gastric pH
              * Affects dissolution of pH-sensitive drugs
            - Splanchnic blood flow:
              * Reduced in shock, heart failure
              * Exercise can reduce GI blood flow
            - P-glycoprotein expression:
              * High P-gp = drug efflux back into lumen
              * Affected by genetics and drug interactions

    Clinical Applications of Ka:

        Drug Development:
            - Ka is measured during formulation development
            - Used to compare bioequivalence of generic vs brand
            - Determines if sustained-release formulation is feasible
            - Required for population PK modeling

        Clinical Practice:
            - Ka rarely calculated directly in routine care
            - Used indirectly when choosing formulation:
              * Acute pain: Choose rapid-release for quick onset
              * Chronic conditions: Choose sustained-release for steady levels
            - Guides timing of doses relative to meals
            - Explains inter-patient variability in drug response

        Dosing Implications:
            - High Ka: Quick onset, potentially high Cmax (toxicity risk)
            - Low Ka: Slow onset (may need loading dose), lower Cmax
            - Sustained-release (low Ka by design): Less dosing frequency,
              better adherence, reduced peak-trough fluctuation

    Clinical Significance:
        - Ka determines onset of action for oral medications
        - Higher Ka = faster onset but potentially higher Cmax
        - Lower Ka = slower onset, more prolonged absorption
        - When Ka < K: "flip-flop" kinetics occurs (see flip_flop_check)
        - Ka variability contributes to inter-individual response differences

    Mathematical Relationships:
        - Absorption half-life: t1/2,abs = 0.693 / Ka
        - Time to peak (Tmax): Tmax = ln(Ka/K) / (Ka - K)
        - At Tmax: Rate of absorption = Rate of elimination

    Notes:
        - Assumes first-order absorption kinetics (rate proportional to amount)
        - Does not account for lag time before absorption begins
        - Does not account for bioavailability (F) - calculate separately
        - For zero-order absorption (some controlled release), different models apply
        - Real absorption may be more complex (multiple absorption sites, saturable transport)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("rate_of_absorption", False)
    b = kwargs.get("amount_awaiting", False)
    c = kwargs.get("Ka", False)

    string, quantity = generic_a_eq_b_x_c(
        a,
        b,
        c,
        ["Rate of Absorption", "Amount Awaiting Absorption", "Absorption Rate Constant (Ka)"],
    )

    return format_output(quantity, string, output_unit, decimals)


def tmax(**kwargs):
    """
    Calculate time to reach maximum plasma concentration after oral administration.

    After oral dosing, plasma concentration rises during the absorption phase, reaches
    a peak (Cmax at time Tmax), then declines during the elimination phase. Tmax is
    the time from administration to peak concentration.

    The physiological basis: At Tmax, the instantaneous rate of drug absorption equals
    the instantaneous rate of drug elimination. Before Tmax, absorption exceeds elimination
    (concentration rises). After Tmax, elimination exceeds absorption (concentration falls).
    This balance point occurs at Tmax regardless of dose or bioavailability.

    Formula: Tmax = ln(Ka/K) / (Ka - K)

    Special case when Ka = K: Tmax = 1/K (derived via L'Hopital's rule)

    Where:
        - Tmax: Time to maximum concentration (e.g., hours)
        - Ka: Absorption rate constant (e.g., 1/hour)
        - K: Elimination rate constant (e.g., 1/hour)
        - ln: Natural logarithm

    This formula assumes:
        - One-compartment model with instantaneous distribution
        - First-order absorption and elimination kinetics
        - No lag time before absorption begins
        - Complete drug availability for absorption (no precipitation)

    Args (provide both parameters):
        Ka (str): Absorption rate constant (e.g., '1.5 1/hour')
        K (str): Elimination rate constant (e.g., '0.2 1/hour')

    Optional kwargs:
        output_unit (str or False): Desired output unit, or False to keep original
        decimals (int): Number of decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        - name (str): 'Time to Maximum Concentration (Tmax)'
        - magnitude (float): Numeric value of Tmax
        - unit_string (str): Unit as string (e.g., 'hour')
        - formatted_string (str): Complete formatted result
        - pint_quantity: Full Quantity object for further calculations

    Examples:
        Calculate Tmax with Ka=1.5/hour and K=0.2/hour:
            >>> tmax(Ka='1.5 1/hour', K='0.2 1/hour')
            ('Time to Maximum Concentration (Tmax)', 1.55, 'hour', '1.55 hour', ...)

        Calculate Tmax when Ka approximately equals K:
            >>> tmax(Ka='0.2 1/hour', K='0.2 1/hour')
            ('Time to Maximum Concentration (Tmax)', 5.0, 'hour', '5.0 hour', ...)
            # Uses special case formula: Tmax = 1/K

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 9: Extravascular Administration
        - Section 9.4: Cmax and Tmax (pages 4864-4868)
        - Section 9.3: Changing the rate of absorption
        - Figure 9.2: Concentration versus time curve showing Tmax (page 4739)
        - Figure 9.3: Effect of changing absorption rate on Tmax (page 4825)

    Typical Tmax Values by Drug Class and Formulation:

        Very Rapid Absorption (Tmax < 30 minutes):
            - Sublingual nitroglycerin: 2-5 minutes
            - Alcohol (fasted): 15-30 minutes
            - Oral solutions (most drugs): 15-30 minutes
            - Acetaminophen liquid: 15-30 minutes

        Rapid Absorption (Tmax 0.5-1.5 hours):
            - Aspirin tablets: 0.5-1 hour
            - Acetaminophen tablets: 0.5-1 hour
            - Ibuprofen: 1-2 hours
            - Most NSAIDs: 1-2 hours
            - Immediate-release tablets (typical): 1-2 hours

        Moderate Absorption (Tmax 1.5-4 hours):
            - Propranolol: 1-2 hours
            - Digoxin: 1-3 hours
            - Theophylline (immediate-release): 1-2 hours
            - Enteric-coated aspirin: 3-4 hours (delayed by coating)
            - Many antibiotics: 1-3 hours

        Slow Absorption (Tmax > 4 hours):
            - Theophylline sustained-release: 4-8 hours
            - Morphine extended-release: 3-6 hours
            - Griseofulvin: 4-8 hours (dissolution-limited)
            - Sustained-release formulations: 4-12 hours
            - Transdermal patches: 12-24+ hours to steady input

    Clinical Drug Examples with Typical Tmax:

        Analgesics:
            - Acetaminophen tablet: 0.5-1 hour (quick onset for pain)
            - Aspirin: 1-2 hours
            - Ibuprofen: 1-2 hours
            - Morphine IR: 1-1.5 hours
            - Morphine ER: 3-6 hours (by design)

        Cardiovascular:
            - Propranolol IR: 1-2 hours
            - Propranolol LA: 6-14 hours
            - Digoxin: 1-3 hours (but effect delayed due to distribution)
            - Metoprolol IR: 1-2 hours
            - Metoprolol ER: 7-8 hours

        Anti-infectives:
            - Amoxicillin: 1-2 hours
            - Ciprofloxacin: 1-2 hours
            - Fluconazole: 1-2 hours

        CNS:
            - Diazepam: 0.5-1 hour (very lipophilic)
            - Phenytoin: 4-12 hours (slow, erratic absorption)

    Formulation Effects on Tmax (Same Drug, Different Forms):

        From Rowe Chapter 9 (Figure 9.3):
            - Rapid-release: Lower Tmax, higher Cmax, earlier decline
            - Slow-release: Higher Tmax, lower Cmax, sustained levels

        Theophylline Example:
            - Solution: Tmax approximately 0.5-1 hour
            - Immediate-release tablet: Tmax approximately 1-2 hours
            - Sustained-release: Tmax approximately 4-8 hours

        Clinical Implications:
            - Headache relief: Want rapid Tmax (quick hit, high peak)
            - Chronic asthma (theophylline): Want delayed Tmax (sustained, moderate levels)
            - Sleep aid: Time Tmax to coincide with desired sleep onset

    Food Effects on Tmax:

        Food Generally DELAYS Tmax:
            - Mechanism: Food delays gastric emptying, keeping drug in stomach
              longer before delivery to absorptive small intestine
            - Acetaminophen: Tmax delayed from 0.5h to 2-3h with fatty meal
            - Aspirin: Tmax delayed 1-2 hours with food
            - Most drugs show 0.5-2 hour delay with standard meal

        Food May Have MINIMAL Effect on Tmax:
            - Some sustained-release formulations designed for consistency
            - Drugs absorbed throughout GI tract
            - Theophylline: Minimal Tmax change with food

        Food May DECREASE Tmax (rare):
            - Highly lipophilic drugs may dissolve faster in fat
            - Some drugs have enhanced gastric emptying with specific meals

        Clinical Guidance:
            - "Take on empty stomach" often means faster Tmax is desired
            - "Take with food" may be for GI tolerance, accepting delayed Tmax
            - Consistency matters: Same conditions each dose for predictability

    Factors Affecting Tmax - Mathematical Relationships:

        Ka/K Ratio Effects:
            - Ka >> K (ratio > 10): Tmax approximately 1/Ka (absorption determines peak time)
            - Ka = K (ratio = 1): Tmax = 1/K = 1/Ka (special case)
            - Ka << K (ratio < 0.1): Tmax approximately 1/Ka (flip-flop, absorption still determines)

        Numerical Examples (K = 0.1/hour fixed):
            - Ka = 2/hour (ratio 20): Tmax = 1.6 hours
            - Ka = 1/hour (ratio 10): Tmax = 2.6 hours
            - Ka = 0.5/hour (ratio 5): Tmax = 4.0 hours
            - Ka = 0.2/hour (ratio 2): Tmax = 9.2 hours
            - Ka = 0.1/hour (ratio 1): Tmax = 10 hours (special case, equals 1/K)
            - Ka = 0.05/hour (ratio 0.5): Tmax = 13.9 hours (flip-flop territory)

    Clinical Interpretation of Tmax:

        Therapeutic Timing:
            - PRN analgesics: Take 30-60 min before anticipated pain (procedure)
            - Sleep medications: Take 30-60 min before desired sleep time
            - Antihypertensives: May time to counter morning BP surge
            - Antibiotics: Less critical (steady-state matters more)

        Bioequivalence Assessment:
            - Tmax is key parameter in generic drug approval
            - 90% CI for Tmax ratio must fall within acceptable range
            - Different Tmax may indicate formulation differences

        Patient Counseling:
            - "You may not feel effect immediately" - explain expected Tmax
            - "Don't take more if first dose hasn't worked yet" - wait for Tmax
            - Food effects: Explain how meals affect timing of drug action

    Relationship to Other Parameters:
        - Independent of dose (first-order kinetics)
        - Independent of bioavailability (F)
        - Dependent ONLY on Ka and K ratio
        - Cmax depends on Tmax (concentration at time Tmax)
        - AUC is independent of Tmax (total exposure unchanged)

    Mathematical Derivation:
        At Tmax, dC/dt = 0 (concentration neither rising nor falling)
        Setting absorption rate = elimination rate:
            Ka * Aa(Tmax) = K * Ab(Tmax)
        Solving for the time when this occurs yields:
            Tmax = ln(Ka/K) / (Ka - K)

    Notes:
        - The formula ln(Ka/K)/(Ka-K) is mathematically undefined when Ka = K
        - When Ka is approximately equal to K (within tolerance), uses Tmax = 1/K
        - Always returns positive Tmax (negative would be non-physical)
        - Assumes instantaneous drug availability for absorption (no lag time)
        - Real absorption may have lag time (add to calculated Tmax for total time)
        - For flip-flop kinetics (Ka < K), formula still works but interpretation changes
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    Ka_qty = Q_(kwargs.get("Ka"))
    K_qty = Q_(kwargs.get("K"))

    # Convert to compatible units for comparison
    # Use Ka's units as the reference
    K_converted = K_qty.to(Ka_qty.units)

    Ka_magnitude = Ka_qty.magnitude
    K_magnitude = K_converted.magnitude

    # Check if Ka approximately equals K (within tolerance)
    if abs(Ka_magnitude - K_magnitude) < _RATE_CONSTANT_TOLERANCE:
        # Special case: Tmax = 1/K when Ka = K
        tmax_qty = 1 / K_qty
    else:
        # General formula: Tmax = ln(Ka/K) / (Ka - K)
        ratio = Ka_magnitude / K_magnitude
        ln_ratio = math.log(ratio)

        # Calculate (Ka - K) in original units
        difference = Ka_qty - K_converted

        # Tmax = ln(Ka/K) / (Ka - K)
        # ln(Ka/K) is dimensionless, so result has units of 1/(1/time) = time
        tmax_qty = ln_ratio / difference

    return format_output(
        tmax_qty,
        "Time to Maximum Concentration (Tmax)",
        output_unit,
        decimals
    )


def cmax(**kwargs):
    """
    Calculate maximum plasma concentration after oral administration.

    Cmax is the peak plasma concentration achieved after oral dosing. It occurs at
    time Tmax and represents the instantaneous balance point where the rate of drug
    absorption equals the rate of drug elimination. Cmax is critical for both efficacy
    (must exceed minimum effective concentration, MEC) and safety (must not exceed
    minimum toxic concentration, MTC).

    The physiological basis: At early times, absorption exceeds elimination and
    concentration rises. At Tmax, the rates are equal (dC/dt = 0). After Tmax,
    elimination exceeds the declining absorption rate. Cmax represents the highest
    concentration achieved and is a key determinant of both therapeutic effect and
    toxicity risk.

    Formula: Cmax = (F x D x Ka) / (V x (Ka - K)) x (e^(-K x Tmax) - e^(-Ka x Tmax))

    Simplified: Cmax = (F x D / V) x [Ka / (Ka - K)] x [e^(-K*Tmax) - e^(-Ka*Tmax)]

    Where:
        - Cmax: Maximum plasma concentration (e.g., mg/L)
        - F: Bioavailability (dimensionless, 0-1, default 1.0)
        - D: Dose administered (e.g., mg)
        - V: Volume of distribution (e.g., L)
        - Ka: Absorption rate constant (e.g., 1/hour)
        - K: Elimination rate constant (e.g., 1/hour)
        - Tmax: Time to maximum concentration (calculated internally from Ka and K)
        - e: Euler's number (base of natural logarithm)

    Args (all required except bioavailability):
        dose (str): Dose administered (e.g., '500 mg')
        volume (str): Volume of distribution (e.g., L)
        Ka (str): Absorption rate constant (e.g., '1.5 1/hour')
        K (str): Elimination rate constant (e.g., '0.2 1/hour')
        bioavailability (str): Fraction absorbed into systemic circulation
            (e.g., '0.8', dimensionless). Default: '1.0' (100% bioavailability)

    Optional kwargs:
        output_unit (str or False): Desired output unit, or False to keep original
        decimals (int): Number of decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        - name (str): 'Maximum Concentration (Cmax)'
        - magnitude (float): Numeric value of Cmax
        - unit_string (str): Unit as string (e.g., 'milligram / liter')
        - formatted_string (str): Complete formatted result
        - pint_quantity: Full Quantity object for further calculations

    Examples:
        Calculate Cmax for 500 mg dose with V=50 L, Ka=1.5/hr, K=0.2/hr:
            >>> cmax(dose='500 mg', volume='50 L', Ka='1.5 1/hour', K='0.2 1/hour')
            ('Maximum Concentration (Cmax)', 7.25, 'milligram / liter', '7.25 mg/L', ...)

        Calculate Cmax with 80% bioavailability:
            >>> cmax(dose='500 mg', volume='50 L', Ka='1.5 1/hour', K='0.2 1/hour',
            ...      bioavailability='0.8')
            ('Maximum Concentration (Cmax)', 5.8, 'milligram / liter', '5.8 mg/L', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 9: Extravascular Administration
        - Section 9.3: Changing the rate of absorption (pages 4802-4863)
        - Section 9.4: Cmax and Tmax (pages 4864-4868)
        - Figure 9.2: Concentration vs time showing Cmax approximately 0.5 mg/L
        - Figure 9.3: Effect of absorption rate on Cmax (page 4825)

    Typical Cmax Values and Therapeutic Ranges:

        Analgesics:
            - Acetaminophen: Cmax 10-20 mcg/mL (therapeutic), >150 mcg/mL (toxic)
            - Aspirin (as salicylate): Cmax 20-100 mcg/mL (anti-inflammatory)
            - Ibuprofen: Cmax 20-50 mcg/mL (typical therapeutic)

        Cardiovascular:
            - Digoxin: Cmax 1-2 ng/mL (therapeutic), >2 ng/mL (toxic risk)
            - Propranolol: Cmax varies widely (individualize to response)
            - Theophylline: Cmax 10-20 mcg/mL (therapeutic), >20 mcg/mL (toxic)

        Antibiotics (concentration-dependent killing):
            - Aminoglycosides: Cmax/MIC ratio > 8-10 desired for efficacy
              Peak target: 5-10 mcg/mL (gentamicin)
            - Fluoroquinolones: Cmax/MIC ratio important for efficacy

    Formulation Effects on Cmax (from Rowe Chapter 9):

        Rapid-Release vs Slow-Release (Same Dose, Same Drug):
            - Rapid-release: Higher Cmax, earlier Tmax
              * Most of dose absorbed quickly before much elimination
              * High peak, then rapid decline
            - Slow-release: Lower Cmax, later Tmax
              * Absorption spread over time, concurrent elimination
              * Moderate sustained levels, less fluctuation

        Clinical Example from Text (Figure 9.3):
            "If you woke up with a dreadful headache, you would probably welcome
            rapid absorption with a quick hit and high peak level of analgesic."
            vs.
            "If you are taking theophylline for your asthma, the pattern of
            sustained, moderate levels with a slow release formulation is far
            preferable. The rapid absorption profile would threaten early toxicity
            and/or a return of symptoms at later times."

        Practical Implications:
            - Immediate-release: Higher Cmax, more peak-related effects (good or bad)
            - Extended-release: Lower Cmax, same AUC, steadier drug levels
            - Choosing formulation trades off Cmax against Tmax and fluctuation

    Bioavailability Effects on Cmax:

        Cmax is Directly Proportional to F:
            - If F = 1.0: Cmax = full calculated value
            - If F = 0.5: Cmax = 50% of full value (all else equal)
            - If F = 0.25: Cmax = 25% of full value

        First-Pass Effect Examples (from Rowe):
            - Propranolol: F approximately 36% (high first-pass)
              Oral Cmax much lower than equivalent IV dose would predict
            - Digoxin: F approximately 75% (moderate first-pass, dissolution-limited)
            - Verapamil, Lidocaine: Very low oral F, must give IV

        Bioequivalence Implications:
            - Generic drugs must show Cmax within 80-125% of reference
            - Different bioavailability = different Cmax = potential therapeutic failure
            - This is why bioequivalence testing is critical

    Food Effects on Cmax:

        Food Generally DECREASES Cmax:
            - Delayed gastric emptying spreads absorption over longer time
            - Peak is lower but more sustained (often similar AUC)
            - Acetaminophen: Cmax reduced 50% with meal
            - Most drugs show some Cmax reduction with food

        Food May INCREASE Cmax (selected drugs):
            - Highly lipophilic drugs: Fat enhances solubilization
            - Griseofulvin: Cmax increased 2-3x with fatty meal
            - Some drugs: Reduced first-pass when taken with food

        Food May Have MINIMAL Effect:
            - Drugs absorbed throughout GI tract
            - Some sustained-release formulations
            - Theophylline: Generally unaffected by food

        Clinical Management:
            - Narrow therapeutic index drugs: Consistent food conditions important
            - If "take on empty stomach": Usually means Cmax-dependent efficacy
            - If "take with food": Often for GI tolerance, accepting lower Cmax

    Factors Affecting Cmax - Detailed Analysis:

        Parameters That INCREASE Cmax:
            - Higher dose (D): Cmax directly proportional to dose
            - Higher bioavailability (F): Cmax directly proportional to F
            - Faster absorption (higher Ka):
              * Less time for elimination during absorption phase
              * More drug accumulates before peak
            - Slower elimination (lower K):
              * Less drug removed during absorption phase
              * Higher peak concentration achieved
            - Smaller volume of distribution (V):
              * Same amount of drug in smaller apparent volume
              * Higher concentration results

        Parameters That DECREASE Cmax:
            - Lower dose
            - Lower bioavailability (first-pass, incomplete absorption)
            - Slower absorption (sustained-release, food effect for some drugs)
            - Faster elimination (enzyme induction, renal function)
            - Larger volume of distribution (lipophilic drugs, obesity effects)

        Ka/K Ratio Effects on Cmax:
            - Very high Ka/K (>>10): Absorption complete before much elimination
              Cmax approaches F*D/V (like IV bolus)
            - Moderate Ka/K (2-10): Significant elimination during absorption
              Cmax notably less than F*D/V
            - Ka/K approximately 1: Maximum "blunting" of Cmax
            - Ka < K (flip-flop): Unusual kinetics, interpretation changes

    Clinical Applications of Cmax:

        Therapeutic Drug Monitoring:
            - Aminoglycosides: Monitor Cmax to ensure adequate peaks
              Target Cmax: 5-10 mcg/mL (gentamicin)
              Too low: Inadequate antibacterial effect
              Too high: Ototoxicity, nephrotoxicity risk
            - Theophylline: Cmax should stay in 10-20 mcg/mL range
              >20 mcg/mL: Seizure risk, cardiac arrhythmias
            - Digoxin: Sample timing critical (distribution phase)

        Dosing Decisions Based on Cmax:
            - Cmax too low (subtherapeutic):
              * Increase dose
              * Switch to faster-absorbing formulation
              * Improve bioavailability (different salt, take fasted)
            - Cmax too high (toxicity risk):
              * Reduce dose
              * Switch to extended-release formulation
              * Divide into more frequent smaller doses

        Bioequivalence Studies:
            - Cmax must fall within 80-125% of reference product
            - Both rate (Tmax, Cmax) and extent (AUC) are evaluated
            - Failed Cmax comparison = not bioequivalent

        Concentration-Dependent Effects:
            - Antibiotics: Higher Cmax = better bactericidal effect
              Cmax/MIC ratio predicts efficacy for aminoglycosides, fluoroquinolones
            - Analgesics: Higher Cmax = faster, stronger pain relief
            - Sedatives: Higher Cmax = deeper, faster sedation
            - Toxicity: Many adverse effects are concentration-dependent

    Relationship to Other PK Parameters:

        Direct Proportionalities:
            - Cmax proportional to Dose (double dose = double Cmax)
            - Cmax proportional to F (halve F = halve Cmax)

        Inverse Proportionality:
            - Cmax inversely proportional to V (double V = halve Cmax)

        Complex Dependencies:
            - Cmax depends on Ka/K ratio, not individually
            - Higher Ka/K ratio = higher Cmax (all else equal)

        Independence:
            - AUC = F * D / Cl (total exposure, independent of Ka)
            - Cmax and AUC are related but distinct parameters
            - Same AUC can have very different Cmax (formulation effect)

    Special Cases and Limit Behavior:

        When Ka = K:
            - Uses L'Hopital's rule: Cmax = (F * D / V) * e^(-1)
            - This is approximately 0.37 * (F * D / V)
            - Maximum "blunting" of peak compared to IV

        When Ka >> K (very fast absorption):
            - Cmax approaches F * D / V (instantaneous absorption limit)
            - Behaves similar to IV bolus
            - Tmax approaches zero

        When Ka << K (flip-flop kinetics):
            - Terminal phase reflects Ka, not K
            - Cmax still calculable but interpretation changes
            - See flip_flop_check function for detection

    Inter-Individual Variability:

        Sources of Cmax Variability:
            - Absorption variability: Food, GI motility, pH, blood flow
            - Bioavailability variability: First-pass metabolism differences
            - Volume variability: Body composition, protein binding
            - Elimination variability: Genetic polymorphisms, organ function

        High-Variability Drugs:
            - Narrow therapeutic index: Small Cmax differences matter
            - Cyclosporine, tacrolimus: High variability, TDM essential
            - Phenytoin: Saturable metabolism adds complexity

        Clinical Management:
            - Consistent administration conditions (food, timing)
            - Therapeutic drug monitoring for critical medications
            - Individualized dosing based on measured concentrations

    Notes:
        - Assumes one-compartment model with first-order kinetics
        - Does not account for nonlinear (saturable) absorption or elimination
        - Does not include distribution phase (two-compartment effects)
        - Food effects can significantly alter Cmax
        - Inter-individual variability in Cmax can be 2-3 fold or more
        - For multiple dosing, steady-state Cmax will be higher (accumulation)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Get parameters
    dose = Q_(kwargs.get("dose"))
    volume = Q_(kwargs.get("volume"))
    Ka_qty = Q_(kwargs.get("Ka"))
    K_qty = Q_(kwargs.get("K"))
    F = Q_(kwargs.get("bioavailability", "1.0"))  # Default 100% bioavailability

    # Convert K to Ka's units for comparison
    K_converted = K_qty.to(Ka_qty.units)

    Ka_magnitude = Ka_qty.magnitude
    K_magnitude = K_converted.magnitude

    # First calculate Tmax
    if abs(Ka_magnitude - K_magnitude) < _RATE_CONSTANT_TOLERANCE:
        # Special case: Tmax = 1/K when Ka = K
        tmax_qty = 1 / K_qty
        tmax_magnitude = tmax_qty.magnitude

        # When Ka = K, use L'Hopital's rule result:
        # Cmax = (F * D * Ka * Tmax * e^(-K*Tmax)) / V
        # Which simplifies to: Cmax = (F * D / V) * e^(-1)
        exp_term = math.exp(-1)
        cmax_qty = (F * dose / volume) * exp_term
    else:
        # General formula for Tmax
        ratio = Ka_magnitude / K_magnitude
        ln_ratio = math.log(ratio)
        difference_qty = Ka_qty - K_converted
        tmax_qty = ln_ratio / difference_qty

        # Get Tmax magnitude in units compatible with K and Ka
        # Tmax has units of time (e.g., hours)
        tmax_magnitude = tmax_qty.to("hour").magnitude
        K_per_hour = K_qty.to("1/hour").magnitude
        Ka_per_hour = Ka_qty.to("1/hour").magnitude

        # Calculate exponential terms
        exp_K_tmax = math.exp(-K_per_hour * tmax_magnitude)
        exp_Ka_tmax = math.exp(-Ka_per_hour * tmax_magnitude)

        # Cmax = (F * D * Ka) / (V * (Ka - K)) * (e^(-K*Tmax) - e^(-Ka*Tmax))
        # Note: Ka - K has units of 1/time
        numerator = F * dose * Ka_qty
        denominator = volume * difference_qty
        exp_difference = exp_K_tmax - exp_Ka_tmax

        cmax_qty = (numerator / denominator) * exp_difference

    return format_output(
        cmax_qty,
        "Maximum Concentration (Cmax)",
        output_unit,
        decimals
    )


def flip_flop_check(Ka: str, K: str) -> bool:
    """
    Check if flip-flop kinetics applies for given absorption and elimination rates.

    Flip-flop kinetics occurs when the absorption rate constant (Ka) is smaller than
    the elimination rate constant (K). In this situation, absorption becomes the
    rate-limiting step, and the terminal phase of the plasma concentration curve
    reflects absorption rather than elimination - the opposite of normal kinetics.

    The term "flip-flop" describes how the roles of Ka and K are reversed in determining
    the shape of the concentration-time curve. In normal kinetics, the terminal slope
    reflects K (elimination). In flip-flop kinetics, the terminal slope reflects Ka
    (absorption), and the "true" elimination rate K is masked.

    The physiological basis: After extravascular administration, concentration rises
    during absorption and falls during elimination. Normally, absorption is fast and
    finishes early, leaving a terminal phase dominated by elimination. When absorption
    is very slow (sustained-release, depot injection), drug continues to enter the
    body even as elimination is occurring. The terminal phase then reflects the
    rate-limiting absorption process, not elimination.

    Condition: Flip-flop occurs when Ka < K

    This is clinically important because:
        - The measured "half-life" from the terminal phase actually reflects absorption
        - Standard pharmacokinetic interpretation may lead to incorrect conclusions
        - Dosing calculations based on apparent elimination may be wrong
        - True elimination may be much faster than the apparent terminal half-life

    Args:
        Ka (str): Absorption rate constant with units (e.g., '0.1 1/hour')
        K (str): Elimination rate constant with units (e.g., '0.5 1/hour')

    Returns:
        bool: True if flip-flop kinetics applies (Ka < K), False otherwise

    Examples:
        Normal kinetics (absorption faster than elimination):
            >>> flip_flop_check(Ka='1.5 1/hour', K='0.2 1/hour')
            False  # Ka > K, normal kinetics

        Flip-flop kinetics (absorption slower than elimination):
            >>> flip_flop_check(Ka='0.1 1/hour', K='0.5 1/hour')
            True  # Ka < K, flip-flop kinetics

        Edge case (approximately equal rates):
            >>> flip_flop_check(Ka='0.2 1/hour', K='0.2 1/hour')
            False  # Ka = K, not flip-flop (boundary case)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 9: Extravascular Administration
        - Section 9.2-9.3: Parallel absorption and elimination processes
        - Pages 4793-4801: Terminal phase interpretation
        - Figure 9.2: Points A, B, C showing absorption/elimination balance
        - Section on terminal linear portion (pages 5136-5161)

    Understanding Normal vs Flip-Flop Kinetics:

        Normal Kinetics (Ka > K, the usual case):
            - Absorption is FASTER than elimination
            - Drug enters body quickly, then is slowly eliminated
            - Concentration-time curve:
              * Rising phase: Dominated by absorption (fast)
              * Falling phase: Dominated by elimination (slow)
            - Terminal slope = -K (elimination rate constant)
            - Half-life from terminal phase = true elimination half-life
            - Example: Most oral tablets (Ka = 1-2/hr, K = 0.1-0.3/hr)

        Flip-Flop Kinetics (Ka < K, special situations):
            - Absorption is SLOWER than elimination
            - Drug enters body slowly from depot/sustained-release
            - Concentration-time curve:
              * Rising phase: Slow, prolonged absorption
              * Falling phase: Appears slow, but actually limited by ongoing absorption
            - Terminal slope = -Ka (absorption rate constant)
            - Apparent half-life = absorption half-life (NOT elimination!)
            - True elimination is FASTER than apparent terminal decline
            - Example: Sustained-release (Ka = 0.1/hr), depot injection (Ka = 0.02/hr)

    When Flip-Flop Kinetics Occurs - Detailed Scenarios:

        Drug Properties Leading to Slow Absorption:
            - Very poorly water-soluble drugs:
              * Griseofulvin: Dissolution-limited absorption
              * Phenytoin: Erratic, slow absorption
              * Carbamazepine: Variable absorption rate
            - Large molecular weight compounds:
              * Peptides (if oral absorption possible)
              * Some biologics (SC/IM administration)
            - Drugs requiring active transport:
              * Saturable mechanisms limit absorption rate
              * At high doses, absorption becomes rate-limiting

        Formulation-Induced Flip-Flop:
            - Sustained-release/Extended-release tablets:
              * Theophylline SR: Ka deliberately reduced by formulation
              * Morphine ER (MS Contin): Slow release over 8-12 hours
              * Metoprolol ER: Once-daily formulation
              * Nifedipine GITS: Zero-order-like slow release
            - Enteric-coated tablets:
              * Delayed release until intestine
              * Once release begins, may still be slow
            - Matrix tablets:
              * Drug embedded in slowly-eroding matrix
              * Release rate controlled by matrix dissolution

        Route-Induced Flip-Flop:
            - Intramuscular depot injections:
              * Haloperidol decanoate: Monthly IM, Ka very low
              * Fluphenazine decanoate: Depot antipsychotic
              * Medroxyprogesterone (Depo-Provera): 3-month depot
              * Penicillin G benzathine: Slow-release IM
            - Subcutaneous administration:
              * Insulin glargine: Slow, prolonged absorption
              * Some vaccines: Depot effect at injection site
            - Transdermal patches:
              * Fentanyl patch: Ka = 0.02-0.05/hr typically
              * Nicotine patch: Continuous slow input
              * Estradiol patch: Steady hormone delivery

    Clinical Examples of Flip-Flop with Drug Details:

        Fentanyl Transdermal (Classic Example):
            - True elimination t1/2: 3-4 hours (from IV data)
            - Apparent t1/2 from patch removal: 17-24 hours
            - Reason: Skin depot continues releasing drug slowly
            - Ka (patch) approximately 0.03/hr, K approximately 0.2/hr
            - Clinical implication: Remove patch, but effect persists for hours

        Morphine Extended-Release:
            - Immediate-release morphine: Ka = 1-2/hr, terminal t1/2 reflects K
            - Extended-release (MS Contin): Ka = 0.1-0.2/hr
            - If K = 0.35/hr (t1/2 = 2 hrs), then Ka < K: flip-flop
            - Apparent oral t1/2 longer than true elimination t1/2

        Procainamide Sustained-Release:
            - Immediate-release: Normal kinetics
            - Sustained-release: Absorption rate deliberately slowed
            - May show flip-flop depending on formulation design

        Depot Antipsychotics:
            - Haloperidol decanoate: Ka approximately 0.01/hr (weeks to absorb)
            - True K much faster (hours to clear from blood)
            - Terminal "half-life" reflects absorption from depot
            - Effect persists for weeks after injection

    Implications for Pharmacokinetic Analysis:

        Misinterpretation Risk:
            - Analyst measures terminal slope, calculates "half-life"
            - In flip-flop: This is absorption half-life, not elimination!
            - Dosing based on this "half-life" will be incorrect

        Example of Misinterpretation:
            - Observed: Terminal t1/2 = 24 hours after patch removal
            - Incorrect conclusion: "Drug has 24-hour elimination half-life"
            - Correct interpretation: "Drug has 24-hour absorption half-life from depot"
            - True elimination might be 4 hours (must check IV data)

        Steady-State Implications:
            - Normally: Time to steady-state = 4-5 elimination half-lives
            - In flip-flop: Time to steady-state = 4-5 absorption half-lives
            - If Ka-determined t1/2 is longer, steady-state takes longer
            - Accumulation predictions require knowing the TRUE rate constants

    Detection Methods for Flip-Flop:

        Gold Standard - Compare Oral and IV Half-Lives:
            - Give same drug IV (no absorption phase)
            - Measure true elimination half-life from IV data
            - Compare to apparent half-life from oral/extravascular data
            - If oral t1/2 > IV t1/2: Flip-flop is occurring

        Single-Dose Analysis Clues:
            - Unusually prolonged terminal phase for drug class
            - Absorption phase longer than expected
            - Model fitting gives Ka < K

        Population Data Analysis:
            - Simultaneous fitting of IV and oral data
            - Deconvolution analysis
            - Compartmental modeling with absorption compartment

    Dosing Implications in Flip-Flop:

        Dosing Interval Selection:
            - In flip-flop, extended dosing intervals are appropriate
            - Based on absorption rate, not elimination rate
            - Example: Fentanyl patch changed every 72 hours
              (Matches slow Ka, not fast K)

        Loading Dose Calculations:
            - Standard loading dose formula assumes known Vd and target Css
            - In flip-flop, don't need "extra" loading based on apparent t1/2
            - True elimination is fast; drug doesn't accumulate as much

        Accumulation Predictions:
            - Accumulation factor = 1 / (1 - e^(-Ka*tau)) approximately
            - Based on SLOWER rate constant (which is Ka in flip-flop)
            - May accumulate more than expected from "apparent" half-life

        Cessation of Therapy:
            - After stopping sustained-release or depot:
            - Drug continues to appear from depot
            - Effect persists longer than "true" half-life would suggest
            - Example: Fentanyl patch removed, but effect lasts 12-24 hours

    Clinical Management Considerations:

        For Sustained-Release Formulations:
            - Understand that measured "half-life" may be formulation-dependent
            - Don't assume oral t1/2 = true elimination t1/2
            - Use consistent formulation for PK interpretation

        For Depot Injections:
            - Effect duration determined by ABSORPTION from depot
            - Stopping depot: Residual effect until depot depleted
            - Switching drugs: Account for residual depot drug

        For Transdermal Patches:
            - Remove patch: Drug continues from skin depot
            - Overdose concerns: Cutting/heating patch increases absorption
            - Disposal: Used patches still contain drug (safety concern)

        Drug Interaction Considerations:
            - Enzyme inducers/inhibitors affect K, not Ka
            - In flip-flop, Ka dominates terminal phase
            - Drug interactions may have less effect on apparent half-life
            - But true elimination (and thus some effects) still affected

    Mathematical Relationships in Flip-Flop:

        Terminal Slope:
            - Normal: slope = -K, t1/2 = 0.693/K
            - Flip-flop: slope = -Ka, apparent t1/2 = 0.693/Ka

        Tmax Formula Still Valid:
            - Tmax = ln(Ka/K) / (Ka - K)
            - When Ka < K: ln(Ka/K) is negative, (Ka - K) is negative
            - Result: Still positive Tmax (mathematically consistent)

        Cmax Calculation:
            - Same formula applies
            - In flip-flop: Cmax tends to be lower, Tmax longer
            - "Blunted" peak compared to same drug immediate-release

    Notes:
        - This function only checks the simple condition Ka < K
        - Does not account for multi-compartment distribution kinetics
        - Does not consider absorption lag time (delay before absorption starts)
        - Real flip-flop determination often requires concentration-time data
        - Some drugs may show partial flip-flop (Ka approximately K)
        - Flip-flop is a CHARACTERISTIC, not a problem to fix
        - Designed flip-flop (sustained-release) is therapeutically beneficial
    """
    Ka_qty = Q_(Ka)
    K_qty = Q_(K)

    # Convert to compatible units for comparison
    K_converted = K_qty.to(Ka_qty.units)

    Ka_magnitude = Ka_qty.magnitude
    K_magnitude = K_converted.magnitude

    # Flip-flop occurs when absorption is slower than elimination
    return Ka_magnitude < K_magnitude
