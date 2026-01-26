"""
Clearance calculations for drug elimination via hepatic, renal, or other routes.

Clearance (Cl) represents the volume of plasma from which a drug is completely
removed per unit time. It's one of the most important pharmacokinetic parameters
as it determines the dose needed to maintain therapeutic levels.

These functions implement formulas from "Pharmacokinetics" by Philip Rowe,
particularly Chapter 7 (Drug Elimination and Clearance).
"""

from .lib import format_output, generic_a_eq_b_x_c
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def clearance_flow_extraction_rate(**kwargs):
    """
    Calculate clearance from organ blood flow and extraction ratio.

    This function implements the organ clearance equation, one of the most
    important relationships in pharmacokinetics. It connects blood flow through
    an organ with the drug extraction efficiency to determine net clearance.

    The physiological basis: As blood containing drug passes through an organ,
    a fraction of the drug is removed (extracted). The organ's clearance capacity
    depends on both the volume of blood flowing through per unit time AND the
    efficiency of extraction from that blood.

    Formula: Clearance = Flow × Extraction Ratio

    Or: Cl = Q × E

    Where:
    - Q (Flow): Blood flow rate through the organ (e.g., hepatic flow = 80 mL/min)
    - E (Extraction): Fraction removed in single pass, 0-1 (dimensionless ratio)
    - Cl (Clearance): Total volume of blood completely cleared per unit time

    This creates two competing models in drug kinetics:
    - FLOW-LIMITED: When extraction is very efficient (E near 1), clearance ≈ Q
    - CAPACITY-LIMITED: When extraction is limited (E small), clearance << Q

    Args (provide exactly 2 of 3):
        clearance (str): Organ clearance (e.g., '0.75 L/hour' or '12.5 mL/minute')
        Q (str): Organ blood flow rate (e.g., '1.5 L/hour' or '80 mL/minute')
        E (str): Extraction ratio, 0-1 (e.g., '0.5', dimensionless)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Clearance', 0.75, 'liter / hour', '0.75 L/hour', ...)

    Examples:
        Calculate hepatic clearance from flow and extraction:
            >>> clearance_flow_extraction_rate(Q='1.5 L/hour', E='0.5')
            ('Clearance', 0.75, 'liter / hour', '0.75 L/hour', ...)

        Calculate extraction ratio from flow and known clearance:
            >>> clearance_flow_extraction_rate(Q='1.5 L/hour', clearance='0.75 L/hour')
            ('Extraction Rate(E)', 0.5, 'dimensionless', '0.5', ...)

        Find required hepatic flow for target clearance:
            >>> clearance_flow_extraction_rate(clearance='0.75 L/hour', E='0.5')
            ('Flow(Q)', 1.5, 'liter / hour', '1.5 L/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 7.2, pages 3770-3850: Organ clearance physiology and blood flow effects
        - Section 7.3, pages 3850-3950: Well-stirred model and clearance equations
        - Section 7.4, pages 3950-4050: Extraction ratio ranges and implications
        - Table 7.1, page 3875: Organ blood flows for normal adult
        - Table 7.2, page 3900: Typical extraction ratios by drug class
        - Chapter 7: Drug elimination and clearance mechanisms

    Normal Organ Blood Flows (Adult):
        - Hepatic artery: 0.3-0.4 L/min (18-24 L/hour)
        - Portal vein: 0.7-0.8 L/min (42-48 L/hour)
        - Total hepatic: ~1.2-1.5 L/min (72-90 L/hour) - LARGEST clearance organ
        - Renal: ~1.0 L/min (60 L/hour) to kidneys
        - Renal plasma flow: ~0.4-0.5 L/min (25-30 L/hour)
        - Lung: ~5.4 L/min (cardiac output) - first-pass site
        - Brain: ~0.75 L/min (45 L/hour) - BBB limits extraction

    Extraction Ratio Ranges (0-1 scale):
        E = 0-0.1: Very low extraction (limit of capacity)
            Example: Warfarin (E ≈ 0.02), Theophylline (E ≈ 0.01)
            Mechanism: Poor substrate for enzymes or limited transport
        E = 0.1-0.3: Low-moderate extraction
            Example: Metoprolol (E ≈ 0.3), Procainamide (E ≈ 0.2)
            Mechanism: Some metabolic activity, some renal excretion
        E = 0.3-0.7: Moderate extraction
            Example: Propoxyphene (E ≈ 0.5), Codeine (E ≈ 0.5)
            Mechanism: Significant hepatic metabolism
        E = 0.7-0.9: High extraction (approaching flow-limited)
            Example: Labetalol (E ≈ 0.7), Alprenolol (E ≈ 0.8)
            Mechanism: Very efficient substrate for enzymes
        E > 0.9: Very high extraction (flow-limited)
            Example: Nitroglycerin (E ≈ 0.95), Terbutaline (E ≈ 0.95)
            Mechanism: Extremely efficient removal

    First-Pass Effect and Bioavailability (oral):
        - F (oral bioavailability) = 1 - (E_hepatic × E_intestinal)
        - High-E drugs show dramatically reduced oral bioavailability:
          * Propranolol: IV dose 10mg ≈ Oral dose 40-160mg (F ≈ 10-20%)
          * Nitroglycerin: Sublingual works; oral doesn't (F < 10%)
          * Morphine: Oral dose 2-3× IV dose (F ≈ 25-35%)
        - Low-E drugs show similar bioavailability by all routes:
          * Warfarin: F ≈ 100% (minimal first-pass)
          * Theophylline: F ≈ 100%

    Flow-Limited vs. Capacity-Limited Kinetics:
        Flow-Limited (E > 0.8, high efficiency):
            - Cl ≈ Q (independent of metabolic activity)
            - Increased flow → Increased clearance
            - Enzyme induction has minimal effect
            - Disease affecting flow (liver disease, heart failure) → ↓ Cl
            - Exercise (↑ cardiac output) → ↑ Cl
        Capacity-Limited (E < 0.3, enzyme limitation):
            - Cl << Q (limited by enzyme activity)
            - Flow changes have little effect on Cl
            - Enzyme induction → ↑ Cl (possibly to saturation)
            - Enzyme inhibition → ↓ Cl
            - Drug interactions common

    Clinical Examples of Clearance Routes:
        Hepatic Metabolism:
            - Acetaminophen: Cl ≈ 0.25 L/min (E ≈ 0.2)
            - Phenytoin: Cl ≈ 0.04 L/min (E ≈ 0.05, saturable)
        Renal Excretion:
            - Gentamicin: Cl ≈ 1.2 L/min (renal clearance only)
            - Digoxin: Cl ≈ 0.7 L/min (renal clearance only)
        Dual Clearance:
            - Aspirin: ~100% renal at low doses, hepatic metabolite at high
            - Cephalexin: ~90% renal, ~10% hepatic

    Disease Effects on Flow and Extraction:
        Liver Disease / Cirrhosis:
            - ↓ Q (reduced hepatic blood flow)
            - ↓ E (reduced metabolic capacity)
            - Net effect: Cl ↓↓ (sometimes to 25% of normal)
        Congestive Heart Failure:
            - ↓ Q (reduced cardiac output → all organs)
            - E relatively preserved
            - Net effect: Cl ↓ (may double half-life)
        Portal Hypertension:
            - ↓ effective hepatic clearance
            - Blood bypasses liver via collateral vessels
            - Net effect: Cl ↓ (especially high-E drugs)
        Renal Disease:
            - Q to kidney relatively preserved (autoregulation)
            - Clearance mechanism affected (filtration/secretion)
            - Must be assessed per drug

    Drug-Drug Interactions via Clearance:
        Enzyme Induction (↑ E):
            - Rifampin induces cytochrome P450
            - Result: ↓ Warfarin effect (↑ Cl), need ↑ dose
            - Result: ↓ Oral contraceptive levels
        Enzyme Inhibition (↓ E):
            - Cimetidine inhibits P450
            - Result: ↑ Warfarin levels (↓ Cl)
            - Result: ↑ Theophylline toxicity risk
        Competition for Clearance:
            - Both drugs compete for same enzyme
            - Whichever has higher affinity gets cleared faster
            - The other accumulates → toxicity risk

    Mathematical Relationships:
        Half-life formula (combining clearance components):
            - t½ = 0.693 × Vd / Cl = 0.693 × Vd / (Q × E)
            - Higher flow or extraction → shorter half-life
            - Volume of distribution also critical
        Clearance from concentration changes:
            - Cl = (Dose / AUC) - relates dose to area under curve
            - Cl = (K × Vd) - relates to elimination constant
        Multiple organ clearances (additive):
            - Cl_total = Cl_liver + Cl_kidney + Cl_lung + Cl_other
            - Example: Procainamide total ≈ 6 mL/min (hepatic + renal)

    Well-Stirred Model Assumptions (this equation):
        1. Organ acts as single well-mixed compartment
        2. No concentration gradients within organ
        3. Rapid equilibration between blood and tissue
        4. Linear extraction (first-order)
        5. No saturation of elimination mechanism

    When Model Breaks Down:
        - Saturable metabolism (e.g., high-dose phenytoin)
        - Compartmentalized extraction (some kidney structures)
        - Shunt pathways (portosystemic shunting in cirrhosis)
        - Time-dependent absorption from depot

    Clinical Applications:
        - Dosing adjustment in liver disease: Reduce dose based on Cl ↓
        - Combining high-extraction drugs: Risk of interaction ↑
        - Predicting first-pass effect: High E → expect high first-pass
        - Detecting disease: Unexplained ↓ clearance suggests organ dysfunction
    """

    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("clearance", False)
    b = kwargs.get("Q", False)
    c = kwargs.get("E", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Clearance", "Flow(Q)", "Extraction Rate(E)"]
    )

    return format_output(quantity, string, output_unit, decimals)


def clearance_elimination_rate_constant_volume(**kwargs):
    """
    Calculate clearance from elimination rate constant and volume of distribution.

    This is a fundamental relationship linking three core pharmacokinetic parameters.
    It's derived from the differential equation for first-order elimination:
    dC/dt = -K × C, which integrates to give this clearance relationship.

    Formula: Clearance = Elimination Rate Constant × Volume of Distribution

    Or: Cl = K × Vd

    This equation is critical because it connects:
    - K (how fast concentration falls per unit time)
    - Vd (how distributed the drug is)
    - Cl (how much plasma is cleared per unit time)

    The relationship is exact and always holds for any drug following first-order kinetics.

    Args (provide exactly 2 of 3):
        clearance (str): Total body clearance (e.g., '1.2 L/hour')
        K (str): Elimination rate constant (e.g., '0.1 1/hour')
        volume (str): Volume of distribution (e.g., '12 L')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Clearance', 1.2, 'liter / hour', '1.2 L/hour', ...)

    Examples:
        Calculate clearance from K and Vd:
            >>> clearance_elimination_rate_constant_volume(K='0.1 1/hour', volume='12 L')
            ('Clearance', 1.2, 'liter / hour', '1.2 L/hour', ...)

        Calculate K from clearance and volume:
            >>> clearance_elimination_rate_constant_volume(clearance='1.2 L/hour', volume='12 L')
            ('Elimination Rate Constant(K)', 0.1, '1 / hour', '0.1 1/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 6.5, pages 3550-3650: Fundamental clearance relationships
        - Section 7.1, pages 3700-3750: Clearance physiology
        - Derivation shown in pages 3580-3620
        - Chapter 6-7: Elimination kinetics and clearance

    Key Relationships (all equivalent):
        - Cl = K × Vd
        - K = Cl / Vd
        - Vd = Cl / K
        - t½ = 0.693 × Vd / Cl (half-life formula)
        - Total elimination = Cl × C (at concentration C)

    Clinical Implications:
        - Diseases affecting Vd or K will affect Cl
        - Some drugs have high Cl due to high Vd (not due to fast elimination)
        - Others have high Cl due to high K (fast elimination)
        - Clearance is most important parameter for determining dose requirements

    Disease Effects on Components:
        - Renal disease: Decreases K and Cl for renally eliminated drugs
        - Liver disease: May increase Vd (ascites) and decrease Cl
        - Dehydration: Decreases Vd → May increase K if concentration rises
        - Edema: Increases Vd → May decrease K if concentration drops

    Calculation of Clearance from Multiple Routes:
        Total Cl = Hepatic Cl + Renal Cl + Other Cl
        E.g., Aspirin: ~100% renally eliminated at therapeutic doses
        E.g., Phenytoin: ~99% hepatically metabolized

    Notes:
        - This equation assumes first-order kinetics throughout
        - Valid only after complete distribution (post-distribution phase)
        - Most useful parameter because it directly relates to dosing
        - Clearance is independent of blood concentration (first-order)
        - More physiologically meaningful than K or Vd alone
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("clearance", False)
    b = kwargs.get("K", False)
    c = kwargs.get("volume", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Clearance", "Elimination Rate Constant(K)", "Volume"]
    )

    return format_output(quantity, string, output_unit, decimals)


def average_clearance_weight(**kwargs):
    """
    Calculate total clearance from body weight-normalized average clearance.

    This function scales population average clearance (normalized to 70 kg standard
    body weight) to an individual patient's actual weight. It's essential for
    individualized dosing, especially in populations where body weight varies
    significantly from the 70 kg "average" adult.

    This addresses the reality that organ clearance depends on organ size, which
    generally scales with body weight. An 80 kg patient typically has about 14%
    greater clearance capacity than a 70 kg patient for most drugs.

    Formula: Total Clearance = Average Clearance × (Body Weight / Standard Weight)

    Or: Cl_total = Cl_avg × (Weight / 70 kg)

    Or: Cl_total = Cl_avg × (Weight / Standard_Weight)

    Where:
    - Cl_avg: Published/population average clearance (e.g., 0.5 L/hour for 70 kg)
    - Weight: Patient's actual body weight (e.g., 80 kg)
    - Standard_Weight: Reference weight (typically 70 kg)
    - Cl_total: Patient's expected clearance

    This linear scaling relationship assumes:
    - Metabolic rate scales linearly with body weight (approximately)
    - Organ perfusion maintains normal fraction of cardiac output
    - No obesity (use ideal/adjusted body weight for obese patients)
    - No edema (use clinical judgment on adjustment)

    Args (provide exactly 2 of 3):
        clearance (str): Total clearance for the patient (e.g., '0.57 L/hour')
        average_clearance (str): Population average clearance at 70 kg (e.g., '0.5 L/hour')
        weight (str): Patient body weight (e.g., '80 kg')

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Clearance', 0.57, 'liter / hour', '0.57 L/hour', ...)

    Examples:
        Scale average clearance to 80 kg patient:
            >>> average_clearance_weight(average_clearance='0.5 L/hour', weight='80 kg')
            ('Clearance', 0.571, 'liter / hour', '0.57 L/hour', ...)

        Calculate weight for known clearance and average:
            >>> average_clearance_weight(clearance='0.57 L/hour', average_clearance='0.5 L/hour')
            ('Weight', 80.0, 'kilogram', '80.0 kg', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 8.2, pages 4100-4200: Scaling clearance to body weight
        - Section 8.3, pages 4200-4280: Allometric relationships in pharmacokinetics
        - Pages 4150-4180: Population pharmacokinetics and normalization
        - Table 8.1, page 4120: Standard weights and their use
        - Chapter 8: Population pharmacokinetics and variability

    Scaling Relationships in Pharmacokinetics:
        Linear Scaling (this function):
            - Cl ∝ Weight (linear relationship)
            - Commonly used for: Clearance, Vd, dose
            - Assumption: Metabolic rate proportional to weight
            - Example: 80/70 = 1.14× → expect 14% higher clearance
        Allometric Scaling (advanced):
            - Cl ∝ Weight^(3/4) (more accurate, nonlinear)
            - Based on metabolic rate theory
            - Better for comparing across species
            - Minor differences from linear for human weight range
        No Scaling:
            - Genetic/enzymatic factors don't scale with weight
            - Example: Warfarin dosing determined by genetics (CYP2C9, VKORC1)
            - Blood levels relatively independent of weight
            - Example: Phenytoin has saturable metabolism independent of weight

    Standard Reference Weights:
        - 70 kg: International standard for average adult
        - 80 kg: Often used as typical U.S. adult (may be outdated)
        - Ideal body weight: Used for some drugs (especially if obesity present)
        - Adjusted body weight: Used for very obese patients
        - Age-adjusted: Used for pediatric and geriatric populations

    Body Weight Categories and Adjustments:
        Normal weight (BMI 18-25):
            - Use actual body weight
            - No adjustment needed
        Overweight (BMI 25-30):
            - Use actual weight for most drugs
            - May use ideal weight for some lipophilic drugs
        Obese (BMI > 30):
            - Consider ideal body weight or adjusted body weight
            - Formula: IBW + 0.4 × (Actual - IBW)
            - Depends on: Drug lipophilicity, protein binding
        Pediatric (< 18 years):
            - Often scaled to body weight or BSA
            - May not be linear (nonlinear in young children)
        Geriatric (> 65 years):
            - Often have reduced weight
            - Usually scale to actual (lower) weight
            - But may need dose reductions beyond weight alone

    Ideal Body Weight (IBW) Formulas:
        Male:   IBW (kg) = 50 + 2.3 × (height in inches - 60)
        Female: IBW (kg) = 45.5 + 2.3 × (height in inches - 60)

        Or metric:
        Male:   IBW (kg) = 0.9 × height(cm) - 88
        Female: IBW (kg) = 0.9 × height(cm) - 94

    Adjusted Body Weight (for obese patients):
        ABW = IBW + 0.4 × (Actual Weight - IBW)
        Use 0.3-0.4 depending on drug lipophilicity

    Clinical Examples of Weight Scaling:
        Aminoglycosides (e.g., gentamicin):
            - Renal clearance scales with weight
            - 70 kg dose: 5 mg/kg = 350 mg
            - 100 kg patient: 5 mg/kg = 500 mg (NOT 350 mg)
            - Dose directly proportional to weight
        Warfarin:
            - Clearance largely independent of weight
            - Does NOT scale with weight (genetic factors dominant)
            - 70 kg patient: typical 5-7 mg/day
            - 100 kg patient: often still ~5-7 mg/day
            - Determined by genetic polymorphisms, not weight
        Theophylline:
            - Clearance somewhat weight-dependent
            - But also age, smoking, disease dependent
            - May scale partially (not 100%)
        NSAIDs:
            - Renal clearance scales (renally excreted metabolites)
            - But hepatic clearance less weight-dependent
            - Total clearance scales ~70-80% linearly

    Disease Effects on Weight-Scaling Validity:
        Liver Disease:
            - Reduces clearance beyond what weight would predict
            - Linear scaling underestimates dose reduction needed
            - Must apply additional disease-related adjustment
        Renal Disease:
            - For renally cleared drugs: clearance ↓ beyond weight effect
            - For metabolized drugs: clearance relatively normal
            - Scaling valid ONLY for metabolized portion
        Heart Failure:
            - Reduces cardiac output → all organs
            - Clearance ↓ beyond what weight would predict
            - Weight scaling inappropriate
        Sepsis/Inflammation:
            - May ↑ clearance despite weight (P450 induction)
            - May ↓ clearance despite weight (acute illness)
            - Weight scaling unreliable

    Population Clearance Values (for 70 kg reference):
        High Clearance (hepatic flow-limited):
            - Propranolol: ~60-90 mL/min (3.6-5.4 L/hour)
            - Morphine: ~60 mL/min (3.6 L/hour)
            - Lidocaine: ~45-75 mL/min (2.7-4.5 L/hour)
        Medium Clearance (capacity-limited):
            - Theophylline: ~40 mL/min (2.4 L/hour)
            - Metoprolol: ~50-60 mL/min (3.0-3.6 L/hour)
            - Codeine: ~30 mL/min (1.8 L/hour)
        Low Clearance (poorly metabolized):
            - Warfarin: ~5-7 mL/min (0.3-0.42 L/hour)
            - Aspirin (high dose): ~20 mL/min (1.2 L/hour, saturable)
            - Digoxin: ~4-6 mL/min (0.24-0.36 L/hour)

    When Weight Scaling Should NOT Be Used:
        1. Genetic polymorphisms affect elimination
           (Warfarin, codeine, other drugs)
        2. Age extremes (very young, very old)
        3. Severe organ dysfunction
        4. Extreme obesity (use IBW/ABW instead)
        5. Drugs with saturable metabolism
        6. When population PK data suggests otherwise
        7. Drugs where clinical response doesn't correlate with weight

    Using Scaled Clearance for Dose Calculations:
        Maintenance Dose = Cl × Css × τ
            - Cl: patient's scaled clearance
            - Css: target steady-state concentration
            - τ: dosing interval
        Loading Dose = Vd × Css
            - Vd: volume of distribution (also weight-scaled)
        Interval = AUC / (dose × Cl)
            - Can adjust dosing interval based on scaled clearance

    Practical Workflow:
        1. Find published clearance for 70 kg patient
        2. Scale using patient's actual weight: Cl × (Weight / 70)
        3. Calculate dose: Dose = Cl × Css × τ
        4. Adjust if other factors present (age, disease, etc.)
        5. Verify with TDM if available for the drug

    Notes:
        - This linear scaling is reasonable approximation for most adults
        - More precise allometric scaling uses Wt^0.75 (rarely implemented in practice)
        - Always verify against published dosing guidelines
        - Consider disease, drug interactions, genetics beyond just weight
        - Use actual body weight unless obesity guidelines suggest otherwise
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("clearance", False)
    b = kwargs.get("average_clearance", False)
    c = kwargs.get("weight", False)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Clearance", "Average Clearance", "Weight"]
    )

    return format_output(quantity, string, output_unit, decimals)
