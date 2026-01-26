"""
Renal clearance calculations for drug elimination via the kidneys.

This module implements renal clearance formulas from "Pharmacokinetics" by Philip Rowe,
particularly Chapter 14 (Renal Clearance and Drug Dosing in Kidney Disease).

The Cockcroft-Gault equation is the most widely used clinical formula for estimating
creatinine clearance (CrCl) as a measure of renal function. It allows dose adjustment
for renally eliminated drugs based on estimated kidney function.
"""

from .lib import format_output
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def cockcroft_gault_male(**kwargs):
    """
    Calculate creatinine clearance for men using the Cockcroft-Gault equation.

    The Cockcroft-Gault equation estimates creatinine clearance (CrCl) from serum
    creatinine, age, and body weight. It remains the standard for drug dosing
    adjustments in renal impairment despite newer equations like CKD-EPI and MDRD.

    The physiological basis: Serum creatinine concentrations depend upon the rates
    of production and clearance of this waste product. Body weight, gender, and age
    largely determine muscle mass, which determines the rate of creatinine production.
    Creatinine production rate and creatinine clearance jointly determine serum
    creatinine concentration.

    Cockcroft & Gault (1976) obtained weight, gender, age, serum creatinine, and
    measured creatinine clearance (via 24-hour urine collection) for a range of
    individuals. They then produced empirical equations that related creatinine
    clearance to the other factors, with constants determined to match calculated
    values as closely as possible to observed values.

    Formula (Male): CrCl = 1.23 x (140 - Age) x Weight / SrCr

    Where:
    - CrCl: Creatinine clearance (mL/min)
    - Age: Patient age (years)
    - Weight: Body weight (kg)
    - SrCr: Serum creatinine (micromol/L)

    The constant 1.23 for males accounts for higher muscle mass and creatinine
    production compared to females. Male patients' bodies generally contain a
    higher proportion of muscle than those of females.

    Note: When using conventional creatinine units (mg/dL instead of micromol/L),
    the formula becomes: CrCl = (140 - Age) x Weight / (72 x SrCr)
    This function uses SI units (micromol/L) as specified in Rowe's textbook.

    Args (provide exactly 3 of 4):
        creatinine_clearance (str): Estimated CrCl (e.g., '80 mL/min')
        age (str): Patient age in years (e.g., '65 year')
        weight (str): Body weight (e.g., '70 kg')
        serum_creatinine (str): Serum creatinine level (e.g., '100 micromol/L')

    Optional:
        output_unit (str): Desired output unit (e.g., 'L/hour')
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the trio provided.
        Example: ('Creatinine Clearance', 80.0, 'milliliter / minute', '80.0 mL/min', ...)

    Examples:
        Calculate CrCl for a 65-year-old male, 70 kg, serum creatinine 100 micromol/L:
            >>> cockcroft_gault_male(age='65 year', weight='70 kg',
            ...                      serum_creatinine='100 micromol/L')
            ('Creatinine Clearance', 64.58, 'milliliter / minute', '64.58 mL/min', ...)

        Calculate age from known CrCl, weight, and serum creatinine:
            >>> cockcroft_gault_male(creatinine_clearance='64.58 mL/min', weight='70 kg',
            ...                      serum_creatinine='100 micromol/L')
            ('Age', 65.0, 'year', '65.0 year', ...)

        Book Example (Rowe, Chapter 14): 55-year-old male, 75 kg, SrCr 110 micromol/L:
            >>> cockcroft_gault_male(age='55 year', weight='75 kg',
            ...                      serum_creatinine='110 micromol/L')
            ('Creatinine Clearance', 71.3, 'milliliter / minute', '71.3 mL/min', ...)

        Convert output to L/hour for dosing calculations:
            >>> cockcroft_gault_male(age='55 year', weight='75 kg',
            ...                      serum_creatinine='110 micromol/L', output_unit='L/hour')
            ('Creatinine Clearance', 4.28, 'liter / hour', '4.28 L/hour', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 14: Creatinine Clearance, pages 125-129
        - Section 14.1: Clearance of creatinine and various drugs
        - Section 14.1.1: Estimation of creatinine clearance
        - Figure 14.1: Interconnections between body weight, gender, age, CrCl and SrCr
        - Practice questions with worked solutions

        Original reference:
        Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum
        creatinine. Nephron. 1976;16(1):31-41.

    Renal Function Categories (CKD Staging based on GFR/CrCl):
        Stage 1 - Normal or high: >= 90 mL/min
            - No dose adjustment needed for most drugs
            - Monitor renal function periodically
        Stage 2 - Mildly decreased: 60-89 mL/min
            - Often no dose adjustment needed
            - May need monitoring for nephrotoxic drugs
            - Watch for progression
        Stage 3a - Mild to moderately decreased: 45-59 mL/min
            - Consider dose reduction for renally cleared drugs
            - Avoid nephrotoxic drugs if possible
            - Monitor drug levels when available
        Stage 3b - Moderate to severely decreased: 30-44 mL/min
            - Dose reduction often required (25-50%)
            - Increased monitoring recommended
            - Consider alternative drugs with non-renal clearance
        Stage 4 - Severely decreased: 15-29 mL/min
            - Significant dose reduction needed (50-75%)
            - Many drugs require individual TDM
            - Prepare for potential dialysis
        Stage 5 - Kidney failure: < 15 mL/min (or dialysis)
            - Major dose reductions (75%+) or drug avoidance
            - Consider dialysis drug removal
            - Specialist consultation recommended

    Drug Dosing Adjustments by Renal Function Category:
        Aminoglycosides (gentamicin, tobramycin, amikacin):
            - CrCl > 80: Standard dosing (5-7 mg/kg once daily)
            - CrCl 60-80: Extend interval to 36 hours
            - CrCl 40-59: Extend interval to 48 hours
            - CrCl 20-39: Consider traditional dosing (1-2 mg/kg q8-12h) with TDM
            - CrCl < 20: Single dose then redose based on levels
            - Dialysis: Dose post-dialysis with level monitoring

        Digoxin:
            - CrCl > 90: Standard dosing (0.125-0.25 mg daily)
            - CrCl 50-89: 0.125 mg daily or 0.25 mg every other day
            - CrCl 30-49: 0.0625-0.125 mg daily
            - CrCl 10-29: 0.0625 mg daily or 0.125 mg every other day
            - CrCl < 10: 0.0625 mg every other day with TDM
            - Always monitor levels (target 0.8-2.0 ng/mL, some prefer 0.5-1.0 for HF)

        Vancomycin:
            - CrCl > 90: 15-20 mg/kg q12h
            - CrCl 50-89: 15-20 mg/kg q24h
            - CrCl 30-49: 15-20 mg/kg q24-48h
            - CrCl 10-29: 15-20 mg/kg then redose based on levels
            - CrCl < 10: Load 15-20 mg/kg, then 500-1000 mg based on levels

        Metformin:
            - CrCl >= 45: No dose adjustment needed
            - CrCl 30-44: Reduce dose by 50%, monitor closely
            - CrCl < 30: Contraindicated (risk of lactic acidosis)

        NSAIDs:
            - CrCl > 50: Use with caution, avoid if possible
            - CrCl 30-50: Avoid if possible, short course only if needed
            - CrCl < 30: Generally contraindicated

        Fluoroquinolones (ciprofloxacin, levofloxacin):
            - CrCl > 50: Standard dosing
            - CrCl 30-50: Reduce dose by 50%
            - CrCl < 30: Reduce dose by 75% or extend interval

    Ideal Body Weight (IBW) Formulas:
        For Obese Patients (BMI > 30), use IBW instead of actual weight:
        Male:   IBW (kg) = 50 + 2.3 x (height in inches - 60)
        Male:   IBW (kg) = 0.9 x height(cm) - 88

        Adjusted Body Weight (ABW) for dosing some drugs:
        ABW = IBW + 0.4 x (Actual Weight - IBW)

    When to Use Actual vs. Ideal Body Weight:
        Use ACTUAL Body Weight:
            - Normal BMI (18.5-24.9)
            - Slightly overweight (BMI 25-29.9) for most drugs
            - When drug distributes into adipose tissue

        Use IDEAL Body Weight:
            - Obesity (BMI >= 30)
            - Morbid obesity (BMI >= 40)
            - When drug has low lipophilicity
            - Aminoglycosides in obese patients
            - Digoxin clearance calculation in obese patients

        Use ADJUSTED Body Weight:
            - Very obese patients for some antimicrobials
            - When drug partially distributes into adipose

    Population-Specific Adjustments:
        Elderly Patients (> 65 years):
            - Bodies tend to lose muscle and gain fat with age
            - CrCl naturally declines with age
            - May have low serum creatinine despite impaired function
            - Consider dose reduction beyond what CrCl suggests
            - Monitor closely for toxicity

        Cachectic/Low Muscle Mass Patients:
            - Equation may significantly overestimate renal function
            - Very low serum creatinine despite poor function
            - Consider measured CrCl or cystatin C-based equations
            - Use clinical judgment in dosing

        Athletes/High Muscle Mass:
            - Equation may underestimate renal function
            - Higher creatinine production from muscle
            - May tolerate higher doses than equation suggests

        Critically Ill Patients:
            - Renal function may fluctuate rapidly
            - Equation may not reflect current status
            - Consider augmented renal clearance in sepsis
            - Monitor drug levels when available

    Alternative Equations and When to Use Them:
        MDRD (Modification of Diet in Renal Disease):
            - Better for CKD staging
            - Normalized to body surface area (mL/min/1.73m2)
            - Not validated for drug dosing (use Cockcroft-Gault)

        CKD-EPI (Chronic Kidney Disease Epidemiology Collaboration):
            - More accurate at higher GFR levels
            - Normalized to BSA
            - FDA still recommends Cockcroft-Gault for drug dosing

        Measured Creatinine Clearance (24-hour urine):
            - Gold standard but impractical
            - Use when equation accuracy is critical
            - Requires complete urine collection

        Cystatin C-based equations:
            - Better for patients with altered muscle mass
            - Not affected by diet or muscle
            - More expensive, less widely available

    Limitations and When NOT to Use the Equation:
        The Cockcroft-Gault equation should NOT be used in:
            1. Acute Kidney Injury (AKI):
               - Serum creatinine not at steady state
               - Function changing rapidly
               - Use measured clearance or clinical judgment

            2. Extremes of Age:
               - Very young (< 18 years): Use pediatric equations
               - Very old (> 90 years): May be less accurate

            3. Extremes of Body Composition:
               - Severe malnutrition
               - Limb amputations
               - Bodybuilders/extreme athletes
               - Pregnancy (physiological changes)

            4. Rapidly Changing Renal Function:
               - Acute illness
               - Post-surgery
               - Nephrotoxic drug exposure

            5. Dietary Extremes:
               - Vegetarian/vegan diets (lower creatinine production)
               - High protein diets (higher creatinine)
               - Creatine supplements

            6. Drugs Affecting Creatinine Secretion:
               - Trimethoprim (inhibits tubular secretion)
               - Cimetidine (inhibits tubular secretion)
               - May falsely elevate serum creatinine

    Clinical Applications:
        - Drug dosing adjustment for renally cleared drugs
        - Estimating digoxin clearance (see digoxin_clearance function)
        - Aminoglycoside dosing calculations
        - Determining if contrast dye is safe
        - Screening for chronic kidney disease
        - Monitoring progression of renal disease

    Notes:
        - The 1.23 constant was empirically derived, not theoretically calculated
        - Most drug package inserts reference Cockcroft-Gault for dosing
        - Convert units: 96 mL/min = 5.76 L/hour (multiply by 60/1000)
        - Always verify against published dosing guidelines
        - Consider therapeutic drug monitoring when available
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Convert all inputs to pint Quantities
    crcl = kwargs.get("creatinine_clearance", False)
    age = kwargs.get("age", False)
    weight = kwargs.get("weight", False)
    srcr = kwargs.get("serum_creatinine", False)

    if crcl:
        crcl = Q_(crcl)
    if age:
        age = Q_(age)
    if weight:
        weight = Q_(weight)
    if srcr:
        srcr = Q_(srcr)

    # Count provided parameters
    provided = sum([bool(crcl), bool(age), bool(weight), bool(srcr)])
    if provided != 3:
        raise ValueError(
            f"cockcroft_gault_male requires exactly 3 of 4 parameters. "
            f"Got {provided}: creatinine_clearance={crcl is not False}, "
            f"age={age is not False}, weight={weight is not False}, "
            f"serum_creatinine={srcr is not False}"
        )

    # Constant for male equation
    MALE_CONSTANT = Q_(1.23, "mL/min * micromol/L / kg / year")

    # Formula: CrCl = 1.23 x (140 - Age) x Weight / SrCr
    # Rearranged forms for solving each variable

    if not crcl:
        # Solve for CrCl: CrCl = 1.23 x (140 - Age) x Weight / SrCr
        age_factor = Q_(140, "year") - age.to("year")
        quantity = MALE_CONSTANT * age_factor * weight.to("kg") / srcr.to("micromol/L")
        string = "Creatinine Clearance"

    elif not age:
        # Solve for Age: (140 - Age) = CrCl x SrCr / (1.23 x Weight)
        # Age = 140 - (CrCl x SrCr / (1.23 x Weight))
        age_factor = crcl.to("mL/min") * srcr.to("micromol/L") / (MALE_CONSTANT * weight.to("kg"))
        quantity = Q_(140, "year") - age_factor
        string = "Age"

    elif not weight:
        # Solve for Weight: Weight = CrCl x SrCr / (1.23 x (140 - Age))
        age_factor = Q_(140, "year") - age.to("year")
        quantity = crcl.to("mL/min") * srcr.to("micromol/L") / (MALE_CONSTANT * age_factor)
        string = "Weight"

    else:  # not srcr
        # Solve for SrCr: SrCr = 1.23 x (140 - Age) x Weight / CrCl
        age_factor = Q_(140, "year") - age.to("year")
        quantity = MALE_CONSTANT * age_factor * weight.to("kg") / crcl.to("mL/min")
        string = "Serum Creatinine"

    return format_output(quantity, string, output_unit, decimals)


def cockcroft_gault_female(**kwargs):
    """
    Calculate creatinine clearance for women using the Cockcroft-Gault equation.

    The female version of the Cockcroft-Gault equation uses a lower constant (1.04)
    compared to males (1.23) to account for lower muscle mass and creatinine
    production in women. Female patients' bodies generally contain a lower proportion
    of muscle than those of males of similar body weight.

    The physiological basis: Serum creatinine concentrations depend upon the rates
    of production and clearance of this waste product. Body weight, gender, and age
    largely determine muscle mass, which determines the rate of creatinine production.
    Creatinine production rate and creatinine clearance jointly determine serum
    creatinine concentration.

    Formula (Female): CrCl = 1.04 x (140 - Age) x Weight / SrCr

    Where:
    - CrCl: Creatinine clearance (mL/min)
    - Age: Patient age (years)
    - Weight: Body weight (kg)
    - SrCr: Serum creatinine (micromol/L)

    The ratio of female to male constants (1.04/1.23 = 0.85) reflects the
    approximately 15% lower muscle mass in women compared to men of similar
    body weight.

    Note: When using conventional creatinine units (mg/dL instead of micromol/L),
    the formula becomes: CrCl = 0.85 x (140 - Age) x Weight / (72 x SrCr)
    This function uses SI units (micromol/L) as specified in Rowe's textbook.

    Args (provide exactly 3 of 4):
        creatinine_clearance (str): Estimated CrCl (e.g., '80 mL/min')
        age (str): Patient age in years (e.g., '65 year')
        weight (str): Body weight (e.g., '60 kg')
        serum_creatinine (str): Serum creatinine level (e.g., '100 micromol/L')

    Optional:
        output_unit (str): Desired output unit (e.g., 'L/hour')
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the trio provided.
        Example: ('Creatinine Clearance', 46.8, 'milliliter / minute', '46.8 mL/min', ...)

    Examples:
        Calculate CrCl for a 65-year-old female, 60 kg, serum creatinine 100 micromol/L:
            >>> cockcroft_gault_female(age='65 year', weight='60 kg',
            ...                        serum_creatinine='100 micromol/L')
            ('Creatinine Clearance', 46.8, 'milliliter / minute', '46.8 mL/min', ...)

        Calculate serum creatinine from known CrCl, age, and weight:
            >>> cockcroft_gault_female(creatinine_clearance='46.8 mL/min', age='65 year',
            ...                        weight='60 kg')
            ('Serum Creatinine', 100.0, 'micromol / L', '100.0 micromol/L', ...)

        Book Practice Question 1 (Rowe, Chapter 14):
        Female patient aged 62, weighs 58 kg, serum creatinine 49 micromol/L:
            >>> cockcroft_gault_female(age='62 year', weight='58 kg',
            ...                        serum_creatinine='49 micromol/L')
            ('Creatinine Clearance', 96.0, 'milliliter / minute', '96.0 mL/min', ...)

        Convert to L/hour (as required in practice question):
            >>> cockcroft_gault_female(age='62 year', weight='58 kg',
            ...                        serum_creatinine='49 micromol/L', output_unit='L/hour')
            ('Creatinine Clearance', 5.76, 'liter / hour', '5.76 L/hour', ...)

        Book Practice Question 2 (Rowe, Chapter 14):
        Female patient aged 44, weighs 54 kg, serum creatinine 127 micromol/L:
            >>> cockcroft_gault_female(age='44 year', weight='54 kg',
            ...                        serum_creatinine='127 micromol/L')
            ('Creatinine Clearance', 42.5, 'milliliter / minute', '42.5 mL/min', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 14: Creatinine Clearance, pages 125-129
        - Section 14.1: Clearance of creatinine and various drugs
        - Section 14.1.1: Estimation of creatinine clearance
        - Section 14.3: Practice questions with worked solutions (pages 129, 146)
        - Figure 14.1: Interconnections between body weight, gender, age, CrCl and SrCr

        Original reference:
        Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum
        creatinine. Nephron. 1976;16(1):31-41.

    Renal Function Categories (CKD Staging based on GFR/CrCl):
        Stage 1 - Normal or high: >= 90 mL/min
            - No dose adjustment needed for most drugs
            - Monitor renal function periodically
        Stage 2 - Mildly decreased: 60-89 mL/min
            - Often no dose adjustment needed
            - May need monitoring for nephrotoxic drugs
            - Watch for progression
        Stage 3a - Mild to moderately decreased: 45-59 mL/min
            - Consider dose reduction for renally cleared drugs
            - Avoid nephrotoxic drugs if possible
            - Monitor drug levels when available
        Stage 3b - Moderate to severely decreased: 30-44 mL/min
            - Dose reduction often required (25-50%)
            - Increased monitoring recommended
            - Consider alternative drugs with non-renal clearance
        Stage 4 - Severely decreased: 15-29 mL/min
            - Significant dose reduction needed (50-75%)
            - Many drugs require individual TDM
            - Prepare for potential dialysis
        Stage 5 - Kidney failure: < 15 mL/min (or dialysis)
            - Major dose reductions (75%+) or drug avoidance
            - Consider dialysis drug removal
            - Specialist consultation recommended

    Drug Dosing Adjustments by Renal Function Category:
        Aminoglycosides (gentamicin, tobramycin, amikacin):
            - CrCl > 80: Standard dosing (5-7 mg/kg once daily)
            - CrCl 60-80: Extend interval to 36 hours
            - CrCl 40-59: Extend interval to 48 hours
            - CrCl 20-39: Consider traditional dosing (1-2 mg/kg q8-12h) with TDM
            - CrCl < 20: Single dose then redose based on levels
            - Dialysis: Dose post-dialysis with level monitoring

        Digoxin:
            - CrCl > 90: Standard dosing (0.125-0.25 mg daily)
            - CrCl 50-89: 0.125 mg daily or 0.25 mg every other day
            - CrCl 30-49: 0.0625-0.125 mg daily
            - CrCl 10-29: 0.0625 mg daily or 0.125 mg every other day
            - CrCl < 10: 0.0625 mg every other day with TDM
            - Always monitor levels (target 0.8-2.0 ng/mL, some prefer 0.5-1.0 for HF)

        ACE Inhibitors (lisinopril, enalapril):
            - CrCl > 30: Start at normal dose, titrate carefully
            - CrCl 10-30: Start at 50% of normal dose
            - CrCl < 10: Start at 25% of normal dose
            - Monitor potassium and creatinine closely

        Gabapentin/Pregabalin:
            - CrCl > 60: Normal dosing
            - CrCl 30-60: Reduce dose by 50%
            - CrCl 15-29: Reduce dose by 75%
            - CrCl < 15: Further reduction with TDM

    Ideal Body Weight (IBW) Formulas:
        For Obese Patients (BMI > 30), use IBW instead of actual weight:
        Female: IBW (kg) = 45.5 + 2.3 x (height in inches - 60)
        Female: IBW (kg) = 0.9 x height(cm) - 94

        Adjusted Body Weight (ABW) for dosing some drugs:
        ABW = IBW + 0.4 x (Actual Weight - IBW)

    Special Considerations for Women:
        Pregnancy:
            - GFR increases by 40-50% during pregnancy
            - Peak increase in second trimester
            - Returns to baseline postpartum
            - May need HIGHER doses of renally cleared drugs
            - Avoid nephrotoxic drugs when possible

        Menopause and Aging:
            - Gradual decline in muscle mass
            - May have lower serum creatinine despite reduced GFR
            - Equation may overestimate function significantly
            - Consider measured CrCl in post-menopausal women

        Oral Contraceptives:
            - Generally no significant effect on renal function
            - Some may cause mild fluid retention
            - No dose adjustment of equation needed

        Hormone Replacement Therapy:
            - Generally no significant effect on renal function
            - Some protective effect on kidney function suggested
            - No dose adjustment of equation needed

        Polycystic Ovary Syndrome (PCOS):
            - May be associated with metabolic syndrome
            - Monitor for diabetes and hypertension effects
            - Standard equation application

    Population-Specific Adjustments:
        Elderly Women (> 65 years):
            - Very common to have low muscle mass
            - Equation often significantly overestimates function
            - May have serum creatinine < 50 micromol/L
            - Consider cystatin C or measured CrCl
            - More conservative dosing often appropriate

        Cachectic/Malnourished Women:
            - Extremely low muscle mass
            - Very low serum creatinine (< 40 micromol/L)
            - Equation grossly overestimates function
            - Use clinical judgment and measured CrCl

        Obese Women (BMI >= 30):
            - Use ideal body weight, not actual weight
            - Excess adipose does not produce creatinine
            - Using actual weight overestimates function

        Female Athletes:
            - Higher muscle mass than typical female
            - May have higher serum creatinine
            - Equation may underestimate function
            - Consider this when dosing

    Alternative Equations for Women:
        The Cockcroft-Gault equation was derived primarily from male subjects.
        Alternative equations may be more accurate in some female populations:

        CKD-EPI (preferred for staging):
            - Race and sex coefficients included
            - More accurate at higher GFR
            - Normalized to BSA (mL/min/1.73m2)

        MDRD (older, less accurate at high GFR):
            - Includes race and sex coefficients
            - Normalized to BSA

        For Drug Dosing:
            - FDA still recommends Cockcroft-Gault
            - Package inserts typically reference CrCl (not GFR)
            - Be aware of limitations in your patient

    Limitations and When NOT to Use the Equation:
        1. Acute Kidney Injury:
           - Serum creatinine not at steady state
           - Wait 48-72 hours or use measured CrCl

        2. Pregnancy:
           - Physiological increase in GFR
           - Lower serum creatinine is normal
           - Consult pregnancy-specific guidelines

        3. Severe Muscle Wasting:
           - Cancer cachexia
           - Chronic illness wasting
           - Anorexia nervosa
           - Use measured CrCl or cystatin C

        4. Rapidly Changing Renal Function:
           - Acute illness, sepsis
           - Post-surgery
           - Nephrotoxic drug exposure

        5. Extremes of Body Weight:
           - Very underweight (BMI < 18.5)
           - Morbid obesity (BMI > 40)

    Clinical Applications:
        - Drug dosing adjustment for renally cleared drugs
        - Estimating digoxin clearance for maintenance dosing
        - Aminoglycoside dosing calculations
        - Contrast dye safety assessment
        - Chronic kidney disease screening
        - Monitoring renal disease progression

    Notes:
        - The 1.04 constant was empirically derived from the original study
        - Female constant is 85% of male constant (1.04/1.23 = 0.85)
        - Most drug package inserts reference Cockcroft-Gault
        - Convert units: 96 mL/min = 5.76 L/hour (multiply by 60/1000)
        - Always verify against published dosing guidelines
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Convert all inputs to pint Quantities
    crcl = kwargs.get("creatinine_clearance", False)
    age = kwargs.get("age", False)
    weight = kwargs.get("weight", False)
    srcr = kwargs.get("serum_creatinine", False)

    if crcl:
        crcl = Q_(crcl)
    if age:
        age = Q_(age)
    if weight:
        weight = Q_(weight)
    if srcr:
        srcr = Q_(srcr)

    # Count provided parameters
    provided = sum([bool(crcl), bool(age), bool(weight), bool(srcr)])
    if provided != 3:
        raise ValueError(
            f"cockcroft_gault_female requires exactly 3 of 4 parameters. "
            f"Got {provided}: creatinine_clearance={crcl is not False}, "
            f"age={age is not False}, weight={weight is not False}, "
            f"serum_creatinine={srcr is not False}"
        )

    # Constant for female equation
    FEMALE_CONSTANT = Q_(1.04, "mL/min * micromol/L / kg / year")

    # Formula: CrCl = 1.04 x (140 - Age) x Weight / SrCr
    # Rearranged forms for solving each variable

    if not crcl:
        # Solve for CrCl: CrCl = 1.04 x (140 - Age) x Weight / SrCr
        age_factor = Q_(140, "year") - age.to("year")
        quantity = FEMALE_CONSTANT * age_factor * weight.to("kg") / srcr.to("micromol/L")
        string = "Creatinine Clearance"

    elif not age:
        # Solve for Age: (140 - Age) = CrCl x SrCr / (1.04 x Weight)
        # Age = 140 - (CrCl x SrCr / (1.04 x Weight))
        age_factor = crcl.to("mL/min") * srcr.to("micromol/L") / (FEMALE_CONSTANT * weight.to("kg"))
        quantity = Q_(140, "year") - age_factor
        string = "Age"

    elif not weight:
        # Solve for Weight: Weight = CrCl x SrCr / (1.04 x (140 - Age))
        age_factor = Q_(140, "year") - age.to("year")
        quantity = crcl.to("mL/min") * srcr.to("micromol/L") / (FEMALE_CONSTANT * age_factor)
        string = "Weight"

    else:  # not srcr
        # Solve for SrCr: SrCr = 1.04 x (140 - Age) x Weight / CrCl
        age_factor = Q_(140, "year") - age.to("year")
        quantity = FEMALE_CONSTANT * age_factor * weight.to("kg") / crcl.to("mL/min")
        string = "Serum Creatinine"

    return format_output(quantity, string, output_unit, decimals)


def cockcroft_gault(gender: str, **kwargs):
    """
    Calculate creatinine clearance using the Cockcroft-Gault equation.

    This is a convenience function that dispatches to either cockcroft_gault_male()
    or cockcroft_gault_female() based on the gender parameter. It implements the
    equations developed by Cockcroft and Gault in 1976, which remain the standard
    for drug dosing adjustments despite the availability of newer equations.

    The Cockcroft-Gault equation allows estimation of creatinine clearance from
    a single blood sample for serum creatinine, combined with patient demographics.
    This is much quicker and easier than the 24-hour urine collections that were
    formerly necessary for measuring actual creatinine clearance.

    The physiological basis: Body weight, gender, and age largely determine muscle
    mass. Muscle mass determines the rate of creatinine production (creatinine is
    a waste product produced in muscles). Creatinine production rate and creatinine
    clearance jointly determine serum creatinine concentration.

    Formula (Male): CrCl = 1.23 x (140 - Age) x Weight / SrCr
    Formula (Female): CrCl = 1.04 x (140 - Age) x Weight / SrCr

    Where:
    - CrCl: Creatinine clearance (mL/min)
    - Age: Patient age (years)
    - Weight: Body weight (kg)
    - SrCr: Serum creatinine (micromol/L)

    The numerical constants (1.23 for males, 1.04 for females) were empirically
    determined by Cockcroft & Gault to match calculated values as closely as
    possible to observed values in their study population.

    Args:
        gender (str): Patient gender - 'male', 'm', 'female', or 'f' (case-insensitive)

    Kwargs (provide exactly 3 of 4):
        creatinine_clearance (str): Estimated CrCl (e.g., '80 mL/min')
        age (str): Patient age in years (e.g., '65 year')
        weight (str): Body weight (e.g., '70 kg')
        serum_creatinine (str): Serum creatinine level (e.g., '100 micromol/L')

    Optional:
        output_unit (str): Desired output unit (e.g., 'L/hour')
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the trio provided.

    Examples:
        Calculate CrCl for a male patient:
            >>> cockcroft_gault('male', age='65 year', weight='70 kg',
            ...                 serum_creatinine='100 micromol/L')
            ('Creatinine Clearance', 64.58, 'milliliter / minute', '64.58 mL/min', ...)

        Calculate CrCl for a female patient:
            >>> cockcroft_gault('F', age='65 year', weight='60 kg',
            ...                 serum_creatinine='100 micromol/L')
            ('Creatinine Clearance', 46.8, 'milliliter / minute', '46.8 mL/min', ...)

        Book Example (Rowe, Chapter 14): Male, 55 years, 75 kg, SrCr 110 micromol/L:
            >>> cockcroft_gault('male', age='55 year', weight='75 kg',
            ...                 serum_creatinine='110 micromol/L')
            ('Creatinine Clearance', 71.3, 'milliliter / minute', '71.3 mL/min', ...)

        Book Practice Question 1: Female, 62 years, 58 kg, SrCr 49 micromol/L:
            >>> cockcroft_gault('female', age='62 year', weight='58 kg',
            ...                 serum_creatinine='49 micromol/L', output_unit='L/hour')
            ('Creatinine Clearance', 5.76, 'liter / hour', '5.76 L/hour', ...)

        Book Practice Question 2: Female, 44 years, 54 kg, SrCr 127 micromol/L:
            >>> cockcroft_gault('female', age='44 year', weight='54 kg',
            ...                 serum_creatinine='127 micromol/L')
            ('Creatinine Clearance', 42.5, 'milliliter / minute', '42.5 mL/min', ...)

    Raises:
        ValueError: If gender is not recognized as male or female

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 14: Creatinine Clearance, pages 125-129
        - Section 14.1: Clearance of creatinine and various drugs
        - Section 14.1.1: Estimation of creatinine clearance
        - Figure 14.1: Interconnections between body weight, gender, age, CrCl and SrCr
        - Section 14.3: Practice questions with worked solutions

        Original reference:
        Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum
        creatinine. Nephron. 1976;16(1):31-41.

    Why Creatinine Clearance Matters for Drug Dosing:
        From Rowe's textbook: "We are mainly concerned with those drugs that are
        simply filtered in the renal glomerulus after which there is no significant
        re-absorption from, or active secretion into, the urine. Because they are
        all handled in a similar manner, these drugs have virtually identical renal
        clearances."

        Typical examples include:
        - Aminoglycoside antibiotics (gentamicin, tobramycin, amikacin)
        - Digoxin (with additional hepatic elimination)
        - Vancomycin
        - Many beta-lactam antibiotics

        Because creatinine is excreted by simple filtration in the kidneys (like
        these drugs), creatinine clearance provides a good estimate for the renal
        clearance of such drugs.

    Gender Considerations:
        Biological Sex Differences:
            - Males have higher muscle mass proportion (constant 1.23)
            - Females have lower muscle mass proportion (constant 1.04)
            - Ratio: 1.04/1.23 = 0.85 (females produce ~15% less creatinine)

        Transgender Patients:
            - Hormone therapy affects muscle mass over time
            - Trans women on estrogen: Consider female equation after 1-2 years
            - Trans men on testosterone: Consider male equation after 1-2 years
            - Consult endocrinology or clinical pharmacy guidelines
            - Individual assessment recommended

        Non-binary Patients:
            - Some clinicians use an intermediate value
            - Average of male and female constants: (1.23 + 1.04) / 2 = 1.135
            - Clinical judgment based on body composition
            - Consider measured CrCl if accuracy critical

    Renal Function Categories (CKD Staging based on GFR/CrCl):
        Stage 1 - Normal: >= 90 mL/min (no dose adjustment)
        Stage 2 - Mild: 60-89 mL/min (monitor, usually no adjustment)
        Stage 3a - Mild-moderate: 45-59 mL/min (consider dose reduction)
        Stage 3b - Moderate-severe: 30-44 mL/min (dose reduction often needed)
        Stage 4 - Severe: 15-29 mL/min (significant dose reduction)
        Stage 5 - Kidney failure: < 15 mL/min (major reduction, dialysis)

    Drugs Requiring Renal Dose Adjustment:
        The following drugs require dose adjustment based on CrCl:

        Aminoglycosides (gentamicin, tobramycin):
            - Almost exclusively eliminated by renal excretion
            - CrCl directly predicts drug clearance
            - Adjust dose AND/OR interval based on CrCl

        Digoxin:
            - Renal excretion is primary route (60-80%)
            - Also has hepatic component (0.33 mL/min/kg)
            - Use digoxin_clearance() function for total clearance

        Vancomycin:
            - Primarily renally eliminated
            - Requires TDM in renal impairment
            - Adjust interval rather than dose typically

        Fluoroquinolones:
            - Variable renal elimination (30-90%)
            - Check specific drug for adjustment

        Beta-lactams:
            - Most are renally eliminated
            - Reduce dose in renal impairment

    Body Weight Considerations:
        Use ACTUAL Body Weight when:
            - Normal BMI (18.5-24.9)
            - Slightly overweight (BMI 25-29.9)
            - From Rowe: "Build is normal - no need to use ideal body weight"

        Use IDEAL Body Weight when:
            - Obese (BMI >= 30)
            - From Rowe: "For a very obese patient, ideal body weight should
              be substituted for actual body weight"

        IBW Formulas:
            Male:   IBW = 50 + 2.3 x (height in inches - 60)
            Female: IBW = 45.5 + 2.3 x (height in inches - 60)

    Clinical Workflow for Drug Dosing:
        1. Determine patient parameters (age, weight, gender, serum creatinine)
        2. Check if weight adjustment needed (use IBW if obese)
        3. Calculate CrCl using Cockcroft-Gault
        4. Look up drug-specific dosing recommendations
        5. Calculate appropriate dose/interval
        6. Monitor drug levels if TDM available
        7. Reassess if renal function changes

    Limitations and Contraindications:
        DO NOT use this equation when:
        - Acute kidney injury (creatinine not at steady state)
        - Rapidly changing renal function
        - Extremes of body composition (cachexia, bodybuilders)
        - Pregnancy (physiological changes in GFR)
        - Pediatric patients (use age-appropriate equations)
        - Drugs affecting creatinine secretion (trimethoprim, cimetidine)

    Alternative Equations:
        CKD-EPI: Better for staging CKD (normalized to BSA)
        MDRD: Older equation for CKD staging
        Schwartz: For pediatric patients
        Measured CrCl: Gold standard (24-hour urine)
        Cystatin C: For patients with altered muscle mass

        Note: FDA drug labeling typically references Cockcroft-Gault,
        not CKD-EPI or MDRD, for dosing adjustments.

    Notes:
        - This convenience function dispatches to gender-specific functions
        - All detailed documentation available in those functions
        - Convert units: mL/min x 60 / 1000 = L/hour
        - Always verify dosing against current drug references
    """
    gender_lower = gender.lower().strip()

    if gender_lower in ("male", "m"):
        return cockcroft_gault_male(**kwargs)
    elif gender_lower in ("female", "f"):
        return cockcroft_gault_female(**kwargs)
    else:
        raise ValueError(
            f"Unrecognized gender '{gender}'. "
            f"Use 'male', 'm', 'female', or 'f' (case-insensitive)."
        )


def digoxin_clearance(**kwargs):
    """
    Calculate total digoxin clearance from creatinine clearance and body weight.

    Digoxin is cleared primarily by the kidneys in patients free of hepatic or
    renal disease, but there is additional hepatic metabolism that needs to be
    taken into account. This function implements the standard equation for
    estimating total digoxin clearance as described in Rowe's Pharmacokinetics.

    The physiological basis: Creatinine and digoxin are handled similarly within
    the kidneys - both are simply filtered in the renal glomerulus with no
    significant re-absorption or active secretion. Because of this similar
    handling, we assume that the renal clearance of digoxin equals that of
    creatinine. The hepatic (non-renal) component is estimated on a simple
    body weight basis.

    Formula: Digoxin Cl = CrCl + 0.33 mL/min/kg x Body Weight

    Where:
    - Digoxin Cl: Total digoxin clearance (mL/min)
    - CrCl: Creatinine clearance (mL/min) - renal component
    - 0.33 mL/min/kg: Non-renal clearance constant (hepatic metabolism)
    - Body Weight: Patient weight (kg)

    The equation reflects that digoxin clearance has two components:
    1. Renal clearance: Approximately equal to CrCl (filtration)
    2. Non-renal clearance: About 0.33 mL/min per kg body weight (hepatic)

    In patients with normal renal function, approximately 60-80% of digoxin
    clearance is renal. In renal failure, non-renal clearance becomes dominant.

    Args (provide exactly 2 of 3):
        digoxin_clearance (str): Total digoxin clearance (e.g., '100 mL/min' or '6 L/hour')
        creatinine_clearance (str): Patient's CrCl (e.g., '80 mL/min')
        weight (str): Body weight (e.g., '70 kg')

    Optional:
        output_unit (str): Desired output unit (e.g., 'L/hour')
        decimals (int): Decimal places for rounding (default: 2)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing from the pair provided.
        Example: ('Digoxin Clearance', 103.1, 'milliliter / minute', '103.1 mL/min', ...)

    Examples:
        Calculate digoxin clearance for patient with CrCl 80 mL/min, weight 70 kg:
            >>> digoxin_clearance(creatinine_clearance='80 mL/min', weight='70 kg')
            ('Digoxin Clearance', 103.1, 'milliliter / minute', '103.1 mL/min', ...)

        Calculate in L/hour:
            >>> digoxin_clearance(creatinine_clearance='80 mL/min', weight='70 kg',
            ...                   output_unit='L/hour')
            ('Digoxin Clearance', 6.19, 'liter / hour', '6.19 L/hour', ...)

        Calculate required CrCl for a specific digoxin clearance:
            >>> digoxin_clearance(digoxin_clearance='103.1 mL/min', weight='70 kg')
            ('Creatinine Clearance', 80.0, 'milliliter / minute', '80.0 mL/min', ...)

        Book Example (Rowe, Chapter 14 - Section 14.2):
        Male patient, 55 years, 75 kg, SrCr 110 micromol/L:
        First calculate CrCl = 71.3 mL/min, then:
            >>> digoxin_clearance(creatinine_clearance='71.3 mL/min', weight='75 kg')
            ('Digoxin Clearance', 96.05, 'milliliter / minute', '96.05 mL/min', ...)
            # Which equals 5.77 L/hour when converted

        Book Practice Question 2 (Rowe, Chapter 14):
        Female patient, 44 years, 54 kg, SrCr 127 micromol/L:
        First calculate CrCl = 42.5 mL/min, then:
            >>> digoxin_clearance(creatinine_clearance='42.5 mL/min', weight='54 kg')
            ('Digoxin Clearance', 60.32, 'milliliter / minute', '60.32 mL/min', ...)
            # Which equals 3.62 L/hour when converted

    Reference:
        Rowe, P. Pharmacokinetics
        - Chapter 14: Creatinine Clearance, pages 125-129
        - Section 14.2: Digoxin Dosing, pages 127-129
        - Worked examples with complete calculations
        - Practice questions with solutions (pages 129, 146)

    Complete Digoxin Dosing Workflow (from the book):
        The book provides a complete workflow for calculating digoxin dose:

        Step 1: Determine patient parameters
            - Gender, Age, Weight, Serum creatinine
            - Example: Male, 55 years, 75 kg, SrCr 110 micromol/L

        Step 2: Calculate CrCl using Cockcroft-Gault
            - For men: CrCl = 1.23 x (140 - Age) x Wt / SrCr
            - Example: CrCl = 1.23 x (140-55) x 75 / 110 = 71.3 mL/min
            - Assume digoxin renal clearance = CrCl

        Step 3: Calculate hepatic (non-renal) clearance
            - Hepatic clearance = 0.33 mL/min/kg x Weight
            - Example: 0.33 x 75 = 24.8 mL/min

        Step 4: Calculate total body clearance
            - Total Cl = Renal Cl + Hepatic Cl
            - Example: 71.3 + 24.8 = 96.1 mL/min

        Step 5: Convert to L/hour for dosing calculations
            - Cl (L/h) = Cl (mL/min) x 60 / 1000
            - Example: 96.1 x 60 / 1000 = 5.77 L/h

        Step 6: Calculate dose using steady-state equation
            - Css,av = F x D / (Cl x tau)
            - Therefore: D = Css,av x Cl x tau / F
            - Target Css,av: 0.8-2.0 microgram/L (use 1.4 as midpoint)
            - Oral bioavailability (F): 70% for tablets
            - Dosing interval (tau): 24 hours for once daily
            - Example: D = 1.4 x 5.77 x 24 / 0.7 = 277 microgram

        Step 7: Round to available tablet strength
            - Digoxin tablets: multiples of 62.5 microgram
            - Available: 62.5, 125, 187.5, 250, 312.5 microgram
            - Example: 277 -> nearest is 4 x 62.5 = 250 microgram

    Digoxin Pharmacokinetics:
        Absorption:
            - Oral bioavailability: 60-80% (tablets), 70-85% (elixir), 90-100% (capsules)
            - Food does not significantly affect absorption
            - Relatively polar drug - absorbed rather inefficiently
            - Significant part of dose may reach colon for degradation by bacteria

        Distribution:
            - Volume of distribution: 5-7 L/kg (large, extensive tissue binding)
            - Highly bound to skeletal muscle
            - Crosses blood-brain barrier (CNS toxicity possible)
            - Minimal protein binding (~25%)

        Elimination:
            - Half-life (normal renal function): 36-48 hours
            - Half-life (ESRD): 3.5-5 days (non-renal only)
            - Time to steady state: 5-7 days (normal), 2-3 weeks (ESRD)
            - Renal: Filtration (no significant reabsorption or secretion)
            - Hepatic: ~20-40% of total clearance

        Therapeutic Monitoring:
            - Therapeutic range: 0.8-2.0 ng/mL (or microgram/L)
            - Heart failure: Some prefer 0.5-1.0 ng/mL
            - Atrial fibrillation: Often tolerate higher levels
            - Toxic levels: > 2.0 ng/mL (increased arrhythmia risk)
            - Sample timing: Trough level (just before next dose)
            - Wait 6+ hours post-dose if not trough

    Renal Dosing Adjustments for Digoxin:
        CrCl > 90 mL/min (Normal):
            - Standard dosing: 0.125-0.25 mg daily
            - No adjustment needed
            - Monitor levels periodically

        CrCl 60-89 mL/min (Mild impairment):
            - May use standard doses with monitoring
            - Consider 0.125 mg daily
            - Check levels in 1-2 weeks

        CrCl 30-59 mL/min (Moderate impairment):
            - 0.0625-0.125 mg daily
            - 25-50% dose reduction
            - More frequent monitoring

        CrCl 15-29 mL/min (Severe impairment):
            - 0.0625 mg daily or 0.125 mg every other day
            - 50-75% dose reduction
            - Weekly monitoring initially

        CrCl < 15 mL/min (ESRD):
            - 0.0625 mg every other day
            - 75%+ dose reduction
            - Non-renal clearance becomes dominant
            - Specialist consultation recommended

        Hemodialysis:
            - Digoxin is NOT significantly removed by dialysis
            - Large Vd prevents efficient removal
            - No supplemental dosing after dialysis
            - Monitor levels to guide dosing

    Drug Interactions Affecting Digoxin:
        Increase Digoxin Levels (reduce digoxin dose):
            - Amiodarone: Reduce digoxin by 50%
            - Verapamil: Reduce digoxin by 25-50%
            - Quinidine: Reduce digoxin by 50%
            - Spironolactone: May increase levels
            - Clarithromycin/Erythromycin: Inhibit P-glycoprotein

        Decrease Digoxin Levels (may need higher dose):
            - Rifampin: Induces P-glycoprotein
            - Antacids: Reduce absorption (separate by 2 hours)
            - Cholestyramine: Binds digoxin in gut
            - St. John's Wort: Induces P-glycoprotein

        Increase Toxicity Risk (monitor closely):
            - Hypokalemia (from diuretics): Increases toxicity
            - Hypomagnesemia: Increases toxicity
            - Hypercalcemia: Increases toxicity
            - Hypothyroidism: Increases sensitivity

    Body Weight Considerations for Digoxin:
        Use ACTUAL Body Weight:
            - Normal build (BMI 18.5-29.9)
            - From the book: "Build is normal - no need to use ideal body weight"
            - Hepatic clearance scales with actual weight

        Use IDEAL Body Weight:
            - Obese patients (BMI >= 30)
            - From the book: "For a very obese patient, ideal body weight
              should be substituted for actual body weight"
            - Digoxin does not distribute well into adipose tissue
            - Using actual weight overestimates clearance

        IBW Formulas:
            Male:   IBW = 50 + 2.3 x (height in inches - 60)
            Female: IBW = 45.5 + 2.3 x (height in inches - 60)

    Digoxin Toxicity:
        Signs and Symptoms:
            - GI: Nausea, vomiting, anorexia (often first signs)
            - Cardiac: Bradycardia, AV block, ventricular arrhythmias
            - CNS: Confusion, visual disturbances (yellow-green halos)
            - Fatigue, weakness

        Risk Factors for Toxicity:
            - Renal impairment (reduced clearance)
            - Electrolyte abnormalities (K+, Mg2+, Ca2+)
            - Drug interactions
            - Advanced age
            - Low body weight
            - Hypothyroidism

        Management:
            - Hold digoxin
            - Correct electrolytes
            - Digoxin-specific antibody fragments (Digibind) for severe toxicity
            - Supportive care

    Disease State Considerations:
        Heart Failure:
            - Lower target levels (0.5-1.0 ng/mL) may be preferred
            - Mortality benefit seen at lower levels
            - Watch for worsening renal function

        Atrial Fibrillation:
            - Rate control indication
            - May tolerate higher levels
            - Often combined with other rate-control agents

        Hepatic Disease:
            - Minimal effect on digoxin clearance (20-40% hepatic)
            - Main concern is volume of distribution changes (ascites)
            - Loading dose may need adjustment

        Thyroid Disease:
            - Hypothyroidism: Increased sensitivity, reduce dose
            - Hyperthyroidism: Decreased sensitivity, may need higher dose
            - Adjust based on thyroid status

    Population-Specific Adjustments:
        Elderly Patients:
            - Reduced renal function (even with normal creatinine)
            - Lower muscle mass may overestimate CrCl
            - Increased sensitivity to toxicity
            - Start low, go slow
            - Target lower end of therapeutic range

        Pediatric Patients:
            - Different PK parameters than adults
            - Weight-based dosing in mg/kg
            - Consult pediatric dosing references

        Critically Ill:
            - Rapidly changing renal function
            - Volume shifts affect Vd
            - More frequent monitoring needed

    Clinical Pearls:
        1. The Cockcroft-Gault equation estimates renal clearance
        2. Add hepatic clearance (0.33 mL/min/kg x weight) for total
        3. Convert to L/hour: multiply mL/min by 60/1000
        4. Target Css,av: 0.8-2.0 microgram/L (use 1.4 as midpoint)
        5. Digoxin tablets: multiples of 62.5 microgram
        6. Round calculated dose to nearest tablet strength
        7. Monitor levels after 5-7 half-lives (7+ days)
        8. Always correct electrolytes before adjusting digoxin

    Notes:
        - The 0.33 mL/min/kg constant is from Rowe's textbook
        - Some sources use slightly different values (0.33-0.57)
        - Excel spreadsheets available from phrData.co.uk for calculations
        - Always verify against current clinical guidelines
        - Therapeutic drug monitoring essential for this narrow TI drug
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    # Convert inputs to pint Quantities
    dig_cl = kwargs.get("digoxin_clearance", False)
    crcl = kwargs.get("creatinine_clearance", False)
    weight = kwargs.get("weight", False)

    if dig_cl:
        dig_cl = Q_(dig_cl)
    if crcl:
        crcl = Q_(crcl)
    if weight:
        weight = Q_(weight)

    # Count provided parameters
    provided = sum([bool(dig_cl), bool(crcl), bool(weight)])
    if provided != 2:
        raise ValueError(
            f"digoxin_clearance requires exactly 2 of 3 parameters. "
            f"Got {provided}: digoxin_clearance={dig_cl is not False}, "
            f"creatinine_clearance={crcl is not False}, weight={weight is not False}"
        )

    # Non-renal clearance constant
    NON_RENAL_CONSTANT = Q_(0.33, "mL/min/kg")

    # Formula: Digoxin Cl = CrCl + 0.33 mL/min/kg x Weight
    # Rearranged forms for solving each variable

    if not dig_cl:
        # Solve for Digoxin Clearance
        non_renal = NON_RENAL_CONSTANT * weight.to("kg")
        quantity = crcl.to("mL/min") + non_renal
        string = "Digoxin Clearance"

    elif not crcl:
        # Solve for CrCl: CrCl = Digoxin Cl - 0.33 x Weight
        non_renal = NON_RENAL_CONSTANT * weight.to("kg")
        quantity = dig_cl.to("mL/min") - non_renal
        string = "Creatinine Clearance"

    else:  # not weight
        # Solve for Weight: Weight = (Digoxin Cl - CrCl) / 0.33
        quantity = (dig_cl.to("mL/min") - crcl.to("mL/min")) / NON_RENAL_CONSTANT
        string = "Weight"

    return format_output(quantity, string, output_unit, decimals)
