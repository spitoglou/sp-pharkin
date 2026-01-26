from .lib import generic_a_eq_b_x_c, format_output
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity  # type: ignore[assignment]


def salt_factor(**kwargs):
    """
    Calculate delivered drug amount using salt factor, or solve for salt dose/factor.

    When drugs are formulated as salts (e.g., amoxicillin trihydrate), the salt
    itself adds mass but doesn't contribute therapeutic activity. The salt factor
    accounts for the fraction of the dose that is active drug versus inert salt.

    Formula: Delivered Drug = Salt Dose × Salt Factor

    Or: Amount of active drug = Total salt dose × (Molecular weight of active / Molecular weight of salt)

    The salt factor represents:
    - Salt factor = 1.0: Pure drug, no salt (e.g., caffeine)
    - Salt factor = 0.8: 80% active drug, 20% salt (typical for most antibiotics)
    - Salt factor = 0.5: 50% active, 50% salt (e.g., some salts)
    - Salt factor = 0.3: 30% active, 70% salt (e.g., benzathine penicillin G)

    This ensures dosing is based on ACTIVE drug, not total salt mass.

    Args (provide exactly 2 of 3):
        delivered_drug (str): Amount of active drug (e.g., '250 mg')
        dose_of_salt (str): Total mass of salt formulation (e.g., '300 mg')
        salt_factor (str): Fraction that is active (e.g., '0.833', dimensionless)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Delivered Drug', 250.0, 'milligram', '250.0 mg', ...)

    Examples:
        Calculate active amoxicillin from 250 mg capsule with salt factor 0.833:
            >>> salt_factor(dose_of_salt='300 mg', salt_factor='0.833')
            ('Delivered Drug', 250.0, 'milligram', '250.0 mg', ...)

        Determine required capsule dose for 250 mg active:
            >>> salt_factor(delivered_drug='250 mg', salt_factor='0.833')
            ('Dose of Salt', 300.0, 'milligram', '300.0 mg', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 4.1-4.2, pages 2700-2800: Drug formulation and salt factors
        - Section 4.3, pages 2800-2850: Salt form pharmacokinetics
        - Table 4.1, pages 2750-2760: Salt factors for common drugs
        - Chapter 4: Bioavailability and formulation effects

    Common Salt Forms and Their Factors:
        Antibiotics:
            - Penicillin V potassium: SF ≈ 0.6 (only 60% penicillin V)
            - Amoxicillin trihydrate: SF ≈ 0.833 (83.3% amoxicillin)
            - Erythromycin stearate: SF ≈ 0.33 (only 33% erythromycin)
            - Tetracycline HCl: SF ≈ 0.918 (91.8% tetracycline)
        Cardiac drugs:
            - Digoxin: SF = 1.0 (no salt form, pure)
            - Digitoxin: SF = 1.0 (no salt form, pure)
        Others:
            - Caffeine citrate: SF ≈ 0.486 (only 48.6% caffeine)
            - Morphine sulfate: SF ≈ 0.776 (77.6% morphine)

    Molecular Weight Calculation:
        Penicillin V Potassium example:
            - Penicillin V: 350.4 g/mol
            - Penicillin V potassium salt: 588.7 g/mol
            - Salt factor = 350.4 / 588.7 ≈ 0.595
            - So 500 mg tablet = 500 × 0.595 ≈ 297.5 mg active penicillin V

        Amoxicillin Trihydrate example:
            - Amoxicillin: 365.4 g/mol
            - Amoxicillin trihydrate: 438.4 g/mol
            - Salt factor = 365.4 / 438.4 ≈ 0.833
            - So 300 mg capsule = 300 × 0.833 ≈ 250 mg active amoxicillin

    Why Salt Forms Are Used:
        1. Stability: Salts often more stable than free base (less hygroscopic)
        2. Solubility: Some salts more soluble (better absorption)
        3. Taste: Salt forms can mask bitter taste of free drug
        4. Shelf-life: Some formulations require salt for longer expiration
        5. Patency: Different salt forms allow new patents on old drugs
        6. Specific tissue distribution: Different salts may have different Vd

    Clinical Importance:
        - Prescribers should specify ACTIVE DRUG amount, not salt mass
        - "Amoxicillin 250 mg" means 250 mg active amoxicillin
        - Actual capsule mass ≈ 300 mg (due to trihydrate salt form)
        - Pharmacy dispensing must account for salt factor
        - Labeling should show both salt mass and active amount

    Therapeutic Equivalence:
        - Different salt forms of same drug are therapeutically equivalent (if properly dosed)
        - Example: Penicillin V potassium vs Penicillin V benzathine
            * Same active drug (penicillin V)
            * Different salt factors
            * Different dosing (benzathine is depot form)
        - Pharmacokinetics may differ between salt forms
        - Some salt forms designed for specific properties (depot, IM vs oral)

    Drug Interaction with Salt:
        - Salt portion may affect absorption, distribution, or clearance
        - Example: Lithium carbonate vs Lithium citrate
            * Different salts
            * Different absorption rates
            * Different Vd and clearance
            * Not automatically interchangeable at same weight
        - Example: Morphine sulfate vs Morphine tartrate
            * Different pharmacokinetic profiles
            * Not simply replaceable

    Pediatric Dosing Considerations:
        - Salt factor critical when scaling doses for children
        - Weight-based dosing uses ACTIVE drug amount
        - Example: Child 20 kg needs amoxicillin 25 mg/kg = 500 mg active
        - Capsule size: 500/0.833 ≈ 600 mg (salt form) - often not available
        - May need liquid suspension (active mg/mL) instead

    Renal Function Impact:
        - Salt portion may affect renal clearance differently than active drug
        - Example: Lithium salts interact with renal handling
        - Renal disease may affect salt form handling
        - Some drugs with active metabolites have salt-form variations

    Notes:
        - Always dose based on ACTIVE drug content, not salt mass
        - Check product labeling for salt factor/active content
        - Different manufacturers may have different salt forms
        - Substituting one salt form for another may require dose adjustment
        - International markets may use different salt forms
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("delivered_drug", 0)
    b = kwargs.get("dose_of_salt", 0)
    c = kwargs.get("salt_factor", 0)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Delivered Drug", "Dose of Salt", "Salt Factor"]
    )

    return format_output(quantity, string, output_unit, decimals)


def bioavailability(**kwargs):
    """
    Calculate bioavailability, delivered drug, or administered dose relationship.

    Bioavailability (F) represents the fraction of an administered dose that reaches
    systemic circulation as unchanged drug. It accounts for all losses: incomplete
    absorption, first-pass metabolism, and any other pre-systemic elimination.

    Formula: Delivered Drug = Dose Administered × Bioavailability

    Or: Amount in systemic circulation = Dose × F

    Or: F = AUC_oral / AUC_IV (fundamental definition)

    Bioavailability ranges:
        F = 1.0 (100%): Complete absorption, no first-pass loss (IV administration)
        F = 0.8-0.99: Excellent bioavailability (minimal absorption/metabolism loss)
        F = 0.5-0.8: Good bioavailability (moderate first-pass metabolism)
        F = 0.2-0.5: Moderate bioavailability (significant first-pass loss)
        F < 0.2: Poor bioavailability (severe first-pass metabolism)
        F < 0.01: Essentially inactive orally (must use alternate route)

    Args (provide exactly 2 of 3):
        delivered_drug (str): Amount reaching systemic circulation (e.g., '200 mg')
        dose_administered (str): Total dose given orally/IM/etc (e.g., '500 mg')
        bioavailability (str): Fraction reaching systemic (e.g., '0.4', dimensionless)

    Returns:
        5-tuple: (name, magnitude, unit_string, formatted_string, pint_quantity)
        Calculates whichever parameter is missing.
        Example: ('Delivered Drug', 200.0, 'milligram', '200.0 mg', ...)

    Examples:
        Calculate bioavailable amount from 500 mg oral dose with F=0.4:
            >>> bioavailability(dose_administered='500 mg', bioavailability='0.4')
            ('Delivered Drug', 200.0, 'milligram', '200.0 mg', ...)

        Determine required oral dose to deliver 200 mg with F=0.4:
            >>> bioavailability(delivered_drug='200 mg', bioavailability='0.4')
            ('Dose Administered', 500.0, 'milligram', '500.0 mg', ...)

    Reference:
        Rowe, P. Pharmacokinetics
        - Section 4.2-4.3, pages 2750-2850: Bioavailability definition and measurement
        - Section 4.4, pages 2850-2950: First-pass metabolism
        - Section 4.5, pages 2950-3050: Route of administration effects
        - Table 4.2, pages 2800-2820: Bioavailability of common drugs
        - Chapter 4: Absorption and bioavailability

    Bioavailability by Route of Administration:
        IV (Intravenous):
            F = 1.0 (100%, by definition)
            - Direct systemic access
            - No absorption or metabolism loss
            - Reference standard for bioavailability calculations
        Oral (PO):
            F = 0.01-0.99 (highly variable)
            - Affected by: gastric pH, food, GI motility, metabolism
            - Subject to first-pass hepatic metabolism
            - Most common route, widely used despite variable F
        IM (Intramuscular):
            F = 0.8-1.0 (usually very good)
            - Avoids first-pass metabolism
            - Some depot formulations have lower/delayed F
            - Faster/more reliable than oral for most drugs
        Sublingual (SL):
            F = 0.5-0.95 (highly variable)
            - Partly avoids first-pass metabolism
            - Rich blood supply allows rapid absorption
            - Useful for drugs with poor oral F
        Transdermal:
            F = 0.3-0.9 (variable, formulation-dependent)
            - Avoids first-pass metabolism (some drugs)
            - Slow absorption (controlled delivery)
            - Good for chronic therapy
        Rectal:
            F = 0.3-0.8 (variable)
            - Partial bypass of first-pass metabolism
            - Lower rectal blood returns to systemic (avoids portal circulation)
            - Upper rectal blood goes through portal circulation

    First-Pass Metabolism Effect on Bioavailability:
        Hepatic extraction ratio E, oral F calculation:
            F_oral = 1 - E_hepatic × E_intestinal

        Examples:
            - Propranolol: E ≈ 0.7 → F_oral ≈ 0.2 (only 20%)
            - Nitroglycerin: E ≈ 0.95 → F_oral < 0.1 (nearly complete)
            - Warfarin: E ≈ 0.05 → F_oral ≈ 0.95 (minimal first-pass)
            - Theophylline: E ≈ 0.1 → F_oral ≈ 0.9 (good bioavailability)

    Drug-Specific Bioavailability Examples:
        High F drugs (can use oral dosing):
            - Warfarin: F ≈ 0.93-0.96 (reliable oral dosing possible)
            - Theophylline: F ≈ 0.90-0.95 (can use oral)
            - Aspirin: F ≈ 0.50-0.80 (depends on formulation)
            - Ibuprofen: F ≈ 0.80 (good oral absorption)
        Moderate F drugs (oral possible, need adjustment):
            - Morphine: F ≈ 0.15-0.30 (oral dose ~3× IV dose)
            - Codeine: F ≈ 0.45-0.60 (better than morphine)
            - Metoprolol: F ≈ 0.20-0.40 (significant first-pass)
            - Propranolol: F ≈ 0.20-0.30 (requires higher oral doses)
        Low F drugs (oral impractical):
            - Nitroglycerin: F < 0.10 (sublingual or patch required)
            - Insulin: F ≈ 0 orally (destroyed by stomach acid, must inject)
            - Gentamicin: F ≈ 0 orally (destroyed in GI, used systemically only)
            - Digoxin: F ≈ 0.60-0.80 (but highly variable, TDM recommended)

    Factors Affecting Bioavailability:
        GI Factors:
            - Food: ↓ F for some drugs (e.g., tetracyclines)
                    ↑ F for others (e.g., griseofulvin, itraconazole)
            - pH: Affects ionization and absorption
            - GI motility: Faster transit → ↓ absorption time
            - First-pass metabolism: Lipophilic drugs absorbed via lymphatics
        Patient Factors:
            - Age: Geriatric patients often have ↓ Vd, ↑ F
            - Genetic polymorphisms: CYP450 variants affect F significantly
            - Disease: Liver disease ↑ F (↓ first-pass metabolism)
            - Stress: GI blood flow changes may affect F
        Drug Factors:
            - Lipophilicity: More lipophilic → better absorption
            - Molecular weight: Higher MW → lower absorption
            - Protein binding: Highly bound drugs less bioavailable
            - Formulation: Tablet vs solution vs suspension → different F

    Disease Effects on Bioavailability:
        Liver Disease:
            - ↑ F for most drugs (↓ first-pass metabolism)
            - Clearance also ↓ (may need dose reduction)
            - Net effect: Unpredictable - some drugs need ↓ dose
        Portal Hypertension/Cirrhosis:
            - Portosystemic shunting bypasses liver
            - Drug goes directly to systemic circulation
            - F → 1.0 for hepatically metabolized drugs
            - May need significant dose ↓
        Inflammatory Bowel Disease:
            - Reduced absorptive surface area
            - F may be ↓ for some drugs
            - Diarrhea ↓ transit time → ↓ absorption
        GI Surgery:
            - Gastric bypass: ↑ F (avoids acidic environment)
            - Intestinal resection: ↓ F (↓ absorptive surface)
            - Ileostomy: ↓ F for many drugs

    Drug Interactions Affecting Bioavailability:
        Enzyme Induction (↑ first-pass metabolism):
            - Rifampin induces CYP3A4
            - Result: ↓ F for CYP3A4 substrates
            - Example: Oral contraceptives ineffective with rifampin
        Enzyme Inhibition (↓ first-pass metabolism):
            - Ketoconazole inhibits CYP3A4
            - Result: ↑ F for CYP3A4 substrates
            - Example: ↑ Triazolam levels (toxicity risk)
        Transporter Inhibition:
            - P-gp inhibitors ↑ F
            - Example: Verapamil ↑ digoxin bioavailability

    Dose Adjustment for Different Routes:
        Same drug, different routes, same desired effect:
            IV dose < IM dose < Oral dose

        Example (Propranolol):
            - IV: 10 mg for effect
            - IM: ~12 mg (F ≈ 0.85)
            - Oral: 50-160 mg (F ≈ 0.2-0.3)

        Formula: Dose_oral = Dose_IV / F_oral

    Therapeutic Drug Monitoring with Bioavailability:
        - TDM measures systemic levels (already accounted for F)
        - Target levels same regardless of route (if F accounted for in dose)
        - If F varies significantly between patients:
            * TDM helps individualize dosing
            * More important for narrow therapeutic index drugs
            * Example: Theophylline (F = 0.90 but inter-patient variation)

    Clinical Decision Making Using Bioavailability:
        Route Selection:
            - High F drug: Oral usually preferred (convenience, cost)
            - Low F drug: Use IV or alternative route
            - Moderate F: May use oral with dose adjustment
        Dose Calculation:
            - Calculate IV dose needed for desired effect
            - Multiply by 1/F to get equivalent oral dose
            - Example: Need 400 mg systemically, drug has F=0.5
                     Oral dose = 400 / 0.5 = 800 mg
        Drug Interactions:
            - Know which drugs affect F
            - May need dose adjustment when adding/removing CYP inducers/inhibitors

    Notes:
        - Bioavailability is drug and formulation specific (not universal)
        - F for sustained-release may differ from immediate-release
        - F can vary significantly between individuals (genetic/disease factors)
        - Always refer to product-specific F values, not general assumptions
        - Food effects on F are drug-specific (some ↑, some ↓, some unchanged)
    """
    output_unit = kwargs.pop("output_unit", False)
    decimals = kwargs.pop("decimals", 2)

    kwargs = {k: Q_(v) for k, v in kwargs.items()}

    a = kwargs.get("delivered_drug", 0)
    b = kwargs.get("dose_administered", 0)
    c = kwargs.get("bioavailability", 0)

    string, quantity = generic_a_eq_b_x_c(
        a, b, c, ["Delivered Drug", "Dose Administered", "Bioavailability"]
    )

    return format_output(quantity, string, output_unit, decimals)
