from .lib import format_output, generic_a_eq_b_x_c
from .reduction_factors import salt_factor, bioavailability
from .functions import *
from .expo import solve_for_c_t, solve_for_c_0, solve_for_k, solve_for_t
from .advanced import (
    loading_dose,
    maintenance_dose,
    time_to_steady_state,
    accumulation_factor,
    infusion_rate,
)

# AUC calculations (Chapter 5, 9)
from .auc import auc_dose_clearance, auc_tail, trapezoidal_auc

# Multiple dosing peak/trough (Chapter 10)
from .multiple_dosing import css_max, css_min, fluctuation

# Renal clearance and Cockcroft-Gault (Chapter 14)
from .renal import (
    cockcroft_gault_male,
    cockcroft_gault_female,
    cockcroft_gault,
    digoxin_clearance,
)

# Oral/extravascular administration (Chapter 9)
from .oral import absorption_rate, tmax, cmax, flip_flop_check

# Two-compartment model (Chapter 7)
from .two_compartment import (
    concentration_two_compartment,
    alpha_beta_from_micro,
    volume_central,
    volume_steady_state,
    distribution_half_life,
    terminal_half_life,
)

# Non-linear kinetics - Michaelis-Menten (Chapter 11)
from .nonlinear import (
    michaelis_menten_rate,
    apparent_km,
    is_linear,
    phenytoin_steady_state,
    time_to_eliminate_nonlinear,
)
