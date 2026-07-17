"""
Plot the saturation curve of a fluid along with the saturated liquid and
vapor density using CoolProp.

Install dependencies first if needed:
    pip install CoolProp matplotlib numpy
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# ---------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------
fluid = "oxygen"          # Any fluid name supported by CoolProp, e.g. "R134a", "Nitrogen", "CO2"
n_points = 200            # Number of points along the saturation curve

# ---------------------------------------------------------------------
# Get critical point and triple point temperatures to bound the curve
# ---------------------------------------------------------------------
T_triple = CP.PropsSI(fluid, "Ttriple")
T_crit = CP.PropsSI(fluid, "Tcrit")

# Stay a hair away from the critical point to avoid numerical issues
T_min = T_triple + 0.5
T_max = T_crit - 0.5

T = np.linspace(T_min, T_max, n_points)

P_sat = np.zeros(n_points)
rho_liq = np.zeros(n_points)
rho_vap = np.zeros(n_points)

for i, Ti in enumerate(T):
    P_sat[i] = CP.PropsSI("P", "T", Ti, "Q", 0, fluid)
    rho_liq[i] = CP.PropsSI("D", "T", Ti, "Q", 0, fluid)
    rho_vap[i] = CP.PropsSI("D", "T", Ti, "Q", 1, fluid)

fig, ax2 = plt.subplots(figsize=(9, 6))

# Density along the saturation curve (liquid and vapor branches)
ax2.plot(T, rho_liq, label="Saturated liquid", color="tab:orange")
ax2.plot(T, rho_vap, label="Saturated vapor", color="tab:purple")
ax2.set_xlabel("Temperature [°K]")
ax2.set_ylabel("Density [kg/m³]")
ax2.set_title(f"{fluid} Saturation Density")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()