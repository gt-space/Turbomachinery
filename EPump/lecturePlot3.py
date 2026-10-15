"""Plot the thermal suppression head (TSH) of oxygen vs. temperature (80 K to 120 K) using CoolProp.

Model (quasi-static energy balance):
    dT   = B * rho_v * h_fg / (rho_l * cp_l)      cooling of the liquid at the cavity
    dp_v = p_sat(T) - p_sat(T - dT)               resulting drop in vapor pressure
    TSH  = dp_v / (rho_l * g)                     expressed as head of liquid

B is the vapor-to-liquid volume ratio: the volume of vapor formed divided by the
volume of liquid that supplies the latent heat (Stepanoff's "B-factor"; written as
V_v/V_l by Moore & Ruggeri, NASA TN D-5292). It is not
a fluid property: it depends on the inducer geometry, speed and cavity size, so the
curves below show the fluid-property trend for a few assumed values of B.

Requires: pip install CoolProp numpy matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

FLUID = "Oxygen"
T_MIN = 80.0    # K
T_MAX = 120.0   # K
N_POINTS = 200
B_VALUES = [0.5, 1.0, 2.0]                     # assumed vapor-to-liquid volume ratios
COLORS = ["#86b6ef", "#2a78d6", "#104281"]     # one hue, light -> dark with increasing B
G = 9.80665     # m/s^2

T = np.linspace(T_MIN, T_MAX, N_POINTS)


def thermal_suppression_head(T, B):
    """Return TSH in metres of liquid oxygen at temperature(s) T (K) for a given vapor-to-liquid volume ratio B."""
    h_fg = PropsSI("H", "T", T, "Q", 1, FLUID) - PropsSI("H", "T", T, "Q", 0, FLUID)  # J/kg
    rho_v = PropsSI("D", "T", T, "Q", 1, FLUID)   # kg/m^3
    rho_l = PropsSI("D", "T", T, "Q", 0, FLUID)   # kg/m^3
    cp_l = PropsSI("C", "T", T, "Q", 0, FLUID)    # J/kg/K
    p_sat = PropsSI("P", "T", T, "Q", 0, FLUID)   # Pa

    dT = B * rho_v * h_fg / (rho_l * cp_l)                    # K
    dp_v = p_sat - PropsSI("P", "T", T - dT, "Q", 0, FLUID)   # Pa
    return dp_v / (rho_l * G)


fig, ax = plt.subplots(figsize=(8, 5))
for B, color in zip(B_VALUES, COLORS):
    tsh = thermal_suppression_head(T, B)
    ax.plot(T, tsh, linewidth=2, color=color, label=f"B = {B:g}")
    # Direct label at the right-hand end of each curve


ax.set_yscale("log")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Thermal suppression head (m)")

ax.set_xlim(T_MIN, T_MAX)
ax.grid(True, which="major", alpha=0.3)
ax.grid(True, which="minor", alpha=0.12)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(title="Vapor-to-liquid volume ratio, B", frameon=False, loc="upper left")
fig.tight_layout()

fig.savefig("oxygen_tsh.png", dpi=150)
plt.show()

# Print a few values for reference
T_ref = np.array([80.0, 90.0, 100.0, 110.0, 120.0])
print("T (K)  " + "  ".join(f"TSH(B={B:g}) m" for B in B_VALUES))
for i, T_i in enumerate(T_ref):
    print(f"{T_i:5.0f}  " + "  ".join(f"{thermal_suppression_head(T_ref, B)[i]:12.3f}" for B in B_VALUES))