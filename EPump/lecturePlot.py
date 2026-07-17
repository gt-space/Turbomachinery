"""
Plot the saturated liquid-to-vapor density ratio vs. temperature to explain
why cavitation is more destructive at COLDER temperatures, even though the
vapor pressure (and therefore NPSH risk of onset) is lower.

Physics:
  When local pressure dips below vapor pressure, a small mass of liquid
  flashes to vapor. The VOLUME of vapor that mass occupies is set by the
  liquid/vapor density ratio (rho_liq / rho_vap). Near the critical point
  this ratio -> 1 (vapor is nearly as dense as liquid), so flashing barely
  changes volume and is self-limiting. Near the triple point (cold water),
  vapor is 10,000-200,000x less dense than liquid, so the SAME mass of
  flashed liquid creates a vastly larger vapor volume -> bigger, more
  violent cavitation bubbles, more collapse energy, more damage -- despite
  a lower absolute vapor pressure and larger nominal NPSH margin.

Install dependencies first if needed:
    pip install CoolProp matplotlib numpy
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# ---------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------
fluid = "Oxygen"
n_points = 300

# ---------------------------------------------------------------------
# Bound temperature range between triple point and critical point
# ---------------------------------------------------------------------
T_triple = CP.PropsSI(fluid, "Ttriple")
T_crit = CP.PropsSI(fluid, "Tcrit")

T_min = T_triple + 0.5
T_max = T_crit - 0.5

T = np.linspace(T_min, T_max, n_points)
T_K = T  # already in Kelvin

rho_liq = np.array([CP.PropsSI("D", "T", Ti, "Q", 0, fluid) for Ti in T])
rho_vap = np.array([CP.PropsSI("D", "T", Ti, "Q", 1, fluid) for Ti in T])
P_sat = np.array([CP.PropsSI("P", "T", Ti, "Q", 0, fluid) for Ti in T]) / 1e5  # bar

ratio = rho_liq / rho_vap  # volume expansion factor when liquid flashes to vapor

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(T_K, ratio, color="tab:red", linewidth=2.5)
ax.set_yscale("log")
ax.set_xlabel("Temperature [K]")
ax.set_ylabel(r"Liquid / Vapor Density Ratio  ($\rho_{liq} / \rho_{vap}$)")
ax.set_title(f"{fluid}: Vapor Bubble Expansion Factor Along Saturation Curve")
ax.grid(True, which="both", alpha=0.3)

# Annotate a few representative operating points
annotate_T = [60, 70, 80, 90.2, 110, 140]
for Tc in annotate_T:
    idx = (np.abs(T_K - Tc)).argmin()
    ax.plot(T_K[idx], ratio[idx], "o", color="black", zorder=5)
    ax.annotate(f"{ratio[idx]:,.0f}x\n@ {Tc}K",
                (T_K[idx], ratio[idx]),
                textcoords="offset points", xytext=(10, 10),
                fontsize=9, ha="left")


plt.tight_layout()
plt.show()

