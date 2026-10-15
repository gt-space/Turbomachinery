"""
  1. Whirl amplitude vs frequency ratio (w/wn) from the full ROSS FE
     model, one curve per ISO 21940-11 balance grade.
  2. Analytical single-DOF (Jeffcott) overlay:
         X = M e w^2 / ((K - w^2 M) + i w C)
  3. Effective stiffness identification:
         K_eff = M * wn^2   
     plus a bearing-stiffness sweep and optional back-solve from a
     measured critical speed.

Run:  python RotorModelV2.py   (writes interactive .html plots)
"""

import os

import numpy as np
import plotly.graph_objects as go
import ross as rs

# Geometry and disks come straight from V1 -- single source of truth.
from buildRotor import (
    rotor,
    shaft_elements,
    disk1, disk2,
    K_RADIAL,
    DOWNLOADS_DIR,
    G_GRADES,
    OMEGA_SERVICE,
    BEARING_CXX,
    MEASURED_CRITICAL_HZ,
    RESPONSE_NODE,
    R_MAX,
    DISK_NODES,
    OPERATING_SPEED_RPM,
    enableDashboard
)

import ModeShapeAnalysis
from analysisDashboard import build_dashboard

# ----------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------



# ----------------------------------------------------------------------
# Rotor (V1 geometry + damped bearings)
# ----------------------------------------------------------------------
def rotorGeometry():
    # --- Rotor geometry plot ---------------------------------------------
    os.makedirs(DOWNLOADS_DIR, exist_ok=True)

    try:
        fig = rotor.plot_rotor()
        out_path = os.path.join(DOWNLOADS_DIR, "rotor_geometry.html")
        fig.write_html(out_path)
        print(f"Saved rotor geometry -> {out_path}")
    except Exception as exc:  # plotting is optional
        print(f"[plot_rotor skipped] {exc}")
def cambellDiagram():
    # --- Campbell diagram ------------------------------------------------
    # Adjust the speed range to your machine's operating range (rad/s).
    speed_range = np.linspace(0, 5000, 41)
    try:
        campbell = rotor.run_campbell(speed_range)
        out_path = os.path.join(DOWNLOADS_DIR, "campbell.html")
        campbell.plot().write_html(out_path)
        print(f"\nSaved Campbell diagram -> {out_path}")
    except Exception as exc:
        print(f"\n[Campbell skipped] {exc}")

def undampedCriticalSpeedMap():
    # --- Undamped critical speed map -------------------------------------
    try:
        ucs = rotor.run_ucs()
        out_path = os.path.join(DOWNLOADS_DIR, "ucs_map.html")
        ucs.plot().write_html(out_path)
        print(f"Saved undamped critical speed map -> {out_path}")
    except Exception as exc:
        print(f"[UCS skipped] {exc}")

os.makedirs(DOWNLOADS_DIR, exist_ok=True)


def first_critical(rtr):
    """First undamped natural frequency [rad/s] at rest."""
    return rtr.run_modal(speed=0.0).wn[0]


def _magnitude_at_node(results, node):
    ndof = results.rotor.number_dof
    return np.abs(np.asarray(results.forced_resp))[ndof * node, :]


# ----------------------------------------------------------------------
# 1 + 2. Whirl amplitude vs frequency ratio (FE + Jeffcott overlay)
# ----------------------------------------------------------------------
def whirl_vs_frequency_ratio(wn, filename="whirl_vs_freq_ratio.html"):
    speed_range = np.linspace(0.02 * wn, R_MAX * wn, 400)
    r = speed_range / wn

    fig = go.Figure()


    # --- FE model response, one curve per balance grade -------------
    m_rotor = rotor.m
    for G in G_GRADES:
        # permissible unbalance sized at the SERVICE speed (ISO 21940-11)
        U_total = m_rotor * (G / 1000.0) / OMEGA_SERVICE      # kg.m
        resp = rotor.run_unbalance_response(
            node=list(DISK_NODES),
            unbalance_magnitude=[U_total * 0.3, U_total * 0.7],
            unbalance_phase=[0.0, np.pi],
            frequency=speed_range,
        )
        amp_um = _magnitude_at_node(resp, RESPONSE_NODE) * 1e6
        fig.add_trace(go.Scatter(
            x=r, y=amp_um, mode="lines",
            name=f"FE model, G{G} (U = {U_total*1e6:.4f} g.mm)"))

    # --- analytical Jeffcott overlay --------------
    # X = M e w^2 / sqrt((K - w^2 M)^2 + (w C)^2), with K = M wn^2
    K_eff = m_rotor * wn**2
    C_eff = 2 * BEARING_CXX                    # two bearings in parallel
    zeta = C_eff / (2 * m_rotor * wn)
    for G in G_GRADES:
        e = (G / 1000.0) / OMEGA_SERVICE       # eccentricity [m]
        w = speed_range
        X = (m_rotor * e * w**2 /
             np.sqrt((K_eff - m_rotor * w**2)**2 + (w * C_eff)**2))
        fig.add_trace(go.Scatter(
            x=r, y=X * 1e6, mode="lines", line=dict(dash="dot"),
            name=f"Jeffcott 1-DOF, G{G}", visible="legendonly"))

    fig.add_vline(x=OMEGA_SERVICE / wn, line_dash="dash",
                  annotation_text=f"operating speed (r = "
                                  f"{OMEGA_SERVICE / wn:.2f})")
    fig.update_xaxes(title="Frequency ratio  w / wn")
    fig.update_yaxes(title="Whirl amplitude X [um, 0-pk]")
    fig.update_layout(
        title=(f"Whirl amplitude vs frequency ratio  "
               f"(wn = {wn/(2*np.pi):.1f} Hz = "
               f"{wn*60/(2*np.pi):.0f} rpm,  zeta = {zeta:.4f})"))
    out_path = os.path.join(DOWNLOADS_DIR, filename)
    fig.write_html(out_path)
    print(f"Saved whirl amplitude plot -> {out_path}")
    return K_eff, zeta


# ----------------------------------------------------------------------
# 3. Stiffness identification
# ----------------------------------------------------------------------
def stiffness_study(filename="stiffness_vs_critical.html"):
    k_values = np.logspace(5, 10, 40)
    n_modes = 4
    wn_table = np.zeros((len(k_values), n_modes))
    for i, k in enumerate(k_values):
        wn_table[i, :] = (rotor.build_rotor(k_radial=k, cxx=BEARING_CXX)
                          .run_modal(speed=0.0).wn[:n_modes]) / (2 * np.pi)

    fig = go.Figure()
    for m in range(n_modes):
        fig.add_trace(go.Scatter(x=k_values, y=wn_table[:, m],
                                 mode="lines", name=f"mode {m + 1}"))
    fig.add_vline(x=K_RADIAL, line_dash="dash",
                  annotation_text=f"current k = {K_RADIAL:.2e} N/m")
    fig.update_xaxes(type="log", title="Bearing radial stiffness [N/m]")
    fig.update_yaxes(type="log", title="Natural frequency [Hz]")
    fig.update_layout(title="Natural frequencies vs bearing stiffness")
    out_path = os.path.join(DOWNLOADS_DIR, filename)
    fig.write_html(out_path)
    print(f"Saved stiffness study -> {out_path}")

    if MEASURED_CRITICAL_HZ is not None:
        f1 = wn_table[:, 0]
        if f1.min() <= MEASURED_CRITICAL_HZ <= f1.max():
            k_id = np.interp(MEASURED_CRITICAL_HZ, f1, k_values)
            print(f"[stiffness ID] {MEASURED_CRITICAL_HZ:.1f} Hz measured "
                  f"-> bearing k_radial ~= {k_id:.3e} N/m "
                  f"(model uses {K_RADIAL:.2e} N/m)")
        else:
            print("[stiffness ID] measured value outside sweep; "
                  "widen k_values.")
    else:
        print("[stiffness ID] set MEASURED_CRITICAL_HZ to back-solve "
              "bearing stiffness from a measured resonance.")


# ----------------------------------------------------------------------
def main():
    print("=" * 60)
    print("Rotor Model Analysis")
    print("=" * 60)
    print()
    print("=" * 60)
    print("Geometry Summary")
    print("=" * 60)
    print(f"Number of shaft elements : {len(shaft_elements)}")
    print(f"Number of nodes          : {rotor.nodes[-1] + 1}")
    print(f"Total length (m)         : {rotor.nodes_pos[-1]:.6f}")
    print(f"Total mass (kg)          : {rotor.m:.6f}")
    print(f"Centre of gravity (m)    : {rotor.CG:.6f}")
    print("Disk 1: m = {:.4e} kg  Id = {:.4e}  Ip = {:.4e}  (kg.m^2)".format(disk1.m, disk1.Id, disk1.Ip))
    print("Disk 2: m = {:.4e} kg  Id = {:.4e}  Ip = {:.4e}  (kg.m^2)".format(disk2.m, disk2.Id, disk2.Ip))
    print()
    rotorGeometry()
    cambellDiagram()
    undampedCriticalSpeedMap()
    print("\n" + "=" * 60)
    print("Critical Speeds and Whirl Amplitude")
    wn = first_critical(rotor)
    print(f"First critical, wn      : {wn:.1f} rad/s = "
          f"{wn/(2*np.pi):.1f} Hz = {wn*60/(2*np.pi):.0f} rpm")
    print(f"Operating speed         : {OPERATING_SPEED_RPM:.0f} rpm "
          f"(r = {OMEGA_SERVICE/wn:.2f})")
    # --- Critical speeds in rpm ---------------------------------------
    RAD_S_TO_RPM = 60.0 / (2 * np.pi)
    try:
        cs = rotor.run_critical_speed(num_modes=8)
        wn_attr = cs.wn
        if callable(wn_attr):
            crit_rpm = np.atleast_1d(wn_attr(frequency_units="RPM"))
        else:
            crit_rpm = np.atleast_1d(wn_attr) * RAD_S_TO_RPM
        label = "synchronous (gyroscopic) critical speeds"
    except Exception:
        crit_rpm = np.atleast_1d(rotor.run_modal(speed=0.0).wn[:8]) \
                   * RAD_S_TO_RPM
        label = "natural frequencies at 0 rpm (no gyroscopics)"

    print(f"\nCritical speeds -- {label}:")
    for i, rpm in enumerate(crit_rpm):
        print(f"  critical {i + 1}: {rpm:10.0f} rpm ({rpm / 60:8.1f} Hz)")
        if OPERATING_SPEED_RPM > 0:
            margin = (rpm - OPERATING_SPEED_RPM) / OPERATING_SPEED_RPM * 100
            print(f"               margin to operating speed: {margin:+.1f} %")
    
    K_eff, zeta = whirl_vs_frequency_ratio(wn)
    print(f"\nEffective 1-DOF stiffness K = M*wn^2 = {K_eff:.4e} N/m")
    print(f"Effective damping ratio zeta        = {zeta:.4f}")

    print("\nPermissible unbalance at service speed:")
    for G in G_GRADES:
        U = rotor.m * (G / 1000.0) / OMEGA_SERVICE
        print(f"  G{G:<4}: U_per = {U*1e6:8.4f} g.mm  "
              f"({U*1e6/2:.4f} g.mm per plane)")

    stiffness_study()

    print("\n" + "=" * 60)
    print("3D Mode Shapes")
    print("=" * 60)
    ModeShapeAnalysis.run()

if enableDashboard == True:
    print("\n" + "=" * 60)
    print("Dashboard")
    print("=" * 60)
    build_dashboard()


if __name__ == "__main__":
    main()