"""
ModeShapeAnalysis.py
========================================================================
Modal analysis (natural frequencies + 2D mode shapes) for the V1 rotor,
computed on a temporarily REFINED shaft mesh for smoother mode-shape
curves. This does not affect RotorModelV1.py or RotorModelV2.py's own
(coarser, 15-element) rotor -- it builds its own rotor here by splitting
each of V1's shaft elements into REFINEMENT pieces.

Run:  python ModeShapeAnalysis.py    (writes mode_shapes.html)
      or import run() from RotorModelV2.py.
"""

import os

import numpy as np
import plotly.graph_objects as go
import ross as rs

from buildRotor import (
    MATERIAL,
    shaft_table,
    disk1, disk2,
    DOWNLOADS_DIR,
    REFINEMENT, NUM_MODES,
    bearing1, bearing2
)

# ----------------------------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Refined shaft: every element split into REFINEMENT equal pieces, so all
# of V1's node numbers (disk/bearing seats) land at (old_node * REFINEMENT)
# in this mesh.
# ----------------------------------------------------------------------
shaft_elements = [
    rs.ShaftElement(
        L=L / REFINEMENT,
        idl=idl,
        odl=odl,
        material=MATERIAL,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for (L, idl, odl) in shaft_table
    for _ in range(REFINEMENT)
]

disk1.n = disk1.n * REFINEMENT
disk2.n = disk2.n * REFINEMENT

disks = [disk1, disk2]

bearing1.n = bearing1.n * REFINEMENT
bearing2.n = bearing2.n * REFINEMENT

bearings = [
    bearing1,
    bearing2
]

rotor = rs.Rotor(shaft_elements, disks, bearings)


def run(num_modes=NUM_MODES, filename="mode_shapes.html"):
    """Run modal analysis on the refined mesh and save a combined 2D mode
    shape plot (all modes overlaid, one colour/legend entry each), plus
    one 3D mode shape plot per mode (kept separate -- overlaying several
    3D orbit/shaft scenes in one plot is unreadable)."""
    print(f"Refined mesh: {len(shaft_elements)} elements "
          f"({len(shaft_table)} x {REFINEMENT})")

    modal = rotor.run_modal(speed=0.0, num_modes=num_modes)
    print("Undamped natural frequencies at 0 speed:")
    for i, wn in enumerate(modal.wn):
        print(f"  mode {i:>2}: {wn:10.2f} rad/s = {wn / (2 * np.pi):8.2f} Hz")

    palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]
    fig_modes = go.Figure()
    for i, wn in enumerate(modal.wn):
        n_before = len(fig_modes.data)
        fig_modes = modal.plot_mode_2d(mode=i, fig=fig_modes, frequency_type="wn")
        color = palette[i % len(palette)]
        for j, trace in enumerate(fig_modes.data[n_before:]):
            trace.update(
                line=dict(color=color),
                legendgroup=f"mode{i}",
                showlegend=(j == 0),
                name=f"Mode {i} ({wn:.1f} rad/s)",
            )
    fig_modes.update_layout(
        title=f"2D mode shapes (undamped, 0 speed, "
              f"{len(shaft_elements)}-element mesh)",
        xaxis_title="Shaft position (m)",
        yaxis_title="Relative displacement (major axis)",
    )
    os.makedirs(DOWNLOADS_DIR, exist_ok=True)

    out_path = os.path.join(DOWNLOADS_DIR, filename)
    fig_modes.write_html(out_path)
    print(f"Saved combined 2D mode shape plot -> {out_path}")

    for i, wn in enumerate(modal.wn):
        fig_3d = modal.plot_mode_3d(mode=i, frequency_type="wn")
        fig_3d.write_html(os.path.join(DOWNLOADS_DIR, f"mode_shape_3d_{i}.html"))
    print(f"Saved {len(modal.wn)} 3D mode shape plots -> {DOWNLOADS_DIR}/mode_shape_3d_<i>.html")

    return modal


if __name__ == "__main__":
    run()
