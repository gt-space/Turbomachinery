"""
From mass calculates disk properites from CAD (FPump V9) impeller and inducer, from geometry calculates disks as hollow cylinders.
"""

import os

import numpy as np
import ross as rs
import bearingCalcs as bc

# Plots are saved to the Downloads folder of whoever runs this script.
DOWNLOADS_DIR = os.path.join(os.path.expanduser("~"), "Downloads")

# ----------------------------------------------------------------------
# Asks if you want to enable a dashboard at the end of the analysis.
# This compiles all the individual HTML plots into one  page.
# ----------------------------------------------------------------------
enableDashboard = True  # <-- set to False to skip dashboard creation

# ----------------------------------------------------------------------
# Parameters for rotor geometry and operation for analyzeRotor.py and ModeShapeAnalysis.py
# ----------------------------------------------------------------------
MATERIAL = rs.materials.steel  # <-- rho = 7810 kg/m^3, E = 211 GPa


OPERATING_SPEED_RPM = 30000       # <-- EDIT: service speed [rpm]
BEARING_CXX = 20.0                # <-- EDIT: bearing damping [N.s/m].
                                  #     Must be > 0 or the resonance
                                  #     peak is infinite (undamped).
G_GRADES = [0.4, 1.0, 2.5, 6.3]   # ISO balance grades [mm/s]
MEASURED_CRITICAL_HZ = None       # <-- EDIT: measured resonance [Hz]
                                  #     to back-solve bearing stiffness
BEARING_NODES = (7, 12)
DISK_NODES = (2, 4)
RESPONSE_NODE = 4                 # where whirl amplitude is reported
                                  # (node 4 = big disk; try 7 for brg).
                                  # x-translation DOF, resolved via
                                  # rotor.number_dof (6 here, not 4).
R_MAX = 2.5                       # sweep up to 2.5 x first critical

OMEGA_SERVICE = OPERATING_SPEED_RPM * 2 * np.pi / 60.0   # rad/s

DISK_METHOD = "from_mass"  # "from_geometry" or "from_mass"

SHAFT_OD_AT_DISKS = 0.0100764  # shaft OD where both disk collars are seated (precise from CAD)

# Don't take the comments next to the table too seriously, I lost track after a few design changes
#            L,            idl,         odl
shaft_table = [
    (0.0053122,   0.0034544,   0.010),      # 0  hollow lead-in
    (0.00494535,  0.0,         0.010),      # 1  -> node 2 (disk 1 centre)
    (0.01025755,  0.0,         0.010),      # 2  -> node 3 (disk boundary)
    (0.0047857,   0.0,         0.010),      # 3  -> node 4 (disk 2 centre)
    (0.0047855,   0.0,         0.010),      # 4  -> node 5 (step to 12 mm)
    (0.0066269,   0.0,         0.012),      # 5  -> node 6 (step to 15 mm)
    (0.0055,      0.0,         0.017),      # 6  -> node 7 (bearing 1 centre)
    (0.0055,      0.0,         0.017),      # 7  -> node 8 (step to 20 mm)
    (0.00635,     0.0,         0.02),       # 8  -> node 9 (first third of 20 mm)
    (0.00635,     0.0,         0.02),       # 9  -> node 10 (second third of 20 mm)
    (0.00635,     0.0,         0.02),       # 10  -> node 11 (final third of 20 mm)
    (0.0055,      0.0,         0.017),      # 11  -> node 12 (bearing 2 centre)
    (0.0055,      0.0,         0.017),      # 12 -> node 13 (step to 10 mm)
    (0.0055684,   0.0,         0.01),       # 13 -> node 14 (free end, 0.0889 m)
    (0.0055684,   0.0,         0.01),       # 14 -> node 15 (free end, 0.0889 m))
]

# Defines disk nodes for modal analysis and Campbell diagram.  The disk nodes are defined in the shaft_table above.
# Disk 1: 0 -> 0.0205151 m, OD = 0.0287528 m  -> node 2
D1_NODE, D1_WIDTH, D1_OD = 2, 0.0205151, 0.0287528
# Disk 2: 0.0205151 -> 0.0300865 m, OD = 0.0506984 m  -> node 4
D2_NODE, D2_WIDTH, D2_OD = 4, (0.0300865 - 0.0205151), 0.0506984

# ----------------------------------------------------------------------
# Parameters for increasing mesh quality for mode shape analysis (ModeShapeAnalysis.py)
# ----------------------------------------------------------------------
REFINEMENT = 2   # <-- split each of V1's 15 shaft elements into this many
                 #     pieces (15 * REFINEMENT = number of elements here)
NUM_MODES = 12   # <-- passed to run_modal(num_modes=...); ARPACK returns
                 #     conjugate pairs, so you get roughly NUM_MODES / 2
                 #     real mode shapes back.




shaft_elements = [
    rs.ShaftElement(
        L=L,
        idl=idl,
        odl=odl,
        material=MATERIAL,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for (L, idl, odl) in shaft_table
]

if DISK_METHOD == "from_geometry":
    disk1 = rs.DiskElement.from_geometry(
        n=D1_NODE,
        material=MATERIAL,
        width=D1_WIDTH,
        i_d=SHAFT_OD_AT_DISKS,
        o_d=D1_OD,
        tag="Disk_1",
    )
    disk2 = rs.DiskElement.from_geometry(
        n=D2_NODE,
        material=MATERIAL,
        width=D2_WIDTH,
    i_d=SHAFT_OD_AT_DISKS,
    o_d=D2_OD,
    tag="Disk_2",
    )
elif DISK_METHOD == "from_mass":
    disk1 = rs.DiskElement(n=2, m=0.02836, Id=1.61671e-6, Ip=1.49207e-6, tag="Inducer")
    disk2 = rs.DiskElement(n=4, m=0.06801, Id=1.0123405e-5, Ip=1.915531e-5, tag="Impeller")

disks = [disk1, disk2]

# ----------------------------------------------------------------------
# Bearing elements
# ----------------------------------------------------------------------
K_RADIAL = bc.radial_Stiffness_Calcs(radial_Force=182, ball_Diameter=0.0056, num_Of_Balls=13, contact_Angle=15)
K_AXIAL = bc.axial_Stiffness_Calcs(preload=300, thrust=682.4, ball_Diameter=0.0056, num_Of_Balls=13, contact_Angle=15)


print(f"Calculated bearing stiffnesses: K_radial = {K_RADIAL:.3e} N/m, K_axial = {K_AXIAL:.3e} N/m")
BEARING_NODES = (7, 12)

bearing1 = rs.BearingElement(
    n=BEARING_NODES[0], kxx=K_RADIAL, kyy=K_RADIAL, kzz=K_AXIAL, cxx=0.0, tag="Bearing_1"
)
bearing2 = rs.BearingElement(
    n=BEARING_NODES[1], kxx=K_RADIAL, kyy=K_RADIAL, kzz=K_AXIAL, cxx=0.0, tag="Bearing_2"
)
bearings = [bearing1, bearing2]

# ----------------------------------------------------------------------
# Assemble the rotor
# ----------------------------------------------------------------------
rotor = rs.Rotor(shaft_elements, disks, bearings)