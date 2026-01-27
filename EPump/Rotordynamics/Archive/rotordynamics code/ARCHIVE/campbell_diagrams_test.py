import ross as rs
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

## Using Example 26 - Isotropic System example from ROSS manual (https://ross.readthedocs.io/en/stable/user_guide/example_26.html)

Q_ = rs.Q_

steel = rs.Material("Steel", E=211e9, G_s=81.1e9, rho=7810)

# create rotor
N = 6
L = 1.5/N
idl = 0
odl = 0.050   # shaft diameter
shaft = [rs.ShaftElement(L=L, idl=idl, odl=odl, material=steel) for i in range(N)]
bearings = [
    rs.BearingElement(n=0, kxx=1e6, kyy=1e6, cxx=100, cyy=100),
    rs.BearingElement(n=6, kxx=1e6, kyy=1e6, cxx=100, cyy=100),
]
disks = [
    rs.DiskElement.from_geometry(
        n=N/3, material=steel, width=0.070, i_d=odl, o_d=0.280, scale_factor="mass"
    ),
    rs.DiskElement.from_geometry(
        n=2*N/3, material=steel, width=0.070, i_d=odl, o_d=0.350, scale_factor="mass"
    )
]
rotor = rs.Rotor(shaft_elements=shaft, disk_elements=disks, bearing_elements=bearings)
rotor.plot_rotor().show()

speed_range = np.linspace(0, 6500, 65)*2*np.pi/60 # rad/s

# plot campbell automatically using run_campbell
campbell = rotor.run_campbell(speed_range)
print(f"Campbell results: {campbell.wd[1,:]}")
fig = campbell.plot()
fig.show()

# plot campbell using modal analysis
damp_freqs = np.ndarray(shape=(len(speed_range),6), dtype=np.double)
for i in range(len(speed_range)):
    speed = speed_range[i]
    modal = rotor.run_modal(speed)
    damp_freqs[i,:] = modal.wd*30/np.pi

fig, ax = plt.subplots()
for i in range(len(damp_freqs[0,:])):
    ax.scatter(speed_range*30/np.pi, damp_freqs[:,i])
plt.show()


