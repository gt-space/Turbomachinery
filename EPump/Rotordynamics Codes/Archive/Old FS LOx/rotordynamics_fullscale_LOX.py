from pathlib import Path
import ross as rs
import numpy as np
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as pgo
import math
from UTILS.annular_seal_calc import childs_seal_calc
import matplotlib.pyplot as plt
pio.renderers.default = "notebook"

N = 20000 # RPM | max shaft speed

#### Setting up shaft & elements ####
## Shaft Elements
print("Setting up shaft elements...")
# INPUTS
num_shaft_elements = 7
num_nodes_per_element = 2
shaft_elem_ODs = [0.67, 0.67, 1.1, 0.94, 0.79, 0.67, 0.39] # in
shaft_elem_lens = [0.66, 0.67, 0.12, 0.98, 0.98, 0.59, 0.65] # in

# convert shaft ODs & lengths from in to m
shaft_elem_lens = [len/39.37 for len in shaft_elem_lens]
shaft_elem_ODs = [OD/39.37 for OD in shaft_elem_ODs]

steel = rs.Material(name="Steel", rho=7810, E=211e9, G_s=81.2e9) # from E and G_s

l = np.array([])
o_d = np.array([])
i_d = np.array([])
for i in range(num_shaft_elements):
    l = np.append(l, np.ones(num_nodes_per_element)*((shaft_elem_lens[i])/num_nodes_per_element))
    o_d = np.append(o_d, np.ones(num_nodes_per_element)*(shaft_elem_ODs[i]))

l = np.array(l)
o_d = np.array(o_d)
i_d = np.array(np.zeros(len(l))) #Internal diameter of the shaft

shaft_elements = []
i_vec = np.arange(0, len(l))
for i in i_vec:
   shaft_elements.append(
        rs.ShaftElement(
            L=l[i],
            idl=i_d[i],
            odl=o_d[i],
            material=steel,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
        )
   )

## Ball Bearings
print("Setting up ball bearing elements...")
# Inputs
# for n_b and d, currently using 7202 bearing dimensions because I dont know the dimensions of 7203
n_b = 10 # number of balls in bearing
d = 6.747 * 1e-3 # (m) diameter of steel balls
f_s = 185 + 1550 # (N) preload + assembly weight
alpha = 15 # (deg) ball contact angle
kb = 13e6 # (N^(2/3)/m^(4/3)) constant given in Friswell pg. 183
n_1 = 2*num_nodes_per_element # bearing 1 placed at end of second element
n_2 = 5*num_nodes_per_element # bearing 2 placed at beginning of sixth element

# calculate bearing stiffness using eqn from Friswell pg. 183
k_xx = kb*(n_b**(2/3))*(d**(1/3))*(f_s**(1/3))*(np.cos(math.radians(alpha))**(5/3))

bearing1 = rs.BearingElement(n=n_1, kxx=k_xx, cxx=0)
bearing2 = rs.BearingElement(n=n_2, kxx=k_xx, cxx=0)

'''
## Annular Seal - deprecated. Now we should calculate annular seal constants at each shaft speed
print("Setting up annular seal element...")
r_seal = 0.8/39.37 # m
r_impeller = 1.33/39.37 # m
L = 0.18/39.37 # m
c_r = 1e-3 # m
omega = N*2*np.pi/60 # rad/s
rho = 1141 # kg/m^3
mu = 6.94e-6 # Pa*s at 90.1 degrees K | from https://www.engineeringtoolbox.com/oxygen-O2-dynamic-kinematic-viscosity-temperature-pressure-d_2081.html?vA=30&degree=K#
deltaP_impeller = 545-150 # psi, from requirements
v_0 = 0
[k_d, k_c, c_d, c_c, m] = childs_seal_calc(r_seal, r_impeller, L, c_r, omega, rho, mu, deltaP_impeller, v_0)
annular_seal = rs.BearingElement(n=0, kxx=k_d, kxy=k_c, kyx=k_c, kyy=k_d, cxx=c_d, cxy=c_c, cyx=c_c, cyy=c_d, m=m)
'''

## Rotor Disk
print("Setting up disk element...")
# Disk Inputs
disk_node_location = 1
disk_mass = 0.56 # lbs
disk_polar_inertia = 0.49 # lb*in^2
disk_diametral_inertia = 0.26 # lb*in^2
# unit conversion
disk_mass = disk_mass/2.205 # lb -> kg
disk_polar_inertia = disk_polar_inertia/(2.205*1550*1550) # lb*in^2 -> kg*m^2
disk_diametral_inertia = disk_diametral_inertia/(2.205*1550*1550) # lb*in^2 -> kg*m^2
disk = [rs.DiskElement(
     n=disk_node_location, #disk nodal location
     m=disk_mass, # (kg) disk mass
     Ip=disk_polar_inertia, #(kg*m2) Polar moment inertia
     Id=disk_diametral_inertia, # #(kg*m2) diametral moment inertia
     tag="Disk"
 )]

mu = disk[0].Id/(disk[0].m * l[0])

## Unbalance
G = 100 #Considering a G40 balancing grade
e_mm = G/(N*2*np.pi/60)
e = e_mm*10**(-3)
n_unb = 0
magnitude = 0.195*e #(kg*m)
phase = 0
frequency_range = np.linspace(0,N*2*np.pi/60,20)

#### Calculating dampened natural frequencies using modal analysis ####
print("Starting dampened natural frequencies calculations...")
# Annular seal inputs
r_seal = 0.8/39.37 # m
r_impeller = 1.33/39.37 # m
L = 0.18/39.37 # m
c_r = 1e-3 # m
rho = 1141 # kg/m^3
mu = 6.94e-6 # Pa*s at 90.1 degrees K | from https://www.engineeringtoolbox.com/oxygen-O2-dynamic-kinematic-viscosity-temperature-pressure-d_2081.html?vA=30&degree=K#
max_deltaP_i = 545-150 # psi, from requirements
v_0 = 0 # swirl, i dont know what this "should" be
speed_range=np.linspace(0.1,N,20) # rpm
damp_freqs = np.ndarray(shape=(len(speed_range),6), dtype=np.double)

# use modal analysis to get number of dampened natural frequencies
print("Getting number of dampened natural frequencies...")
omega = speed_range[1]*2*np.pi/60
deltaP = max_deltaP_i*(speed_range[1]/N)**2
[k_d, k_c, c_d, c_c, m] = childs_seal_calc(r_seal, r_impeller, L, c_r, omega, rho, mu, deltaP, v_0)
annular_seal = rs.BearingElement(n=0, kxx=k_d, kxy=k_c, kyx=k_c, kyy=k_d, cxx=c_d, cxy=c_c, cyx=c_c, cyy=c_d, m=m)
bearings = [bearing1, bearing2, annular_seal]
rotor = rs.Rotor(shaft_elements, disk, bearings)
num_freqs = len(rotor.run_modal(omega).wd)

# calculate dampened frequencies
damp_freqs = np.ndarray(shape=(len(speed_range),6),dtype=np.double)
for i in range(len(speed_range)):
    omega = speed_range[i]*2*np.pi/60
    print(f"Calculating dampened natural frequencies for w={omega} rads/s")
    deltaP = max_deltaP_i*(speed_range[i]/N)**2
    [k_d, k_c, c_d, c_c, m] = childs_seal_calc(r_seal, r_impeller, L, c_r, omega, rho, mu, deltaP, v_0)
    annular_seal = rs.BearingElement(n=0, kxx=k_d, kxy=k_c, kyx=k_c, kyy=k_d, cxx=c_d, cxy=c_c, cyx=c_c, cyy=c_d, m=m)
    bearings = [bearing1, bearing2, annular_seal]
    rotor = rs.Rotor(shaft_elements, disk, bearings)
    modal = rotor.run_modal(omega)
    damp_freqs[i,:] = modal.wd

# NOTE: rotor after this loop will be the configuration at maximum speed (N)

print("Getting system responses & plots...")
#Modal response of the rotor
modal = rotor.run_modal(N*2*np.pi/60)

# NOTE: these unbalance results will only be valid for maximum frequency (N) since that is what the rotor annular seal is set up for
results_unbalance = rotor.run_unbalance_response(n_unb, magnitude, phase, frequency_range)

#### Plots and results ####
#plot the rotor
rotor.plot_rotor().show(renderer='browser')

# plot campbell using manual function
# fig = px.line(x=speed_range, y=[speed_range, *[damp_freqs[:,i]*30/np.pi for i in range(len(damp_freqs[0,:]))]], title='Campbell Dampened Frequencies using Manual Iteration')
fig = px.line(title='Campbell Dampened Frequencies using Manual Iteration')
fig['layout']['xaxis']['title'] = {'text': 'Rotor Speed (rpm)'}
fig['layout']['yaxis']['title'] = {'text': 'Natural Frequencies (rpm)'}
fig.update_layout(xaxis_range=[0,N])
fig.update_layout(yaxis_range=[0,500000])
fig.add_trace(pgo.Line(x=speed_range, y=speed_range, mode='lines', line=dict(color='blue', dash='dash')))
for i in range(len(damp_freqs[0,:])):
    if i % 2 == 0:
        symbol = "triangle-down"
    else:
        symbol = "triangle-up"
    fig.add_trace(pgo.Scatter(x=speed_range, y=damp_freqs[:,i]*30/np.pi, mode='markers', marker_size=5, marker_symbol=symbol, marker_color="red"))
print(fig)
# formatting
fig.show(renderer='browser')

# plot campbell using automatic function
campbell = rotor.run_campbell(speed_range*np.pi/30)
campbell.plot(frequency_units="rpm", title="Campbell Dampened Frequencies using Built-in ROSS Function").show(renderer='browser')
#Modal responde for the first mode
modal.plot_mode_3d(1).show(renderer='browser')
#Forced response
results_unbalance.plot_deflected_shape(speed=N*2*np.pi/60).show(renderer='browser')

print("="*36)
print(f"Young's Modulus: {steel.E}")
print(f"Shear Modulus:    {steel.G_s}")
print("="*36)
print("Rotor total mass = ", np.round(rotor.m, 2))
print("Rotor center of gravity =", np.round(rotor.CG, 2))
print("="*36)
print("mu = ",mu)
print("="*36)
print("Undamped natural frequencies:\n", modal.wn)
