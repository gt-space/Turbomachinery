from pathlib import Path
import ross as rs
import numpy as np
import plotly.io as pio
import math
pio.renderers.default = "notebook"

# # Shaft Inputs
elements_1 = 2
elements_2 = 2
elements_3 = 2
elements_4 = 2
elements_5 = 2

steel = rs.Material(name="Steel", rho=7810, E=211e9, G_s=81.2e9) # from E and G_s
l_first = np.ones(elements_1)
l_first = np.array(l_first)*((5.25e-3)/len(l_first))
od_first = np.ones(len(l_first))*8e-3

l_second = np.ones(elements_2)
l_second = np.array(l_second)*((18.45e-3)/len(l_second))
od_second = np.ones(len(l_second))*15.1e-3

l_third= np.ones(elements_3)
l_third = np.array(l_third)*((18.5e-3)/len(l_third))
od_third = np.ones(len(l_third))*18e-3

l_fourth = np.ones(elements_4)
l_fourth = np.array(l_fourth)*((18.54e-3)/len(l_fourth))
od_fourth = np.ones(len(l_fourth))*15.1e-3

l_fifth = np.ones(elements_5)
l_fifth = np.array(l_fifth)*((42.46e-3)/len(l_fifth))
od_fifth = np.ones(len(l_fifth))*10e-3

l = np.concatenate((l_first,l_second,l_third,l_fourth,l_fifth))
o_d = np.concatenate((od_first,od_second,od_third,od_fourth,od_fifth))
i_d = np.array(np.zeros(len(l))) #Internal diameter of the shaft

n = len(l)
N = 25000 #(RPM) shaft speed

#Bearing inputs

n_b = 13 #number of steel balls in the ball bearings
d = 4.762e-3 # (meters) diameter of the steel balls
f_s =  200 + 1.81*9.81 # (Newtons) preload + assembly weight
alpha = 15 #(deg) ball contact angle
kb = 13e6 #N^(2/3) m^(-4/3)
n_1 = len(l_first) + len(l_second)  #Bearing 1 nodal location
n_2 = len(l_first) + len(l_second) + len(l_third) #Bearing 2 nodal location
#Get the radial stiffness constant (Friswell, Dynamics of Rotating Machinery, page 183)
k_xx = kb*(n_b**(2/3))*(d**(1/3))*(f_s**(1/3))*(np.cos(math.radians(alpha))**(5/3))

# Disk Inputs
disk = [rs.DiskElement(
     n=0, #disk nodal location
     m=47.43e-3, # (kg) disk mass
     Ip=12765.3e-9, #(kg*m2) Polar moment inertia
     Id=6619.35e-9, # #(kg*m2) diametral moment inertia
     tag="Disk"
 )]

#Unbalance
G = 100 #Considering a G40 balancing grade
e_mm = G/(N*2*np.pi/60)
e = e_mm*10**(-3)
n_unb = 0
magnitude = 0.195*e #(kg*m)
phase = 0
frequency_range=np.linspace(0*2*np.pi/60,25000*2*np.pi/60,100)



# Cylindrical shaft elements

i_vec = np.arange(0, n)

shaft_elements = []
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


#Bearings


bearing1 = rs.BearingElement(n=n_1,kxx = k_xx,cxx=0)
bearing2 = rs.BearingElement(n=n_2,kxx = k_xx,cxx=0)
#bearing3 = rs.BearingElement(n=5,kxx = 5.4e6,cxx=0)


bearings = [bearing1,bearing2]


mu = disk[0].Id/(disk[0].m * l[0])



#Initialize the rotor
rotor1 = rs.Rotor(shaft_elements, disk, bearings)


#Modal response of the rotor
modal = rotor1.run_modal(25000*2*np.pi/60)

#Forced response
# G = 40 #Considering a G40 balancing grade
# e_mm = G/(N*2*np.pi/60)
# e = e_mm*10**(-3)
# F_unb =  0.195*0.0153e-3


results_unbalance = rotor1.run_unbalance_response(n_unb, magnitude, phase, frequency_range)



#get the campbell diagram

samples = 200
speed_range = np.linspace(0, 30000*2*np.pi/60, samples)
campbell = rotor1.run_campbell(speed_range)

#Plots and results

#plot the rotor
rotor1.plot_rotor().show(renderer='browser')
#Campbell diagram
campbell.plot(frequency_units="rpm").show(renderer='browser')
#Modal responde for the first mode
modal.plot_mode_3d(1).show(renderer='browser')
#Forced response
results_unbalance.plot_deflected_shape(speed=25000*2*np.pi/60).show(renderer='browser')

print("="*36)
print(f"Young's Modulus: {steel.E}")
print(f"Shear Modulus:    {steel.G_s}")
print("="*36)
print("Rotor total mass = ", np.round(rotor1.m, 2))
print("Rotor center of gravity =", np.round(rotor1.CG, 2))
print("="*36)
print("mu = ",mu)
print("="*36)
print("Undamped natural frequencies:\n", modal.wn)