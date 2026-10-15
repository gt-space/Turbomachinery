# Calcs for bearings, radial and axial stiffnesses.
# Source: Dynamics of Rotating Machinery, Friswell, Penny, Garvey & Lees, Chapter 5, Section 5.5, pp. 183, eq. (5.90).
import numpy as np

# Function to calculate bearing stiffness
def radial_Stiffness_Calcs(radial_Force, ball_Diameter, num_Of_Balls, contact_Angle):
    radial_Stiffness =(
        1.3e7 
        * (num_Of_Balls ** (2/3)) 
        * (ball_Diameter ** (1/3)) 
        * (radial_Force ** (1/3)) 
        * ((np.cos(np.deg2rad(contact_Angle))) ** (5/3))
        )
    return radial_Stiffness

# Derived from source equation, 
def axial_Stiffness_Calcs(preload, thrust, ball_Diameter, num_Of_Balls, contact_Angle):
    total_Axial_Force = preload + thrust
    axial_Stiffness = (
        (1.5 / 4.36e-8)
        * (num_Of_Balls ** (2/3))
        * (ball_Diameter ** (1/3))
        * (total_Axial_Force ** (1/3))
        * (np.sin(np.deg2rad(contact_Angle)) ** (5/3))
        )
    return axial_Stiffness