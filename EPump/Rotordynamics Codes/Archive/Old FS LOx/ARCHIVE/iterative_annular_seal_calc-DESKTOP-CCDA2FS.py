from numpy import pi

######################## INPUTS #############################
c_r = 1/1000 # m | clearance between shaft & casing
L = 0.1 # m | length of annular seal
radius_shaft = 16/1000 # m
radius_casing = 17/1000 # m
radius_impeller = 100/1000 # m
N = 1000 # rad/s | angular velocity of impeller
deltaH_impeller = 100 # Pa 
#############################################################

### SEAL CONSTANTS
g = 9.81 # m/s^2 | gravity constant
rho = 1141 # kg/m^3 | density of liquid oxygen
xi = 0.5 # given constant by NPTEL course 
mu = 0.6 # friction coeff between shaft & fluid (steel & liquid oxygen) # TODO: update with actual value
R = (radius_casing + radius_shaft)/2 # m | radius of centerline of seal
nu = 0.5 # poisson ratio of liquid oxygen | incompressible fluids nu=~0.5

# Calculate delta-P across seal using eqn 3.7.3 from Gulich
u_2 = 2*pi*radius_impeller*N # m/s | circumferential velocity at impeller outlet
Re_u2 = u_2*radius_impeller/nu
k = 0.9*((Re_u2**0.3)*(c_r*2*radius_shaft/((2*radius_casing)**2))*(c_r/L)**0.5)**0.087 # rotation factor | Gulich eqn. 3.7.2
deltaH_sp = deltaH_impeller - (k**2)*((u_2**2)/(2*g))*(1-(radius_shaft/radius_impeller)**2) # delta-Head over seal | Gulich eqn. 3.7.3
u_c_sp = (radius_impeller/radius_shaft)*u_2