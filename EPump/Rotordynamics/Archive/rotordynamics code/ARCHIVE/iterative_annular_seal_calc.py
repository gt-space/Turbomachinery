from numpy import pi, abs

######################## INPUTS #############################
L = 0.01 # m | length of annular seal
radius_shaft = 16/1000 # m
radius_casing = 17/1000 # m
radius_impeller = 100/1000 # m
N = 1000 # rad/s | angular velocity of impeller
deltaP_pump = 395 # psi
#############################################################

### SEAL CONSTANTS
g = 9.81 # m/s^2 | gravity constant
rho = 1141 # kg/m^3 | density of liquid oxygen
xi = 0.5 # given constant by NPTEL course
c_r = radius_casing-radius_shaft # m | clearance between shaft & casing
R = (radius_casing + radius_shaft)/2 # m | radius of centerline of seal
nu = 1 # Dynamic viscosity of liquid oxygen | TODO: find this value
deltaH_impeller = deltaP_pump*6895/(g*rho) # m

# Calculate delta-P across seal using eqn 3.7.3 from Gulich
u_2 = radius_impeller*N # m/s | circumferential velocity at impeller outlet
Re_u2 = u_2*radius_impeller/nu
k = 0.9*((Re_u2**0.3)*(c_r*2*radius_shaft/((2*radius_impeller)**2))*(c_r/L)**0.5)**0.087 # rotation factor | Gulich eqn. 3.7.2
deltaH_sp = deltaH_impeller - (k**2)*((u_2**2)/(2*g))*(1-(radius_shaft/radius_impeller)**2) # delta-Head over seal | Gulich eqn. 3.7.3
deltaP_sp = deltaH_sp*rho*g
print(deltaP_sp)
u_c_sp = (radius_impeller/radius_shaft)*u_2
Re_c = c_r*u_c_sp/nu

### ITERATIVE ANALYSIS
# performed using algorithm described in NPTEL course fig. 3.38
sigma = 0
e = 1
while e >= 1E-4:
    sigma += 0.001
    u_a = ((2*deltaP_sp)/(rho*(1+xi+2*sigma)))**(1/2)
    Re_a = c_r*u_a/nu
    b = Re_a/Re_c
    lambd_prime = sigma*c_r/L
    lambd = 0.066*(Re_a**(-1/4))*(1+(2*b)**-2)**0.375
    e = abs(lambd_prime - lambd)

def ann_seal_calc(L, radius_shaft, radius_impeller, radius_casing, deltaP_pump):
    ### SEAL CONSTANTS
    g = 9.81 # m/s^2 | gravity constant
    rho = 1141 # kg/m^3 | density of liquid oxygen
    xi = 0.5 # given constant by NPTEL course
    c_r = radius_casing-radius_shaft # m | clearance between shaft & casing
    R = (radius_casing + radius_shaft)/2 # m | radius of centerline of seal
    nu = 1 # Dynamic viscosity of liquid oxygen | TODO: find this value
    deltaH_impeller = deltaP_pump*6895/(g*rho) # m

    # Calculate delta-P across seal using eqn 3.7.3 from Gulich
    u_2 = radius_impeller*N # m/s | circumferential velocity at impeller outlet
    Re_u2 = u_2*radius_impeller/nu
    k = 0.9*((Re_u2**0.3)*(c_r*2*radius_shaft/((2*radius_impeller)**2))*(c_r/L)**0.5)**0.087 # rotation factor | Gulich eqn. 3.7.2
    deltaH_sp = deltaH_impeller - (k**2)*((u_2**2)/(2*g))*(1-(radius_shaft/radius_impeller)**2) # delta-Head over seal | Gulich eqn. 3.7.3
    deltaP_sp = deltaH_sp*rho*g
    print(deltaP_sp)
    u_c_sp = (radius_impeller/radius_shaft)*u_2
    Re_c = c_r*u_c_sp/nu

    ### ITERATIVE ANALYSIS
    # performed using algorithm described in NPTEL course fig. 3.38
    sigma = 0
    e = 1
    while e >= 1E-4:
        sigma += 0.001
        u_a = ((2*deltaP_sp)/(rho*(1+xi+2*sigma)))**(1/2)
        Re_a = c_r*u_a/nu
        b = Re_a/Re_c
        lambd_prime = sigma*c_r/L
        lambd = 0.066*(Re_a**(-1/4))*(1+(2*b)**-2)**0.375
        e = abs(lambd_prime - lambd)