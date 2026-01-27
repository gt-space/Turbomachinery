from numpy import abs, pi, exp

def NPTEL_seal_calc(L, r_shaft, r_casing, r_impeller, deltaP_pump, N_range):
    '''
    Calculates the stiffness and dampening coefficients of the annular seal
    between the impeller and the inlet casing using the algorithm described in
    NPTEL Figure 3.38 (https://archive.nptel.ac.in/courses/112/103/112103024/#)

    INPUTS
    L_mm - length of seal in m
    r_shaft - radius of shaft at seal location in m
    r_casing - radius of casing at seal location in m
    r_impeller - radius of impeller in m
    deltaP_pump - deltaP of pump in psi
    N_range - range of speeds we want to analyze in rpm

    OUTPUTS
    k_d - direct stiffness coefficient
    k_c - cross-coupled stiffness coefficient
    c_d - direct damping coefficient
    c_c - cross-coupled damping coefficient
    m_d - direct inertia coefficient
    '''

    ### SEAL CONSTANTS
    g = 9.81 # m/s^2
    rho = 995 # kg/m^3
    xi = 0.5 # given constant by NPTEL course
    c_r = r_casing - r_shaft # m
    R = (r_casing + r_shaft) / 2 # m
    nu = 8e-7 # kinematic viscosity of liquid oxygen

    # create empty lists for results to be appended to
    k_d = []
    k_c = []
    c_d = []
    c_c = []
    m_d = []

    for i in range(len(N_range)):
        N = N_range[i]*2*pi/60
        print("calculating for N = " + str(N))

        # calculate pressure drop across seal & circumferential velocity at seal, commented out for now
        # [deltaP_sp, u_c_sp] = gulich_seal_calcs(r_shaft, r_impeller, L, c_r, N, nu, rho, deltaP_pump)

        # test values from NPTEL, delete later after verifying results
        deltaP_sp = (2-1)*100000 # Pa
        R = 22/1000
        L = 11/1000
        c_r = 0.2/1000
        u_c_sp = N*R

        Re_c = c_r*u_c_sp/nu


        ### ITERATIVE ANALYSIS
        sigma = 0
        e = 1
        while e >= 1E-4:
            sigma += 0.001
            u_a = ((2*deltaP_sp)/(rho*(1+xi+2*sigma)))**(1/2)
            Re_a = c_r*u_a/nu
            b = Re_c/Re_a
            lambd_prime = sigma*c_r/L
            lambd = 0.066*(Re_a**(-1/4))*(1+(2*b)**-2)**0.375
            e = abs(lambd_prime - lambd)


        A = (pi*sigma)/(1+xi+2*sigma)
        B = (1+7*b**2)/(1+4*b**2)
        E = (1+xi)/(2*(1+xi+B*sigma))
        V = (2*deltaP_sp/(rho*(1+xi+2*sigma)))**(1/2)
        T = L/V
        a_0 = 2.5*A*E
        a_1 = 2*A*((E/sigma)+(B/2)*(E+(1/6)))
        a_2 = (A/sigma)*(E+(1/6))
        k_star = deltaP_sp*L*R/c_r
        c_star = k_star*T
        m_star = k_star*T**2

        k_d.append((a_0-0.25*a_2*(N*T)**2)*k_star)
        k_c.append(0.5*a_1*N*T*k_star)
        c_d.append(a_1*c_star)
        c_c.append(a_2*N*T*c_star)
        m_d.append(a_2*m_star)

    return [k_d, k_c, c_d, c_c, m_d]

def childs_seal_calc(r_seal, r_impeller, L, c_r, omega, rho, mu, deltaP_impeller, v_0):
    '''
    Calculate the dynamic coefficients of the seal using the equations described in Childs' "Dynamic Analysis of Turbulent Annular Seals Based on Hirs'Lubrication Equation" (1983)
    
    INPUTS:
    r_seal - radius of the shaft at seal in m
    r_impeller - outer radius of the impeller in m
    L - length of the seal in m
    c_r - clearance of the seal in m
    omega - angular speed of shaft/impeller in rad/s
    rho - density of working fluid in kg/m^3
    mu - dynamic viscosity of working fluid in Pa*s 
    deltaP_impeller - pump pressure increase in psi
    v_0 - swirl at the entrace to the seal 

    RETURNS:
    k_d - direct stiffness coefficient (shown as K in Childs)
    k_c - cross-coupled stiffness coefficient (shown as k in Childs)
    c_d - direct damping coefficient (shown as C in childs)
    c_c - cross-coupled damping coefficient (shown as c in childs)
    m - inertia coefficient 
    '''

    # constants used for calculations, defined in Childs
    m_0 = -0.25
    n_0 = 0.066

    # print("STARTING SEAL CALCS...")

    # calculate pressure drop across seal & circumferential velocity at seal, commented out for now
    # print("Calculating seal constants using Gulich equations...")
    [deltaP_seal, u_c] = gulich_seal_calcs(r_seal, r_impeller, L, c_r, omega, mu/rho, rho, deltaP_impeller)
    Re_c = rho*u_c*c_r/mu
    '''
    # TEST VALUES USING EXAMPLE FROM CHILDS
    r_seal = (7.98E-2)/2
    L = 4.32E-2
    c_r = 1.397E-4
    omega = 37360*2*pi/60
    rho = 70.78
    mu = 1.160E-5
    deltaP_seal = 1.492E7
    xi = 0.1
    u_c = r_seal*omega
    Re_c = rho*u_c*c_r/mu
    v_0 = -0.5
    '''

    # typical value according to NPTEL
    xi = 0.5

    # iterative analysis to find approximate value of sigma
    sigma = 0
    e = 1
    # print("Performing iterative analysis to get sigma...")
    while e > 1E-4:
        sigma += 0.0001
        u_a = ((2*deltaP_seal)/(rho*(1+xi+2*sigma)))**(1/2)
        Re_a = rho*c_r*u_a/mu
        b = Re_a/Re_c
        lambd_prime = sigma*c_r/L
        lambd = n_0*(Re_a**m_0)*(1+(1/(4*b**2)))**((1+m_0)/2)
        e = abs(lambd_prime - lambd)
        # print(f"Sigma = {sigma} | error = {e}")
        if (e > 1.0):
            raise ValueError(f"Error is growing indefinitely! current speed = {omega}")

    # print("Sigma successfully calculated...")
    beta = 1/(1+4*(b**2))
    a = sigma*(1+beta*(1+m_0))
    B = 1+4*(b**2)*beta*(1+m_0)
    E = (1+xi)/(2*(1+xi+B*sigma))
    omegaT = L/(r_seal*b)
    T = L/u_a

    # calculate dimensionless stiffness coefficients (shown with tildas in Childs)
    k_d = ((2*sigma**2)/(1+xi+2*sigma))*(E*(1-m_0)-((omegaT**2)/(4*sigma))*(0.5*((1/6)+E)+(2*v_0/a)*((E+(a**-2))*(1-exp(-1*a))-(0.5+(1/a))*exp(-1*a))))
    k_c = (((sigma**2)*omegaT)/(1+xi+2*sigma))*((E/sigma)+(B/2)*((1/6)+E)+(2*v_0/a)*(E*B+((1/sigma)-(B/a))*((1-exp(-1*a))*(E+0.5+(1/a))-1)))
    c_d = ((2*sigma**2)/(1+xi+2*sigma))*((E/sigma)+(B/2)*((1/6)+E))
    c_c = ((2*sigma*omegaT)/(1+xi+2*sigma))*(0.5*((1/6)+E)+(v_0/a)*((1-exp(-1*a))*(E+0.5+(1/(a**2)))-(0.5+((exp(-1*a))/a))))
    m = ((sigma*((1/6)+E))/(1+xi+2*sigma))

    # return stiffness coefficients with dimensions
    k_d = (pi*r_seal*deltaP_seal/lambd)*k_d
    k_c = (pi*r_seal*deltaP_seal/lambd)*k_c
    c_d = T*(pi*r_seal*deltaP_seal/lambd)*c_d
    c_c = T*(pi*r_seal*deltaP_seal/lambd)*c_c
    m = (T**2)*(pi*r_seal*deltaP_seal/lambd)*m
    # print("Returning seal constants...")
    return [k_d, k_c, c_d, c_c, m]

def gulich_seal_calcs(r_shaft, r_impeller, L, c_r, omega, nu, rho, deltaP_impeller):
    '''
    Calculate flow and delta-head across seal using equations described in Gulich Table 3.7 (1)
    INPUTS:
    r_seal - radius of the shaft at seal in m
    r_impeller - outer radius of the impeller in m
    L - length of the seal in m
    c_r - clearance of the seal in m
    omega - angular speed of shaft/impeller in rad/s
    nu - kinematic viscosity of working fluid in m^2/s
    rho - density of working fluid in kg/m^3
    deltaP_impeller - pressure increase across impeller in psi

    OUTPUTS:
    deltaP_sp - deltaP over seal
    u_c_sp - circumferential velocity at seal entrance
    '''

    g = 9.81

    # convert impeller deltaP in psi to delta Head in m
    deltaH_impeller = (deltaP_impeller*6894.76)/(rho*g)

    # Calculate delta-P across seal using eqn 3.7.3 from Gulich
    u_2 = r_impeller*omega # m/s | circumferential velocity at impeller outlet    
    Re_u2 = u_2*r_impeller/nu
    k = 0.9*((Re_u2**0.3)*(c_r*2*r_shaft/((2*r_impeller)**2))*(c_r/L)**0.5)**0.087 # rotation factor | Gulich eqn. 3.7.2
    deltaH_sp = deltaH_impeller - (k**2)*((u_2**2)/(2*g))*(1-(r_shaft/r_impeller)**2) # delta-Head over seal | Gulich eqn. 3.7.3
    deltaP_sp = deltaH_sp*rho*g

    # calculate flow circumferential speed at seal using conservation of angular momentum
    u_c_sp = (r_impeller/r_shaft)*u_2
    
    return [deltaP_sp, u_c_sp]

if __name__ == "__main__":
    [k_d, k_c, c_d, c_c, m] = childs_seal_calc(0, 0, 0, 0, 0, 0, 0, 0, 0)
    print("k_d = " + str(k_d))
    print("k_c = " + str(k_c))
    print("c_d = " + str(c_d))
    print("c_c = " + str(c_c))
    print("m = " + str(m))
