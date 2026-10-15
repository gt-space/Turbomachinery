from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve

def sizePump(
        fluid = "oxygen",
        Nss = 500, #RPM m3 s^-.5 m^-.75
        N = 30000, #RPM
        NPSH_marg = 0.5, #(NPSH_operating - NPSH_bd)/NPSH_bd
        P_0 = 50, #psi
        T_0 = 91, #k
        mdot = 3, #kg/s
        deltaP = 400, #psi
        voluteCP = .5,
        r_ind_hub = .01
        ):

        #conversions
        deltaP = deltaP*6895
        P_0 = P_0 * 6895
        omega = N*np.pi/30

        # ========= INDUCER INLET ===========   

        if fluid == "oxygen":
                rho_0 = PropsSI("D", "T", T_0, "P", P_0, fluid)
                p_vap = PropsSI("P", "T", T_0, "Q", 0, fluid)
        else:
                rho_0 = 805
                p_vap = 3447.38 #Pa

        Q = mdot/rho_0
        
        # Calculate Necessary Nss
        NPSH_op = (N * (Q **.5) / Nss) ** (4/3)
        NPSH_bd = NPSH_op/(1+NPSH_marg)
        Nss_bd = N * (Q**.5) / (NPSH_bd**.75)

        # inlet conditions
        P_1 = P_0 + (deltaP * .1) #initial guess
        rho_1 = rho_0

        # solve for beta
        beta_b1 = np.atan(1/(1/phi)))

        #solve for tip radius from stripling curves
        stripling = lambda eps
        # define min radius, correct if below
        epsilon_min = .15
        if epsilon_1 < epsilon_min:
                epsilon_1 = epsilon_min

        #create velocity triangles
        
        #update Fluid States


        # ========= INDUCER OUTLET =========
        LoverT = 1.5 #length to pitch ratio

        
        # populate stage results
        pump = []
        inducer_in = {"Name": "Inducer Inlet", }
        pump.append(inducer_in)


        return pump
