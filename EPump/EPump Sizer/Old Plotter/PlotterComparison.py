import numpy as np
import math
from defineImpeller import impellerClass
from defineBearing import bearingClass
from defineSeal import sealClass
import matplotlib.pyplot as plt
from Plotter1 import pumpPlot

# current parameters: ("rp1",460,6.17,2.23,285,135,0) #tamb should stay above 293 or the bearing lubrication goes a little iffy
def pumpPlot(prop, deltaP, mdot, MR,Tamb,p_tank):
    plt.close("all")
    #inputs assumed to be same across pumps
    d_H = .023 #m, from key sizing
    e_Rs = .002*1.25
    deltaP = deltaP * 6894.76
    d_w = .016 #from shaft/key sizing
    d_D = d_w #for axial force
    
    # Fluid properties
    if prop == "rp1":
        rho = 804.59  # kg/m^3
        Q = mdot/((1+MR)*rho)
        Cp = 2050
    elif prop == "lox":
        rho =  1140 # kg/m^3
        Q = mdot*MR/((1+MR)*rho)
    else:
        raise ValueError("Unknown propellant")
    H = deltaP/(rho*9.81)

    def pumpPower(n):
        deltaT = 1
        eta_V = .95111944218 #only used in axial force calcs, for now an estimate, later from ratio of flowrates


        ## Define first pump in order to get bearing circulation rate
        #units here are very weird. Impeller and shaft in m, seal and bearing in mm
        #deltat in kelvin, pressures are inputted in psi, math in pa, math for bearing done in mpa.
        #1) size Impeller, do beam calcs to find radial forces on bearings
        impeller1 = impellerClass(Q,n,H,rho,e_Rs,d_D,d_H,eta_V)
        #impeller1.summary()

        #2) select bearings based on rpm and forces. For now most calculations skipped because of selection complications
        # assume design has same shaft dimensioning as model as of 9/13
        f_rlower = impeller1.f_r * 23/74
        f_rupper = f_rlower + impeller1.f_r
        upperBearing = bearingClass("AC",Tamb,deltaT)     
        lowerBearing = bearingClass("DG",0,0)
        upperBearing.heating("AC",n,impeller1.f_ax,f_rupper)
        lowerBearing.heating("DG",n,0,f_rlower)
        #upperBearing.bearingSummary("AC",impeller1.f_ax,f_rupper)
        #lowerBearing.bearingSummary("DG",0,f_rlower)
        #lowerBearing.bearingSummary()
        #3) find heating on bearings, remember to incorporate deltaT's effect on viscosity.

        #4) size seals and find heating
        seal = sealClass(lowerBearing.d1,deltaP,p_tank)
        seal.powerLoss(n)
        #seal.sealSummary()
        #6) find required florwate
        Mcooling = upperBearing.p/(Cp*deltaT) #temp
        Qcooling = Mcooling/rho
        #5)rerun impeller sizing with new mdot
        Qnew = Q + Qcooling
        eta_Vnew = Q/((Qnew)*1.03) # assuming 3% leak rate
        impeller2 = impellerClass(Qnew,n,H,rho,e_Rs,d_D,d_H,eta_Vnew)
        #impeller2.summary()
        p_imp = impeller2.p #converting to kW, apply hydraulic efficiency and leak rate
        #6) apply all power losses and efficiencies
        p_draw1 = ((p_imp/(impeller2.eta_H*1.03)) + seal.p + upperBearing.p + lowerBearing.p)
        p_draw2 = p_imp/(impeller2.eta_H*1.03) * 1.2
        p_draw3 = p_imp/(impeller2.eta_opt*1.03) * 1.2
        # Penalties/constraints
        err = 0
        if impeller2.d_2 > .08:
            err = 1
        if lowerBearing.d1 > 24: #done in mm because that's what the bearing heating calcs are in
            err = 1 # breaks the seal code rn. Also generally good, dont want seal face speeds to get too high.
        return p_draw1, p_draw2, p_draw3, err
    # plotting section
    nrange = np.linspace(15000,50000,100)
    y1, y2, y3, errs = zip(*(pumpPower(xi) for xi in nrange))
    y1 = np.array(y1)
    y2 = np.array(y2)
    y3 = np.array(y3)

    plt.plot(nrange, y1, label='Sizer')
    plt.plot(nrange, y2, label='Shifted Hydraulic Gulich')
    plt.plot(nrange, y3, label='Totall Gulich')

    plt.xlabel('n(RPM)')
    plt.ylabel('Power (kW)')
    plt.legend()
    plt.grid(True)
    plt.show()
    