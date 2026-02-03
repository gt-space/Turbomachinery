import numpy as np
import math
from defineImpeller import impellerClass
from defineBearing import bearingClass
from defineSeal import sealClass
from CoolProp.CoolProp import PropsSI


def pumpPower(x,vec):
        n = x[0]
        deltaT = x[1]
        p_all = (n/50000)*40000 # in kW
        eta_V = .95 #only used in axial force calcs, for now an estimate, later from ratio of flowrates
        
        # Assumptions
        # Change for running different cases/if overall design changes.
        d_H = .023 #m, from key sizing, assumption
        e_Rs = .002*1.25
        d_w = .016 #from shaft/key sizing
        d_D = d_w #for axial force
        
        # unpack vec
        prop = vec[0]
        deltaP = vec[1]
        deltaP = deltaP * 6894.76
        mdot = vec[2]
        MR = vec[3]
        p_tank = vec[4]
        p_tank = p_tank * 6894.76
        
        # Setup
        if prop == "rp1":
            rho = 804.59  # kg/m^3
            Q = mdot/((1+MR)*rho)
            Cp = 2050 #J/kg-K
            visc_1 = 2.1*10**-6 #https://gtvault.sharepoint.com/:x:/r/sites/AE-YellowJacketSpaceProgram/_layouts/15/Doc.aspx?sourcedoc=%7B5C9FA83D-3C55-4725-AD59-68FDFCDB4DA2%7D&file=RP-1%20Viscosity%20Curve.xlsx&action=default&mobileredirect=true
            #print(visc_1)
        elif prop == "lox":
            rho = PropsSI('D','T',90,'P|liquid',(p_tank+deltaP),'Oxygen') 
            Q = mdot*MR/((1+MR)*rho)
            Cp = PropsSI('C','T',90,'P|liquid',(p_tank+deltaP),'Oxygen')
            visc_1 = PropsSI('V','T',90,'P|liquid',(p_tank+deltaP),'Oxygen')/PropsSI("D","T",90,'P|liquid',(p_tank+deltaP),'Oxygen')
            #print(visc_1)
        else:
            raise ValueError("Unknown propellant")
        H = deltaP/(rho*9.81)
        
        ## Define first pump in order to get bearing circulation rate
        #units here are very weird. Impeller and shaft in m, seal and bearing in mm
        #deltat in kelvin, pressures are inputted in psi, math in pa, math for bearing done in mpa.
        #1) size Impeller, do beam calcs to find radial forces on bearings
        impeller1 = impellerClass(Q,n,H,rho,e_Rs,d_D,d_H,eta_V,visc_1)
        impeller2 = impellerClass(Q,n,H/(impeller1.eta_H+.15),rho,e_Rs,d_D,d_H,eta_V,visc_1)
        #2) select bearings based on rpm and forces. For now most calculations skipped because of selection complications
        # assume design has same shaft dimensioning as model as of 9/13
        f_rlower = impeller1.f_r * 23/74
        f_rupper = f_rlower + impeller1.f_r
        if prop == "rp1":
            upperBearing = bearingClass("AC",300,deltaT,prop,p_tank,deltaP)     
            lowerBearing = bearingClass("DG",0,0,prop,p_tank,deltaP)
        elif prop == "lox":
            upperBearing = bearingClass("AC",90,deltaT,prop,p_tank,deltaP)     
            lowerBearing = bearingClass("DG",0,0,prop,p_tank,deltaP)
        #upperBearing.bearingSummary("AC",impeller1.f_ax,f_rupper)
        #lowerBearing.bearingSummary("DG",0,f_rlower)
        #lowerBearing.bearingSummary()
        #3) find heating on bearings, remember to incorporate deltaT's effect on viscosity.
        upperBearing.heating("AC",n,impeller1.f_ax,f_rupper)
        lowerBearing.heating("DG",n,0,f_rlower)
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
        impeller3 = impellerClass(Qnew,n,H/(impeller2.eta_H+.15),rho,e_Rs,d_D,d_H,eta_Vnew,visc_1)
        #6) apply all power losses and efficiencies
        p_draw = ((impeller3.p/(impeller3.eta_H)) + seal.p + upperBearing.p + lowerBearing.p)
        # Penalties/constraints
        err = 0
        if impeller3.d_2 > .08:
            err = 1
        if p_draw > p_all:
            err = 2
        if lowerBearing.d1 > 24: #done in mm because that's what the bearing heating calcs are in
            err = 3 # breaks the seal code rn. Also generally good, dont want seal face speeds to get too high.
        powerInfo = {
            'Power Draw (W)': p_draw,
            'Gulich\'s power(W)': impeller3.p/impeller3.eta_opt,
            'Available Power (W)': p_all,
            'hydraulic Power (W)': impeller3.p/impeller3.eta_H,
            'Hydraulic Efficiency': impeller3.eta_H,
            'Mechanical Parasitic Power (W)': seal.p + upperBearing.p + lowerBearing.p,
            'Disk Friction Losses (W)': impeller3.Prr,
        }
        sizeInfo = {    
            'RPM': n,
            'flow rate (m^3/s)': Qnew,
            'DeltaT (K)': deltaT,
            'Head (m)': H,
            'Specific Speed': impeller3.n_q,
            'impeller Diameter (m)': impeller3.d_2,
            'Reynold\'s number': impeller3.re,
            #'Shaft Length (m)': 1000000,
        }
        powerInfo = {k: float(v) for k, v in powerInfo.items()}
        sizeInfo = {k: float(v) for k, v in sizeInfo.items()}
        return p_draw, err, sizeInfo, powerInfo


def pumpPowerNotebook(x,vec):
        n = x[0]
        deltaT = x[1]
        p_all = (n/50000)*40000 # in kW
        eta_V = .95 #only used in axial force calcs, for now an estimate, later from ratio of flowrates
        
        # Assumptions
        # Change for running different cases/if overall design changes.
        d_H = .023 #m, from key sizing, assumption
        e_Rs = .002*1.25
        d_w = .016 #from shaft/key sizing
        d_D = d_w #for axial force
        
        # unpack vec
        prop = vec[0]
        deltaP = vec[1]
        deltaP = deltaP * 6894.76
        mdot = vec[2]
        MR = vec[3]
        p_tank = vec[4]
        p_tank = p_tank * 6894.76
        
        # Setup
        if prop == "rp1":
            rho = 804.59  # kg/m^3
            Q = mdot/((1+MR)*rho)
            Cp = 2050 #J/kg-K
            visc_1 = 2.1*10**-6 #https://gtvault.sharepoint.com/:x:/r/sites/AE-YellowJacketSpaceProgram/_layouts/15/Doc.aspx?sourcedoc=%7B5C9FA83D-3C55-4725-AD59-68FDFCDB4DA2%7D&file=RP-1%20Viscosity%20Curve.xlsx&action=default&mobileredirect=true
            #print(visc_1)
        elif prop == "lox":
            rho = PropsSI('D','T',90,'P|liquid',(p_tank+deltaP),'Oxygen') 
            Q = mdot*MR/((1+MR)*rho)
            Cp = PropsSI('C','T',90,'P|liquid',(p_tank+deltaP),'Oxygen')
            visc_1 = PropsSI('V','T',90,'P|liquid',(p_tank+deltaP),'Oxygen')/PropsSI("D","T",90,'P|liquid',(p_tank+deltaP),'Oxygen')
            #print(visc_1)
        else:
            raise ValueError("Unknown propellant")
        H = deltaP/(rho*9.81)
        
        ## Define first pump in order to get bearing circulation rate
        #units here are very weird. Impeller and shaft in m, seal and bearing in mm
        #deltat in kelvin, pressures are inputted in psi, math in pa, math for bearing done in mpa.
        #1) size Impeller, do beam calcs to find radial forces on bearings
        impeller1 = impellerClass(Q,n,H,rho,e_Rs,d_D,d_H,eta_V,visc_1)
        impeller2 = impellerClass(Q,n,H/(impeller1.eta_H+.15),rho,e_Rs,d_D,d_H,eta_V,visc_1)
        #2) select bearings based on rpm and forces. For now most calculations skipped because of selection complications
        # assume design has same shaft dimensioning as model as of 9/13
        f_rlower = impeller1.f_r * 23/74
        f_rupper = f_rlower + impeller1.f_r
        if prop == "rp1":
            upperBearing = bearingClass("AC",300,deltaT,prop,p_tank,deltaP)     
            lowerBearing = bearingClass("DG",0,0,prop,p_tank,deltaP)
        elif prop == "lox":
            upperBearing = bearingClass("AC",90,deltaT,prop,p_tank,deltaP)     
            lowerBearing = bearingClass("DG",0,0,prop,p_tank,deltaP)
        #upperBearing.bearingSummary("AC",impeller1.f_ax,f_rupper)
        #lowerBearing.bearingSummary("DG",0,f_rlower)
        #lowerBearing.bearingSummary()
        #3) find heating on bearings, remember to incorporate deltaT's effect on viscosity.
        upperBearing.heating("AC",n,impeller1.f_ax,f_rupper)
        lowerBearing.heating("DG",n,0,f_rlower)
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
        impeller3 = impellerClass(Qnew,n,H/(impeller1.eta_H+.15),rho,e_Rs,d_D,d_H,eta_Vnew,visc_1)
        #6) apply all power losses and efficiencies
        p_draw = ((impeller3.p/(impeller3.eta_H)) + seal.p + upperBearing.p + lowerBearing.p)
        # Penalties/constraints
        err = 0
        if impeller3.d_2 > .08:
            err = 1
        if p_draw > p_all:
            err = 2
        if lowerBearing.d1 > 24: #done in mm because that's what the bearing heating calcs are in
            err = 3 # breaks the seal code rn. Also generally good, dont want seal face speeds to get too high.
        powerInfo = {
            'Power Draw (W)': p_draw,
            'Gulich\'s power(W)': impeller3.p/impeller3.eta_opt,
            'Available Power (W)': p_all,
            'hydraulic Power (W)': impeller3.p/impeller3.eta_H,
            'Hydraulic Efficiency': impeller3.eta_H,
            'Mechanical Parasitic Power (W)': seal.p + upperBearing.p + lowerBearing.p,
            'Disk Friction Losses (W)': impeller3.Prr,
        }
        sizeInfo = {    
            'RPM': n,
            'flow rate (m^3/s)': Qnew,
            'DeltaT (K)': deltaT,
            'Head (m)': H,
            'Specific Speed': impeller3.n_q,
            'impeller Diameter (m)': impeller3.d_2,
            'Reynold\'s number': impeller3.re,
            #'Shaft Length (m)': 1000000,
        }
        powerInfo = {k: float(v) for k, v in powerInfo.items()}
        sizeInfo = {k: float(v) for k, v in sizeInfo.items()}
        return impeller1,impeller2,impeller3, Qcooling, p_draw