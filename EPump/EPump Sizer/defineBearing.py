import math
class bearingClass:
    def __init__(self,bearingType,t,dt):
        if bearingType == "AC":
            self.f = 1
            self.d1 = 20
            self.d2 = 47
            self.tb = t+dt
            self.v1 = y = 8*10**10*self.tb**-4.282

        elif bearingType == "DG":
            self.f = 1
            self.d1 = 20
        else:
            raise ValueError("Unknown bearing type")


    def heating(self,bearingType,n,fax,fr):
        if bearingType == "AC":
            #inputs
            d_m = .5*(self.d1+self.d2) # once bearing selection is included in the function will need to make some of this parametric as the constants could change
            kz = 4.4
            R1 = 4.33*(10**-7)
            R2 = 2.02
            R3 = 2.44 * (10**-12)
            S1 = .0182
            S2 = .71
            S3 = 2.44*(10**(-12))
            self.fr = fr/2 #split in two between both bearings
            B = 14
            H = 47
            if H/d_m > 1.2:
                V_m = .00125
            else:
                V_m = (H/d_m)*(.00125/1.2)
            i_rw = 1
            # Rolling Moment
            phi_ish = 1/(1+((1.84*10**(-9))*((n*d_m)**1.28)*(self.v1**0.64)))
            K_rs = 1*(10**(-8))
            phi_rs = 1 / (math.exp  (K_rs * self.v1 * n * (self.d1 + self.d2) * math.sqrt (kz /(2*(self.d2 - self.d1)) )   )  )
            fg = R3 * (d_m**4)*(n**2)
            G_rr = R1*(d_m**1.97)*((self.fr+fg+(R2*fax))**.54)
            self.M_rr = phi_ish*phi_rs*G_rr*((n*self.v1)**.6) 
            # Sliding Moment
            phi_bl = 1 / math.exp(2.6*(10**(-8)) * d_m * ((n*self.v1)**1.4))
            mu_bl = .12
            mu_ehl = .05
            mu_sl = phi_bl*mu_bl + (1-phi_bl)*mu_ehl
            G_sl = S1 * (d_m**.26) * (((self.fr+fg)**(4/3))+(S2*(fax**(4/3))))
            self.M_sl = G_sl*mu_sl
            # Drag moment
            if H/d_m > 1.2:
                t = 2*math.pi
            else:
                t = 2 * math.acos((.6*d_m - H)/(.6*d_m))
            if t < math.pi:
                ft = math.sin(.5*t)
            else:
                ft = 1
            k_ball = (10**(-12))*i_rw*kz*(self.d2+self.d1)/(self.d2-self.d1)
            fa = 0.05*(kz*(self.d2+self.d1))/(self.d2-self.d1)
            Rs = 0.36*(d_m**2)*(t-math.sin(t))*fa
            int1 = 0.4*k_ball*V_m*(d_m**5)*(n**2)
            int2 = 1.093*(10**-7)*(n**2)*(d_m**3)
            int3 = Rs*((n*(d_m**2)*ft/self.v1)**(-1.379))
            self.M_drag = int1 + (int2*int3)
            self.p1 = ((1.05*(10**(-4))))*(self.M_rr + self.M_sl + self.M_drag)*n

            # Second, unloaded bearing. All else same except no f_ax
            fax = 0
            phi_ish = 1/(1+((1.84*10**(-9))*((n*d_m)**1.28)*(self.v1**0.64)))
            K_rs = 1*(10**(-8))
            phi_rs = 1 / (math.exp  (K_rs * self.v1 * n * (self.d1 + self.d2) * math.sqrt (kz /(2*(self.d2 - self.d1)) )   )  )
            fg = R3 * (d_m**4)*(n**2)
            G_rr = R1*(d_m**1.97)*((self.fr+fg+(R2*fax))**.54)
            self.M_rr = phi_ish*phi_rs*G_rr*((n*self.v1)**.6) 
            # Sliding Moment
            phi_bl = 1 / math.exp(2.6*(10**(-8)) * d_m * ((n*self.v1)**1.4))
            mu_bl = .12
            mu_ehl = .05
            mu_sl = phi_bl*mu_bl + (1-phi_bl)*mu_ehl
            G_sl = S1 * (d_m**.26) * (((self.fr+fg)**(4/3))+(S2*(fax**(4/3))))
            self.M_sl = G_sl*mu_sl
            # Drag moment
            if H/d_m > 1.2:
                t = 2*math.pi
            else:
                t = 2 * math.acos((.6*d_m - H)/(.6*d_m))
            if t < math.pi:
                ft = math.sin(.5*t)
            else:
                ft = 1
            k_ball = (10**(-12))*i_rw*kz*(self.d2+self.d1)/(self.d2-self.d1)
            fa = 0.05*(kz*(self.d2+self.d1))/(self.d2-self.d1)
            Rs = 0.36*(d_m**2)*(t-math.sin(t))*fa
            int1 = 0.4*k_ball*V_m*(d_m**5)*(n**2)
            int2 = 1.093*(10**-7)*(n**2)*(d_m**3)
            int3 = Rs*((n*(d_m**2)*ft/self.v1)**(-1.379))
            self.M_drag = int1 + (int2*int3)
            self.p2 = ((1.05*(10**(-4))))*(self.M_rr + self.M_sl + self.M_drag)*n
            self.p = self.p1 + self.p2
        elif bearingType == "DG":
            self.p = 50 #W, conservative estimate, do actual math later


    def bearingSummary(self,bearingType,fax,fr):
        if bearingType == "AC":
            print("=== Upper Bearings ===")
            print(f"Radial Force f_r: {fr:.4f} N")
            print(f"Radial Force Per Bearing fr: {self.fr:.4f} N")
            print(f"Axial Force f_a: {fax:.4f} N")
            print(f"Bearing Temp Tb: {self.tb:.4f} K")
            print(f"Viscosity V: {self.v1:.4f} mm^2/s")
            print(f"Loaded Power Loss: {self.p1:.4f} W")
            print(f"Unloaded Power Loss: {self.p2:.4f} W")
            print(f"Total Power Loss: {self.p:.4f} W")
        else:
            print("=== Lower Bearing ===")
            print(f"Radial Force f_r: {fr:.4f} kM")
            print(f"Power P: {self.p:.4f} W")