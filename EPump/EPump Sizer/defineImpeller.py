import math
class impellerClass:
    def __init__(self,vdot,n,H,rho,e_Rs,d_D,d_H,eta_V):
        self.Q = vdot
        self.psi = 1.09
        self.n_q = n*(math.sqrt(vdot))*(H**(-3/4))
        self.d_2 = (84.6 / n) * math.sqrt(H / self.psi)
        #efficiencies
        if vdot <= 1:
            a = 1 
        elif vdot > 1:
            a = .5
        m = 0.1 * a * ((1/vdot)**(.15)) * ((45/self.n_q)**(0.06))
        eta_opt = 1 - (0.095*((1/vdot)**m)) - (.3*((.35 - math.log10(self.n_q/23))**2)*((1/vdot)**(.05)))
        m = 0.08 * a * ((1/vdot)**(.15)) * ((45/self.n_q)**(0.06))
        self.eta_H =  1 - (0.055*((1/vdot)**m)) - .2*((.26 - math.log10(self.n_q/25))**2)*((1/vdot)**(.1))
        #  casing gap below impeller
        s_ax = self.d_2*.035
        # inlet diameter
        f_d1 = 1.15
        d_Hstar = d_H/self.d_2
        self.d_1 = self.d_2 * f_d1 * math.sqrt((d_Hstar **2) + 0.00148 * self.psi * self.n_q ** ( 4 / 3 ) / ( (eta_V) ** .67 ))
        # axial force
        d_sp = 2*(.001) + self.d_1 + 2*(e_Rs)
        self.f_ax = .9* rho * 9.81 * H * (math.pi/4) * (d_sp**2 - d_D**2)
        self.p = rho*9.81*H*vdot
        # Radial forces
        eps_sp = 180
        b_2 = b_2star = .017 + self.n_q*(.262/100) - .08*((self.n_q/100)**2) + .0093*((self.n_q/100)**3)
        b_2 = b_2star*self.d_2
        b_2tot = b_2 +2*e_Rs
        F_DSp = 1.75 - 0.0083*(eps_sp)
        k_R0 = -5*(10**(-5))*(self.n_q**2) + 9.2*(10**-3)*(self.n_q)+1.91*(10**-1)
        k_RDSp = F_DSp * k_R0
        F_Rst = k_RDSp * rho * 9.81 * H * self.d_2 * b_2tot
        F_Rdyn = .12 * rho * 9.81 * H * self.d_2 * b_2tot 
        self.f_r = F_Rdyn + F_Rst 
    def summary(self):
        print("=== Impeller ===")
        print(f"Flow rate Q      : {self.Q:.4f} m³/s")
        print(f"Specific speed n_q: {self.n_q:.3f}")
        print(f"Diameter d_2     : {self.d_2:.4f} m")
        print(f"Inlet diameter d_1: {self.d_1:.4f} m")
        print(f"Hydraulic efficiency: {self.eta_H:.3f}")
        print(f"Hydraulic power   : {self.p:.2f} W")
        print(f"Axial force f_ax  : {self.f_ax:.2f} N")
        print(f"Radial force f_r  : {self.f_r:.2f} N")