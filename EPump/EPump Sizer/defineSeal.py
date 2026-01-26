import math
class sealClass: 
    def __init__(self,db,dp,p_tank):
        self.d1 = db
        if self.d1 % 2: #stopgap solution for converting from bearing size to seal size, should probably be fixed later with spreadsheet.
            self.d1 += 1 
        elif self.d1 < 14:
            self.d1 = 14 #minimum seal size
        self.d2 = db + 4
        self.d3 = self.d2 + 14 # don't need to account for anything bigger because of the constraint we set
        self.f = .07
        self.k = .5
        self.dp = ((p_tank + (dp/6894.76) - 95)/1000000)*6894.76 # same vibey estimate for pressure loss before bearings, need to validate later. horrible unit conversions
    def powerLoss(self,n):
        self.fa = math.pi*((self.d3**2) - (self.d1**2))/4
        self.br = (self.d3**2 - self.d2**2)/(self.d3**2 - self.d1**2)
        self.tfp = self.dp*1.35*(self.br-self.k)
        self.mfd = (self.d3 + self.d1)/2
        self.rt = self.mfd*self.tfp*self.f*self.fa/2000
        self.p = (self.rt*n)/9548
    def sealSummary(self):
        print("=== Seal ===")
        print(f"Seal diameter D1 : {self.d1} mm")
        print(f"Power loss       : {(self.p*1000):.4f} W")