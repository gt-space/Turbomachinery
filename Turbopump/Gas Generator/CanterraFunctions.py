import cantera as ct
import numpy as np
import yaml
from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI

psi2pa = 6894.75729

with open("canterraSetup.yaml", "r") as f:
    data = yaml.safe_load(f)
mdotf = data["FuelSide"]["mdot"]
mdoto = data["OxSide"]["mdot"]
T_ox = data["OxSide"]["T_in"]
P_ox = data["OxSide"]["P_in"]*psi2pa


h_RP1 = -1.7698e6  # J/kg, liquid RP-1 at 298.15 K (CEA value, formation basis)
h_lox = (PropsSI("H", "T", T_ox, "P", P_ox, "Oxygen")
         - PropsSI("H", "T", 298.15, "P", 101325, "Oxygen"))     # CoolProp uses its own reference state, so shift it so O2 gas at 298.15 K, 1 atm = 0




@dataclass
class CombustionResult:
    # Equilibrium products
    T: float             # adiabatic flame temperature [K]
    rho: float           # density [kg/m^3]
    h: float             # specific enthalpy [J/kg]
    s: float             # specific entropy [J/kg/K]
    cp: float            # [J/kg/K]
    cv: float            # [J/kg/K]
    gamma: float         # cp/cv (frozen) [-]
    MW: float            # mean molecular weight [kg/kmol]
    R: float             # specific gas constant [J/kg/K]

def combustionProps(MR=0.5, Pcc=60, soot=False, gas=None, graphite=None):
    """
    Gas-generator outlet properties at chamber equilibrium.

    Args:
        MR: mdot_ox / mdot_fuel
        Pcc: chamber pressure [psi]
        soot: if True, allow solid carbon to form (gas-only properties are still returned)
    """
    if gas is None:
        gas = ct.Solution("nDodecane_Reitz.yaml")
    P = Pcc * psi2pa
    h_in = (h_RP1 + MR * h_lox) / (1 + MR)

    gas.TPY = 300, P, {"c12h26": 1, "o2": MR}
    gas.equilibrate("TP")                        
    gas.HP = h_in, P
    gas.equilibrate("HP")

    if soot:
        if graphite is None:
            graphite = ct.Solution("graphite.yaml")
        mix = ct.Mixture([(gas, 1.0), (graphite, 0.0)])
        mix.equilibrate("HP")

    gamma = gas.cp_mass / gas.cv_mass
    R = ct.gas_constant / gas.mean_molecular_weight
    return CombustionResult(
        T=gas.T, rho=gas.density, h=gas.enthalpy_mass, s=gas.entropy_mass,
        cp=gas.cp_mass, cv=gas.cv_mass, gamma=gamma,
        MW=gas.mean_molecular_weight, R=R,
    )

def combustionPropsMethalox(MR=0.5, Pcc=60, gas=None, h_LCH4 = -5.59e6):
    P = Pcc * psi2pa
    h_in = (h_LCH4 + MR * h_lox) / (1 + MR)
 
    gas.TPY = 300, P, {"CH4": 1, "O2": MR}
    gas.equilibrate("TP")
    gas.HP = h_in, P
    gas.equilibrate("HP")
 
    return CombustionResult(
        T=gas.T, rho=gas.density, h=gas.enthalpy_mass, s=gas.entropy_mass,
        cp=gas.cp_mass, cv=gas.cv_mass, gamma=gas.cp_mass / gas.cv_mass,
        MW=gas.mean_molecular_weight, R=ct.gas_constant / gas.mean_molecular_weight,
    )
