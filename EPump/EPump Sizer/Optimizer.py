import numpy as np
from prometheus_client import Info
from scipy.optimize import minimize
from powerDraw import pumpPower

def epumpOpt(prop, deltaP, mdot, MR, p_in):

    # Optimization
    x0 = [40000, 2] 
    bounds = [(15000, 50000), (1, 5)] 
    vec = [prop, deltaP, mdot, MR, p_in]

    # Objective function
    def objective(x):
        power, err, _, _ = pumpPower(x, vec)
        return power + 1e6 * err

    # Run optimization
    optVals = minimize(
        objective,
        x0,
        method='SLSQP',
        bounds=bounds,)

    # Display results
    power, err, sizeInfo, powerInfo = pumpPower(optVals.x, vec)
    if err == 1:
        print("Warning: Impeller too large.")
    elif err == 2:
        print("Warning: Pump power exceeds available power.")
    elif err == 3:
        print("Warning: Bearing size exceeds limit.")
    return sizeInfo,powerInfo

def printInfo(Info):
    if 'Power Draw (W)' in Info:
        print("\nPump Power Results")
        print("-" * 45)
    else:
        print("\nPump Sizing Results")
        print("-" * 45)
    for key, value in Info.items():
        if isinstance(value, float):
            print(f"{key:<30}: {value:>12.4g}")
        else:
            print(f"{key:<30}: {value}")

def run_epump(prop, deltaP, mdot, MR, p_in):
    results = epumpOpt(
        prop=prop,
        deltaP=deltaP,
        mdot=mdot,
        MR=MR,
        p_in=p_in
    )
    printInfo(results[0])
    printInfo(results[1])
    return results
