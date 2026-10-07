from scipy.optimize import minimize, NonlinearConstraint
import numpy as np
def minimizeTurbRad(rpm = 30000,
                    R = 0,
                    t0in = 1000,
                    p0in = 482633,
                    pout = 101325,
                    work = 30000,
                    mdotMax = .5
                    tipspeedMax = 350
                    ):


    # sizing function
    def radius(x):
        phi, psi = x
        U = np.sqrt(dh_s / psi)

        turbine = []
        return turbine.r

    #objective function
    def objective(x):
        return radius(x)   # minimize radius

    # functions used for constraints
    def mdot_calc(x):
        phi, psi = x
        U = np.sqrt(dh_s / psi)
        D = 2 * radius(x)
        return phi * rho_in * U * D**2

    def U_calc(x):
        phi, psi = x
        return np.sqrt(dh_s / psi)



    # Initialization
    cons = [
        NonlinearConstraint(mdot_calc, -np.inf, mdot_max),
        NonlinearConstraint(U_calc,    -np.inf, U_max),
    ]

    bounds = [(0.01, 1.5),   # phi typical range
            (0.5, 3.0)]    # psi typical range

    x0 = [0.3, 1.0]

    # Run Optimizer
    result = minimize(objective, x0, method='SLSQP',
                    bounds=bounds, constraints=cons,
                    options={'disp': True})


    # sizing function
    def getTurbine(x):

        turbine = []
        return turbine
    turbine = getTurbine(result.x)
    return turbine
