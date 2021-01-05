import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt
from math import atan, tan, sin, cos, sqrt

global Mr
Mr = 2  # Mach number from duct
global Pr
Pr = 2e5  # Pressure from duct
Pa = 1e5  # Ambient pressure
PaPrratio = Pa / Pr
global gamma
gamma = 1.4  # Heat capacity ratio of gas
global N
N = 4  # Number of lines in expansion fan


def nu(M):
    nu = np.sqrt((gamma + 1)/(gamma - 1)) * np.arctan(np.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - np.arctan(np.sqrt(M ** 2 - 1))
    return nu

def mu(M):
    mu = np.arctan(1 / (M ** 2 - 1))
    return mu

def findPt(P, M):
    Pt = P * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))
    return Pt

def findMach(Pt, P):
    M = np.sqrt(((Pt/P) ** ((gamma - 1) / gamma) - 1) * 2 / (gamma - 1))
    return M

def isentropic_exp(M1, PaPrratio):
    # Find pressure ratio assuming homentropic flow

    ptpr = (1 + M1 ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
    ptpa = ptpr / PaPrratio

    M = np.sqrt(((ptpa) ** ((gamma - 1) / gamma) - 1) * 2 / (gamma - 1))

    return M


def propcalc(alpha, M_a):
    M = np.zeros(N)
    phi = np.zeros(N)

    def f(M, alpha, M_a, i):
        zero = nu(M) - mu(M) - alpha[i] - nu(M_a)
        return zero

    for i in range(N):
        M[i] = optimize.newton(f, 2, args=[alpha, M_a, i])
        phi[i] = alpha[i] + mu(M[i])

    return phi, M




if  __name__ == "__main__":
    # Find mach number assosiated with ambient pressure
    M_ACD = isentropic_exp(Mr, PaPrratio)

    # Find boundaries of ABC
    alpha_AB = -mu(Mr)
    phi_AC = nu(M_ACD) - nu(Mr)  # V+
    alpha_AC = phi_AC - mu(M_ACD)

    #Find total pressure of exhaust
    Pt = Pr * (1 + (gamma - 1) / 2 * Mr ** 2) ** (gamma / (gamma - 1))
    # Discretize expansion fan ABC
    P_ABC = np.linspace(Pr, Pa, N)
    #Find machnumbers in fan ABC
    M_ABC = findMach(Pt, P_ABC)
    print("Mach",M_ABC)
    #Find mach angle
    mu_ABC = mu(M_ABC)
    print("mu", mu_ABC)
    #Find pradtl angle
    nu_ABC = nu(M_ABC)
    print("nu", nu_ABC)
    #Initialize angle alpha for the shockwaves leaving the edge
    phi_ABC = np.zeros(N)
    dydx_ABC = np.zeros(N)
    dydx_ABC[0] = np.tan(phi_ABC[0] - mu_ABC[0])
    for i in range(1, N):
        phi_ABC[i] = nu_ABC[i - 1] - nu_ABC[i] + phi_ABC[i - 1]
        dydx_ABC[i] = np.tan(phi_ABC[i] - mu_ABC[i])

    print("phi",phi_ABC)
    #dy/dx for Gamma minus line


    # x-coordinates
    y_A = 0.5
    y_BC = np.zeros(N)
    print(dydx_ABC)
    x_ABC = - y_A / dydx_ABC
    print(x_ABC)

    # plt.set_ylim(0, 1.1)
    # plt.set_xlim(0, 3)
    fig, axs = plt.subplots(2, 1)
    for i in range(N):
        axs[0].plot([0, x_ABC[i]], [0.5, 0], 'k')

    axs[0].set_title('Waves', fontsize=10)
    axs[0].set_xlim(0, 3)
    axs[0].set_ylim(0, 1.1)

    axs[1].plot()
    axs[1].axis('equal')
    axs[1].set_title('', fontsize=10)

    fig.tight_layout()

    plt.show()

    plt.plot()
