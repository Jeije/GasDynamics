import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

Mr = 2  # Mach number from duct
Pr = 2e5  # Pressure from duct
Pa = 1e5  # Ambient pressure
PaPrratio = Pa / Pr
global gamma
gamma = 1.4  # Heat capacity ratio of gas
global N
N = 4  # Number of lines in expansion fan


def main():
    print("Hellow")


def nu(M):
    nu = np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(np.sqrt((gamma + 1) / (gamma - 1) * (Mr ** 2 - 1))) - np.arctan(
         np.sqrt(M ** 2 - 1))
    return nu


def mu(M):
    mu = np.arctan(1 / M)
    return mu


def discr(alpha_a, alpha_b):
    alpha = np.zeros(N)
    for i in range(len(alpha)):
        alpha[i] = alpha_a + i * (alpha_a - alpha_b) / N
    print(alpha_a)
    print(alpha_b)
    print(alpha)
    return alpha


def fanprop(alpha, M_a, phi_a, M_b, phi_b):
    #a is incomming mach number and direction and b is outcomming. alpha is hte spacing
    phi = np.zeros(len(alpha))
    M = np.zeros(len(alpha))
    phi[0] = 0
    M[0] = Mr
    nu_a = nu(M_a)
    nu_b = nu(M_b)

    for i in range(1, len(alpha)):
        def f(M, nu_a, phi_a, alpha):
            zero = alpha + nu(M) - mu(M) -1 * (nu_a - phi_a)
            return zero
        M[i] = optimize.newton(f, 2, args=[nu_a, phi_a, alpha[i]])
        phi[i] = mu(M[i]) - alpha[i]
    print(phi)
    print(M)
    return phi, M


def isentropic_exp(M1, PaPrratio):
    # Find pressure ratio assuming homentropic flow
    P0Prat = (1 + M1 ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
    P0Prat_out = P0Prat / PaPrratio
    def f(M, P0Prat_out):
        P0Pratio = (1 + M ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
        zero = P0Pratio - P0Prat_out
        return zero

    M2 = optimize.newton(f, M1, args=[P0Prat_out])
    return M2


if __name__ == "__main__":
    M_ACD = isentropic_exp(Mr, PaPrratio)
    print(M_ACD)
    phi_ACD = nu(M_ACD) - nu(Mr) + 0  # V+
    alpha_A = mu(Mr) - 0
    alpha_B = phi_ACD - mu(M_ACD)
    alpha_ABC = discr(alpha_A, alpha_B)
    print(alpha_ABC)
    phi_ABC, M_ABC = fanprop(alpha_ABC, Mr, 0, M_ACD, phi_ACD)  # phi and M from alpha M and phi
    print(phi_ABC)
    print(M_ABC)

    for i in range(len(alpha_ABC)):
        plt.plot([0, 0 + np.cos(alpha_ABC[i])], [1, 1 - np.sin(alpha_ABC[i])])
    plt.plot([0, 0 + np.cos(phi_ACD)], [1, 1 + np.sin(phi_ACD)], "-.")
    plt.ylim(0, 1.5)
    plt.xlim(0,5)
    plt.axis("equal")
    plt.show()
    # x=np.linspace(-3,10,1000)
    # y=np.zeros(len(x))
    # for i in range(len(x)):
    #     y[i]=nu(x[i])
    # y1=mu(x)
    # plt.plot(x,y)
    # plt.plot(x,y1)
    # plt.show()
