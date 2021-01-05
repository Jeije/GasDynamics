import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt
from math import atan, tan

global Mr
Mr = 2  # Mach number from duct
Pr = 2e5  # Pressure from duct
Pa = 1e5  # Ambient pressure
PaPrratio = Pa / Pr
global gamma
gamma = 1.4  # Heat capacity ratio of gas
global N
N = 4  # Number of lines in expansion fan


def main():
    number_compl_waves = 3
    alpha = x = y = phi = M = np.zeros((number_compl_waves, N, N))
    # Find mach number assosiated with ambient pressure
    M_ACD = isentropic_exp(Mr, PaPrratio)

    # Find boundaries of ABC
    alpha_AB = -mu(Mr)
    phi_AC = nu(M_ACD) - nu(Mr)  # V+
    alpha_AC = phi_AC - mu(M_ACD)

    # Discretize expansion fan
    alpha_ABC = phi_ABC = M_ABC = np.zeros((N,N))
    alpha_ABC[-1, :] = np.linspace(alpha_AB, alpha_AC, N)
    phi_ABC[-1, :], M_ABC[-1, :] = propcalc(alpha_ABC[-1, :], Mr)

    # Set start point wave
    x_A = np.zeros([N, N])
    y_A = np.ones([N, N])


    # Calculate initial wave
    alpha[0], x[0], y[0], phi[0], M[0] = waveinteraction(alpha_ABC, phi_ABC, M_ABC, x_A, y_A)

    # Calculate all waves after
    for n in range(1,number_compl_waves):
        alpha[n], x[n], y[n], phi[n], M[n] = waveinteraction(alpha[n-1], phi[n-1], M[n-1], x[n-1], y[n-1])



def nu(M):
    nu = ((gamma + 1) / (gamma - 1)) ** 0.5 * atan(((gamma + 1) / (gamma - 1) * (M ** 2 - 1)) ** 0.5) - atan(
        (M ** 2 - 1) ** 0.5)
    return nu


def invnu(findnu):
    #M = optimize.fixed_point(nu, findnu)
    M = 2.1
    return M


def mu(M):
    mu = atan(1 / M)
    return mu


def invmu(mu):
    M = tan(mu) ** -1
    return M


def isentropic_exp(M1, PaPrratio):
    # Find pressure ratio assuming homentropic flow
    ptpr = (1 + M1 ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
    ptpa = ptpr / PaPrratio

    def f(M, ptpa):
        ptpa_guess = (1 + M ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
        zero = ptpa - ptpa_guess
        return zero

    M2 = optimize.newton(f, 2, args=[ptpa])
    return M2


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


def waveinteraction(alpha_left, phi_left, M_left, x_left, y_left):
    alpha = np.zeros([N, N])
    alpha[:, 0] = -alpha_left[-1, :]
    x = np.zeros([N, N])
    y = np.zeros([N, N])
    phi = np.zeros([N, N])
    phi[:, 0] = phi_left[-1, :]
    M = np.zeros([N, N])
    M[:, 0] = M_left[-1, :]
    matnu = np.zeros([N,N])
    for n in range(len(M_left[-1, :])):
        matnu[n, 0] = nu(M_left[-1, n])

    for j in range(1, N):
        x[j, j] = y_left[j, j] / tan(alpha_left[-1, j])
        y[j, j] = 0
        alpha[j, j] = -alpha_left[-1, j]
        phi[j, j] = phi_left[-1, 0]
        M[j, j] = M_left[-1, 0]
        matnu[j, j] = nu(M[j, j])
        for i in range(j + 1, N):
            x[i, j] = 0
            y[i, j] = 0
            alpha[i, j] = -alpha_left[-1, j]
            phi[i, j] = 0.5 * (phi[i - 1, j] + phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) - nu(M[i, j - 1]))
            print(i)
            print(j)
            print(M)
            print(phi)
            print(matnu)
            print("nu:")
            print(nu(M[i, j-1]))
            M[i, j] = invnu(0.5 * (phi[i - 1, j] - phi[i, j - 1]) + 0.5 * (matnu[i - 1, j] + matnu[i, j - 1]))
    return alpha, x, y, phi, M


if __name__ == "":
    main()

if __name__ == "__main__":
    # Find mach number assosiated with ambient pressure
    M_ACD = isentropic_exp(Mr, PaPrratio)

    # Find boundaries of ABC
    alpha_AB = -mu(Mr)
    phi_AC = nu(M_ACD) - nu(Mr)  # V+
    alpha_AC = phi_AC - mu(M_ACD)

    # Discretize expansion fan
    alpha_ABC = phi_ABC = M_ABC = np.zeros((N,N))
    alpha_ABC[-1, :] = np.linspace(alpha_AB, alpha_AC, N)
    phi_ABC[-1, :], M_ABC[-1, :] = propcalc(alpha_ABC[-1, :], Mr)
    # phi_ABC[-1, :] = np.linspace(0, phi_AC, N)  # FIX wrong!
    # M_ABC[-1, :] = np.linspace(Mr, M_ACD, N)  # FIX wrong!

    # Set start point wave
    x_A = np.zeros([N, N])
    y_A = np.ones([N, N])

    # Calculate wave interaction BCE
    alpha_BCE, x_BCE, y_BCE, phi_BCE, M_BCE = waveinteraction(alpha_ABC, phi_ABC, M_ABC, x_A, y_A)

    # Get the simple wave region CDEF from BCE array
    alpha_CDEF, phi_CDEF, M_CDEF = alpha_BCE[-1, :], phi_BCE[-1, :], M_BCE[-1, :]

    # Calculate wave interaction DFG
    alpha_DFG, x_DFG, y_DFG, phi_DFG, M_DFG = waveinteraction(alpha_CDEF, phi_CDEF, M_CDEF, x_BCE, y_BCE)

    print(np.transpose(M_DFG))

