import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt
from math import atan, tan, sin, cos, sqrt

global Mr
Mr = 2  # Mach number from duct
Pr = 2e5  # Pressure from duct
Pa = 1e5  # Ambient pressure
PaPrratio = Pa / Pr
global gamma
gamma = 1.4  # Heat capacity ratio of gas
global N
N = 4  # Number of lines in expansion fan


def nu(M):
    nu = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma + 1) / (gamma - 1) * (M ** 2 - 1))) - atan(sqrt(M ** 2 - 1))
    return nu


def invnu(findnu):
    print("find " + str(findnu))

    def f(M):
        zero = findnu - nu(M)
        return zero

    M = optimize.newton(f, 2.1)
    print("found " + str(M) + " with " + str(nu(M)))
    return M


def mu(M):
    mu = np.arctan(1 / (M ** 2 - 1))
    return mu


def invmu(mu):
    M = 1 / sin(mu)
    return M


# def isentropic_exp(M1, PaPrratio):
#     # Find pressure ratio assuming homentropic flow
#     ptpr = (1 + M1 ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
#     ptpa = ptpr / PaPrratio
#
#     def f(M, ptpa):
#         ptpa_guess = (1 + M ** 2 * (gamma - 1) / 2) ** (gamma / (gamma - 1))
#         zero = ptpa - ptpa_guess
#         return zero
#
#     M2 = optimize.newton(f, 2, args=[ptpa])
#     return M2


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


def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0] * p2[1] - p2[0] * p1[1])
    return A, B, -C


def intersection(L1, L2):
    D = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x, y
    else:
        return False


def findintersect(alpha_A_or_phi, alpha_B, x_A, x_B, y_A, y_B):
    LA = line([x_A, y_A], [x_A + cos(alpha_A_or_phi), y_A + sin(alpha_A_or_phi)])
    LB = line([x_B, y_B], [x_B + cos(alpha_B), y_B + sin(alpha_B)])
    x, y = intersection(LA, LB)
    return x, -y


def waveinteractionmirror(alpha_left, phi_left, M_left, x_left, y_left):
    x = np.zeros([N, N])
    y = np.zeros([N, N])
    phi = np.zeros([N, N])
    phi[:, 0] = phi_left
    M = np.zeros([N, N])
    M[:, 0] = M_left
    alpha_min = np.zeros([N, N])
    alpha_plus = np.zeros([N, N])
    alpha_min[:, 0] = alpha_left
    alpha_plus[:, 0] = alpha_left + 2 * mu(M_left)

    for j in range(1, N):
        x[j, j] = - y_left[j, j] / tan(alpha_left[j])
        y[j, j] = 0
        phi[j, j] = phi_left[0]
        M[j, j] = M_left[0]
        alpha_min[j, j] = phi[j, j] - mu(M[j, j])
        alpha_plus[j, j] = phi[j, j] + mu(M[j, j])
        x[j, j], y[j, j] = findintersect(alpha_left[j], 0, x_left[-1, j], 0, y_left[-1, j], 0)
        for i in range(j + 1, N):
            # print("Finding x and y")
            # x[i, j], y[i, j] = findintersect(alpha[i, 0], alpha[0, j], x_left[-1, i], x[i, i], y_left[-1, i], y[i, i])
            phi[i, j] = 0.5 * (phi[i - 1, j] + phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) - nu(M[i, j - 1]))
            M[i, j] = invnu(0.5 * (phi[i - 1, j] - phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) + nu(M[i, j - 1])))
            print(M[i, j])
            alpha_min[i, j] = phi[i, j] - mu(M[i, j])
            alpha_plus[i, j] = phi[i, j] + mu(M[i, j])

    # x[0, 0], y[0, 0] = findintersect(0.5*(alpha_min[0, 0]+alpha_left[0]), alpha_plus[0, 0], x_left[-1, 0], x[0, 0], y_left[-1, 0],
    #                                 y[0, 0])
    # x[0, 0] = -y_left[j, j] / tan(alpha_left[j])
    # y[0, 0] = 0
    # for j in range(N):
    #     for i in range(j + 1, N):
    #         x[i, j], y[i, j] = findintersect(0.5 * (alpha_min[i, j] + alpha_left[i]), alpha_plus[-1, j], x_left[-1, i],
    #                                          x[i, i], y_left[-1, i],
    #                                          y[i, i])

    return alpha_plus, x, y, phi, M


# def waveinteractionjet(alpha_left, phi_left, M_left, x_left, y_left):
#     alpha = np.zeros([N, N])
#     alpha[:, 0] = -alpha_left
#     x = np.zeros([N, N])
#     y = np.zeros([N, N])
#     phi = np.zeros([N, N])
#     phi[:, 0] = phi_left
#     M = np.zeros([N, N])
#     M[:, 0] = M_left
#
#     for j in range(1, N):
#         phi[j, j] = phi_left[0]
#         M[j, j] = M_left[0]
#         alpha[j, j] = alpha_left[j] - 2 * mu(M[j, j])
#         x[j, j], y[j, j] = findintersect(phi[j, j], alpha[-1, j], x_left[-1, j], x[j, j], y_left[-1, j],
#                                          y[j, j])
#         for i in range(j + 1, N):
#             alpha[i, j] = -alpha_left[j]
#             phi[i, j] = 0.5 * (phi[i - 1, j] + phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) - nu(M[i, j - 1]))
#             M[i, j] = invnu(0.5 * (phi[i - 1, j] - phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) + nu(M[i, j - 1])))
#
#     x[0, 0], y[0, 0] = findintersect(alpha_left[0], alpha[-1, 0], x_left[-1, 0], x[0, 0], y_left[-1, 0],
#                                      y[0, 0])
#
#     # for j in range(N):
#     #     for i in range(j + 1, N):
#     #         x[i, j], y[i, j] = findintersect(0.5*(alpha_min[i, j]+alpha_left[i]), alpha_plus[-1, j], x_left[-1, i], x[i, i], y_left[-1, i],
#     #                                          y[i, i])
#     return alpha, x, y, phi, M

def waveinteractionjet(alpha_left, phi_left, M_left, x_left, y_left, x_border, y_border, phi_border):
    alpha = np.zeros([N, N])
    x = np.zeros([N, N])
    y = np.zeros([N, N])
    phi = np.zeros([N, N])
    phi[:, 0] = phi_left
    M = np.zeros([N, N])
    M[:, 0] = M_left

    for j in range(1, N):
        alpha[j, j] = -(alpha_left[j] - 2 * phi_left[0])
        phi[j, j] = phi_left[j]
        M[j, j] = M_left[j]
        print(x_left)
        print(y_left)
        print(phi_border, alpha_left[j], x_left[j], x_border, y_left[j], y_border)
        x[j, j], y[j, j] = findintersect(phi_border, alpha_left[j], x_left[j], x_border, y_left[j], y_border)
        for i in range(j + 1, N):
            alpha[i, j] = -alpha_left[j]
            phi[i, j] = 0.5 * (phi[i - 1, j] + phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) - nu(M[i, j - 1]))
            M[i, j] = invnu(0.5 * (phi[i - 1, j] - phi[i, j - 1]) + 0.5 * (nu(M[i - 1, j]) + nu(M[i, j - 1])))
    return alpha, x, y, phi, M


if __name__ == "__main__":
    # Find mach number assosiated with ambient pressure
    M_ACD = isentropic_exp(Mr, PaPrratio)

    # Find boundaries of ABC
    alpha_AB = -mu(Mr)
    phi_AC = nu(M_ACD) - nu(Mr)  # V+
    alpha_AC = phi_AC - mu(M_ACD)

    # Discretize expansion fan
    alpha_ABC = np.linspace(alpha_AB, alpha_AC, N)
    phi_ABC, M_ABC = propcalc(alpha_ABC, Mr)

    # Set start point wave
    x_A = np.zeros([N, N])
    y_A = np.ones([N, N])

    # Calculate wave interaction BCE
    print("Calculate wave interaction BCE")
    alpha_BCE, x_BCE, y_BCE, phi_BCE, M_BCE = waveinteractionmirror(alpha_ABC, phi_ABC, M_ABC, x_A, y_A)

    # Get the simple wave region CDEF from BCE array
    print("Get the simple wave region CDEF from BCE array")
    alpha_CDEF, phi_CDEF, M_CDEF, x_CDEF, y_CDEF = alpha_BCE[-1, :], phi_BCE[-1, :], M_BCE[-1, :], x_BCE[-1, :], y_BCE[-1, :]

    # Calculate wave interaction DFG
    print("Calculate wave interaction DFG")
    alpha_DFG, x_DFG, y_DFG, phi_DFG, M_DFG = waveinteractionjet(alpha_CDEF, phi_CDEF, M_CDEF, x_CDEF, y_CDEF, x_border=0,
                                                                 y_border=1, phi_border=phi_AC)

    # Get the simple wave region FGHI from DFG array
    alpha_FGHI, phi_FGHI, M_FGHI = alpha_DFG[-1, :], phi_DFG[-1, :], M_DFG[-1, :]

    # Calculate wave interaction HIJ
    alpha_HIJ, x_HIJ, y_HIJ, phi_HIJ, M_HIJ = waveinteractionmirror(alpha_FGHI, phi_FGHI, M_FGHI, x_DFG, y_DFG)

    # Get the simple wave region IJKL from HIJ array
    alpha_IJKL, phi_IJKL, M_IJKL = alpha_HIJ[-1, :], phi_HIJ[-1, :], M_HIJ[-1, :]

    # Plotting sequence
    print("Start plotting")
    fig, axs = plt.subplots(2, 1)

    for i in range(N):
        for j in range(N):
            axs[0].plot(x_A[i, j], y_A[i, j], 'o')

    for i in range(N):
        for j in range(N):
            print("Plotting x= " + str(x_BCE[i, j]) + " y= " + str(y_BCE[i, j]))
            axs[0].plot(x_BCE[i, j], y_BCE[i, j], 'o')

    for i in range(N):
        axs[0].plot([x_A[i, i], x_BCE[i, i]], [y_A[i, i], y_BCE[i, i]], 'k')

    for i in range(N):
        axs[0].plot([x_BCE[i, i], x_DFG[i, i]], [y_BCE[i, i], y_DFG[i, i]], 'k')

    for i in range(N):
        axs[0].plot([x_DFG[i, i], x_HIJ[i, i]], [y_DFG[i, i], y_HIJ[i, i]], 'k')

    axs[0].set_title('Waves', fontsize=10)
    axs[0].set_xlim(0, 3)
    axs[0].set_ylim(0, 1.1)

    axs[1].plot()
    axs[1].axis('equal')
    axs[1].set_title('', fontsize=10)

    fig.tight_layout()

    plt.show()

    plt.plot()
