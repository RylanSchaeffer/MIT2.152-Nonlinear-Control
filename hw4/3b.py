# Slotine and Li, Example 9.3
import matplotlib.pyplot as plt
import numpy as np

# problem coefficients
a_1 = 1 + 1 * (0.5 ** 2) + 0.25 + 2 * (0.6 ** 2) + 2
a_2 = 0.25 + 2 * (0.6 ** 2)
a_3 = 2 * 0.6 * np.cos(.5)
a_4 = 2 * 0.6 * np.sin(0.5)


def calc_H(q):
    H = np.array([
        [a_1 + 2 * a_3 * np.cos(q[1]) + 2 * a_4 * np.sin(q[1]),
         a_2 + a_3 * np.cos(q[1]) + a_4 * np.sin(q[1])],
        [a_2 + a_3 * np.cos(q[1]) + a_4 * np.sin(q[1]),
         a_2]
    ])
    return H


def calc_C(q, qdot):
    h = a_3 * np.sin(q[1]) - a_4 * np.cos(q[1])
    C = np.array([
        [-h * qdot[1], -h * (qdot[0] + qdot[1])],
        [h * qdot[0], 0]
    ])
    return C


def calc_Y(q, qdot, qdot_r, qdotdot_r):
    Y_13 = (2 * qdotdot_r[0] + qdotdot_r[1])*np.cos(q[1]) - \
           (qdot[1]*qdot_r[0] + qdot[0]*qdot_r[1] + qdot[1]*qdot_r[1])*np.sin(q[1])
    Y_14 = (2 * qdotdot_r[0] + qdotdot_r[1])*np.sin(q[1]) + \
           (qdot[1]*qdot_r[0] + qdot[0]*qdot_r[1] + qdot[1]*qdot_r[1])*np.cos(q[1])
    Y = np.array([
        [qdotdot_r[0],
         qdotdot_r[1],
         Y_13,
         Y_14],
        [0,
         qdotdot_r[0] + qdotdot_r[1],
         qdotdot_r[0] * np.cos(q[1]) + qdot[0]*qdot_r[0]*np.sin(q[1]),
         qdotdot_r[1] * np.sin(q[1]) - qdot[0]*qdot_r[0]*np.cos(q[1])],
    ])
    return Y


# simulation parameters
num_dts = 25000
control_history = np.zeros(shape=(3, num_dts, 2))  # three gammas
error_history = np.zeros(shape=(3, num_dts, 2))

dt = 0.0001
K_p, K_d = 100*np.eye(2), 20*np.eye(2)  # controller positive-definite matrices
lambd = 20*np.eye(2)
gamma_1 = np.diag([0.03, 0.05, 0.1, 0.3])
gammas = [gamma_1, 200 * gamma_1, 0.1 * gamma_1]


for i, gamma in enumerate(gammas):
    q = np.zeros(2)  # initial angles, radians
    qdot = np.zeros(2)  # initial velocity
    qdotdot = np.zeros(2)  # initial acceleration
    ahat = np.zeros(4)

    for j in range(num_dts):

        # determine t
        t = dt * j / 100

        # calculate desired terms
        q_desired = np.array([1 - np.exp(-t), 2 - 2*np.exp(-t)])
        qdot_desired = np.array([np.exp(-t), 2.*np.exp(-t)])
        qdotdot_desired = -1. * qdot_desired

        # calculate tilde terms
        qtilde = q - q_desired
        qtildedot = qdot - qdotdot_desired

        # calculate r terms
        qdot_r = qdot_desired - np.matmul(lambd, qtilde)
        qdotdot_r = qdotdot_desired - np.matmul(lambd, qtildedot)

        # calculate s
        s = qdot - qdot_desired + np.matmul(lambd, qtilde)

        # calculate Y
        Y = calc_Y(q=q, qdot=qdot, qdot_r=qdot_r, qdotdot_r=qdotdot_r)

        # update adaptation law
        ahatdot = -np.matmul(gamma, np.matmul(Y.T, s))
        ahat += dt * ahatdot

        C = calc_C(q=q, qdot=qdot)
        H = calc_H(q=q)
        Hinv = np.linalg.inv(H)

        # calculate control law
        tau = np.matmul(Y, ahat) - np.matmul(K_d, s)

        # update system
        qdotdot = np.matmul(Hinv, tau - np.matmul(C, qdot))
        qdot += dt * qdotdot
        q += dt * qdot

        # record system
        control_history[i, j][:] = tau
        error_history[i, j][:] = q - q_desired


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True)
axes[0].set_ylabel(f'Error (q - q_d)')
for i in range(2):
    ax = axes[i]
    ax.set_title(f'Error in q_{i + 1} by Time (s)')
    ax.set_xlabel('Time (s)')
    for j, gamma in enumerate(gammas):
        ax.plot(
            dt * np.arange(num_dts) / 100,
            error_history[j, :, i],
            label=f'Gamma {j+1}')
    ax.legend()
plt.savefig('3b_error.jpg')


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True)
axes[0].set_ylabel(f'Torque')
for i in range(2):
    ax = axes[i]
    ax.set_title(f'Torque by Time (s)')
    ax.set_xlabel('Time (s)')
    for j, gamma in enumerate(gammas):
        ax.plot(
            dt * np.arange(num_dts) / 100,
            error_history[j, :, i],
            label=f'Gamma {j+1}')
    ax.legend()
plt.savefig('3b_torque.jpg')
