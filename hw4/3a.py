# Slotine and Li, Example 9.1
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


def calc_tau(q, qdot):
    tau = -np.matmul(K_p, q - q_desired) - np.matmul(K_d, qdot)
    return tau


# simulation parameters
num_dts = 2000
control_history = np.zeros(shape=(num_dts, 2))
error_history = np.zeros(shape=(num_dts, 2))

dt = 0.004
q_desired = np.array([1, 2])  # desired angles, radians
q = np.zeros(2)  # initial angles, radians
qdot = np.zeros(2)  # initial velocity
qdotdot = np.zeros(2)  # initial acceleration
K_p, K_d = 100*np.eye(2), 20*np.eye(2)  # controller positive-definite matrices

for i in range(num_dts):
    C = calc_C(q=q, qdot=qdot)
    H = calc_H(q=q)
    Hinv = np.linalg.inv(H)
    tau = calc_tau(q=q, qdot=qdot)
    qdotdot = np.matmul(Hinv, tau - np.matmul(C, qdot))
    qdot += dt * qdotdot
    q += dt * qdot
    control_history[i][:] = tau
    error_history[i][:] = q - q_desired


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True)
axes[0].set_ylabel(f'Error (q - q_d)')
for i in range(2):
    ax = axes[i]
    ax.set_title(f'Error in q_{i+1} by Time (s)')
    ax.set_xlabel('Time (s)')
    ax.plot(
        dt * np.arange(num_dts),
        error_history[:, i])
plt.savefig('3a_error.jpg')


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True)
axes[0].set_ylabel(f'Control Torque')
for i in range(2):
    ax = axes[i]
    ax.set_title(f'Control Torque by Time (s)')
    ax.set_xlabel('Time (s)')
    ax.plot(
        dt * np.arange(num_dts),
        control_history[:, i])
plt.savefig('3a_torque.jpg')
