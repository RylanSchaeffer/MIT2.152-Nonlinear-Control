import matplotlib.pyplot as plt
import numpy as np


# switching window
# unknown true parameters
alpha_1 = 5.8
alpha_1_hat = 7
alpha_2 = 3.4
alpha_2_hat = 3.
d = 0.31
eta = 1.1
lamb = 5
dt = 0.01
num_dt = 10000

# create desired trajectory (constant velocity)
time = np.cumsum(np.full(shape=num_dt, fill_value=dt))
xdot_desired = np.full(shape=num_dt, fill_value=0.1)
xdotdot_desired = np.zeros_like(xdot_desired)
x_desired = np.cumsum(dt*xdot_desired)


def calculate_xdotdot(x, xdot, xdotdot_old):
    return f(x=x, xdot=xdot) + u(x=x, xdot=xdot)


def u(x, xdot, xdotdot):
    xtilde_dot =
    return -fhat(x=x, xdot=xdot) + xdotdot_desired - lamb*


def f(x, xdot):
    u = clacul
    prefactor = alpha_1 + alpha_2 * np.square(np.cos(x))
    return -prefactor * np.abs(xdot) * xdot + d


def fhat(x, xdot):
    prefactor = alpha_1_hat + alpha_2_hat * np.square(np.cos(x))
    return -prefactor * np.abs(xdot) * xdot



# create control law, calculate
x, xdot, xdotdot = 0., 0., 0.
xs, xdots = np.zeros(num_dt), np.zeros(num_dt)
for i in range(num_dt):
    xs[i] = x
    xdots[i] = xdot
    xdotdot = calculate_xdotdot(x=x, xdot=xdot)



fig, axes = plt.subplots(
    nrows=1,
    ncols=2)

axes[0].plot(time, x_desired, label='Desired')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('x')
axes[0].legend(numpoints=1, loc='best')

plt.show()