import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from scipy.integrate import quad


x, y = np.linspace(-2, 2, 100), np.linspace(-2, 2, 100)
xx, yy = np.meshgrid(x, y)
x_flat = xx.flatten()
y_flat = yy.flatten()

ke = 0.5 * np.power(y_flat, 2)


def integrand_d(x):
    return np.power(x, 5) - np.power(x, 3) * np.power(np.sin(x * np.pi / 2), 3)


def integrand_e(x):
    return x - np.sin(x * np.pi / 2)


pe = []
for x in x_flat:
    integral, error = quad(
        integrand_d,
        0,
        x)
    pe.append(integral)
pe = np.array(pe)
zz = ke + pe
zz = zz.reshape(xx.shape)


fig = go.Figure(data=[go.Surface(z=zz, x=xx, y=yy)])
fig.show()
