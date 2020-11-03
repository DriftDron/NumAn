import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.integrate import odeint

def f(x, y):
    return math.cos(y) / (1 + x) - 0.5 * y

def model(y,t):
    dydt = math.cos(y) / (1 + t) - 0.5 * y
    return dydt

# initial condition
y0 = 0
# time points
t = np.linspace(0,0.5, 11)
# solve ODE
y = odeint(model,y0,t)
y_math = y[:,0]
plt.plot(t,y)
plt.xlabel('x')
plt.ylabel('y(x)')

h = 0.05
a = 0
b = 0.5
n = (b - a) / h
y = 0
x = 0
y_h = np.array([])
x_m = np.array([])
for i in range(int(n) + 1):
    x_m = np.append(x_m, x)
    y_h = np.append(y_h, y)
    ym = y + h * f(x + h/2, y + h/2 * f(x, y))
    x += h
    y = ym
plt.plot(x_m, y_h, color = 'red')

y = 0
x = 0
h = 0.05 / 2
y_h2 = np.array([])
x_m = np.array([])
for i in range(int(2 * n) + 1):
    x_m = np.append(x_m, x)
    y_h2 = np.append(y_h2, y)
    ym = y + h * f(x + h/2, y + h/2 * f(x, y))
    x += h
    y = ym

plt.plot(x_m, y_h2, color = 'green')

y_h2 = y_h2[::2]
r = (y_h2 - y_h) / (2 ** 2 - 1)
y_rev = y_h2 + r

x_m = x_m[::2]
plt.plot(x_m, y_rev, color = 'yellow')

plt.plot(x_m, y_rev - y_math, color = 'brown')
d = {"x ": x_m, "y_math ": y_math, "y_h  ": y_h, "y_h2  ": y_h2, "y_rev ": y_rev, "y_rev-y_math": y_rev - y_math   }
df = pd.DataFrame(d)
print(df)

y_RK = np.array([])
x = 0
y = 0
h = 0.05
x_RK = np.array([])
for i in range(int(n) + 1):
    y_RK = np.append(y_RK, y)
    x_RK = np.append(x_RK, x)
    k1 = h * f(x, y)
    k2 = h * f(x + h / 2, y + k1 / 2)
    k3 = h * f(x + h / 2, y + k2 / 2)
    k4 = h * f(x + h, y + k3)
    ynext = y + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    x += h
    y = ynext


plt.plot(x_RK, y_RK, color = 'pink')
plt.show()

y_Ad_ex = np.zeros([11])
h = 0.05
x_Ad_ex = np.zeros([11])
for i in range(5):
    y_Ad_ex[i] = y_math[i]
    x_Ad_ex[i] = x_m[i]
y = y_Ad_ex[4]
x = x_Ad_ex[4]

for i in range(5, int(n) + 1):
    q = h * f(x_Ad_ex[i-1], y_Ad_ex[i-1])
    q1 = h * f(x_Ad_ex[i-2], y_Ad_ex[i-2])
    q2 = h * f(x_Ad_ex[i-3], y_Ad_ex[i-3])
    q3 = h * f(x_Ad_ex[i-4], y_Ad_ex[i-4])
    q4 = h * f(x_Ad_ex[i-5], y_Ad_ex[i-5])
    ynext = y + (1901 * q - 2774 * q1 + 2616 * q2 - 1274 * q3 + 251 * q4) / 720
    x += h
    y = ynext
    y_Ad_ex[i] = y
    x_Ad_ex[i] = x

y_Ad_in = np.zeros([11])
x_Ad_in = np.zeros([11])
for i in range(5):
    y_Ad_in[i] = y_math[i]
    x_Ad_in[i] = x_m[i]
y = y_Ad_in[4]
x = x_Ad_in[4]
q0 = h * f(x_Ad_ex[4], y_Ad_ex[4])

for i in range(5, int(n) + 1):
    q = h * f(x_Ad_in[i-1], y_Ad_in[i-1])
    q1 = h * f(x_Ad_in[i-2], y_Ad_in[i-2])
    q2 = h * f(x_Ad_in[i-3], y_Ad_in[i-3])
    q3 = h * f(x_Ad_in[i-4], y_Ad_in[i-4])
    ynext = y + (251 * q0 + 646 * q - 264 * q1 + 106 * q2 - 19 * q3) / 720
    x += h
    y = ynext
    y_Ad_in[i] = y
    x_Ad_in[i] = x
    q0 = h * f(x_Ad_in[i], y_Ad_in[i])


d2 = {"y_math": y_math, "y_math-y_RK": y_math-y_RK, "y_math-y_Ad_ex": y_math-y_Ad_ex, "y_math-y_Ad_in": y_math-y_Ad_in}
df2 = pd.DataFrame(d2)
print(df2)

