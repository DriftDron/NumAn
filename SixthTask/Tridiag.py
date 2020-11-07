import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def p(x):
    return 1/(x+3)


def q(x):
    return -x


def r(x):
    return np.log(2 + x)


def A(x, h):
    return - p(x) / h ** 2 - q(x) / (2 * h)


def B(x, h):
    return - 2 * p(x) / h ** 2 - r(x)


def C(x, h):
    return - p(x) / h ** 2 + q(x) / (2 * h)


def G(x):
    return 1 - x/2

print('O(h), n = 20')
a1, a2, b1, b2, a, b = 0, -1, 1/2, 1, 0, 0
h = 0.1

x = np.arange(-1, 1.1, h)
s = np.empty(x.shape[0])
s[0] = (-a2 / h) / (-a1 - a2 / h)
t = np.empty(x.shape[0])
t[0] = 0
y = np.empty(x.shape[0])
for i in range(1, x.shape[0]-1):
    s[i] = C(x[i], h) / (B(x[i], h) - A(x[i], h) * s[i - 1])
    t[i] = (A(x[i], h) * t[i-1] - G(x[i])) / (B(x[i], h) - A(x[i], h) * s[i-1])
s[20] = C(x[20], h) / (-b1 - b2 / h - (-b2 / h) * s[20 - 1])
t[20] = ((-b2 / h) * t[20 - 1] - b) / ((-b1 - b2 / h) - (-b2 / h) * s[20 - 1])
y[y.shape[0]-1] = t[t.shape[0]-1]

for i in range(x.shape[0] - 2, -1, -1):
    y[i] = s[i] * y[i+1] + t[i]

y_20 = np.array([2.9275, 2.9124, 2.8671, 2.7849, 2.6881, 2.5395, 2.4048, 2.2330, 2.0476,
                 1.8545, 1.6598, 1.5127, 1.2893, 1.1239, 0.9767, 0.8650, 0.7444, 0.6655,
                 0.5944, 0.5467, 0.5142])

plt.plot(x, y_20)
plt.plot(x, y, color = 'pink')

df = pd.DataFrame({'i': np.arange(0, y_20.shape[0], 1),'x': x,
                   'A': np.array([A(x[i], h) for i in range(y_20.shape[0])]),
                   'B': np.array([B(x[i], h) for i in range(y_20.shape[0])]),
                   'C': np.array([C(x[i], h) for i in range(y_20.shape[0])]),
                   'G': np.array([G(x[i]) for i in range(y_20.shape[0])]),
                   's': s,'t': t,'y_ex': y_20,'y': y, '|y-y_ex|': np.absolute(y-y_20)})
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
print(df)

print('O(h), n = 10')
h = 0.2
x = np.arange(-1, 1.1, h)
s = np.empty(x.shape[0])
s[0] = (-a2 / h) / (-a1 - a2 / h)
t = np.empty(x.shape[0])
t[0] = 0
y = np.empty(x.shape[0])
for i in range(1, x.shape[0]-1):
    s[i] = C(x[i], h) / (B(x[i], h) - A(x[i], h) * s[i - 1])
    t[i] = (A(x[i], h) * t[i-1] - G(x[i])) / (B(x[i], h) - A(x[i], h) * s[i-1])
s[10] = C(x[10], h) / (-b1 - b2 / h - (-b2 / h) * s[10 - 1])
t[10] = ((-b2 / h) * t[10 - 1] - b) / ((-b1 - b2 / h) - (-b2 / h) * s[10 - 1])
y[y.shape[0]-1] = t[t.shape[0]-1]

for i in range(x.shape[0] - 2, -1, -1):
    y[i] = s[i] * y[i+1] + t[i]

df = pd.DataFrame({'i': np.arange(0, y.shape[0], 1),'x': x,
                   'A': np.array([A(x[i], h) for i in range(y.shape[0])]),
                   'B': np.array([B(x[i], h) for i in range(y.shape[0])]),
                   'C': np.array([C(x[i], h) for i in range(y.shape[0])]),
                   'G': np.array([G(x[i]) for i in range(y.shape[0])]),
                   's': s,'t': t, 'y': y,'y_ex': y_20[::2],'|y-y_ex|': np.absolute(y-y_20[::2])})

print(df)
plt.plot(x, y, color = 'brown')
plt.show()