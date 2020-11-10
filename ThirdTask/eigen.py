import numpy as np

a = np.array([[-1.53698, -0.19907, 0.95855],
              [-0.19907, 1.17742, 0.06992],
              [0.95855, 0.06992, -1.55151]])
ik, jk = np.unravel_index(np.argmax(np.absolute(np.triu(a, 1))), a.shape)
amax = np.max(np.absolute(np.triu(a, 1)))
eps = 1e-6
k = 0
x = np.identity(3)
while amax > eps and k != 10:
    d = ((a[ik, ik] - a[jk, jk]) ** 2 + 4 * a[ik, jk] ** 2) ** (1 / 2)
    c = (1 / 2 * (1 + abs(a[ik, ik] - a[jk, jk]) / d)) ** (1 / 2)
    part_val = np.abs(a[ik][ik] - a[jk][jk]) / d
    s = np.sign(a[ik, jk] * (a[ik, ik] - a[jk, jk])) * (0.5 * (1 - abs(a[ik, ik] - a[jk, jk]) / d)) ** (1 / 2)
    v = np.identity(3)
    v[ik, ik] = c
    v[jk, jk] = c
    v[ik, jk] = -s
    v[jk, ik] = s
    x = np.dot(x, v)
    anext = np.copy(a)
    for i in range(a.shape[0]-1):
        for j in range(a.shape[0]-1):
            if i != ik and i != jk and j != jk and j != ik:
                anext[i, j] = a[i, j]
            elif i != ik and i != jk:
                anext[i, ik] = c * a[i, ik] + s * a[i, jk]
                anext[i, jk] = - s * a[i, ik] + c * a[i, jk]
    anext[ik, ik] = c ** 2 * a[ik, ik] + 2 * c * s * a[ik, jk] + s ** 2 * a[jk, jk]
    anext[jk, jk] = s ** 2 * a[ik, ik] - 2 * c * s * a[ik, jk] + c ** 2 * a[jk, jk]
    anext[jk, ik] = anext[ik, jk] = 0
    a = np.copy(anext)
    ik, jk = np.unravel_index(np.argmax(np.absolute(np.triu(a, 1))), a.shape)
    amax = np.max(np.absolute(np.triu(a, 1)))
    k += 1
print('Собственные числа методом Якоби:', np.diag(a))
print('Матрица собственных векторов:\n', x)

a = np.array([[-1.53698, -0.19907, 0.95855],
              [-0.19907, 1.17742, 0.06992],
              [0.95855, 0.06992, -1.55151]])
lmax = np.max(np.absolute(np.linalg.eigh(a)[0]))
v = [1, 0, 0]
q = [v]
eps = 10 ** (-3)
lam1, lam2 = 1, 0
while np.abs(lam1 - lam2) > eps:
    q.append(np.dot(a, q[i - 1]))
    lam2 = np.copy(lam1)
    lam1 = q[i][0] / q[i - 1][0]
    i += 1
w_max = lam1
print("Модуль с.ч. степенным методом", abs(w_max))
print("Модуль макс с.ч. встроенным методом", lmax)
print("Абсолютная погрешность", abs(abs(w_max) - lmax))
q[i - 1] = q[i - 1] / np.linalg.norm(q[i - 1])
print("Соответствуйющий собственный вектор", q[i - 1])

v = [1, 0, 0]
q = [v]
lambdas = [v]
i = 1
eps = 10 ** (-3)
lam1, lam2 = 1, 0
while np.abs(lam1 - lam2) > eps:
    q.append(np.dot(a, q[i - 1]))
    lam2 = np.copy(lam1)
    lam1 = np.dot(q[i], q[i-1])/np.dot(q[i-1], q[i-1])
    i += 1

w_max = lam1
print("Модуль макс с.ч. методом скалярных произведений", abs(w_max))
print("Модуль макс с.ч. встроенным методом", lmax)
print("Абсолютная погрешность", abs(abs(w_max) - lmax))








