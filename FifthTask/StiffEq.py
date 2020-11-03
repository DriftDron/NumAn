import numpy as np
import pandas as pd
import math


def y1(x):
    return 0.00407478 * math.e ** (-247.204 * x) + 0.995925 * math.e ** (-0.795942 * x)


def y2(x):
    return -0.00404184 * math.e ** (-247.204 * x) + 1.00404 * math.e ** (-0.795942 * x)


print('Шаг h = 0.001')
x = 0.1
h = 0.1
xv = np.array([0])
y_1 = np.array([1.0])
y_2 = np.array([1.0])
for i in range(5):
    y_1 = np.append(y_1, y1(x))
    y_2 = np.append(y_2, y2(x))
    xv = np.append(xv, x)
    x += h
Y_math = np.array([y_1, y_2])
h = 0.001

t = 0.001
y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])
a = np.array([[-125, 123.2],
             [123.2, -123]])
w = np.identity(2) + a * h
l1 = np.linalg.eigh(w)[0][0]
l2 = np.linalg.eigh(w)[0][1]
if abs(l1) + abs(l2) < 2:
    print('|',l1,'|<1 and |', l2, '|<1 => Явный метод Эйлера устойчив')
else:
    print('Явный метод Эйлера неустойчив')
for i in range(499):
    t += h
    Ynext = np.dot(w, Y)
    Y = Ynext
    if (t % 0.1  < 10e-4):
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])
Y_eu1 = np.array([y_1, y_2])

w = np.linalg.inv(np.identity(2) - a * h)
l1 = np.linalg.eigh(w)[0][0]
l2 = np.linalg.eigh(w)[0][1]
print('Собственные числа обр. м. Эйлера: ', l1, l2)
t = 0.001
y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])
for i in range(499):
    t += h
    Ynext = np.dot(w, Y)
    Y = Ynext
    if (t % 0.1  < 10e-4):
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])
Y_eu2 = np.array([y_1, y_2])

t = 0.001
y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])
for i in range(499):
    t += h
    Ynext = np.dot(np.dot((np.identity(2) + h / 2 * a), np.linalg.inv((np.identity(2) - h / 2 * a))), Y)
    Y = Ynext
    if (t % 0.1  < 10e-4):
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])
Y_Ad = np.array([y_1, y_2])

d = {'t': xv, 'Точное-y*1': Y_math[0], 'Точное-y*2': Y_math[1],
     'М. Эйлера y*1-y1': Y_math[0]-Y_eu1[0], 'М. Эйлера y*2-y2': Y_math[1]-Y_eu1[1],
     'Обр. м. Эйлера y*1-y1': Y_math[0]-Y_eu2[0], 'Обр. м. Эйлера y*2-y2': Y_math[1]-Y_eu2[1],
     'М. Адамса y*1-y1': Y_math[0]-Y_Ad[0], 'М. Адамса y*2-y2': Y_math[1]-Y_Ad[1]}
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
df = pd.DataFrame(d)
print(df)
#------------------------
#------------------------
#------------------------
print('Шаг h = 0.05')
h = 0.05
y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])
a = np.array([[-125, 123.2],
             [123.2, -123]])
w = np.identity(2) + a * h
l1 = np.linalg.eigh(w)[0][0]
l2 = np.linalg.eigh(w)[0][1]
if abs(l1) + abs(l2) < 2:
    print('|',l1,'|<1 and |', l2, '|<1 => Явный метод Эйлера устойчив')
else:
    print(l1,' and ', l2, '=> Явный метод Эйлера неустойчив')
for i in range(10):
    Ynext = np.dot(w, Y)
    Y = Ynext
    if i % 2 == 1:
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])
Y_eu1 = np.array([y_1, y_2])

y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])

w = np.linalg.inv(np.identity(2) - a * h)
l1 = np.linalg.eigh(w)[0][0]
l2 = np.linalg.eigh(w)[0][1]
print('Собственные числа обр. м. Эйлера: ', l1, l2)
for i in range(10):
    Ynext = np.dot(w, Y)
    Y = Ynext
    if i % 2 == 1:
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])

Y_eu2 = np.array([y_1, y_2])

y_1 = np.array([1.0])
y_2 = np.array([1.0])
Y = np.array([y_1, y_2])
for i in range(10):
    t += h
    Ynext = np.dot(np.dot((np.identity(2) + h / 2 * a), np.linalg.inv((np.identity(2) - h / 2 * a))), Y)
    Y = Ynext
    if i % 2 == 1:
        y_1 = np.append(y_1, Y[0])
        y_2 = np.append(y_2, Y[1])
Y_Ad = np.array([y_1, y_2])


d = {'t': xv, 'Точное-y*1': Y_math[0], 'Точное-y*2': Y_math[1],
     'М. Эйлера y*1-y1': Y_math[0]-Y_eu1[0], 'М. Эйлера y*2-y2': Y_math[1]-Y_eu1[1],
     'Обр. м. Эйлера y*1-y1': Y_math[0]-Y_eu2[0], 'Обр. м. Эйлера y*2-y2': Y_math[1]-Y_eu2[1],
     'М. Адамса y*1-y1': Y_math[0]-Y_Ad[0], 'М. Адамса y*2-y2': Y_math[1]-Y_Ad[1]}
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
df = pd.DataFrame(d)
print(df)
