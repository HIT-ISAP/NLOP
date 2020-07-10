import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

name = 'LevenbergMarquardt'

def f(x, y):
    return (1 - x) ** 2 + 100 * (y - x * x) ** 2

n = 256

x = np.linspace(-1, 1.1, n)
y = np.linspace(-0.1, 1.1, n)

X, Y = np.meshgrid(x, y)

plt.figure()

plt.contourf(X, Y, f(X, Y), 5, alpha=0, cmap=plt.cm.hot)

C = plt.contour(X, Y, f(X, Y), 8, locator=ticker.LogLocator(), colors='black', linewidth=0.01)

plt.clabel(C, inline=True, fontsize=10)
# ---------------------

x1 = []
x2 = []
with open('../data/' + name + '.txt') as f:
    for line in f.readlines():
        temp = line.split()
        x1.append(float(temp[0]))
        x2.append(float(temp[1]))



plt.plot(x1[0], x2[0], marker='o', color='green', markersize=10)
plt.plot(x1[len(x1)-1], x2[len(x2)-1], marker='o', color='red', markersize=10)

plt.plot(x1, x2, color='orange')

plt.xlabel('x1')
plt.ylabel('x2')
plt.title(name)
#plt.legend()
plt.savefig("./picture/" + name)
plt.show()

