import numpy as np
from surface.Surface import Surface


v = 10
delta = 0.1
sigma = 0.1
size = (50, 50)


# y(x+h) = y(x) + h/6 * (k1+2k2+2k3+k4)
# x is vector [h(t), v(t)]
def runge_kutta(f, x, h):
    k_1 = k1(f, x, h)
    k_2 = k2(f, x, h, k_1)
    k_3 = k3(f, x, h, k_2)
    k_4 = k4(f, x, h, k_3)
    return x + (h / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)


def k1(f, x, h):
    return f(x)


def k2(f, x, h, k_1):
    return f(x + h * k_1 / 2)


def k3(f, x, h, k_2):
    return f(x + h * k_2 / 2)


def k4(f, x, h, k_3):
    return f(x + h * k_3)


# f([h, v]) = [v, h'']
def f(x):
    up = x[0]
    down = x[1]
    return np.array([down, derivative(up)])


# h'' = Ldh
def derivative(heights):
    der_heights = np.zeros(size, dtype=np.float32)

    for i in range(0, size[0]):
        for k in range(0, size[1]):
            left  = heights[i][(k - 1 + size[1]) % size[1]]
            right = heights[i][(k + 1) % size[1]]
            up    = heights[(i - 1 + size[0]) % size[0]][k]
            down  = heights[(i + 1) % size[0]][k]
            this  = heights[i][k]

            der_heights[i][k] = ((v ** 2) * (sigma ** 2) / (delta ** 2)) * (left + right + up + down - 4 * this)
    return der_heights


def getHeights(surface, h, h_der):
    # first time
    if h is None:
        # return surface.height(0)
        h = np.ones(size, dtype=np.float32) * 0.2
        h[size[0] // 2, size[1] // 2] = 1
        h[size[0] // 2 + 1, size[1] // 2] = 1
        h[size[0] // 2, size[1] // 2 + 1] = 1
        h[size[0] // 2 + 1, size[1] // 2 + 1] = 1

        h_der = np.zeros(size, dtype=np.float32)

        return np.array([h, h_der])
    # other times
    else:
        return runge_kutta(f, np.array([h, h_der]), sigma)


def getNormal(heights):
    normal = np.zeros((size[0], size[1], 2), dtype=np.float32)
    for i in range(0, size[0]):
        for j in range(0, size[1]):
            left = heights[i][(j - 1 + size[1]) % size[1]]
            right = heights[i][(j + 1) % size[1]]
            up = heights[(i - 1 + size[0]) % size[0]][j]
            down = heights[(i + 1) % size[0]][j]

            normal[i][j][0] = (left + right) / (2 * delta)
            normal[i][j][1] = (up + down) / (2 * delta)
    return normal