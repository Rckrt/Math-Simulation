import random

import numpy as np

from water_math.render.methods.generator import HeightsGenerator


class Base:
    def __init__(self, method, v=1, delta=1, sigma=1000, size=(50, 50), max_height=0.1, min_height=-0.1, borders=False,
                 is_shallow=False):
        self.method = method
        self.v = v
        self.delta = delta
        self.sigma = sigma
        self.size = size
        self.max_height = max_height
        self.min_height = min_height
        self.borders = borders
        self.is_shallow = is_shallow
        self.g = 9.81
        self.generator = HeightsGenerator(size)

        self.methods = {
            "peak": self.generator.peak,
            "bubble": self.generator.bubble,
            "vertical": self.generator.vertical,
            "shallow": self.generator.shallow
        }

    def init(self, part=4):
        params = {
            "max_height": self.max_height,
            "min_height": self.min_height,
            "part": part
        }
        print(self.max_height)
        return self.methods[self.method](params)

    # f([h, v]) = [v, h'']
    def f(self, x):
        return np.array([x[1], self.derivative(x[0])])

    # h'' = Ldh
    # !!! indices arrays: 4 arrays
    # No cycles. Indices rolling

    def derivative(self, heights):
        der_heights = np.zeros(self.size, dtype=np.float32)

        if not self.borders:
            for i in range(0, self.size[0]):
                for k in range(0, self.size[1]):
                    left = heights[i][(k - 1 + self.size[1]) % self.size[1]]
                    right = heights[i][(k + 1) % self.size[1]]
                    up = heights[(i - 1 + self.size[0]) % self.size[0]][k]
                    down = heights[(i + 1) % self.size[0]][k]
                    this = heights[i][k]
                    der_heights[i][k] = ((self.v ** 2) * (self.sigma ** 2) / (self.delta ** 2)) * (
                    left + right + up + down - 4 * this)

        else:
            der_heights[1:-1, 1:-1] += ((self.v ** 2) * (self.sigma ** 2) / (self.delta ** 2)) * (
            heights[2:, 1:-1] + heights[0:-2, 1:-1] + heights[1:-1, 2:] + heights[1:-1, 0:-2] - 4 * heights[1:-1, 1:-1])

        return der_heights

    def get_heights(self, h_desc):
        pass

    def get_normal(self, heights):
        normal = np.zeros((self.size[0], self.size[1], 2), dtype=np.float32)
        for i in range(0, self.size[0]):
            for j in range(0, self.size[1]):
                left = heights[i][(j - 1 + self.size[1]) % self.size[1]]
                right = heights[i][(j + 1) % self.size[1]]
                up = heights[(i - 1 + self.size[0]) % self.size[0]][j]
                down = heights[(i + 1) % self.size[0]][j]

                normal[i][j][0] = (left + right) / (2 * self.delta)
                normal[i][j][1] = (up + down) / (2 * self.delta)

        return normal

    def UVg(self, U, h):
        res = np.power(U, 2) / h + self.g * np.power(h, 2) / 2
        return res

    def UVh(self, U, V, h):
        res = np.divide(np.multiply(U, V), h)
        return res

    def f_shallow(self, x):
        h = np.clip(x[0], self.min_height, self.max_height)

        U = np.clip(x[1], 0, self.v)
        V = np.clip(x[2], 0, self.v)

        h_d = np.roll(h, 1, axis=0)
        h_u = np.roll(h, -1, axis=0)
        h_l = np.roll(h, -2, axis=1)
        h_r = np.roll(h, 2, axis=1)

        U_d = np.roll(U, 1, axis=0)
        U_u = np.roll(U, -1, axis=0)
        U_l = np.roll(U, -2, axis=1)
        U_r = np.roll(U, 2, axis=1)

        V_d = np.roll(V, 1, axis=0)
        V_u = np.roll(V, -1, axis=0)
        V_l = np.roll(V, -2, axis=1)
        V_r = np.roll(V, 2, axis=1)

        der_U = 4*U - U_d - U_u - U_l - U_r
        der_V = 4*V - V_d - V_u - V_l - V_r


        der_Ug = 4*self.UVg(U, h) - self.UVg(U_r, h_r) - self.UVg(U_l, h_l) - self.UVg(U_u, h_u) - self.UVg(U_d, h_d)
        der_Vg = 4*self.UVg(V, h) - self.UVg(V_r, h_r) - self.UVg(V_l, h_l) - self.UVg(V_u, h_u) - self.UVg(V_d, h_d)
        der_Uh = 2*self.UVh(U, V, h) - self.UVh(U_d, V_d, h_d) - self.UVh(U_u, V_u, h_u)
        der_Vh = 2*self.UVh(U, V, h) - self.UVh(U_r, V_r, h_r) - self.UVh(U_l, V_l, h_l)

        return np.array([-(der_U + der_V),
                         -(der_Ug + der_Uh),
                         -(der_Vg + der_Vh)]) / self.delta
