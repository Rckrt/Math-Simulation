
import numpy as np

from water_math.Surface import Surface


class Bed(object):
    def __init__(self):
        super().__init__()
        self.random_wave_surface = Surface((50, 50), 5, 0.3).normal(0)[:, :, 0]

    def new_random_surface(self):
        self.random_wave_surface = Surface((50, 50), 5, 0.3).normal(0)[:, :, 0]

    def bed_depths(self, shape):
        depths = np.ndarray([50, 50], dtype=np.float32)
        if shape == "straight":
            for i in range(depths.shape[0]):
                for j in range(depths.shape[1]):
                    if i <= 19:
                        depths[i][j] = 1
                    elif i >= 80:
                        depths[i][j] = 5
                    else:
                        depths[i][j] = (i - 20) * 4 / 58 + 1

        if shape == "random":
            depths = self.random_wave_surface

        if shape == "linspace":
            depths = np.linspace(-1.2, -0.6, 2500, dtype=np.float32).reshape([50, 50])

        if shape == "beach":
            depths = np.linspace(-0.3, 0.3, 2500, dtype=np.float32).reshape([50, 50])

        if shape == "normal":
            depths = np.ones((50, 50), dtype=np.float32) * -0.5

        if shape == "sky":
            for i in range(depths.shape[0]):
                for j in range(depths.shape[1]):
                    print((2 - ((i - 50) / 50) ** 2 - ((j - 50) / 50) ** 2))
                    depths[i][j] = ((20 - ((i - 50) / 50) ** 2 - (
                    (j - 50) / 50) ** 2) ** 0.5 + 0.0000000000001) - 4.36
            return depths

        return depths