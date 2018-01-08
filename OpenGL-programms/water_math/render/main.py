
from vispy import app

from water_math.render.Canvas import Canvas
from water_math.surface.CircularWavesSurface import CircularWavesSurface

if __name__ == '__main__':
    c = Canvas(CircularWavesSurface(max_height=0.05), size=(700, 700))
    app.run()
