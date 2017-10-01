
from vispy import app
from vispy import gloo

import shaders
from Surface import Surface

class Canvas(app.Canvas):

    def __init__(self):
        self.width = 600
        self.height = 600

        app.Canvas.__init__(self, size=(self.width, self.height), title='Water surface simulator')

        gloo.set_state(clear_color=(0, 0, 0, 1), depth_test=False, blend=False)

        self.surface = Surface()

        self.program = gloo.Program(shaders.vert_shader, shaders.frag_shader)
        self.program['a_position'] = self.surface.position()

        self.activate_zoom()
        self.show()

    def activate_zoom(self):
        self.width, self.height = self.size
        gloo.set_viewport(0, 0, *self.physical_size)

    def on_draw(self, event):
        gloo.clear()
        self.program.draw('points')
