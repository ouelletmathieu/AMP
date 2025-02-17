


class Parameter_simulation():

    def __init__(self):
        self.minima = 0.05
        self.bounding_energy = 0.1
        self.K = (1,1)
        self.damping = 1
        self.max_time = 20000
        self.n_step = 4000000
        self.n_plot = 200
        self.n_molecule = 15
        self.box_side = 100
        self.ratio_distance_conformation = 1
        self.cut_off = self.minima*50