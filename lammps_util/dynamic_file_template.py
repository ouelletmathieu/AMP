import pickle

class Save_file_full_sim():

    def __init__(self, key,max_time, n_step, n_plot, n_molecule,box_side, n_node,n_out  ):
        self.key = key
        self.minima = None
        self.bounding_energy = None
        self.K = None
        self.damping = None
        self.max_time = max_time
        self.n_step = n_step
        self.n_plot = n_plot
        self.n_molecule = n_molecule
        self.box_side = box_side
        self.n_node = n_node
        self.n_out = n_out
    

        #default value
        self.ratio_distance_conformation = 1
        self.cut_off = None
        
        #computed value
        self.temp_list_H = []
        self.conn_list_H = []
        self.bounded_list_H = [] 
        self.count_list_H = []
        
        self.temp_list = []
        self.conn_list = []
        self.bounded_list = []
        self.count_list = []

        self.bind_energy_rc = []
        self.bind_energy_energy = []

    def set_physics_param(self, minima, bounding_energy, K, damping):
        self.minima = minima
        self.bounding_energy = bounding_energy
        self.K = K
        self.damping = damping
        self.cut_off = minima*5
        
    def add_healthy_sim(self, temp_list_H, conn_list_H ,bounded_list_H ,count_list_H):
        self.temp_list_H.extend(temp_list_H)
        self.conn_list_H.extend(conn_list_H)
        self.bounded_list_H.extend(bounded_list_H)
        self.count_list_H.extend(count_list_H)
    
    def set_binding_energy(self, rc, ener):
        self.bind_energy_rc = rc
        self.bind_energy_energy = ener

    def add_binded_sim(self, temp_list, conn_list, bounded_list, count_list):
        self.temp_list.extend(temp_list)
        self.conn_list.extend(conn_list)
        self.bounded_list.extend(bounded_list)
        self.count_list.extend(count_list)

    def save(self, name=""):
        
        with open(f'{self.key}_{name}.data', 'wb') as save_file:
            pickle.dump(self, save_file)

    def load(path):
        with open(path, 'rb') as save_file:
            return pickle.load(save_file)