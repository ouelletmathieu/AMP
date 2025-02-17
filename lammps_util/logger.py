import itertools
import math
from .util import Util
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from .sim_param import Parameter_simulation
import pickle

class Logger:
    """class that allows to store information about simulation from lammps 
    """
    
    def __init__(self, my_lammps, struct_to_compare_pos = [], id_struct_to_compare = [], is_2d = True, nb_molecule = 1, n_node = None ):
            
        if len(struct_to_compare_pos)!=len(id_struct_to_compare):
            print("id struct to compare should match")
            raise ValueError("id struct to compare should match")
            
        #lammps instance
        self.L = my_lammps
        self.is_2d = is_2d
        #store atom position in the format [ [[x1(t1), x2(t1), ..., xn(t1)], ..., x1(tN), x2(tN), ..., xn(tN)], 
        #                                    [[y1(t1), y2(t1), ..., yn(t1)], ..., y1(tN), y2(tN), ..., yn(tN)]] 
        self.struct_list = [[], []]
        if not is_2d:
            self.struct_list = [[], [], []]
        #time for each entry stored 
        self.time_list = []
        #store the list of angle for each time, format [ [a1(t1),..., an(t1)], ...., [a1(tN),..., an(tN)]]
        self.angle_list = [] 
        #store the list of potential energy for each time, format:  [ e(t1), ...., e(tN) ]
        self.energy_list = []
        #store the list of kinetic energy for each time, format:  [ e(t1), ...., e(tN) ]
        self.energyKin_list = []
        #identifiant of the structure to compare
        self.id_struct_to_compare = id_struct_to_compare
        #list of expected angle for each of structure to compare 
        self.angle_list_to_compare = {id_struct: Util.get_angles(pos_struct, len(pos_struct)) for pos_struct, id_struct in zip(struct_to_compare_pos,id_struct_to_compare)}
        #nb of molecule in the simulation 
        self.nb_molecule = nb_molecule
        #nb of node per molecule
        self.n_node = n_node
        #parameter used for the simulation
        self.parameter_simulation = None 

        #result for the analysis
        self.dyn_connection = None
        self.dyn_binding = None
        self.dyn_count = None

        #connection to outside node 
        self.conn_out = None
        #connection to inside node
        self.conn_in = None

        if nb_molecule!=1 and n_node is None:
            print("if more than one molecule you need to set the number of node (n_node) per molecule")
            raise Exception("n_node not set")

    def set_connection(self, conn_in, conn_out):
        self.conn_in=conn_in
        self.conn_out=conn_out
        
    def append_energy(self):
        self.energy_list.append(self.L.getEnergy())
        self.energyKin_list.append(self.L.getEnergyKin())
         
    def append_time(self, time):
        self.time_list.append(time)
        
    def append_position(self):

        if self.nb_molecule==1:
            #set up the position

            self.active_x_pos_list = [self.L.lammps().atoms[i].position[0] for i in range(len(self.L.lammps().atoms))]
            self.active_y_pos_list = [self.L.lammps().atoms[i].position[1] for i in range(len(self.L.lammps().atoms))]

            
            self.struct_list[0].append(self.active_x_pos_list)
            self.struct_list[1].append(self.active_y_pos_list)

            if not self.is_2d: 
                self.active_z_pos_list = [self.L.lammps().atoms[i].position[2] for i in range(len(self.L.lammps().atoms))]
                self.struct_list[2].append(self.active_z_pos_list)
        else:
            atom_pos_lammps = self.L.lammps().atoms
            self.active_x_pos_list = []
            self.active_y_pos_list = []

            for mol_id in range(self.nb_molecule):
                pos_x = [ atom_pos_lammps[n_id+self.n_node*mol_id].position[0] for n_id in range(self.n_node) ]
                pos_y = [atom_pos_lammps[n_id+self.n_node*mol_id].position[1] for n_id in range(self.n_node) ]
                self.active_x_pos_list.append(pos_x)
                self.active_y_pos_list.append(pos_y)
            
            self.struct_list[0].append(self.active_x_pos_list)
            self.struct_list[1].append(self.active_y_pos_list)
   
    def append_angle(self):
        active_angle = Util.get_angles_lammps(self.L.lammps().atoms, len(self.L.lammps().atoms)) 
        self.angle_list.append(active_angle)
                  
    def log(self, time, is_angle = True):
        self.append_time(time)
        self.append_energy()
        self.append_position()
        if is_angle:
            self.append_angle()
          
    def get_last_distance(self, index, id_to_compare = None):
        return sum(
            (angle1 - angle2) ** 2
            for angle1, angle2 in zip(
                self.angle_list_to_compare[id_to_compare], self.angle_list[index]
            )
        )

        
    ###############################################
    #             PLOTTING  (1 molecule)          #
    ###############################################
    
    
    def plot_positions(self, time_list = None, index_list = None):
        
        if time_list is None:
            time_list = []
        if index_list is None:
            index_list = []
        if len(time_list)!=0:
            print("not implemented yet")
            raise NotImplementedError("not implemented yet")

        if self.is_2d:
            if len(index_list)!=0:

                for i in range(len(index_list)):
                    plt.scatter(self.struct_list[0][index_list[i]],self.struct_list[1][index_list[i]])

                plt.show()
        else:

            fig = plt.figure()
            ax = plt.axes(projection='3d')

            for i in range(len(index_list)):
                ax.scatter3D(self.struct_list[0][index_list[i]],self.struct_list[1][index_list[i]], self.struct_list[2][index_list[i]])

            plt.show()

    def plot_angle_difference(self, id_to_compare = None ):
        """ plot all angle difference with a given conformation specified by id_to_compare
            the conformation need to be known in advance by the logger object 
        """
        angle_dif_formatted = self._get_angle_difference(id_to_compare)

        for i in range(len(angle_dif_formatted)):
            plt.plot(self.time_list, angle_dif_formatted[i])

        if id_to_compare is None:
            plt.title("angle difference with t=0" , fontsize=17)
        else:
            plt.title(f"angle difference conf= {str(id_to_compare)}", fontsize=17)
        plt.xlabel('time (time)', fontsize=15)
        plt.ylabel('rad difference', fontsize=15)
        plt.show()    
        
    def plot_angular_sum_distance(self, id_to_compare_list = None):
        """compute and display the mean square angular displacement
        """
        
        if id_to_compare_list is None:
            id_to_compare_list = []
        list_of_list_of_distance = {id_:[] for id_ in id_to_compare_list}

        for id_tcp in id_to_compare_list:

            lst = list_of_list_of_distance[id_tcp]

            angle_formatted = self._get_angle_difference( id_to_compare = id_tcp )

            for t in range(len(self.time_list) ):

                total_dist = sum(
                    (angle_formatted[i][t]) ** 2
                    for i in range(len(angle_formatted))
                )


                lst.append(total_dist)

            plt.plot(self.time_list, lst, label=str(id_tcp))

        plt.title("mean square angular displacement" , fontsize=17)
        plt.xlabel('time (time)', fontsize=15)
        plt.ylabel('sqared rad', fontsize=15)
        plt.show()  
        
    def plot_energy(self, index_list = None):
        
        if index_list is None:
            index_list = []
        if len(index_list)==0:
            plt.plot(self.time_list, self.energy_list,  label='potential')
            plt.plot(self.time_list, self.energyKin_list,  label='kinetic')
            plt.plot(self.time_list, [ a+b for a,b in zip(self.energy_list,self.energyKin_list)],  label='all')
        else:
            plt.plot([self.time_list[i] for i in index_list], [ self.energy_list[i]   for i in index_list],  label='potential')
            plt.plot([self.time_list[i] for i in index_list], [ self.energyKin_list[i]   for i in index_list], label='kinetic')

        plt.title("energy per bead" , fontsize=17)
        plt.xlabel('time (time)', fontsize=15)
        plt.ylabel('energy (yes)', fontsize=15)
        plt.legend(loc="upper left")
        plt.show() 

    #########################################
    #           dynamic property            #
    #########################################

    def log_simulation_param(self, parameter_simulation: Parameter_simulation):
        self.parameter_simulation = parameter_simulation

    def log_dynamic(self):

        if self.parameter_simulation is None:
            print("simulation parameter need to be set")
            raise Exception("simulation parameter need to be set")

        # storing connection
        out_conn, in_conn = [], []
        # stroting data about binding
        self_bonded, outer_bounded = [], []
        binding_dict = []
        specific_binded_number = []
        # storing data about conformation
        count_healthy, count_prion, count_other = [], [], []
        count_healthy_even, count_prion_even, count_other_even = [], [], []
        conf_mol_list = []
        distance_list = []

        # for each frame
        for frame in range(len(self.struct_list[0])):

            # compute position and data related to the frame
            frame_in_conn, frame_out_conn = self.get_connection_position_molecule(frame)
            #frame_pos_list, frame_pos_per_mol_list, frame_mol_id_list = main_logger.get_position_molecule(frame)
            out_conn.append(frame_out_conn)
            in_conn.append(frame_in_conn)

            # compute binding
            binding = self.compute_frame_binding(frame, self.parameter_simulation.minima)
            specific_binded_number.append(binding[0])
            binding_dict.append(binding[1])
            self_bonded.append(binding[2])
            outer_bounded.append(binding[3])

            # compute angles
            frame_count = self.compute_frame_angle(frame, self.parameter_simulation.ratio_distance_conformation)
            count_healthy.append(frame_count[0])
            count_prion.append(frame_count[1])
            count_other.append(frame_count[2])
            count_healthy_even.append(frame_count[3])
            count_prion_even.append(frame_count[4])
            count_other_even.append(frame_count[5])
            conf_mol_list.append(frame_count[6])
            distance_list.append(frame_count[7])

            self.dyn_connection =  (out_conn, in_conn)
            self.dyn_binding =  (specific_binded_number,binding_dict, self_bonded, outer_bounded) 
            self.dyn_count = (count_healthy, count_prion, count_other, count_healthy_even, count_prion_even, count_other_even, conf_mol_list, distance_list)

    def get_bond_energy(self, protein_pair, length_list):
        

        self.bond_energy = []

        for frame in range(len(self.struct_list[0])):
            mol_bond_length = []
            
            for mol_id in range(self.nb_molecule):
                conn_bond_length = []
                pos_mol = [ (self.struct_list[0][frame][mol_id][k], self.struct_list[1][frame][mol_id][k] ) for k in range(self.n_node)]

                for conn, init_l in  zip(protein_pair.connection,length_list):
                
                    xy1 = pos_mol[conn[0]]
                    xy2 = pos_mol[conn[1]]
                    l = math.sqrt((xy1[0]-xy2[0])**2 + (xy1[1]-xy2[1])**2)
                    conn_bond_length.append((conn,(init_l-l)**2))

                mol_bond_length.append(conn_bond_length)
        
            self.bond_energy.append(mol_bond_length)

        return self.bond_energy

    def get_dynamic(self):
        return self.dyn_connection, self.dyn_binding, self.dyn_count 

    #########################################
    #                 save                  #
    #########################################


    def save_for_blender(self, path, protein_pair, length_list):

        dict_data = {"n_molecule": self.nb_molecule, "n_node":self.n_node, 
        "conn_out":self.conn_out, "conn_in":self.conn_in, "struct_list":self.struct_list}

        self.get_bond_energy(protein_pair, length_list)
        dict_data["bond_energy"] = self.bond_energy
        dict_data["distance_list"] = self.dyn_count[7]
        dict_data["binding_list"] = self.dyn_binding[1]

        with open(f'{path}.pickle', 'wb') as handle:
            pickle.dump(dict_data, handle, protocol=4)
     

    #########################################
    #                 getter                #
    #########################################
    
    def _get_angle_difference(self, id_to_compare = None ):
        
        angle_n = len(self.angle_list[0])
        angle_formatted = [[] for _ in range(angle_n)]

        for i, j in itertools.product(range(len(self.angle_list)), range(angle_n)):
            if id_to_compare is None:
                angle_formatted[j].append(self.angle_list[i][j] - self.angle_list[0][j] )
            else:
                angle_formatted[j].append(self.angle_list[i][j] - self.angle_list_to_compare[id_to_compare][j] )

        return angle_formatted

    def get_RMS_angle_difference_all_time(self):
            
        angle_n = len(self.angle_list[0])
        angle_sum = 0
        nb_angle = 0

        for i, j in itertools.product(range(len(self.angle_list)), range(angle_n)):
            nb_angle+=1
            angle_sum += (self.angle_list[i][j] - self.angle_list[0][j])**2 

        return math.sqrt(angle_sum/nb_angle)

    def check_multiple_molecule(self):
        if self.nb_molecule ==1 :
            print("this method was not coded for the case n=1")
            raise Exception("TODO code for n=1")

    def get_connection_position_molecule(self, index):

        if self.conn_in is None or self.conn_out is None:
            print("conn_in and conn_out need to be set before calling this method")
            raise Exception("conn_in and conn_out is not set")

        self.check_multiple_molecule()
   
        frame_in_conn, frame_out_conn = [], []

        # for each molecule
        for mol_id in range(self.nb_molecule):
            mol_frame_out_conn, mol_frame_in_conn = [], []
            pos_mol = [ (self.struct_list[0][index][mol_id][k], self.struct_list[1][index][mol_id][k] ) for k in range(self.n_node)]
    
            for c, t in itertools.chain(itertools.product(self.conn_out, {"out"}),  itertools.product(self.conn_in, {"in"})):
                x1, x2 = pos_mol[c[0]][0], pos_mol[c[1]][0]
                y1, y2 = pos_mol[c[0]][1], pos_mol[c[1]][1]
                xl, yl = [x1, x2], [y1, y2]
                if t == "out":
                    mol_frame_out_conn.append([xl, yl])
                else:
                    mol_frame_in_conn.append([xl, yl])

            frame_out_conn.append(mol_frame_out_conn)
            frame_in_conn.append(mol_frame_in_conn)
        
        return frame_in_conn, frame_out_conn

    def get_position_molecule(self, index):

        self.check_multiple_molecule()

        # used to compute the angles difference and the conformation
        frame_pos_list, frame_pos_per_mol_list, frame_mol_id_list = [],[], []

        # for each molecule
        for mol_id in range(self.nb_molecule):
            # get all position
            mol_frame_pos_list = []

            pos_mol = [ (self.struct_list[0][index][mol_id][k], self.struct_list[1][index][mol_id][k] ) for k in range(self.n_node)]
            for n_id in range(self.n_node):
                pos = pos_mol[n_id]
                mol_frame_pos_list.append(pos)
                frame_pos_list.append(pos)
                frame_mol_id_list.append(mol_id)

            frame_pos_per_mol_list.append(mol_frame_pos_list)
        return  frame_pos_list, frame_pos_per_mol_list, frame_mol_id_list

    def compute_frame_binding(self, index, minima):

        self.check_multiple_molecule()

        frame_pos_list, frame_pos_per_mol_list, frame_mol_id_list = self.get_position_molecule(index)

        frame_self_bonded, frame_outer_bounded = 0, 0
        frame_binding_dict = {}
        frame_specific_binding = [0]*(self.n_node**2)

        neigh = NearestNeighbors(radius=1.5*minima)
        neigh.fit(frame_pos_list)
        nbrs = neigh.radius_neighbors(return_distance=False)

        for i, neigh in enumerate(nbrs):
            if len(neigh) > 0:
                for j in neigh:
                    if frame_mol_id_list[i] != frame_mol_id_list[j]:
                        frame_outer_bounded += 1
                    else:
                        frame_self_bonded += 1
                    if frame_mol_id_list[i] < frame_mol_id_list[j]:
                        if (frame_mol_id_list[i], frame_mol_id_list[j]) in frame_binding_dict:
                            frame_binding_dict[(
                                frame_mol_id_list[i], frame_mol_id_list[j])] += 1
                        else:
                            frame_binding_dict[(
                                frame_mol_id_list[i], frame_mol_id_list[j])] = 1

        for _, nb_binded in frame_binding_dict.items():
            frame_specific_binding[nb_binded] += 1

        return frame_specific_binding, frame_binding_dict, frame_self_bonded, frame_outer_bounded

    def compute_frame_angle(self, index, ratio_distance_conformation):
        # sourcery skip: merge-duplicate-blocks, remove-redundant-if

        angles_healthy = self.angle_list_to_compare["healty"]
        angles_prion = self.angle_list_to_compare["prion"]

        _, frame_pos_per_mol_list, _ = self.get_position_molecule(index) 

        allowed_distance = ratio_distance_conformation * \
            Util.get_avg_angle_distance(angles_healthy, angles_prion)

        frame_count_healthy, frame_count_prion, frame_count_other = 0, 0, 0
        frame_count_healthy_even, frame_count_prion_even, frame_count_othe_even = 0, 0, 0
        conf_list = []
        angle_distance_list = []

        for mol_id in range(self.nb_molecule):
            
            pos_active = frame_pos_per_mol_list[mol_id]
            angle_active = Util.get_angles(pos_active)

            dist_to_healthy = Util.get_avg_angle_distance(
                angles_healthy, angle_active)
            dist_to_prion = Util.get_avg_angle_distance(
                angles_prion, angle_active)

            angle_distance_list.append( (dist_to_healthy,dist_to_prion,allowed_distance) )
            
            if dist_to_healthy<dist_to_prion:
                if dist_to_healthy<allowed_distance:
                    frame_count_healthy += 1
                    conf_list.append(0)
                    if mol_id%2==0:
                        frame_count_healthy_even+=1
                else:
                    frame_count_other += 1
                    conf_list.append(2)
                    if mol_id%2==0:
                        frame_count_othe_even+=1

            elif dist_to_prion<allowed_distance:
                frame_count_prion += 1
                conf_list.append(1)
                if mol_id%2==0:
                    frame_count_prion_even+=1
            else:
                frame_count_other += 1
                conf_list.append(2)
                if mol_id%2==0:
                    frame_count_othe_even+=1


        return frame_count_healthy, frame_count_prion, frame_count_other, frame_count_healthy_even, frame_count_prion_even, frame_count_othe_even, conf_list, angle_distance_list
        
    