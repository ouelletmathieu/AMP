import math
from .util import Util
import matplotlib.pyplot as plt


class Logger:
    """class that allows to store information about simulation from lammps 
    """
    
    def __init__(self, my_lammps, struct_to_compare_pos = [], id_struct_to_compare = [], is_2d=True ):
            
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

        
    def append_energy(self):
        self.energy_list.append(self.L.getEnergy())
        self.energyKin_list.append(self.L.getEnergyKin())
        
        
    def append_time(self, time):
        self.time_list.append(time)
        
        
    def append_position(self):
        #set up the position
        self.active_x_pos_list = [self.L.lammps().atoms[i].position[0] for i in range(len(self.L.lammps().atoms))]
        self.active_y_pos_list = [self.L.lammps().atoms[i].position[1] for i in range(len(self.L.lammps().atoms))]
        self.struct_list[0].append(self.active_x_pos_list)
        self.struct_list[1].append(self.active_y_pos_list)

        if not self.is_2d: 
            self.active_z_pos_list = [self.L.lammps().atoms[i].position[2] for i in range(len(self.L.lammps().atoms))]
            self.struct_list[2].append(self.active_z_pos_list)


        
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
        
  

    ###################################
    #             PLOTTING            #
    ###################################
    
    
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



    
    def _get_angle_difference(self, id_to_compare = None ):
        
        angle_n = len(self.angle_list[0])
        angle_formatted = [[] for _ in range(angle_n)]

        for i in range(len(self.angle_list)):
            for j in range(angle_n):
                
                if id_to_compare is None:
                    angle_formatted[j].append(self.angle_list[i][j] - self.angle_list[0][j] )
                else:
                    angle_formatted[j].append(self.angle_list[i][j] - self.angle_list_to_compare[id_to_compare][j] )

        return angle_formatted


    def get_RMS_angle_difference_all_time(self):
            
            angle_n = len(self.angle_list[0])
            angle_sum = 0
            nb_angle = 0

            for i in range(len(self.angle_list)):
                for j in range(angle_n):
                    nb_angle+=1
                    angle_sum += (self.angle_list[i][j] - self.angle_list[0][j])**2 
            
            return math.sqrt(angle_sum/nb_angle)
            
