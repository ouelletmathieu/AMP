import random
from re import X
from tkinter import Y
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import pandas as pd

class Count:

    def __init__(self):
         self.count = 1

    def next(self):
        self.count +=1
        return self.count - 1

class Network:
    """class that allows to store information about simulation from lammps 
    """
    
    def __init__(self):

        #id, x,y,z,
        self.network_point = []
        #connection matrix
        self.matrix_connection = []
        #distance matrix
        self.matrix_distance = []       
        #dictionarry of edge point  (from, to):[id,x,y,z]
        self.edge_points = dict()
        #counter used for id of atoms
        self.counter = Count()

    def sample_spherical(self, npoints, ndim=3):
        vec = np.random.randn(ndim, npoints)
        vec /= np.linalg.norm(vec, axis=0)
        return vec

    def create_network_from_file_RAS(self, file_path, nb_edge, exp_param, edge_point_density):
        df = pd.read_csv(file_path, sep=",")
        col_x  = df['R'].tolist()
        col_y = df['A'].tolist()
        col_z = df['S'].tolist()
        n = len(col_y)

        for i in range(n):
            x,y,z = col_x[i]/100, col_y[i]/100, col_z[i]/100
            self.network_point.append((self.counter.next(), x,y,z) )

        self.matrix_connection = np.zeros((n,n))
        self.matrix_distance = np.zeros((n,n))

        n_tot = 0

        while n_tot < nb_edge:
            i = random.randint(0,n-1)
            j = random.randint(0,n-1)

            if i!=j and self.matrix_connection[i,j]==0:
                
                p = self.network_point
                dx = (p[j][1]-p[i][1])**2
                dy = (p[j][2]-p[i][2])**2
                dz = (p[j][3]-p[i][3])**2
                dist = math.sqrt(dx + dy + dz)

                if random.random() < math.exp(-exp_param*dist):
                    n_tot+=1
                    self.matrix_distance[i,j] = dist
                    self.matrix_connection[i,j] = 1
                    self.edge_points[(i,j)]=[]



        for i,j in self.edge_points:
            
            nb_point = int(round(self.matrix_distance[i,j]*edge_point_density))
            if nb_point<1:
               nb_point=1 
            p = self.network_point
            x1,y1,z1 = p[i][1], p[i][2], p[i][3]
            x2,y2,z2 = p[j][1], p[j][2], p[j][3]
            #print("edges: " + str((i,j)))
            #X-----o----o----o----o----X
            for n in range(1,nb_point+1):
                x_point = (x1-x2)*n/(nb_point+1) + x2
                y_point = (y1-y2)*n/(nb_point+1) + y2
                z_point = (z1-z2)*n/(nb_point+1) + z2
                self.edge_points[(i,j)].append((self.counter.next(), x_point, y_point, z_point))
        

    def create_random_network(self, n=5, density_edge = 0.5, edge_point_density = 5, dist_min = 0, rand = True):
        print("nb edges: " + str(density_edge*n*(n-1)))
        if rand:
            for i in range(n):
                while True:
                    rvec = self.sample_spherical(1, ndim=3)
                    x, y, z  = rvec[0][0],rvec[1][0],rvec[2][0]

                    is_okay = True
                    for pt in self.network_point:
                        dx = (x-pt[1])**2
                        dy = (y-pt[2])**2
                        dz = (z-pt[3])**2
                        dist = math.sqrt(dx+dy+dz)
                        if dist < dist_min:
                            is_okay=False
                            break
                    if is_okay:
                        self.network_point.append( (self.counter.next(), x,y,z) )
                        break
        else:
            self.network_point.append( (self.counter.next(), 0.0,0.0,1.0) )
            self.network_point.append( (self.counter.next(), 1.0,0.0,1.0) )
            self.network_point.append( (self.counter.next(), 1.0,1.0,1.0) )
            self.network_point.append( (self.counter.next(), 0.0,1.0,1.0) )
            self.network_point.append( (self.counter.next(), 0.0,0.0,0.0) )
            self.network_point.append( (self.counter.next(), 1.0,0.0,0.0) )
            self.network_point.append( (self.counter.next(), 1.0,1.0,0.0) )
            self.network_point.append( (self.counter.next(), 0.0,1.0,0.0) )


        self.matrix_connection = np.zeros((n,n))
        self.matrix_distance = np.zeros((n,n))

        if rand:
            for j in range(n):
                for i in range(0,j):
                    if random.random()<density_edge:
                        
                        
                        self.matrix_connection[i,j] = 1
                        self.edge_points[(i,j)]=[]
                        p = self.network_point
                        dx = (p[j][1]-p[i][1])**2
                        dy = (p[j][2]-p[i][2])**2
                        dz = (p[j][3]-p[i][3])**2
                        self.matrix_distance[i,j] = math.sqrt(dx + dy + dz)



        else:
            lcon = ((0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4))
            for i,j in lcon:
                #print((i,j))
                self.matrix_connection[i,j] = 1
                self.edge_points[(i,j)]=[]
                self.matrix_distance[i,j] = 1.0
        
        #print("adding edges")
        for i,j in self.edge_points:
            
            nb_point = int(round(self.matrix_distance[i,j]*edge_point_density))
            if nb_point<1:
               nb_point=1 
            p = self.network_point
            x1,y1,z1 = p[i][1], p[i][2], p[i][3]
            x2,y2,z2 = p[j][1], p[j][2], p[j][3]
            #print("edges: " + str((i,j)))
            #X-----o----o----o----o----X
            for n in range(1,nb_point+1):
                x_point = (x1-x2)*n/(nb_point+1) + x2
                y_point = (y1-y2)*n/(nb_point+1) + y2
                z_point = (z1-z2)*n/(nb_point+1) + z2
                self.edge_points[(i,j)].append((self.counter.next(), x_point, y_point, z_point))
        
    

    def plot(self):

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        for main_point in self.network_point:
            ax.scatter3D(main_point[1], main_point[2], main_point[3], s=40, c='g')
        
        for key in self.edge_points:
            for point in self.edge_points[key]:
                ax.scatter3D(point[1], point[2], point[3], s=10, c='r')

        return fig

    def get_atom_number(self):
        total_n = len(self.network_point)
        for key in self.edge_points:
            total_n+=len(self.edge_points[key])
        return len(self.network_point), total_n

    def get_atom_list(self):
        #type (1 for fixed, 2 for edge), id, x,y,z
        list_atom = []
        for pt in  self.network_point:
            list_atom.append((1,pt[0],pt[1],pt[2],pt[3]))
        
        for key in self.edge_points:
            for pt in self.edge_points[key]:
                list_atom.append((2,pt[0],pt[1],pt[2],pt[3]))

        return list_atom

    #id1, id2, distance, type1, type2
    def get_connection(self):
        connection_list = []
        for key in self.edge_points:
            from_pt_id = key[0]
            to_pt_id = key[1]

            for n, pt in enumerate(self.edge_points[key]):
                if n==0:
                    dx = (pt[1]-self.network_point[from_pt_id][1])**2
                    dy = (pt[2]-self.network_point[from_pt_id][2])**2
                    dz = (pt[3]-self.network_point[from_pt_id][3])**2
                    dist = math.sqrt(dx+dy+dz)
                    connection_list.append( (from_pt_id+1, pt[0], dist, 1, 2))
                else:
                    dx = (pt[1]-self.edge_points[key][n-1][1])**2
                    dy = (pt[2]-self.edge_points[key][n-1][2])**2
                    dz = (pt[3]-self.edge_points[key][n-1][3])**2
                    dist = math.sqrt(dx+dy+dz)
                    connection_list.append( (self.edge_points[key][n-1][0], pt[0], dist, 2,2))



            pt = self.edge_points[key][-1]
            dx = (pt[1]-self.network_point[to_pt_id][1])**2
            dy = (pt[2]-self.network_point[to_pt_id][2])**2
            dz = (pt[3]-self.network_point[to_pt_id][3])**2
            dist = math.sqrt(dx+dy+dz)
            connection_list.append( (pt[0],to_pt_id+1, dist, 2, 1 ))

        return connection_list


    def get_connection_number(self):
        total_n = 0
        #X-----o----o----o----o----X
        for key in self.edge_points:
            total_n+= (len(self.edge_points[key])+1)
        return total_n  

    def create_Lammps_file(self, path : str,  K, mass= 1,  cell= [-6,6,-6,6,-6,6]):  
        """create the lamps file used for simulation 
        cell = simulation cell [minx, maxx, miny, maxy, minz, maxz ]
        name of the file to save
        """     

        
        #remove the file if it already exist
        try:
            os.remove(path)
            print("overwritten file")
        except:
            print("new file")
        
        #starting the new file
        sep = "     "
        f = open(path, "a")

        f.write("# LAMMPS data file for rigid bodies \n \n")

        #declare objects number
        n_main, n_total = self.get_atom_number()
        n_conntection = self.get_connection_number()
        f.write(sep +str(n_total)+sep+"atoms"+"\n")
        f.write(sep +str(n_conntection)+sep+"bonds"+"\n")
        f.write(sep +str(0)+sep+"angles"+"\n")
        f.write(sep +str(0)+sep+"dihedrals"+"\n")
        f.write("\n")

        #declare types number
        f.write(sep +"2"+ sep +"atom types \n")
        f.write(sep +str(n_conntection)+ sep +"bond types \n")
        f.write(sep +"0     angle types"+ "\n")
        f.write(sep +"0     dihedral types"+ "\n")
        f.write("\n")
        f.write(sep +str(cell[0])+ sep +  str(cell[1]) + " xlo xhi"+ "\n")
        f.write(sep +str(cell[2])+ sep +  str(cell[3]) + " ylo yhi"+ "\n")
        f.write(sep +str(cell[4])+ sep +  str(cell[5]) + " zlo zhi"+ "\n")
        f.write("\n")

        #declare masses
        f.write("Masses"+ "\n")
        f.write("\n")
        f.write(sep+"1"+sep+str(1000000000000000000) + "\n")
        f.write(sep+"2"+sep+str(mass) + "\n")
        f.write("\n")

        atom_list = self.get_atom_list()

        #declare atom positions 
        f.write("Atoms"+ "\n")
        f.write("\n")
        #type (1 for fixed, 2 for edge), id, x,y,z
        for i in range(len(atom_list)):
            if atom_list[i][0] == 1:
                f.write(sep +str(atom_list[i][1]) +sep+ "1" +sep+ "1" +sep+ str(atom_list[i][2]) + sep + str(atom_list[i][3]) + sep + str(atom_list[i][4])+ "\n" ) 
            else: 
                f.write(sep +str(atom_list[i][1]) +sep+ "2" +sep+ "2" +sep+ str(atom_list[i][2]) + sep + str(atom_list[i][3]) + sep + str(atom_list[i][4])+ "\n" ) 
        f.write("\n")

        type_list = []
        lenght_list = []
        
        connection_list = self.get_connection()

        #declare bonds length and atom linked  
        f.write("Bonds" + "\n")
        f.write("\n")
        for i in range(len(connection_list)):
            atom1 = connection_list[i][0]
            atom2 = connection_list[i][1]
            lenght_list.append(connection_list[i][2])
            f.write(sep +str(i+1) + sep + str(i+1) + sep + str(atom1) + sep + str(atom2)  + "\n") 
        f.write("\n")
        
        #declare bonds coefficient
        f.write("Bond Coeffs" + "\n")
        f.write("\n")
        for i in range(len(connection_list)):
            f.write(sep +str(i+1) + sep + str(K) + sep + str(lenght_list[i]) +  "\n") 
        f.close()
