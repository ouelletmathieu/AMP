from typing import TypeVar, Union, Literal
import math
import os


class Protein_Template:
    """This is the template to construct a protein. It contain both the definition for the healthy 
        one and the prion one. It is able to instantiate them in simulation. 
        Save them to file and load them from file.
    """

    def __init__(self):
        self.healthy_position = None 
        self.prion_position = None
        self.connection = None
    
    def set_healthy_structure(self, position, connection = None):
        """set the structure to be considered healthy 

        Args:
            position (np.array([a,x,y])): atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside
            connection ([[id1,id2],[...], .. ], optional):  if defined the previous connection pattern is over written. Else it is kept the same
    #            format: [[id1,id2],[...], .. ] id of connected node
        """
        self.healthy_position = position
        self.connection  = connection if connection  is not None else self.connection 
     
    def set_prion_structure(self, position, connection = None):
        """set the structure to be considered the prion 

        Args:
            position (np.array([a,x,y])): atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside
            connection ([[id1,id2],[...], .. ], optional):  if defined the previous connection pattern is over written. Else it is kept the same
    #            format: [[id1,id2],[...], .. ] id of connected node
        """
        self.prion_position = position
        self.connection  = connection if connection  is not None else self.connection    
    
    def set_interaction(self, mass : tuple[float,float] = (1,1) , K : tuple[float,float] = (1,1), 
                        lj_param : Union[list, tuple[float,float]] = []):
        """set the interaction to be considered 

        Args:
            mass (tuple[float,float], optional): mass of the nodes  mass[0] = mass of normal atom mass[1]= mass of the clamping atom. Defaults to (1,1).
            K (tuple[float,float], optional): K[0] = boundary spring constant K[1] = inside spring constant. Defaults to (1,1).
            lj_param (Union[list, tuple[float,float]], optional):  #lj_param: parameter for the Lennard Jones interaction
              format: lj_param = [] if no lennard-jones potential, 
                      lj_param = [d_min, emin] (distance of minimum and absolute energy minimum). Defaults to [].
        """
        self.mass = mass
        self.K = K
        self.lj_param = lj_param
        
    def read_from_file(self, path):
        print("todo not coded")
        
    def save_template(self, path):
        print("todo not coded")
    

    def create_Lammps_file(self, path : str, type_mol : Literal["prion", "healthy"], cell = [-6,6,-6,6,-6,6]):  
        """create the lamps file used for simulation 
        cell = simulation cell [minx, maxx, miny, maxy, minz, maxz ]
        name of the file to save
        """     
        connections = self.connection
        mass = self.mass
        K = self.K
        lj_param = self.lj_param
        
        #set the right positions list
        if type_mol=="prion":
            atom_list = self.prion_position
        elif type_mol=="healthy":
            atom_list = self.healthy_position
        else:
            print("type is not valid ")
            raise ValueError("type is not valid either prion or healthy")
        
        #remove the file if it already exist
        try:
            os.remove(path)
            print("overwritten file")
        except:
            print("new file")
        
        #starting the new file
        sep = "     "
        f = open(path, "a")

        f.write("# LAMMPS data file for rigid bodies \n \n ")

        #declare objects number
        n = len(atom_list)
        f.write(sep +str(n)+sep+"atoms"+"\n")
        f.write(sep +str(len(connections))+sep+"bonds"+"\n")
        f.write(sep +str(0)+sep+"angles"+"\n")
        f.write(sep +str(0)+sep+"dihedrals"+"\n")
        f.write("\n")

        #declare types number
        f.write(sep +"2"+ sep +"atom types \n")
        f.write(sep +str(len(connections))+ sep +"bond types \n")
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
        f.write(sep+"1"+sep+str(mass[0]) + "\n")
        f.write(sep+"2"+sep+str(mass[1]) + "\n")
        f.write("\n")
        
        #declare atom positions 
        f.write("Atoms"+ "\n")
        f.write("\n")
        for i in range(n):
            if i ==0:
                f.write(sep +str(i+1) +sep+ "1" +sep+ "2" +sep+ str(atom_list[i,1]) + sep + str(atom_list[i,2]) + sep + "0.0" + "\n" ) 
            else: 
                f.write(sep +str(i+1) +sep+ "1" +sep+ "1" +sep+ str(atom_list[i,1]) + sep + str(atom_list[i,2]) + sep + "0.0" + "\n" ) 
        f.write("\n")

        type_list = []
        lenght_list = []

        #declare bonds length and atom linked  
        f.write("Bonds" + "\n")
        f.write("\n")
        for i in range(connections.shape[0]):
            atom1 = connections[i,0]
            atom2 = connections[i,1]

            atom1_type = atom_list[atom1,0]
            atom1_pos =  (atom_list[atom1,1],atom_list[atom1,2])

            atom2_type = atom_list[atom2,0]
            atom2_pos =  (atom_list[atom2,1],atom_list[atom2,2])        

            if atom1_type!=atom2_type:
                type_list.append(1)
            else:
                type_list.append(0)
            #length are computed here
            lenght = math.sqrt((atom1_pos[0]-atom2_pos[0])**2 + (atom1_pos[1]-atom2_pos[1])**2)
            lenght_list.append(lenght)
            f.write(sep +str(i+1) + sep + str(i+1) + sep + str(atom1+1) + sep + str(atom2+1)  + "\n") 
        f.write("\n")
        
        #declare bonds coefficient
        f.write("Bond Coeffs" + "\n")
        f.write("\n")
        for i in range(connections.shape[0]):
            if type_list[i]==0:
                f.write(sep +str(i+1) + sep + str(K[0]) + sep + str(lenght_list[i]) +  "\n")
            else:
                f.write(sep +str(i+1) + sep + str(K[1]) + sep + str(lenght_list[i]) +  "\n")    

        if len(lj_param)!=0:
            print("please pass the coeff in the simulation constructor")
            raise ValueError("lj parameter are now passed through the simulator")
        f.close()



    def get_bond_length(self, type_mol : Literal["prion", "healthy"]):

        connections = self.connection
        K = self.K
        lenght_list, type_list = [], []
       

        #set the right positions list
        if type_mol=="prion":
            atom_list = self.prion_position
        elif type_mol=="healthy":
            atom_list = self.healthy_position
        else:
            print("type is not valid ")
            raise ValueError("type is not valid either prion or healthy")       

        for i in range(connections.shape[0]):
            atom1 = connections[i,0]
            atom2 = connections[i,1]

            atom1_type = atom_list[atom1,0]
            atom1_pos =  (atom_list[atom1,1],atom_list[atom1,2])

            atom2_type = atom_list[atom2,0]
            atom2_pos =  (atom_list[atom2,1],atom_list[atom2,2])        
            type = 0
            if atom1_type!=atom2_type:
                type = 1
            type_list.append(type)
            #length are computed here
            lenght = math.sqrt((atom1_pos[0]-atom2_pos[0])**2 + (atom1_pos[1]-atom2_pos[1])**2)
            lenght_list.append(lenght)
        
        return lenght_list, type_list


    