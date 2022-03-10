
import os
import matplotlib.pyplot as plt

class Nudged_Elastic_Band:
    
    
    #method that compute the 
    @staticmethod
    def create_from_protein_template(protein_template):
        print("print not implemented yet")
           
            
    #Create the final atom position file needed to do a NEB analysis
    @staticmethod
    def create_final_file( path, atom_list ):
        try:
            os.remove(path)
        except:
            print("new file")
        sep = "   "    
        f = open(path, "a")
        f.write(str(len(atom_list)) + "\n")
        for i in range(len(atom_list)):
            f.write(str(i+1) + sep + str(atom_list[i,1]) + sep + str(atom_list[i,2]) + sep + "0.0 \n" ) 

        f.close()
        
    #Create the task  file needed to do a NEB analysis
    @staticmethod
    def create_neb_task_file(path, init_path, final_path, path_dir_dump, timestep, spring_inter_replica = 1.0, MaxE = 0.0, nbMin = 2000, nbClimb = 2000, maxF = 0.0, prnt = -1 ):
        try:
            os.remove(path)
        except:
            print("new file")

        sep = "     "
        f = open(path, "a")

        f.write("# LAMMPS task file for neb \n \n ")

        f.write(sep + "dimension   2"+"\n")
        f.write(sep + "atom_style bond" + "\n")
        f.write(sep + "boundary   f f p"+"\n")
        f.write(sep + "bond_style harmonic"+"\n")
        f.write(sep + "atom_modify	map array sort 0 0.0"+"\n")
        f.write("\n")

        f.write(sep + "variable  u uloop 20\n")
        f.write(sep +"read_data  " + init_path + "\n")
        f.write(sep +"timestep  "+ str(timestep) + "\n")
        f.write(sep +"fix 1 all neb " + str(spring_inter_replica) + " parallel ideal"+ "\n")
        f.write("\n")
        f.write(sep + "thermo 100"+ "\n")
        f.write( sep +"fix 2 all enforce2d" + "\n")
        f.write(sep + "dump 1 all atom 10 " + path_dir_dump+"dump.neb.$u"+ "\n")
        f.write(sep + " min_style quickmin"+ "\n")
     

        if MaxE==0.0 and maxF==0.0:
            MaxE = 0.00000001
        if prnt == -1:
            prnt = int(nbMin/50)


        f.write( sep +"neb " +  str(MaxE) +" "+ str(maxF) +" "+ str(nbMin) +" "+ str(nbClimb) +" "+ str(prnt) +" final "+ final_path + "\n")

        f.close()

    @staticmethod
    def restructure_output(neb_final_py_path, dump_dir, output_file = "dump.restructured.final"):
        import subprocess
        
        str_to_run = "python " + str(neb_final_py_path) + " -o "+ str(output_file) + " -r " + dump_dir +"dump.neb.*"
        test = subprocess.Popen([str_to_run], stdout=subprocess.PIPE)
        output = test.communicate()[0]  
        print(output)
        
        return str_to_run

    @staticmethod
    def plot_position(output_file_path, index_list=[], marker = False):
                
        list_to_plot = Nudged_Elastic_Band.get_atom_position(output_file_path)

        for i in range(len(index_list)):
            xpos = list_to_plot[index_list[i]][0]
            ypos = list_to_plot[index_list[i]][1]
            plt.scatter(xpos,ypos)
            
            if marker:
                for j in range(len(xpos)):
                    plt.text(xpos[j], ypos[j], str(j))
        
        plt.show()

    @staticmethod
    def get_atom_position(output_file_path):
                
        file1 = open(output_file_path, 'r')
        Lines = file1.readlines()

        timestep_id = [] 

        for i, line in enumerate(Lines):
            if 'TIMESTEP' in line:
                timestep_id.append(i)


        atom_per_sim = int(Lines[timestep_id[0]+3])

        pos_string_list_list = [[] for _ in range(len(timestep_id))]
        atom_pos = [[[],[]] for _ in range(len(timestep_id))]

        for i in range(len(timestep_id)):

            for j in range(atom_per_sim):

                string = Lines[timestep_id[i]+j+9]

                splitted = string.split()
                map_object = map(float, splitted)
                list_of_float= list(map_object)

                pos_string_list_list[i].append( (list_of_float[2],list_of_float[3] ) ) 
                atom_pos[i][0].append(list_of_float[2])
                atom_pos[i][1].append(list_of_float[3])

        return atom_pos

    @staticmethod    
    def plot_energy(log_file_path):
        
        file1 = open(log_file_path, 'r')
        Lines = file1.readlines()

        timestep_id = [] 
        last_line = Lines[-1]
        splitted = last_line.split()
        map_object = map(float, splitted)
        list_of_float = list(map_object)

        list_of_float = list_of_float[9:]
        n_replica = int(len(list_of_float)/2)

        rc_list = []
        energy_list = []

        for i in range(len(list_of_float)):
            if i%2==0:
                rc_list.append(list_of_float[i])
            else:
                energy_list.append(list_of_float[i])

        plt.plot(rc_list, energy_list)
        plt.show()
