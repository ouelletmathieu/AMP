import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from test_file import Fake_file
from protein_template import Protein_Template
from util import Util
from my_lammps import ID_counter, MyLammps
from lammps import PyLammps
from logger import Logger



def logger_to_flat(logger_struct_list ):
    list_vec = []
    for i in range(len(logger_struct_list[0])):
        new_vec = []
        for j in range(len(logger_struct_list[0][i])):
            new_vec.append(logger_struct_list[0][i][j])
            new_vec.append(logger_struct_list[1][i][j])
        list_vec.append(new_vec)
    return list_vec


# Defining main function
def main():

    cwd = os.getcwd()
    file_path_prion = Util.remove_dropbox_path(str(cwd) + "/src/test/" +"test_prion_1_newclass.lj")
    file_path_healthy = Util.remove_dropbox_path(str(cwd) + "/src/test/" +"test_healthy_1_newclass.lj")

    ffile = Fake_file()    
    protein = Protein_Template()
    protein.set_healthy_structure(ffile.healthy_position, ffile.connection)
    protein.set_interaction( (1,1) , (1,1), [])
    protein.set_prion_structure(ffile.prion_position)
    protein.create_Lammps_file(file_path_healthy, type_mol = "healthy")
    protein.create_Lammps_file(file_path_prion, type_mol = "prion")

    file_path_1 = file_path_prion
    file_path_2 = file_path_healthy

    temperature = 0.002
    damping = 5
    max_time = 20
    n_step = 200000
    n_plot = 200

    L_main = PyLammps()
    id_counter = ID_counter()
    L = MyLammps(L_main,id_counter )

    vec_list_all = []
    energy_list_all = []
    energyKin_list_all = []

    for i in range(40):
        
        temperature = 0.001 + 0.0005*(i/4)
        
        L.command("clear")
        
        if i%2==0:
            file_path = file_path_1
        else:
            file_path = file_path_2
        
        L.create_molecule_2d(file_path)


        main_logger = Logger(L, [ffile.healthy_position, ffile.prion_position], id_struct_to_compare=["healty", "prion"])
        main_logger.log(0)

        L.run_brownian(temperature, damping, max_time, n_step, n_plot, main_logger)

        main_logger.plot_energy()
        main_logger.plot_angular_sum_distance(id_to_compare_list = ["healty", "prion"] )

        vec_list_all.extend(logger_to_flat(main_logger.struct_list))
        energy_list_all.extend(main_logger.energy_list)
        energyKin_list_all.extend(main_logger.energyKin_list)



# __name__
if __name__=="__main__":
    main()

