import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from test_file import Fake_file
from analysis import get_protein_max_stability_temp
from protein_template import Protein_Template
from util import Util


# Defining main function
def main():
    """Compute the max temperature a protein can be considered stable ()
    """
    cwd = os.getcwd()
    file_path_prion = Util.remove_dropbox_path(str(cwd) + "/src/test/" +"test_prion_1_newclass.lj")
    file_path_healthy = Util.remove_dropbox_path(str(cwd) + "/src/test/" +"test_healthy_1_newclass.lj")


    print("\n " + file_path_prion + " \n")

    ffile = Fake_file()    
    protein = Protein_Template()
    protein.set_healthy_structure(ffile.healthy_position, ffile.connection)
    protein.set_interaction( (1,1) , (1,1), [])
    protein.set_prion_structure(ffile.prion_position)
    protein.create_Lammps_file(file_path_healthy, type_mol = "healthy")
    protein.create_Lammps_file(file_path_prion, type_mol = "prion")


    max_temp_totest = 0.004
    square_angle_sum = 25

    max_stable_temp_healthy =   get_protein_max_stability_temp(file_path_healthy, protein.healthy_position, protein.prion_position, 'healthy', max_temp_totest, square_angle_sum, 
                                max_time= 30, n_plot = 30, nb_needed_for_stability = 10,  max_try = 5)
    max_stable_temp_prion =     get_protein_max_stability_temp(file_path_prion, protein.healthy_position, protein.prion_position,'prion', max_temp_totest, square_angle_sum, 
                                max_time= 30, n_plot = 30, nb_needed_for_stability = 10,  max_try = 5)

    print("max temperature healthy protein: " + str(max_stable_temp_healthy))
    print("max temperature prion protein: " + str(max_stable_temp_prion))

  

# __name__
if __name__=="__main__":
    main()

