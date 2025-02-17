import unittest
import sys
import os
import pathlib

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from analysis import get_max_stable_temp
import random
from test_file import Fake_file
from analysis import get_protein_max_stability_temp
from protein_template import Protein_Template
from util import Util

class Test_Protein_Template(unittest.TestCase):

    
    def test_temperature_stability_fake_experiment(self):
        

        for i in range(10):
            fake_exp = lambda temp : 40 if temp<0.002 else 60*random.random()
            max_time_test = 40
            max_temp_test = 0.004
            min_temp_test = 0
            max_test_test = 3
            max_try = 20

            max_temp_ok = False
            max_temp = -1
            while not max_temp_ok:
                try :
                    max_temp = get_max_stable_temp(min_temp_test, max_temp_test, fake_exp, max_test_test, max_time_test, 0.0001, max_try)
                    if max_temp is not None :
                        max_temp_ok=True 
                    else:
                        print('did not converge, restarting', max_test_test+1, sep=" : ")
                        max_test_test = max_test_test+1
                except ValueError:
                    max_temp_test = max_temp_test*2
                    print('increasing max temperature', max_temp_test*2, sep=" : ")

        print(max_temp)

        acceptable_range =  (max_temp > 0.00196875 - 0.0001) and (max_temp < 0.00196875 + 0.0001)
        self.assertTrue(acceptable_range) 


    def test_temperature_stability_fake_experiment(self):

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

        max_stable_temp_healthy =  get_protein_max_stability_temp(file_path_healthy, protein.healthy_position, protein.prion_position, 'healthy', max_temp_totest, square_angle_sum)
        max_stable_temp_prion =    get_protein_max_stability_temp(file_path_prion, protein.healthy_position, protein.prion_position,'prion', max_temp_totest, square_angle_sum)

        print(max_stable_temp_healthy,max_stable_temp_prion)
        
    
if __name__ == '__main__':
    unittest.main()






