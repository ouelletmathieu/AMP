import unittest
import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from protein_template import Protein_Template
from test_file import Fake_file


class Test_Protein_Template(unittest.TestCase):

    
    def test_main(self):

        #testing Protein_Template
        structure2 = Protein_Template()

        #loading the fake file 
        ffile = Fake_file()
        structure2.set_healty_structure(ffile.healty_position, ffile.connection)
        structure2.set_interaction( (1,1) , (1,1), [])
        structure2.set_prion_structure(ffile.prion_position)


        structure2.create_Lammps_file("test_healty_1_newclass.lj", type_mol = "healthy")
        structure2.create_Lammps_file("test_prion_1_newclass.lj", type_mol = "prion")

        self.assertTrue(True)  


    
if __name__ == '__main__':
    unittest.main()
