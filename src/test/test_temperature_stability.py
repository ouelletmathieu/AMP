import unittest
import sys
import os

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from analysis import get_max_stable_temp
import random


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


    
    
if __name__ == '__main__':
    unittest.main()






