import numpy as np
import itertools
import math


class Util:
    
    #
    def get_angles_lammps(list_atom : list[tuple[float,float,float]], n_atoms: int) -> list[float]:
        """return a list of angle between all 3 group of atom. (highly inneficient way to characterize a confomarmation)
            expect atoms in the Lammps format i.e. ((x1,y1,z1),,...,(xn,yn,zn))

        Args:
            list_atom (list[tuple[float,float,float]]): atoms in the Lammps format i.e. ((x1,y1,z1),,...,(xn,yn,zn))
            n_atoms (int): number of atoms in total

        Returns:
            list[float]: list of angle in the molecule
        """
        list_angle = []
        
        for nangle in itertools.combinations(list(range(n_atoms)), 3):
            a = np.array([list_atom[nangle[0]].position[0],list_atom[nangle[0]].position[1]])
            b = np.array([list_atom[nangle[1]].position[0],list_atom[nangle[1]].position[1]])
            c = np.array([list_atom[nangle[2]].position[0],list_atom[nangle[2]].position[1]])
            ba = a - b
            bc = c - b
            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle = np.arccos(cosine_angle)
            list_angle.append(angle)

        return list_angle

    #
    def get_angles(list_atom : list[tuple[int, float,float]], n_atoms:int ) -> list[float]:
        """return a list of angle between all 3 group of atom. (highly inneficient way to characterize a confomarmation) 
           expect atoms in the Lammps format i.e. ((type1,x1,y1),,...,(typen,x1,yn)) 

        Args:
            list_atom (list[tuple[int, float,float]]): atoms in the Lammps format i.e. ((type1,x1,y1),,...,(typen,x1,yn)) 
            n_atoms (int): number of atoms

        Returns:
            list[float]:  list of angle in the molecule
        """

        list_angle = []
        
        for nangle in itertools.combinations(list(range(n_atoms)), 3):
            a = np.array([list_atom[nangle[0],1],list_atom[nangle[0],2]])
            b = np.array([list_atom[nangle[1],1],list_atom[nangle[1],2]])
            c = np.array([list_atom[nangle[2],1],list_atom[nangle[2],2]])
            ba = a - b
            bc = c - b
            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle = np.arccos(cosine_angle)
            list_angle.append(angle)

        return list_angle

    def get_avg_angle_distance(list_angle_A : list[float], list_angle_B : list[float]) -> float:
        """get distance form two list of angles 

        Args:
            list_angle_A (list[float]): list of angle in the molecule
            list_angle_B (list[float]): list of angle in the molecule

        Returns:
            float: distnace (euclideen)
        """
        total_dist = 0
        for i in range(len(list_angle_A)):
            total_dist += (list_angle_A[i]-list_angle_B[i])**2

        return math.sqrt(total_dist)


    