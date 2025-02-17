import numpy as np

class Fake_file:
    """This class contain a fake file so we do not have to load a real file for testing 
    """
    
    def __init__(self):
        # 0:boundary    1:inside 
        self.healthy_position =       \
            np.array([                \
                [0, 0.5,     0.1340], \
                [0, -0.5,    0.1340], \
                [0, 0.3660,  0.6340], \
                [0, -0.3660, 0.6340], \
                [0, 0,       0],      \
                [0, 0,       1.0000], \
                [1, 0.0,    0.3330],  \
                [1, 0.0096, 0.1710],  \
                [1, 0.0001, 0.9497] ])

        self.connection = np.array([ [2, 5],\
                                       [2, 4],\
                                       [4, 6],\
                                       [6, 3],\
                                       [3, 1],\
                                       [1, 5],\
                                       [2, 7],\
                                       [4, 7],\
                                       [3, 7],\
                                       [1, 7],\
                                       [4, 8],\
                                       [3, 8],\
                                       [5, 8],\
                                       [1, 9],\
                                       [2, 9],\
                                       [6, 9],\
                                        ])

        self.connection = self.connection - 1

        self.prion_position =        \
            np.array([                \
                [0, 0.4314,  0.1964], \
                [0, -0.4314, 0.1964], \
                [0, 0.4314,  0.7140], \
                [0, -0.4314, 0.7140], \
                [0, 0,       0.4824], \
                [0, 0,       1.0000], \
                [1, 0.0,    0.5180],  \
                [1, 0.0081, 0.3113],  \
                [1, 0.0002, 1.0503] ])

        self.mass = [1,1]

        self.K = [10000000,1]