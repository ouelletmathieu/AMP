import lammps_util.dump as dump
import numpy as np
import dill
import math

def run_dill_encoded(payload):
    args = dill.loads(payload)
    return get_conformation_from_DUMP(args)



def get_conformation_from_DUMP(data_for_task):

    path_dump, reac_coor_funcs, other = data_for_task
    print(path_dump)

    d = dump.dump(path_dump) 
    (func_translate, func_translate_binded), (endpoint_pca, endpoint_pca_bind) = reac_coor_funcs
    healthy_pc1, prion_pc1 = endpoint_pca[0][0], endpoint_pca[1][0]
    direction = 1
    line33, line66 = 0.3333*(prion_pc1-healthy_pc1)+healthy_pc1, 0.6666*(prion_pc1-healthy_pc1)+healthy_pc1
    mol_conf_list = []

    nsnaps = len(d.snaps)
    natoms = d.snaps[0].natoms
    data_array_x = np.zeros((nsnaps,natoms))
    data_array_y = np.zeros((nsnaps,natoms))

    columns = []
    for name in ['x','y']:
        columns.append(d.names[name])

    id = d.names["id"]


    for i, snap in enumerate(d.snaps):
        atoms = snap.atoms
        for j in range(natoms):
            id_atom =  int(atoms[j][0])-1
            data_array_x[i,id_atom] = atoms[j][columns[0]]
            data_array_y[i,id_atom] = atoms[j][columns[1]]


    for kn in range(10000000):
        conf_list = []
        try:
            xlist = np.array([data_array_x[:,i-1].tolist() for i in range(1+kn*7,8+kn*7)]).T
            ylist = np.array([data_array_y[:,i-1].tolist() for i in range(1+kn*7,8+kn*7)]).T
            frame_pos_list = [ [xlist[l,:], ylist[l,:] ] for l in range(len(ylist))]
        except :
            break

        pos_translated = func_translate(frame_pos_list)
        #if prion is lower than healthy on manifold reverse it
        if line66<line33:
            direction=-1
        #set first conformation
        if  abs(pos_translated[0][0] - healthy_pc1) < abs(pos_translated[0][0] - prion_pc1):
            current_conf=0
        else:
            current_conf=1

        for frame in range(len(frame_pos_list)):
            if direction==1 :
                if pos_translated[frame][0] > line66 and current_conf==0:
                    current_conf=1
                elif pos_translated[frame][0] < line33 and current_conf==1:
                    current_conf=0
            elif direction==-1:
                if pos_translated[frame][0] > line33 and current_conf==1:
                    current_conf=0
                elif pos_translated[frame][0] < line66 and current_conf==0:
                    current_conf=1
            conf_list.append(current_conf)

        mol_conf_list.append(conf_list)

    return mol_conf_list