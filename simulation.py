import database
import os
import math
import lammps_util.protein_template
import numpy as np
import lammps_util.dynamic_file_template as dynamic_file_template
import lammps_util.protein_template as protein_template
from lammps import PyLammps
import lammps_util.my_lammps as mylammps
import lammps_util.sim_param as sim_param
import copy
import polygon_util as pu
from shapely.geometry import Polygon, Point
import lammps_util.Nudged_Elastic_Band as NEB
from os.path import exists
import lammps_util.logger as logger
from gif_util import create_gif_multiple_molecule
from scipy.stats import qmc
import random

PRE = "AMP/src/polygon/"
NB_TEMP = 5
MIN_TEMP = 1E-4
MAX_TEMP = 1E-3
MAX_TIME = 144000*5#48000#12000 
N_PER_TIME =300
N_STEP = int(MAX_TIME*N_PER_TIME)
N_DUMP = 4*N_PER_TIME
N_PLOT = 600
PRE_DUMP = "/Volumes/Lab_math/sim_prion_2/" #"/Users/mathieuouellet/Desktop/AMP/" #

E_MIN = 0.1
MINIMA = 0.05
K_SPRING = (1,1)


def get_protein_pairs_binded_divided_bypass(pos_h, pos_p, translation, inside_node, conn, hp_vec, pp_vec, n_out, translated = True ):


    best_sol_hp, binded_hp =  hp_vec
    best_sol_pp, binded_pp =  pp_vec

    # get back the binding PH and the binding PP
    healthy_poly, prion_poly = Polygon(pos_h), Polygon(pos_p)

    pos_h_2, pos_p_2 = copy.copy(pos_h), copy.copy(pos_p)


    # get position with inside node for prion / healthy
    all_h, all_p = NEB.add_center_node(pos_h_2, pos_p_2, inside_node)
    full_healthy_hp, full_prion_pp_p2 = np.copy(all_h), np.copy(all_p)

    # get binded node from the prion that do not move
    healthy_binded, prion_binded = [p[0] for p in binded_hp], [p[1] for p in binded_hp]
    init_prion_3 = [list(prion_poly.exterior.coords)[i] for i in prion_binded]

    # transform the inside point by refitting the subset of bonded point for HP
    sol_hp_3 = [list(best_sol_hp.exterior.coords)[i]for i in prion_binded]
    full_prion_hp, trans_hp = NEB.transform_all_point_back(
        sol_hp_3, init_prion_3, best_sol_hp, all_p, False)
        
    # transform the inside point by refitting the subset of bonded point for PP
    sol_pp_3 = [list(best_sol_pp.exterior.coords)[i]
                for i in prion_binded]
    full_prion_pp_p1, trans_pp = NEB.transform_all_point_back(
        sol_pp_3, init_prion_3, best_sol_pp, all_p, False)

    #remove the closed polygon 
    full_healthy_hp = np.delete(full_healthy_hp, obj=n_out, axis=0 )
    full_prion_hp  = np.delete(full_prion_hp, obj=n_out, axis=0 )
    full_prion_pp_p1 = np.delete(full_prion_pp_p1, obj=n_out, axis=0 )
    full_prion_pp_p2 = np.delete(full_prion_pp_p2, obj=n_out, axis=0 )

    #add the type of atom
    full_healthy_hp = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_healthy_hp))]) , full_healthy_hp ]      
    full_prion_hp = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_hp))]) , full_prion_hp ]
    full_prion_pp_p1 = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_pp_p1))]) , full_prion_pp_p1 ]
    full_prion_pp_p2 = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_pp_p2))]) , full_prion_pp_p2 ]

    full_healthy_hp_trans = np.zeros( ( len(full_healthy_hp), 3) )
    full_prion_pp_p1_trans = np.zeros( ( len(full_healthy_hp), 3) )


    if translated:
        translation = copy.copy(translation)
        for i in range(n_out, len(full_healthy_hp)):
            translation[i]=i

        for k_trans in  translation.keys():
            full_healthy_hp_trans[k_trans] = full_healthy_hp[translation[k_trans]]
            full_prion_pp_p1_trans[k_trans] = full_prion_pp_p1[translation[k_trans]]
        return full_healthy_hp_trans, full_prion_hp, full_prion_pp_p1_trans, full_prion_pp_p2

    else:
        return full_healthy_hp, full_prion_hp, full_prion_pp_p1, full_prion_pp_p2

def get_protein_pairs_binded_divided(pos_h, pos_p, binded_hp, inside_node, translation, target_previous_pp_fitness, target_previous_hp_fitness, min_r_squared_default, n_point_max, n_out, translated = True):

    # correct connection (fit the indices and not matlab)
    #conn = [(c[0]-1, c[1]) for c in conn]
    #print(conn)

    # get back the binding PH and the binding PP
    healthy_poly, prion_poly = Polygon(pos_h), Polygon(pos_p)

    pos_h_2, pos_p_2 = copy.copy(pos_h), copy.copy(pos_p)



    fit_ordered_hp, sol_ordered_hp, output_index_hp, _ = pu.estimate_stacking_icp(
        healthy_poly, prion_poly, avg_dist_tol=0.3, n_point_max=n_point_max,  multiple_bind=False, output_index=True)
    hp_fit_max = pu.rescale_fit(fit_ordered_hp[-1], min_r_squared_default)
    fit_ordered_pp, sol_ordered_pp, output_index_pp, _ = pu.estimate_stacking_icp(
        prion_poly, prion_poly, avg_dist_tol=0.3, n_point_max=n_point_max,  multiple_bind=False, output_index=True)
    pp_fit_max = pu.rescale_fit(fit_ordered_pp[-1], min_r_squared_default)

    # todo assert if they are kinda the same
    assert (hp_fit_max-target_previous_hp_fitness) < 0.01, "there is a big difference in the fitness" + \
        str(hp_fit_max-target_previous_hp_fitness)
    assert (pp_fit_max-target_previous_pp_fitness) < 0.01, "there is a big difference in the fitness" + \
        str(pp_fit_max-target_previous_pp_fitness)

    # get position with inside node for prion / healthy
    all_h, all_p = NEB.add_center_node(pos_h_2, pos_p_2, inside_node)
    full_healthy_hp, full_prion_pp_p2 = np.copy(all_h), np.copy(all_p)

    # get binded node from the prion that do not move
    healthy_binded, prion_binded = [p[0] for p in binded_hp], [p[1] for p in binded_hp]
    init_prion_3 = [list(prion_poly.exterior.coords)[i] for i in prion_binded]

    # transform the inside point by refitting the subset of bonded point for HP
    sol_hp_3 = [list(sol_ordered_hp[-1].exterior.coords)[i]
                for i in prion_binded]
    full_prion_hp, trans_hp = NEB.transform_all_point_back(
        sol_hp_3, init_prion_3, sol_ordered_hp[-1], all_p, False)
        
    # transform the inside point by refitting the subset of bonded point for PP
    sol_pp_3 = [list(sol_ordered_pp[-1].exterior.coords)[i]
                for i in prion_binded]
    full_prion_pp_p1, trans_pp = NEB.transform_all_point_back(
        sol_pp_3, init_prion_3, sol_ordered_pp[-1], all_p, False)

    #remove the closed polygon 
    full_healthy_hp = np.delete(full_healthy_hp, obj=n_out, axis=0 )
    full_prion_hp  = np.delete(full_prion_hp, obj=n_out, axis=0 )
    full_prion_pp_p1 = np.delete(full_prion_pp_p1, obj=n_out, axis=0 )
    full_prion_pp_p2 = np.delete(full_prion_pp_p2, obj=n_out, axis=0 )

    #add the type of atom
    full_healthy_hp = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_healthy_hp))]) , full_healthy_hp ]      
    full_prion_hp = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_hp))]) , full_prion_hp ]
    full_prion_pp_p1 = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_pp_p1))]) , full_prion_pp_p1 ]
    full_prion_pp_p2 = np.c_[ np.array([ ( 0 if i < n_out else 1) for i in range(len(full_prion_pp_p2))]) , full_prion_pp_p2 ]

    full_healthy_hp_trans = np.zeros( ( len(full_healthy_hp), 3) )
    full_prion_pp_p1_trans = np.zeros( ( len(full_healthy_hp), 3) )

    if translated:
        translation = copy.copy(translation)
        for i in range(n_out, len(full_healthy_hp)):
            translation[i]=i

        for k_trans in  translation.keys():
            full_healthy_hp_trans[k_trans] = full_healthy_hp[translation[k_trans]]
            full_prion_pp_p1_trans[k_trans] = full_prion_pp_p1[translation[k_trans]]
        return full_healthy_hp_trans, full_prion_hp, full_prion_pp_p1_trans, full_prion_pp_p2

    else:
        return full_healthy_hp, full_prion_hp, full_prion_pp_p1, full_prion_pp_p2


def set_up_molecule_simulation(L: mylammps.MyLammps, path_molecules, length_list,  box_side,  bounding_energy, minima, cut_off, K):
    
    L.L.units("lj")
    L.L.command("dimension 2")
    L.L.command("atom_style full")
    L.L.command("boundary p p p")
    L.L.command("neighbor 3 bin")
    #L.L.command("neigh_modify every 2 delay 10 ")

    L.L.command("bond_style harmonic")
    L.L.command(
        f"region mybox block {-1*box_side} {box_side} {-1*box_side} {box_side} {-1*box_side} {box_side}")
    # create_box N_atoms region-ID keyword value ...
    L.L.command(
        f"create_box 2 mybox bond/types {len(length_list)} extra/bond/per/atom {len(length_list)} extra/special/per/atom  {len(length_list)} ")

    if isinstance(path_molecules, list):
        for i in range(len(path_molecules)):
            L.L.command(f"molecule test_mol_{i} {path_molecules[i]}")
    else:
        L.L.command(f"molecule test_mol {path_molecules}")

    for i in range(len(length_list)):
        L.L.command(
            f"bond_coeff    {str(i+1)}     {str(K[0])}     {str(length_list[i])}")

    L.L.command("mass * 1.0")
    L.L.command(f"pair_style lj/cut {cut_off}")
    L.L.command(f"pair_coeff 1 1 {bounding_energy} {minima}")
    L.L.command(f"pair_coeff 1 2 0.0 {minima}")
    L.L.command(f"pair_coeff 2 2 0.0 {minima}")


def insert_avoiding_molecule(L: mylammps.MyLammps, n_molecule,  box_side, n_type = 1):
    # create spaced molecule
    sampler = qmc.Halton(d=2, scramble=False)
    sample = sampler.random(n=n_molecule)

    #fake_point = [(-15,-15),(-15,15),(15,-15),(15,15)]
    for i in range(n_molecule):
        x_pos, y_pos, z_pos = 2*(sample[i, 0]-0.5)*0.33 * \
            box_side, 2*(sample[i, 1]-0.5)*0.33*box_side, 0

        #x_pos, y_pos, z_pos = fake_point[i][0], fake_point[i][1], 0 

        pos_str = f'{str(x_pos)} {str(y_pos)} {z_pos} '
        #select type if more than one possible
        typ = "" if n_type == 1 else f'_{i%n_type}'
        L.L.command(
            f"create_atoms 0 single {pos_str} mol test_mol{typ} 1 units box ")
        L.L.command(f"group mol_{i} molecule {i+1}")

def insert_avoiding_binded_molecule(L: mylammps.MyLammps, n_molecule,  box_side, delta_vec):
    # create spaced molecule
    sampler = qmc.Halton(d=2, scramble=False)
    sample = sampler.random(n=n_molecule)

    for i in range(n_molecule):
        x_pos, y_pos, z_pos = 2*(sample[i, 0]-0.5)*0.33* \
            box_side, 2*(sample[i, 1]-0.5)*0.33*box_side, 0


        pos_str_0 = f'{str(x_pos)} {str(y_pos)} {z_pos} '
        pos_str_1 = f'{str(x_pos+delta_vec[0])} {str(y_pos+delta_vec[1])} {z_pos+delta_vec[2]} '
        #select type if more than one possible
        L.L.command(
            f"create_atoms 0 single {pos_str_0} mol test_mol_0 1 rotate 0 0 0 1 units box ")
        L.L.command(f"group mol_{2*i} molecule {2*i+1}")
        L.L.command(
            f"create_atoms 0 single {pos_str_1} mol test_mol_1 2 rotate 0 0 0 1 units box ")
            #   rotate 0 0 0 0
        L.L.command(f"group mol_{2*i+1} molecule {2*i+2}")



def simulate_HP(db, keys, is_binded, temperature, minima_, bounding_energy_, K_, is_gif = False, is_print= False, sim_param = None, is_logger = False, is_PP = False, pre_path = ""):

    if sim_param is None:
        sim_param_ = lammps_util.sim_param.Parameter_simulation()
    else:
        sim_param_=sim_param
        
    #TODO move those they should not be needed
    n_node = 7
    n_out = 5
    sim_param_.minima = minima_
    sim_param_.bounding_energy = bounding_energy_
    sim_param_.K = K_

    k_h, k_p, k_s = keys
    pos_h = db.db_pairs[k_h]['struct_healthy']
    pos_p = db.db_pairs[k_h][k_p]['struct']
    translation = db.db_pairs[k_h][k_p]['translation']
    binded_hp = db.db_pairs[k_h][k_p]['binded_hp']
    inside_node = db.db_pairs[k_h][k_p]['sols'][k_s]['inside_node']
    conn = db.db_pairs[k_h][k_p]['sols'][k_s]['conn']
    n_out = len(pos_h)-1
    # load properties
    target_previous_pp_fitness = db.db_pairs[k_h][k_p]['pp_fitness']
    target_previous_hp_fitness = db.db_pairs[k_h][k_p]['hp_fitness']
    n_point_max = db.db_parameter["n_point_max"]
    min_r_squared_default = db.min_r_squared_default

    if  is_binded:

        if not 'hp_fit_dic' in  db.db_pairs[k_h][k_p].keys():
            full_healthy_hp, full_prion_hp, prion_binded, full_prion_pp_p1 = get_protein_pairs_binded_divided(pos_h, pos_p, binded_hp, inside_node, translation, target_previous_pp_fitness, target_previous_hp_fitness, min_r_squared_default, n_point_max, n_out)

        else:
            print("bypass in angle")
            hp_fit_dic = db.db_pairs[k_h][k_p]['hp_fit_dic']
            pp_fit_dic = db.db_pairs[k_h][k_p]['pp_fit_dic']
            hp_vec = hp_fit_dic['sol_ordered_hp'][-1], hp_fit_dic['binded_hp']
            pp_vec = pp_fit_dic['sol_ordered_pp'][-1], pp_fit_dic['binded_pp']

            print( 'get_protein_pairs_binded_divided_bypass' )
            full_healthy_hp, full_prion_hp, prion_binded, full_prion_pp_p1 =   get_protein_pairs_binded_divided_bypass(pos_h, pos_p, translation, inside_node, conn, hp_vec, pp_vec, n_out )
            print( 'get_protein_pairs_binded_divided_bypass out' )


        if is_PP:
            full_healthy_hp, full_prion_hp = prion_binded, full_prion_pp_p1
        
        all_conn = [(p[0]-1,p[1]-1) for p in conn] + [(i,i+1) for i in range(n_out-1)] + [(n_out-1,0)]
        all_conn = np.array(all_conn)

        protein_pair = protein_template.Protein_Template()
        protein_pair.set_prion_structure(full_prion_hp, connection = all_conn)
        protein_pair.set_healthy_structure(full_healthy_hp)

        h_c_x = sum(full_healthy_hp[i, 1] for i in range(len(full_healthy_hp))) / len(full_healthy_hp)
        h_c_y = sum(full_healthy_hp[i, 2] for i in range(len(full_healthy_hp))) / len(full_healthy_hp)
        p_c_x = sum(full_prion_hp[i,1] for i in range(len(full_prion_hp)))/len(full_prion_hp)
        p_c_y = sum(full_prion_hp[i,2] for i in range(len(full_prion_hp)))/len(full_prion_hp)

        delta = (p_c_x - h_c_x+0.1, p_c_y - h_c_y+0.1, 0)

    else:
        protein_pair = db.load_protein_pair(keys[0], keys[1], keys[2])

    #protein_pair.plot()
    #print(protein_pair.connection)
    #protein_pair_unbinded = db.load_protein_pair(keys[0], keys[1], keys[2])

    L_main = PyLammps()
    id_counter = mylammps.ID_counter()
    L = mylammps.MyLammps(L_main, id_counter)

    protein_pair.mass = (1, 1)
    protein_pair.K = sim_param_.K
    protein_pair.lj_param = []

    worker_id = os.getpid()
    # create a molecule file
    path_molecule_healthy = f"molecule_healthy_{worker_id}.lj"
    length_list_healthy = protein_pair.create_Lammps_molecule_file(path_molecule_healthy, type_mol="healthy", n_out=n_out)
    path_molecule_prion = f"molecule_prion_{worker_id}.lj"
    length_list_prion = protein_pair.create_Lammps_molecule_file(path_molecule_prion, type_mol="prion", n_out=n_out)

    # set up the simulation box
    set_up_molecule_simulation(L, [path_molecule_healthy, path_molecule_prion], length_list_prion, sim_param_.box_side, sim_param_.bounding_energy, sim_param_.minima,  sim_param_.cut_off, sim_param_.K)
    # add the molecule
    if is_binded:
        insert_avoiding_binded_molecule(L, sim_param_.n_molecule,  sim_param_.box_side, delta)
    else:
        insert_avoiding_molecule(L, sim_param_.n_molecule,  sim_param_.box_side, n_type=2)

    # set up logger 
    main_logger = logger.Logger(L, [protein_pair.healthy_position,
                                protein_pair.prion_position], id_struct_to_compare=["healty", "prion"], nb_molecule= sim_param_.n_molecule, n_node = n_node)
    main_logger.log(0, is_angle=False)


    ranint = random.randint(1,100000)
    name_data = f"_{keys[0]}_{keys[1]}_{keys[2]}_{is_binded}_{temperature}_{ranint}_"
    #run simulation
    end_time = L.run_brownian(temperature, sim_param_.damping, sim_param_.max_time, sim_param_.n_step,
                              sim_param_.n_plot, main_logger, is_angle=False, name_log = name_data, n_dump = N_DUMP, pre_path=pre_path)
    # get the connection
    conn = db.db_pairs[keys[0]][keys[1]]['sols'][keys[2]]['conn']
    conn_in = [[x[0]-1, x[1]-1] for x in conn]# get the connection from the inside
    conn_out = [(i, (i+1) % n_out) for i in range(n_out)] # get connection for outside node

    main_logger.set_connection(conn_in, conn_out)
    main_logger.log_simulation_param(sim_param_)
    main_logger.log_dynamic()
    conn_val, bounded_val, count_val = main_logger.get_dynamic()


    # CREATE THE GIF
    if is_gif:
        path_gif = f'gif_out/{sim_param_.n_molecule}_{is_binded}_{sim_param_.bounding_energy}_{temperature}_{sim_param_.minima}_{hash(keys)}_{is_PP}.gif'
        
        while exists(path_gif):
            path_gif=path_gif[0:-4]+"_"+str(1)+".gif"

        conn_to_plot = conn_val
        count_to_plot = count_val[0], count_val[1], count_val[2]
        bound_to_plot = bounded_val[2],  bounded_val[3], bounded_val[0]
        create_gif_multiple_molecule(
            path_gif, bound_to_plot,  count_to_plot, conn_to_plot, sim_param_.n_molecule,  sim_param_.box_side)

    #   Return data if no gif needed
    if is_print:
        return conn_val, bounded_val, count_val, 

    if is_logger:
        path_blender = f'{sim_param_.n_molecule}_{is_binded}_{sim_param_.bounding_energy}_{temperature}_{sim_param_.minima}_{hash(keys)}'
        main_logger.save_for_blender(path_blender, protein_pair, length_list_prion)
        return main_logger


def simulate( db: database.Database, worker_id, key_selected,
    num =  NB_TEMP, temperature_min = MIN_TEMP, temperature_max = MAX_TEMP,  PRE = PRE, is_PP = False, temperature_list = None, sim_param=None ):

    if sim_param is None:
        sim_param = lammps_util.sim_param.Parameter_simulation()
        sim_param.n_step = N_STEP
        sim_param.max_time = MAX_TIME
        sim_param.n_plot = N_PLOT

    if temperature_list is None:
        temperature_range = np.logspace( math.log(temperature_min), math.log(temperature_max), num=num, base=math.e)
    else:
        temperature_range = temperature_list

    temp_list, conn_list, bounded_list, count_list= [], [], [], []
    for temp in temperature_range:
        print("here1")
        conn_val, bounded_val, count_val = simulate_HP(db, key_selected, True,  temp, MINIMA, E_MIN, K_SPRING, sim_param= sim_param, is_gif= False, is_print= True, pre_path = PRE_DUMP, is_PP = is_PP )
        print("here2")
        temp_list.append(temp)
        conn_list.append(conn_val)
        bounded_list.append(bounded_val)
        count_list.append(count_val)

    temp_list_H, conn_list_H, bounded_list_H, count_list_H= [], [], [], []
    for temp in temperature_range:
        conn_val, bounded_val, count_val = simulate_HP(db, key_selected, False, temp, MINIMA, E_MIN, K_SPRING, sim_param=sim_param, is_gif= False, is_print= True, pre_path = PRE_DUMP, is_PP = is_PP)
        temp_list_H.append(temp)
        conn_list_H.append(conn_val)
        bounded_list_H.append(bounded_val)
        count_list_H.append(count_val)

    binded_energy = db.db_pairs[key_selected[0]][key_selected[1]]['sols'][key_selected[2]]['NEB_HP_PP_energy']

    svf = dynamic_file_template.Save_file_full_sim(key = key_selected,max_time= sim_param.max_time ,
    n_step=sim_param.n_step, n_plot= sim_param.n_plot, n_molecule= 15,box_side= 15, n_node = 7,n_out = 5 )   
    svf.set_physics_param( minima=MINIMA, bounding_energy=E_MIN, K=K_SPRING, damping = 1)
    svf.set_binding_energy([], binded_energy)
    svf.add_healthy_sim( temp_list_H, conn_list_H ,bounded_list_H ,count_list_H)
    svf.add_binded_sim( temp_list, conn_list, bounded_list, count_list)
    svf.save(name = str(worker_id))

    
            

if __name__ == "__main__":

    path = "/Users/mathieuouellet/Desktop/n=5_5deg_2_fig4b"
    db = database.Database(path)
    

    key_ps_pair = None # simulate all solution 
    list_temp = None # use the standard range defined in the header of this file 
    sim_param = None #use the standard param


    for npair, (key_h, key_p,key_s) in enumerate(key_ps_pair):
        print(key_h,key_p,key_s)
        worker_id = os.getpid()
        simulate( db, worker_id, (key_h,key_p,key_s), num =  NB_TEMP, temperature_min = MIN_TEMP, temperature_max = MAX_TEMP,  PRE = PRE, is_PP = False, temperature_list=  list_temp, sim_param=sim_param)


