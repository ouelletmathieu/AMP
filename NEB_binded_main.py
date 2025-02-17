import random
import subprocess
import database
import lammps_util.Nudged_Elastic_Band as NEB
import copy
import os
from lammps_util.util import Util

TIMESTEP = 0.001#1e-2
PRE = "AMP/src/polygon/"
PATH_SCR = "lammps_util/script/"
PATH_INIT = PRE + PATH_SCR + "test_healty.lj"
PATH_FINAL = PRE + PATH_SCR + "prion_neb_pos.lj"
PATH_NEB = PRE + PATH_SCR + "neb_task.lj"
PATH_DUMP = f'{PRE}temp/'
NB_NODE = 12
SIM_BOX_SIDE = 12
SPRING_INTER_REPLICA = 10000#1.0

def create_executable(script_path, log_path, plog_file, pscreen_file, in_file):

    with open(script_path, 'w') as f:
        f.write("#!/bin/bash \n")
        f.write("\n")
        f.write("export TMPDIR=/tmp \n")

        script = (
            f"mpirun -np {str(NB_NODE)} --oversubscribe lmp_mpi  -partition "
            + str(NB_NODE)
            + "x1 "
        )

        script += f" -log {log_path}"
        script += f"  -plog {plog_file}"
        script += f"  -pscreen {pscreen_file}"
        script += f"  -in {in_file}"
        f.write(script)

def check_valid_binding(binded_hp):
    if NEB.count_binding(binded_hp) != 1:
        print("binding not okay")
        print(binded_hp)
        return False
    return True

def load_value(db, id, key_healthy, key_prion, key_sol):
    return {"mass": db.db_parameter['mass'], "K": db.db_parameter['K'], "n_point_max": db.db_parameter["n_point_max"],
            "min_r_squared_default": db.min_r_squared_default, "pos_h": db.db_pairs[key_healthy]['struct_healthy'],
            "pos_p": db.db_pairs[key_healthy][key_prion]['struct'], "translation": db.db_pairs[key_healthy][key_prion]['translation'],
            "binded_hp": db.db_pairs[key_healthy][key_prion]['binded_hp'], "inside_node": db.db_pairs[key_healthy][key_prion]['sols'][key_sol]['inside_node'],
            "conn": db.db_pairs[key_healthy][key_prion]['sols'][key_sol]['conn'], "target_previous_pp_fitness": db.db_pairs[key_healthy][key_prion]['pp_fitness'],
            "target_previous_hp_fitness": db.db_pairs[key_healthy][key_prion]['hp_fitness']}

def get_NEB(db, id, key_healthy, key_prion, key_sol, PRE=PRE, PATH_DUMP=PATH_DUMP, PATH_SCR=PATH_SCR, PATH_NEB=PATH_NEB, PATH_FINAL=PATH_FINAL):

    value_dict = load_value(db, id, key_healthy, key_prion, key_sol)

    n_out = len(value_dict["pos_h"])-1

    # check binding
    if not check_valid_binding(value_dict["binded_hp"]):
        return False

    # load previous energy found
    target_previous_pp_fitness = db.db_pairs[key_healthy][key_prion]['pp_fitness']
    target_previous_hp_fitness = db.db_pairs[key_healthy][key_prion]['hp_fitness']

    # load protein
    protein_pair_unbinded = db.load_protein_pair(
        key_healthy, key_prion, key_sol)
    protein_pair_unbinded.set_interaction(
        value_dict["mass"], value_dict["K"], lj_param=[])
    fpath_binded = f'{PATH_DUMP}hp_{str(key_healthy)}_{str(key_prion)}_{id}.txt'


    if not 'hp_fit_dic' in  db.db_pairs[key_healthy][key_prion].keys():
        print("no bypass in fit")
        protein_pair_binded, dict_pos_in_protein, edges = NEB.get_protein_pairs_binded(value_dict["pos_h"], value_dict["pos_p"], value_dict["translation"], binded_hp=value_dict["binded_hp"], inside_node=value_dict["inside_node"], conn=value_dict["conn"], target_previous_pp_fitness=target_previous_pp_fitness,
                                                                            target_previous_hp_fitness=target_previous_hp_fitness, min_r_squared_default=value_dict["min_r_squared_default"], n_point_max=value_dict["n_point_max"])
    
    else:
        print("bypass in fit")
        hp_fit_dic = db.db_pairs[key_healthy][key_prion]['hp_fit_dic']
        pp_fit_dic = db.db_pairs[key_healthy][key_prion]['pp_fit_dic']
        hp_vec = hp_fit_dic['sol_ordered_hp'][-1], hp_fit_dic['binded_hp']
        pp_vec = pp_fit_dic['sol_ordered_pp'][-1], pp_fit_dic['binded_pp']


        protein_pair_binded, dict_pos_in_protein, edges = NEB.get_protein_pairs_binded_bypass(value_dict["pos_h"], value_dict["pos_p"], value_dict["translation"], inside_node=value_dict["inside_node"], conn=value_dict["conn"], hp_vec=hp_vec, pp_vec=pp_vec)
    
    
    """
    def get_protein_pairs_binded_bypass(pos_h, pos_p, translation, inside_node, conn, hp_vec, pp_vec):
    
    n_out = len(pos_h)-1
    best_sol_hp, binded_hp =  hp_vec
    best_sol_pp, binded_pp =  pp_vec

    """
    
    
    protein_pair_binded.set_interaction(
        value_dict["mass"], value_dict["K"], lj_param=[])

    # create files
    protein_pair_binded.create_Lammps_file(fpath_binded, type_mol="healthy")
    NEB.create_final_file(PATH_FINAL+id, protein_pair_binded.prion_position)
    NEB.create_neb_task_file(PATH_NEB+id, fpath_binded,
                             PATH_FINAL+id, PATH_DUMP+id, TIMESTEP,  minimize=True, spring_inter_replica= SPRING_INTER_REPLICA)

    # create executable
    path_script = PRE + PATH_SCR + "script_neb_binded_" + id
    path_log = f'{PATH_DUMP}log.lammps_binded_{id}'
    path_screen = f'{PATH_DUMP}screen_binded_{id}'
    create_executable(path_script, path_log, path_log,
                      path_screen, PATH_NEB+id)

    
    # run the script
    result = subprocess.run(
        ["sh", path_script], capture_output=True, text=True
    )

    print(result.stderr)

    # gather the reuslt
    reaction_coordinate_list, energy_list = NEB.get_energy(path_log)

    # gather the position
    atom_list_list = [
        NEB.get_atom_position(
            PATH_DUMP + id + "dump.neb." + str(i)
        )[-1]
        for i in range(1, NB_NODE + 1)
    ]

    return reaction_coordinate_list, energy_list, atom_list_list, protein_pair_binded, edges

def get_angle_difference(db, key_healthy, key_prion, key_sol, path_pos, NB_NODE=NB_NODE):

    value_dict = load_value(db, None, key_healthy, key_prion, key_sol)

    n_out = len(value_dict["pos_h"])-1
    # load the protein
    protein_pair_unbinded = db.load_protein_pair(
        key_healthy, key_prion, key_sol)
        
    if not 'hp_fit_dic' in  db.db_pairs[key_healthy][key_prion].keys():
        protein_pair_binded, dict_pos_in_protein, (in_edge,out_edge)  = NEB.get_protein_pairs_binded(value_dict["pos_h"], value_dict["pos_p"], value_dict["translation"], binded_hp=value_dict["binded_hp"], inside_node=value_dict["inside_node"], conn=value_dict["conn"], target_previous_pp_fitness=value_dict["target_previous_pp_fitness"],
                                                                            target_previous_hp_fitness=value_dict["target_previous_hp_fitness"], min_r_squared_default=value_dict["min_r_squared_default"], n_point_max=value_dict["n_point_max"])
    else:
        print("bypass in angle")
        hp_fit_dic = db.db_pairs[key_healthy][key_prion]['hp_fit_dic']
        pp_fit_dic = db.db_pairs[key_healthy][key_prion]['pp_fit_dic']
        hp_vec = hp_fit_dic['sol_ordered_hp'][-1], hp_fit_dic['binded_hp']
        pp_vec = pp_fit_dic['sol_ordered_pp'][-1], pp_fit_dic['binded_pp']
        protein_pair_binded, dict_pos_in_protein, edges = NEB.get_protein_pairs_binded_bypass(value_dict["pos_h"], value_dict["pos_p"], value_dict["translation"], inside_node=value_dict["inside_node"], conn=value_dict["conn"], hp_vec=hp_vec, pp_vec=pp_vec)
    
    
    # get translation
    binded_translated = {
        value_dict["translation"][p[0]]: p[1] for p in value_dict["binded_hp"]}

    pos_unbind = [[], []]
    for i in range(len(protein_pair_unbinded.healthy_position)):
        pos = protein_pair_unbinded.healthy_position[i]
        pos_unbind[0].append([pos[1], pos[2]])
        pos = protein_pair_unbinded.prion_position[i]
        pos_unbind[1].append([pos[1], pos[2]])

    # get the postiion for each reaction coordinate
    list_angle_bottom,  list_top_top = [], []
    list_pos_bottom,  list_pos_top = [], []

    for i in range(1, NB_NODE+1):
        all_pos = NEB.get_atom_position(path_pos+str(i))
        all_pos = [[all_pos[-1][0][i], all_pos[-1][1][i]]
                   for i in range(len(all_pos[-1][0]))]
        pos_bottom, pos_top = NEB.get_position_healthy_prion(
            all_pos, dict_pos_in_protein, binded_translated, len(protein_pair_unbinded.prion_position), n_out)
        bot_top_bind_angle = list(
            map(Util.get_angles, [pos_bottom, pos_top]))
        list_pos_bottom.append(pos_bottom)
        list_pos_top.append(pos_top)
        list_angle_bottom.append(bot_top_bind_angle[0])
        list_top_top.append(bot_top_bind_angle[1])

    list_angle_diff = []
    for i in range(NB_NODE):
        active_angle = [list_angle_bottom[i],
                        list_angle_bottom[i], list_top_top[i], list_top_top[i]]
        distance = list(map(Util.get_avg_angle_distance, [
                        list_angle_bottom[0], list_top_top[0], list_angle_bottom[-1], list_top_top[-1]], active_angle))
        list_angle_diff.append(distance)

    return list_pos_bottom, list_pos_top, list_angle_diff

def neb_fill_database(db, worker_id, list_key_dict = None, overwrite = False):
    # sourcery skip: avoid-builtin-shadow

    total_done = 0
    worker_id = str(worker_id)

    # shuffle list to minimize worker all checking the same file at the same time
    healthy_key_list = list(db.db_pairs.keys())
    random.shuffle(healthy_key_list)

    if list_key_dict is not None: #if ovveride the keys
        list_key_healthy = list_key_dict.keys()

    # for each healthy
    for i_h, key_healthy in enumerate(list_key_healthy):

        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)):
            db.lock(int(key_healthy))  # lock the healhy

            prion_key_dict = None
            if list_key_dict is not None: #if overide the keys
                prion_key_dict = list_key_dict[key_healthy]
                prion_key_list = prion_key_dict.keys()
            else:
                prion_key_list = copy.copy(list(db.db_pairs[key_healthy].keys()))
                prion_key_list.remove('struct_healthy')


            # for each potential prion
            for i_p, key_prion in enumerate(prion_key_list):
                if (not prion_key_dict is None) and prion_key_dict[key_prion] is not None: #if overide the keys
                    sol_key_list = prion_key_dict[key_prion]
                else:
                    sol_key_list,  sol_value = db.get_all_sol(
                        key_healthy, key_prion)

                # for each solution
                for i_s, key_sol in enumerate(sol_key_list):
                    

                    if not 'NEB_HP_PP_energy' in db.db_pairs[key_healthy][key_prion]['sols'][key_sol] or overwrite:

                        reaction_coordinate_list, energy_list, atom_list_list,protein_pair_binded, edges = get_NEB(db, worker_id, key_healthy, key_prion, key_sol)
                        path_pos = PATH_DUMP + str(worker_id) + "dump.neb." 
                        #value_dict = load_value(db, worker_id, key_healthy, key_prion, key_sol)
                        list_pos_bottom, list_pos_top, list_angle_diff = get_angle_difference(
                            db, key_healthy, key_prion, key_sol, path_pos)
                        db.add_NEB_HP_PP(key_healthy, key_prion, key_sol,
                                        reaction_coordinate_list, energy_list, (list_pos_bottom, list_pos_top), list_angle_diff, atom_list_list=atom_list_list)
                        total_done += 1
                        print(i_h, i_p, i_s, total_done)

                    else:
                        print("already done")

                directory = os.listdir(PATH_DUMP)
                for fname in directory:
                    if (
                        os.path.isfile(PATH_DUMP + os.sep + fname)
                        and worker_id in fname
                    ):
                        os.remove(PATH_DUMP + os.sep + fname)

            db.unlock(int(key_healthy))  # unlock the healthy

if __name__ == "__main__":

    
    list_key_dict = None
    #list_key_dict = {1737806584325749351: {-8771402139612191995:[-4281357692739835561]}}

    path = "src/polygon/data/n=5_5deg_2"

    db = database.Database(path)
    worker_id = os.getpid()
    neb_fill_database(db, worker_id, list_key_dict=list_key_dict, overwrite=True)
