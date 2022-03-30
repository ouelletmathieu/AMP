import random
import subprocess
import database
import lammps_util.Nudged_Elastic_Band as NEB
import copy
import os
import math
import lammps_util.protein_template
from lammps_util.util import Util

TIMESTEP = 1e-1
PRE = "src/polygon/"
PATH_SCR = "lammps_util/script/"
PATH_INIT = PRE + PATH_SCR + "test_healty.lj"
PATH_FINAL = PRE + PATH_SCR + "prion_neb_pos.lj"
PATH_NEB = PRE + PATH_SCR + "neb_task.lj"
PATH_DUMP = f'{PRE}temp/'
NB_NODE = 12
SIM_BOX_SIDE = 12


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


def neb_fill_database(db: database.Database, worker_id):
    # sourcery skip: avoid-builtin-shadow

    total_done = 0
    id = str(worker_id)
    healthy_key_list = list(db.db_pairs.keys())
    random.shuffle(healthy_key_list)

    for i_h, key_healthy in enumerate(healthy_key_list):
        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)):
            db.lock(int(key_healthy))  # lock the healhy
            list_prion_key = copy.copy(list(db.db_pairs[key_healthy].keys()))
            list_prion_key.remove('struct_healthy')

            for i_p, key_prion in enumerate(list_prion_key):

                list_sol_key,  sol_value = db.get_all_sol(
                    key_healthy, key_prion)

                for i_s, key_sol in enumerate(list_sol_key):

                    # param for physical simulation
                    mass, K = db.db_parameter['mass'], db.db_parameter['K']
                    pos_h = db.db_pairs[key_healthy]['struct_healthy']
                    pos_p = db.db_pairs[key_healthy][key_prion]['struct']
                    translation = db.db_pairs[key_healthy][key_prion]['translation']
                    binded_hp = db.db_pairs[key_healthy][key_prion]['binded_hp']
                    inside_node = db.db_pairs[key_healthy][key_prion]['sols'][key_sol]['inside_node']
                    conn = db.db_pairs[key_healthy][key_prion]['sols'][key_sol]['conn']
                    n_out = len(pos_h)-1

                    if NEB.count_binding(binded_hp) != 1:
                        print("binding not okay")
                        print(binded_hp)
                        continue

                    # load properties
                    target_previous_pp_fitness = db.db_pairs[key_healthy][key_prion]['pp_fitness']
                    target_previous_hp_fitness = db.db_pairs[key_healthy][key_prion]['hp_fitness']
                    n_point_max = db.db_parameter["n_point_max"]
                    min_r_squared_default = db.min_r_squared_default

                    # load protein
                    protein_pair_unbinded = db.load_protein_pair(
                        key_healthy, key_prion, key_sol)
                    protein_pair_unbinded.set_interaction(mass, K, lj_param=[])

                    fpath_binded = f'{PATH_DUMP}hp_{str(key_healthy)}_{str(key_prion)}_{id}.txt'
                    protein_pair_binded, dict_pos_in_protein = NEB.get_protein_pairs_binded(pos_h, pos_p, translation, binded_hp=binded_hp, inside_node=inside_node, conn=conn, target_previous_pp_fitness=target_previous_pp_fitness,
                                                                                            target_previous_hp_fitness=target_previous_hp_fitness, min_r_squared_default=min_r_squared_default, n_point_max=n_point_max)
                    protein_pair_binded.set_interaction(mass, K, lj_param=[])
                    protein_pair_binded.create_Lammps_file(
                        fpath_binded, type_mol="healthy")
                    NEB.create_final_file(
                        PATH_FINAL+id, protein_pair_binded.prion_position)
                    NEB.create_neb_task_file(
                        PATH_NEB+id, fpath_binded, PATH_FINAL+id, PATH_DUMP+id, TIMESTEP, minimize=True)

                    path_script = PRE + PATH_SCR + "script_neb_binded_" + id
                    path_log = f'{PATH_DUMP}log.lammps_binded_{id}'
                    path_screen = f'{PATH_DUMP}screen_binded_{id}'

                    create_executable(path_script, path_log,
                                      path_log, path_screen, PATH_NEB+id)

                    # run the script
                    result = subprocess.run(
                        ["sh", path_script], capture_output=True, text=True
                    )
                    print(result.stderr)
                    # gather the reuslt
                    reaction_coordinate_list, energy_list = NEB.get_energy(
                        path_log)
                    print(i_h, i_p, i_s, total_done)
                    for i in range(len(reaction_coordinate_list)):
                        print(
                            f'{str(reaction_coordinate_list[i])} , {str(energy_list[i])}')

                    #pos_unbind = [ [p[1],p[2]] for p in protein_pair_unbinded.healthy_position], [ [p[1],p[2]] for p in protein_pair_unbinded.prion_position]

                    # def get_atom_position(output_file_path):
                    binded_translated = {
                        translation[p[0]]: p[1] for p in binded_hp}

                    pos_unbind = [[], []]
                    for i in range(len(protein_pair_unbinded.healthy_position)):
                        pos = protein_pair_unbinded.healthy_position[i]
                        pos_unbind[0].append([pos[1], pos[2]])
                        pos = protein_pair_unbinded.prion_position[i]
                        pos_unbind[1].append([pos[1], pos[2]])

                    ang_un = list(map(Util.get_angles, pos_unbind))

                    # get the postiion for each reaction coordinate
                    list_angle_bottom,  list_top_top = [], []
                    list_pos_bottom,  list_pos_top = [], []
                    for i in range(1, 13):
                        path_pos = PATH_DUMP + id + "dump.neb." + str(i)
                        all_pos = NEB.get_atom_position(path_pos)

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


                    for i in range(12):
                        active_angle = [
                            list_angle_bottom[i], list_angle_bottom[i], list_top_top[i], list_top_top[i]]

                        distance = list(map(Util.get_avg_angle_distance, [
                                        list_angle_bottom[0], list_top_top[0], list_angle_bottom[-1], list_top_top[-1]], active_angle))
                        list_angle_diff.append(distance)

                    db.add_NEB_HP_PP(key_healthy, key_prion, key_sol,
                                     reaction_coordinate_list, energy_list, (list_pos_bottom, list_pos_top), list_angle_diff)



                directory = os.listdir(PATH_DUMP)
                for fname in directory:
                    if (
                        os.path.isfile(PATH_DUMP + os.sep + fname)
                        and id in fname
                    ):
                        os.remove(PATH_DUMP + os.sep + fname)

            db.unlock(int(key_healthy))  # unlock the healthy


if __name__ == "__main__":
    db = database.Database("/Users/mathieuouellet/Desktop/n=5_5deg_2")
    worker_id = os.getpid()
    neb_fill_database(db, worker_id)
