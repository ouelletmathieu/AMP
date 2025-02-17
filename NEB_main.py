import subprocess
import database
import lammps_util.Nudged_Elastic_Band as NEB
import copy
import os
import math
import lammps_util.protein_template
import random

TIMESTEP = 0.001#1e-1
PRE = "AMP/src/polygon/"
PATH_SCR = "lammps_util/script/"
PATH_INIT = PRE + PATH_SCR + "test_healty.lj"
PATH_FINAL = PRE + PATH_SCR + "prion_neb_pos.lj"
PATH_NEB = PRE + PATH_SCR + "neb_task.lj"
PATH_DUMP = f'{PRE}temp/'
NB_NODE = 12
SIM_BOX_SIDE = 12
SPRING_INTER_REPLICA = 10000 #1.0


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

def _get_energy(protein_pair: lammps_util.protein_template.Protein_Template, initial_lenght_list, initial_type_list, new_atom_list, K):

    connections = protein_pair.connection
    energy = 0
    for i in range(connections.shape[0]):
        atom1 = connections[i, 0]
        atom2 = connections[i, 1]

        atom1_pos = (new_atom_list[0][atom1], new_atom_list[1][atom1])
        atom2_pos = (new_atom_list[0][atom2], new_atom_list[1][atom2])

        lenght = math.sqrt(
            (atom1_pos[0]-atom2_pos[0])**2 + (atom1_pos[1]-atom2_pos[1])**2)

        if initial_type_list[i] == 0:
            energy += K[0]*(SIM_BOX_SIDE*lenght - initial_lenght_list[i])**2
        else:
            energy += K[1]*(SIM_BOX_SIDE*lenght - initial_lenght_list[i])**2

    return energy

def get_energy(id, protein_pair, K):

    initial_lenght_list, initial_type_list = protein_pair.get_bond_length(
        type_mol="healthy")

    energy_list = []
    for i in range(1, NB_NODE+1):
        new_atom_list = NEB.Nudged_Elastic_Band.get_atom_position(
            PATH_DUMP+id+"dump.neb."+str(i))
        energy = _get_energy(protein_pair, initial_lenght_list,
                             initial_type_list, new_atom_list[-1], K)
        energy_list.append(energy)

    return energy_list

def get_NEB( db: database.Database,id, key_healthy, key_prion, key_sol, PRE = PRE, PATH_DUMP = PATH_DUMP, PATH_SCR=PATH_SCR, PATH_NEB= PATH_NEB, PATH_FINAL=PATH_FINAL ):
    # param for physical simulation
    mass = db.db_parameter['mass']
    K = db.db_parameter['K']

    # load protein
    protein_pair = db.load_protein_pair(key_healthy, key_prion, key_sol)
    protein_pair.set_interaction(mass, K, lj_param=[])

    fpath_h = f'{PATH_DUMP}h_{str(key_healthy)}_{str(key_prion)}_{id}.txt'
    #fpath_p = PATH_DUMP+"p_" + str(key_healthy)+"_"+str(key_prion)+"_"+id+".txt"

    protein_pair.create_Lammps_file(fpath_h, type_mol="healthy")
    #protein_pair.create_Lammps_file(fpath_p, type_mol="prion")

    # protein_pair.create_Lammps_file(p_init,"healthy")
    NEB.create_final_file(PATH_FINAL+id, protein_pair.prion_position)

    print(PATH_FINAL+id)

    NEB.create_neb_task_file(PATH_NEB+id, fpath_h,
                             PATH_FINAL+id, PATH_DUMP+id, TIMESTEP, minimize=True, spring_inter_replica= SPRING_INTER_REPLICA)
    print(PATH_NEB+id)

    path_script = PRE + PATH_SCR + "script_neb_" + id
    path_log = f'{PATH_DUMP}log.lammps_{id}'
    path_screen = f'{PATH_DUMP}screen_{id}'


    create_executable(path_script, path_log, path_log,
                      path_screen, PATH_NEB+id)

    print("exec  created")
    print(path_script)
    # run the script
    result = subprocess.run(
        ["sh", path_script], capture_output=True, text=True
    )
    print(result.stdout)
    print("result.stderr")
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

    return reaction_coordinate_list, energy_list, atom_list_list

def neb_fill_database(db: database.Database, worker_id, list_key_dict = None, overwrite = False):

    total_done = 0
    worker_id = str(worker_id)
    list_key_healthy = list(db.db_pairs.keys())
    random.shuffle(list_key_healthy)

    if list_key_dict is not None: #if ovveride the keys
        list_key_healthy = list_key_dict.keys()


    #for each healthy
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

            #for each potential prion
            for i_p, key_prion in enumerate(prion_key_list):

                if (not prion_key_dict is None) and prion_key_dict[key_prion] is not None: #if overide the keys
                    sol_key_list = prion_key_dict[key_prion]
                else:
                    sol_key_list,  sol_value = db.get_all_sol(
                        key_healthy, key_prion)

                #for each solution
                for i_s, key_sol in enumerate(sol_key_list):
                    
                    if not 'NEB_HP_energy' in db.db_pairs[key_healthy][key_prion]['sols'][key_sol] or overwrite:


                        reaction_coordinate_list, energy_list, atom_list_list = get_NEB( db,worker_id,
                            key_healthy, key_prion, key_sol)

                        db.add_NEB_HP(key_healthy, key_prion, key_sol,
                                    reaction_coordinate_list, energy_list, atom_list_list = atom_list_list)

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
    
    #can optionally gives key so that it does not compute for all the pairs in the databse but only some
    #list_key_dict = {1737806584325749351: {-8771402139612191995:[-4281357692739835561]}}

    path = "src/polygon/data/n=5_5deg_2"

    db = database.Database(path)
    worker_id = os.getpid()

    neb_fill_database(db, worker_id, list_key_dict=list_key_dict, overwrite=True)

    