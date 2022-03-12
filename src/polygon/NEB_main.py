import subprocess
import database 
import lammps_util.Nudged_Elastic_Band as NEB
import copy
import os
import math
import lammps_util.protein_template

TIMESTEP   = 1e-1
PRE        = "src/polygon/"
PATH_SCR   = "lammps_util/script/"
PATH_INIT  = PRE + PATH_SCR  +  "test_healty.lj"
PATH_FINAL = PRE + PATH_SCR  +  "prion_neb_pos.lj"
PATH_NEB   = PRE + PATH_SCR  +  "neb_task.lj"
PATH_DUMP  = PRE + "temp/"
NB_NODE = 12
SIM_BOX_SIDE = 12

def create_executable(script_path, log_path, plog_file, pscreen_file, in_file):

    with open(script_path,'w') as f:
        f.write("#!/bin/bash \n")
        f.write("\n")
        f.write("export TMPDIR=/tmp \n")

        script = "mpirun -np "+str(NB_NODE)+" --oversubscribe lmp_mpi  -partition " +str(NB_NODE)+"x1 "
        script+= " -log " + log_path
        script+= "  -plog " + plog_file
        script+= "  -pscreen " + pscreen_file
        script+= "  -in " + in_file
        f.write(script)


def _get_energy(protein_pair:lammps_util.protein_template.Protein_Template, initial_lenght_list, initial_type_list, new_atom_list, K):

        connections = protein_pair.connection
        energy = 0
        for i in range(connections.shape[0]):
            atom1 = connections[i,0]
            atom2 = connections[i,1]

            atom1_pos =  (new_atom_list[0][atom1],new_atom_list[1][atom1])
            atom2_pos =  (new_atom_list[0][atom2],new_atom_list[1][atom2])        
            
            lenght = math.sqrt((atom1_pos[0]-atom2_pos[0])**2 + (atom1_pos[1]-atom2_pos[1])**2)
            
            if initial_type_list[i]==0:
                energy+= K[0]*(SIM_BOX_SIDE*lenght - initial_lenght_list[i])**2
            else:
                energy+= K[1]*(SIM_BOX_SIDE*lenght -initial_lenght_list[i])**2
        
        return energy



def get_energy(id, protein_pair, K):

    initial_lenght_list, initial_type_list = protein_pair.get_bond_length(type_mol = "healthy")

    energy_list = []
    for i in range(1,NB_NODE+1):
        new_atom_list = NEB.Nudged_Elastic_Band.get_atom_position(PATH_DUMP+id+"dump.neb."+str(i))
        energy = _get_energy(protein_pair, initial_lenght_list, initial_type_list, new_atom_list[-1], K)
        energy_list.append(energy)
    
    return energy_list
        

def neb_fill_database(db: database.Database, worker_id):

    total_done = 0
    id = str(worker_id)

    for i_h, key_healthy in enumerate(db.db_pairs.keys()):
        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)):
            db.lock(int(key_healthy))  # lock the healhy
            list_prion_key = copy.copy(list(db.db_pairs[key_healthy].keys()))
            list_prion_key.remove('struct_healthy')

            for i_p, key_prion in enumerate(list_prion_key):

                list_sol_key,  sol_value = db.get_all_sol(key_healthy, key_prion)
                for i_s, key_sol in enumerate(list_sol_key):
                    if 'RMS_hp' in db.db_pairs[key_healthy][key_prion]['sols'][key_sol]:
                        #param for physical simulation 
                        mass = db.db_parameter['mass']
                        K = db.db_parameter['K']

                        # load protein
                        protein_pair = db.load_protein_pair(key_healthy, key_prion, key_sol)
                        protein_pair.set_interaction(mass, K, lj_param=[])

                        fpath_h = PATH_DUMP+"h_" + str(key_healthy)+"_"+str(key_prion)+"_"+id+".txt"
                        #fpath_p = PATH_DUMP+"p_" + str(key_healthy)+"_"+str(key_prion)+"_"+id+".txt"

                        protein_pair.create_Lammps_file(fpath_h, type_mol="healthy")
                        #protein_pair.create_Lammps_file(fpath_p, type_mol="prion")
                        
                        #protein_pair.create_Lammps_file(p_init,"healthy")
                        NEB.Nudged_Elastic_Band.create_final_file(PATH_FINAL+id, protein_pair.prion_position)
                        NEB.Nudged_Elastic_Band.create_neb_task_file(PATH_NEB+id, fpath_h, PATH_FINAL+id, PATH_DUMP+id, TIMESTEP)

                        path_script =  PRE + PATH_SCR +  "script_neb_" + id
                        path_log    =  PATH_DUMP      +  "log.lammps_" + id
                        path_screen =  PATH_DUMP      +  "screen_" + id
                        
                        create_executable(path_script, path_log, path_log, path_screen, PATH_NEB+id)

                        #run the script 
                        result = subprocess.run(
                            ["sh", path_script], capture_output=True, text=True
                        )
                        print(result.stderr)
                        #gather the reuslt 
                        reaction_coordinate_list, energy_list = NEB.Nudged_Elastic_Band.get_energy(path_log)

                        db.add_NEB_HP(key_healthy, key_prion, key_sol, reaction_coordinate_list, energy_list)

                        #
                        #energy_recomputed = get_energy(id, protein_pair, K)
                        #for x,a,b in zip(reaction_coordinate_list,energy_list,energy_recomputed):
                        #    print(x,a,b)
                        #assert False
                        print(i_h, i_p, i_s, total_done)
                        total_done+=1

                    """
                    directory = os.listdir(PATH_DUMP)
                    for fname in directory:
                        if os.path.isfile(PATH_DUMP + os.sep + fname):
                            # Full path
                            if id in fname:
                                os.remove(PATH_DUMP + os.sep + fname)
                    """    


            db.unlock(int(key_healthy))  # unlock the healthy



if __name__ == "__main__":
    db = database.Database("src/polygon/data/n=5_5deg_2")
    worker_id = os.getpid()
    neb_fill_database(db, worker_id)

    