import subprocess
import database 
import lammps_util.Nudged_Elastic_Band as NEB


def main():
    #should be called from terminal from the AMP directory (the one containing src)

    ##########################
    #          PARAM         #
    ##########################
    dict_name = "n=5_5deg_2"
    #test param
    key_pairs = (-3087218339716488137, -1917219715052742202)
    

    #fixed param
    timestep = 1e-1
    pre = "src/polygon/"
    path_scr = "lammps_util/script/"
    p_init = pre+path_scr  + "test_healty.lj"
    p_final = pre+path_scr + "prion_neb_pos.lj"
    p_neb = pre+path_scr   + "neb_task.lj"
    path_dump = pre+ "temp/"

    #load protein pair
    db = database.Database(pre+"data/"+dict_name)
    protein_pair = db.load_protein_pair(key_pairs[0], key_pairs[1])
    #TODO move to the database the parameters
    protein_pair.set_interaction( mass=(1,1) , K=(1,1), lj_param=[])

    #protein_pair.create_Lammps_file(p_init,"healthy")
    NEB.Nudged_Elastic_Band.create_final_file(p_final, protein_pair.prion_position)
    NEB.Nudged_Elastic_Band.create_neb_task_file(p_neb, p_init, p_final, path_dump, timestep )

    result = subprocess.run(
        ["sh", pre + path_scr + "script_neb"], capture_output=True, text=True
    )

    print(result)


    print("done")

if __name__ == "__main__":
    main()

    