import math
import numpy as np
import database
import copy
import lammps_util.my_lammps as mylammps
from lammps_util.util import Util
import lammps_util.logger as logger
from lammps import PyLammps
import lammps
import scipy.stats
import contextlib


SIM_PER_TEMP_SECOND_RANGE = 3
STABILITY_DIST_RATIO = 0.1
MIN_POINT = 6
DAMPING = 1



def run_time(fpath_molecule, id_compare_to, protein_pair, temperature, maxtime, square_angle_sum):
    
    with  contextlib.redirect_stdout(None):

        lmp = lammps.lammps(cmdargs=[ "-log", "none"])
        L_main = PyLammps(ptr=lmp)
        id_counter = mylammps.ID_counter()
        L = mylammps.MyLammps(L_main, id_counter)
        L.create_molecule_2d(fpath_molecule)
        main_logger = logger.Logger(L, [protein_pair.healthy_position,
                                    protein_pair.prion_position], id_struct_to_compare=["healthy", "prion"])
        #print(id_counter.next(), protein_pair, temperature, 1, maxtime, 10000, maxtime, main_logger,id_compare_to, square_angle_sum)
        time = L.run_brownian(temperature, DAMPING, maxtime, 1000, maxtime, main_logger, stop_change_conf=True,
                            id_conf=id_compare_to, square_angle_sum=square_angle_sum, random_seed=True)
        
    return time


def simulate_stability(fpath_molecule, id_compare_to, protein_pair, start_temp, max_temp, n_test, maxtime, square_angle_sum):

    list_temp_all, list_time_all = [], []
    active_temp = start_temp
    # set up the first range to find the limit of fully stable, aka minimum temperature
    delta_1 = math.pow(max_temp, 1/(n_test-1)) / \
        math.pow(start_temp, 1/(n_test-1))
    list_temp_1 = [start_temp*(math.pow(delta_1, i)) for i in range(n_test)]
    counter_max_time_reached = 0  # count the number of time we reach the max
    # run for the first range
    for temperature in reversed(list_temp_1):
        time = run_time(fpath_molecule, id_compare_to,
                        protein_pair, temperature, maxtime, square_angle_sum)
        list_time_all.append(time)
        list_temp_all.append(temperature)
        # if we reach a point where it is possible to pass the maxtime without changing conformation
        if time == maxtime:
            counter_max_time_reached += 1
        if counter_max_time_reached == 1:
            break
        active_temp = temperature
    # set up the second range
    delta_2 = math.pow(max_temp, 1/(n_test-1)) / \
        math.pow(active_temp, 1/(n_test-1))
    list_temp_2 = [active_temp*(math.pow(delta_2, i)) for i in range(n_test)]
    for temp in reversed(list_temp_2):
        for i in range(SIM_PER_TEMP_SECOND_RANGE):
            time = run_time(fpath_molecule, id_compare_to,
                            protein_pair, temperature, maxtime, square_angle_sum)
            list_time_all.append(time)
            list_temp_all.append(temp)

    return list_temp_all, list_time_all


def get_fit(list_temp_all, list_time_all, maxtime):
    print("starting fit")
    if len(list_temp_all) >= MIN_POINT:
        new_t = []
        new_dat = []
        for i in range(len(list_temp_all)):
            if list_time_all[i]<maxtime:
               new_t.append(list_temp_all[i]) 
               new_dat.append(list_time_all[i])

        # put in a dict to get the average
        """
        dct_temp = dict()
        for i, t in enumerate(list_temp_all):
            # filter the one that went over equal to maxtime
            if list_time_all[i] != 0 and list_time_all[i] < maxtime:
                if 1/t in dct_temp:
                    dct_temp[1/t].append(list_time_all[i])
                else:
                    dct_temp[1/t] = [list_time_all[i]]
        new_t = [temp for temp in dct_temp.keys()]
        new_dat = [sum(dct_temp[temp])/len(dct_temp[temp])
                   for temp in dct_temp.keys()]
        """
        m1, b1, r_value, p_value, std_err = scipy.stats.linregress(
            new_t, new_dat)

        print(m1, b1, r_value**2)
        for i in range(len(list_temp_all)):
            print(list_temp_all[i], list_time_all[i])

        assert False
        return m1, b1, r_value**2
    else:
        print("not enough point")
        return 0, 0, 0


def fill_stability_time_to_out(db: database.Database):

    for i_h, key_healthy in enumerate(db.db_pairs.keys()):
        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)):
            db.lock(int(key_healthy))  # lock the healhy
            list_prion_key = copy.copy(list(db.db_pairs[key_healthy].keys()))
            list_prion_key.remove('struct_healthy')

            for i_p, key_prion in enumerate(list_prion_key):

                list_sol_key = db.get_all_sol_key(key_healthy, key_prion)
                for i_s, key_sol in enumerate(list_sol_key):

                    # TODO add finder of iteration that works
                    # load simulation parameter
                    mass = db.db_parameter['mass']
                    K = db.db_parameter['K']
                    # load protein
                    protein_pair = db.load_protein_pair(
                        key_healthy, key_prion, key_sol)
                    protein_pair.set_interaction(mass, K, lj_param=[])
                    # create file path for lammps to read the file
                    fpath_dir = "src/polygon/temp/"
                    fpath_h = fpath_dir+"h_" + \
                        str(key_healthy)+"_"+str(key_prion)+".txt"
                    fpath_p = fpath_dir+"p_" + \
                        str(key_healthy)+"_"+str(key_prion)+".txt"
                    # save the molecule files
                    protein_pair.create_Lammps_file(
                        fpath_h, type_mol="healthy")
                    protein_pair.create_Lammps_file(fpath_p, type_mol="prion")
                    # compute the distance needed to consider we escaped
                    angles_healthy = Util.get_angles(
                        protein_pair.healthy_position, len(protein_pair.healthy_position))
                    angles_prion = Util.get_angles(
                        protein_pair.prion_position, len(protein_pair.healthy_position))
                    avg_distance = Util.get_avg_angle_distance(
                        angles_healthy, angles_prion)
                    square_angle_sum = STABILITY_DIST_RATIO * \
                        (avg_distance**2)  # put the ratio in the db
                    # parameter for the simulation and range of temperature
                    n_test = db.db_parameter['n_test_temperature']
                    start_temp = db.db_parameter['start_temp']
                    max_temp = db.db_parameter['max_temp']
                    maxtime = db.db_parameter['maxtime']
                    # get list of temperature and max stable time for both
                    list_temp_all_healthy, list_time_all_healthy = simulate_stability(
                        fpath_h, "healthy", protein_pair, start_temp, max_temp, n_test, maxtime, square_angle_sum)
                    list_temp_all_prion, list_time_all_prion = simulate_stability(
                        fpath_p, "prion", protein_pair, start_temp, max_temp, n_test, maxtime, square_angle_sum)
                    # get the best fit
                    m1_h, b1_h, rs2_h = get_fit(
                        list_temp_all_healthy, list_time_all_healthy, maxtime)
                    m1_v, b1_v, rs2_v = get_fit(
                        list_temp_all_prion, list_time_all_prion, maxtime)

                    db.add_stability(key_healthy, key_prion, key_sol,
                                     (m1_h, b1_h, rs2_h), (m1_v, b1_v, rs2_v))
                    print(i_h, i_p, i_s)

            db.unlock(int(key_healthy))  # unlock the healthy


def run_time_angle(fpath_molecule, protein_pair, temperature, maxtime):
    
    with  contextlib.redirect_stdout(None):

        lmp = lammps.lammps(cmdargs=[ "-log", "none"])
        L_main = PyLammps(ptr=lmp)
        id_counter = mylammps.ID_counter()
        L = mylammps.MyLammps(L_main, id_counter)
        L.create_molecule_2d(fpath_molecule)
        main_logger = logger.Logger(L, [protein_pair.healthy_position,
                                    protein_pair.prion_position], id_struct_to_compare=["healthy", "prion"])
        #print(id_counter.next(), protein_pair, temperature, 1, maxtime, 10000, maxtime, main_logger,id_compare_to, square_angle_sum)
        time = L.run_brownian(temperature, DAMPING, maxtime, 1000, maxtime, main_logger, stop_change_conf=False, random_seed=True)
        avg_angle = main_logger.get_RMS_angle_difference_all_time()
        
    return  avg_angle


def fill_stability_average_distance(db: database.Database):

    for i_h, key_healthy in enumerate(db.db_pairs.keys()):
        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)):
            db.lock(int(key_healthy))  # lock the healhy
            list_prion_key = copy.copy(list(db.db_pairs[key_healthy].keys()))
            list_prion_key.remove('struct_healthy')

            for i_p, key_prion in enumerate(list_prion_key):

                list_sol_key,  sol_value = db.get_all_sol(key_healthy, key_prion)
                for i_s, key_sol in enumerate(list_sol_key):

                    # TODO add finder of iteration that works
                    # load simulation parameter
                    mass = db.db_parameter['mass']
                    K = db.db_parameter['K']
                    # load protein
                    protein_pair = db.load_protein_pair(
                        key_healthy, key_prion, key_sol)
                    protein_pair.set_interaction(mass, K, lj_param=[])
                    # create file path for lammps to read the file
                    fpath_dir = "src/polygon/temp/"
                    fpath_h = fpath_dir+"h_" + \
                        str(key_healthy)+"_"+str(key_prion)+".txt"
                    fpath_p = fpath_dir+"p_" + \
                        str(key_healthy)+"_"+str(key_prion)+".txt"
                    # save the molecule files
                    protein_pair.create_Lammps_file(
                        fpath_h, type_mol="healthy")
                    protein_pair.create_Lammps_file(fpath_p, type_mol="prion")
                    # compute the distance needed to consider we escaped
                    angles_healthy = Util.get_angles(
                        protein_pair.healthy_position, len(protein_pair.healthy_position))
                    angles_prion = Util.get_angles(
                        protein_pair.prion_position, len(protein_pair.healthy_position))
                    avg_distance = Util.get_avg_angle_distance(
                        angles_healthy, angles_prion)

                    # parameter for the simulation and range of temperature
                    n_test = db.db_parameter['n_test_temperature']
                    start_temp = db.db_parameter['start_temp']
                    max_temp = db.db_parameter['max_temp']
                    maxtime = db.db_parameter['maxtime']

                    delta_1 = math.pow(max_temp, 1/(n_test-1)) / math.pow(start_temp, 1/(n_test-1))
                    list_temp_1 = [start_temp*(math.pow(delta_1, i)) for i in range(n_test)]
                    healthy_avg, prion_avg = [], []
                    for temp in list_temp_1:
                        sum_avg_angle_prion = 0
                        sum_avg_angle_healthy = 0
                        for i in range(3):
                            sum_avg_angle_healthy+=run_time_angle(fpath_h, protein_pair, temp, maxtime)
                            sum_avg_angle_prion+=run_time_angle(fpath_p, protein_pair, temp, maxtime)

                        healthy_avg.append(sum_avg_angle_healthy/3)
                        prion_avg.append(sum_avg_angle_prion/3)
                    
                    db.add_RMS_stability(key_healthy, key_prion, key_sol, list_temp_1,  healthy_avg, prion_avg, avg_distance)
                    #db.add_stability(key_healthy, key_prion, key_sol,(m1_h, b1_h, rs2_h), (m1_v, b1_v, rs2_v))
                    print(i_h, i_p, i_s)

            db.unlock(int(key_healthy))  # unlock the healthy

if __name__ == "__main__":
    db = database.Database("src/polygon/data/n=5_5deg_2")
    #fill_stability_time_to_out(db)
    fill_stability_average_distance(db)