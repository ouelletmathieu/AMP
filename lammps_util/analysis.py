from lammps import PyLammps
from .my_lammps import MyLammps, ID_counter
from .logger import Logger
from .error import ConvergenceError


def find_appropriate_time_step(max_exponent, file_path, temperature, damping, max_time, healthy_positions, prion_positions, square_angle_sum=25):
    """find the time step needed to have at least a simulation that converge up to max_time 

    Args:
        max_exponent ([int]): the maximum step that can be considered such that max_step = 10**max_exponent
        file_path ([string]): file_path 
        temperature ([float]): temperature for the langevin simulation
        damping ([float]): damping for the langevin simulation 
        max_time ([float]): maximum time 
        healthy_positions (np.array([a , x, y ])): position of the atoms in the healthy protein. atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside
        prion_positions (np.array([a , x, y ])): position of the atoms in the prion protein. atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside

    Returns:
        [int]: number of step needed
    """
    # TODO move this function to somewhere else
    max_temp = temperature
    n_step = 0

    for i in range(3, max_exponent):

        converged = True

        for nbtry in range(3):
            n_step = 10**i
            n_plot = 30

            L_main = PyLammps()
            id_counter = ID_counter()
            L = MyLammps(L_main, id_counter)
            print("\n\n " + file_path + "\n\n ")
            L.create_molecule_2d(file_path)
            main_logger = Logger(
                L, [healthy_positions, prion_positions], id_struct_to_compare=["healthy", "prion"])
            main_logger.log(0)
            time = L.run_brownian(temperature, damping, max_time, n_step, n_plot, main_logger,
                                  stop_change_conf=False, id_conf=[], square_angle_sum=square_angle_sum, random_seed=True)

            if time == -1:
                converged = False

        if converged:
            return n_step

def get_max_stable_temp(min_temp, max_temp, simulate, max_test, max_time, target_precision, max_try):
    """this function find the maximum temperature for which the molecule is stable for at least the time max_time
        at a given precision given by target_precision. A molecule is considered stable if we can repeat the experiment 
        max_test number of time. Is the code get stuck in an unstable area max_time need to be lowered

    Args:
        min_temp ([type]): minimum temperature to consider (need to be stable at that temperature)
        max_temp ([type]): maximum temperature to consider  (need to be unstable at that temperature)
        simulate ([type]): function that depend on the temperature and return time 
        max_test ([type]): number of simulation needed to be considered stable
        max_time ([type]): time used to simulate the molecule. A molecule in the same conformation after max_time 
                            is considered stable 
        target_precision ([float]): Max uncertainty accepted on the temperature.
        max_try ([int]): Maximum number of time max_test can be increased

    Raises:
        ValueError: If max temperature is stable increase max temp

    Returns:
        [type]: highest stable temperature
    """
    # check if unstable at max temp
    max_unstable = False
    n_test = 0
    while n_test < max_test and not max_unstable:
        time = simulate(max_temp)
        if time < max_time:
            max_unstable = True
        n_test += 1
    if not max_unstable:
        print("max temperature is stable increase max temp")
        raise ValueError("max temperature is stable increase max temp")

    # try to find the hisghest stable temp
    try:
        highest_stable_temp = _go_lower(
            min_temp, max_temp, simulate, max_test, max_time, target_precision, max_try, 0)
    except ValueError:
        return None

    return highest_stable_temp

def _go_lower(highest_stable, lowest_unstable, simulate, max_test, max_time, target_precision, max_try, try_nb):
    """method that divide the search interval in two to find the when the stability change happen
    """
    # check if we are not stuck in a false spot
    try_nb += 1
    if try_nb > max_try:
        raise ValueError("did not converge increase the max_test")

    # test if stable up to our criterion
    unstable = False
    n_test = 0
    new_temp = (highest_stable+lowest_unstable)/2
    while n_test < max_test and not unstable:
        time = simulate(new_temp)
        if time == -1:
            raise ConvergenceError(
                "did not converge increase the number of steps")

        if time < max_time:
            unstable = True
        n_test += 1

    # compute precision
    precision = lowest_unstable - highest_stable

    # go lower if needed (unstable or not enough precision)
    if unstable:
        print("unstable: [ ", highest_stable, new_temp, " ]", sep=",")
        return _go_lower(highest_stable, new_temp, simulate, max_test, max_time, target_precision, max_try, try_nb)
    else:
        if precision <= target_precision:
            return new_temp
        else:
            print("precision: [ ", new_temp, lowest_unstable, " ]", sep=",")
            return _go_lower(new_temp, lowest_unstable, simulate, max_test, max_time, target_precision, max_try, try_nb)

def get_protein_max_stability_temp(file_path, healthy_positions, prion_positions, to_check, max_temp_totest, square_angle_sum,  max_time=30, n_plot=30, nb_needed_for_stability=10,  max_try=5):
    """find the maximal temperature for the healthy and prion protein given a maximal distance (square_angle_sum)

    Args:
        file_path ([string]): file_path 
        healthy_positions (np.array([a , x, y ])): position of the atoms in the healthy protein. atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside
        prion_positions (np.array([a , x, y ])): position of the atoms in the prion protein. atom position list in the form np.array([a , x, y ]) where a = 0 if outside and a=1 if inside

        max_temp_totest ([type]): max temperature to test for stability (need to be unstable)
        square_angle_sum ([type]): the distance needed to consider that the conformer is changed
        max_time (int, optional): Time to simulate. Defaults to 30.
        n_plot (int, optional): Number of time the conformer is checked . Defaults to 30.
        nb_needed_for_stability (int, optional): Number of full simulation needed to consider the compound to be stable Defaults to 10.
        max_try (int, optional): maximum number of try of the full process. Defaults to 5.

    Returns:
        [type]: [description]
    """
    if to_check not in ['prion', 'healthy']:
        raise ValueError("to_check is either prion or healthy")

    damping = 1
    min_temp_totest = 0
    target_precision = 0.01*max_temp_totest
    n_step = find_appropriate_time_step(10, file_path,  max_temp_totest, damping, max_time,  healthy_positions, prion_positions)
    n_step = n_step

    nb_try = 0

    # define the simulation function
    # TODO should not start Lammps each time and rebuild everything
    def simulate(temperature):
        L_main = PyLammps()
        id_counter = ID_counter()
        L = MyLammps(L_main, id_counter)
        L.create_molecule_2d(file_path)
        main_logger = Logger(
            L, [healthy_positions, prion_positions], id_struct_to_compare=["healthy", "prion"])
        main_logger.log(0)
        time = L.run_brownian(temperature, damping, max_time, n_step, n_plot, main_logger,
                              stop_change_conf=True, id_conf=to_check, square_angle_sum=square_angle_sum, random_seed=True)
        
        return time

    max_temp_ok = False
    max_temp = -1

    while not max_temp_ok:
        if nb_try >= max_try:
            print('cannot converge even with stability equal to',
                  nb_needed_for_stability, sep=" : ")
            return -1
        try:
            max_temp = get_max_stable_temp(
                min_temp_totest, max_temp_totest, simulate, nb_needed_for_stability, max_time, target_precision, max_try)
            if max_temp is not None:
                max_temp_ok = True
            else:
                print('did not converge, restarting',
                      nb_needed_for_stability+5, sep=" : ")
                nb_needed_for_stability = nb_needed_for_stability+5
                nb_try += 1
        except ValueError:
            max_temp_test = max_temp_test*2
            print('increasing max temperature', max_temp_test*2, sep=" : ")

    return max_temp
