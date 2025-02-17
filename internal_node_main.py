import math
import numpy as np
import random
from shapely.geometry import Polygon, Point
import database
import copy
import matlab.engine
import itertools


def start_matlab_engine():
    eng = matlab.engine.start_matlab()
    s = eng.genpath('matlab/')
    eng.addpath(s, nargout=0)
    return eng


def get_random_conn(n_out, nedge, npoint):
    """
    Generate a random connection between the internal and external nodes.

    Args:
        n_out (int): Number of nodes (External).
        nedge (int): Total number of edges.
        npoint (int): Number of internal nodes

    Returns:
        list: A list of connections where each sublist contains nodes that are connected.
    """
    found_connected_to = False # Flag to check if all output nodes are connected

    while not found_connected_to:

        node_list = list(range(n_out))
        found = False # Flag to check if a valid partition is found

        while not found:
            okay = True
            P = nedge
            I = npoint

            # Generate random split points to divide edges into groups
            mylist = range(1, P + 1)
            split_points = np.random.choice(P - 2, I - 1, replace=False) + 1
            split_points.sort()
            result = np.split(mylist, split_points)
            for block in result:
                if len(block) <= 1:
                    okay = False
            if okay:
                found = True

        connected_to = [random.sample(node_list, len(result[i]))
                        for i in range(len(result))]
        
        # Ensure no group is empty or has only one edge
        node_used = set()
        for block in connected_to:
            for b in block:
                node_used.add(b)

        if node_used == set(node_list):
            found_connected_to = True

    return connected_to


def compare_connection(connection1, connection2):
    #ex: [[4, 6], [3, 6], [2, 6], [3, 7], [1, 7], [5, 7]]

    list_inside_node1, list_inside_node2 = list(
        {p[1] for p in connection1}), list({p[1] for p in connection2})
    if len(list_inside_node1) != len(list_inside_node2):
        return False, [], []

    dict_conn1, dict_conn2 = {p: set() for p in list_inside_node1}, {
        p: set() for p in list_inside_node2}
    for pair in connection1:
        dict_conn1[pair[1]].add(pair[0])
    for pair in connection2:
        dict_conn2[pair[1]].add(pair[0])
    # this is slow to check all permutation but it is easy to code and the list will contains max 5 element (120 ways)
    for perm in itertools.permutations([i for i in range(len(list_inside_node1))]):

        found = True
        for i, in_node_1 in enumerate(list_inside_node1):
            in_perm_node_2 = list_inside_node2[perm[i]]
            if dict_conn1[in_node_1] != dict_conn2[in_perm_node_2]:
                found = False
                break

        if found:
            translation = {node1: list_inside_node2[perm[i]]
                           for i, node1 in enumerate(list_inside_node1)}
            return True, translation, perm

    return False, [], []


def compare_distance(inside_node_1, inside_node_2, perm):
    #list_inside_node1, list_inside_node2 = list({p[1] for p in connection1}), list({p[1] for p in connection2})
    in_healthy_x_1, in_healthy_y_1 = inside_node_1[0], inside_node_1[1]
    in_prion_x_1, in_prion_y_1 = inside_node_1[2], inside_node_1[3]
    in_healthy_x_2, in_healthy_y_2 = [inside_node_2[0][i]
                                      for i in perm], [inside_node_2[1][i] for i in perm]
    in_prion_x_2, in_prion_y_2 = [inside_node_2[2][i]
                                  for i in perm], [inside_node_2[3][i] for i in perm]

    dist_list = [math.sqrt((x1-x2)**2 + (y1-y2)**2) for x1, y1, x2, y2 in zip(
        in_healthy_x_1, in_healthy_y_1, in_healthy_x_2, in_healthy_y_2)]
    dist_list.extend([math.sqrt((x1-x2)**2 + (y1-y2)**2) for x1, y1, x2,
                     y2 in zip(in_prion_x_1, in_prion_y_1, in_prion_x_2, in_prion_y_2)])
    avg_dist = sum(dist_list)/len(dist_list)

    return avg_dist


def check_if_solution_already_found(sol_list, inside_node, connection, avg_dist_node_threshold):
    # try to find a solution that match
    #sol_all_dict = self.db_pairs[key_healthy][key_prion]['sols']
    found = False
    for sol_inside_node, sol_conn in sol_list:
        same_conn, _, perm = compare_connection(connection, sol_conn)
        if same_conn:
            avg_dist = compare_distance(inside_node, sol_inside_node, perm)
            if avg_dist < avg_dist_node_threshold:
                found = True
                break
    return found


def find_solution(db: database.Database, eng, list_key_dict= None, nb_solution_wanted = None, max_nb_try = None, min_avg_distance_sol = None, nb_node_in = None, conn_func = None):
    lock = False
    
    #to specify the key to check list_key_dict {h_key:[p_key,...],...}
    #shuffle the key it helps with avoiding two worker working on the same element without solution one after the other
    list_key_healthy = list(db.db_pairs.keys())
    random.shuffle(list_key_healthy)

    if list_key_dict is not None: #if ovveride the keys
        list_key_healthy = list_key_dict.keys()

    for i_h, key_healthy in enumerate(list_key_healthy):

        # if healthy not locked by someone we can work on it
        if not db.is_locked(int(key_healthy)) or not lock:
            if lock:
                db.lock(int(key_healthy))  # lock the healhy

            n_poly = db.db_parameter["n_node"]

            list_prion_key = copy.copy(list(db.db_pairs[key_healthy].keys()))
            list_prion_key.remove('struct_healthy')
            

            if list_key_dict is not None: #if overide the keys
                list_prion_key = list_key_dict[key_healthy]
            if max_nb_try is  None:
                max_nb_try = db.db_parameter['max_nb_try']
            if nb_solution_wanted is None: #if not override
                nb_solution_wanted = db.db_parameter['nb_solution_wanted']
            if min_avg_distance_sol is None:
                min_avg_distance_sol = db.db_parameter['min_avg_distance_sol']
            if nb_node_in is None:
                nb_node_in = db.db_parameter['nb_node_in']

            print(f"nb_node_in:  {nb_node_in}")
            print(f"min_avg_distance_sol:  {min_avg_distance_sol}")
            print(f"nb_solution_wanted:  {nb_solution_wanted}")
            print(f"max_nb_try:  {max_nb_try}")

            for i_p, key_prion in enumerate(list_prion_key):

                healthy_pos = db.db_pairs[key_healthy]['struct_healthy']
                poly_healthy = Polygon(healthy_pos)
                pos_healty_no_loop = healthy_pos[:-1]

                # try it if there is less than nb_solution_wanted
                # if not store nb_solution_wanted-already found
                already_found_sol = []
                if 'sols' in db.db_pairs[key_healthy][key_prion]:
                    list_sol_key,  sol_value = db.get_all_sol(key_healthy, key_prion)
                    for key, items in zip(list_sol_key, sol_value):
                        already_found_sol.append(
                            (items["inside_node"], items["conn"]))

                nb_already_found = len(already_found_sol)
                nb_already_try=0
                if 'sols' in db.db_pairs[key_healthy][key_prion] and 'nb_try' in db.db_pairs[key_healthy][key_prion]['sols']:
                    nb_already_try = db.db_pairs[key_healthy][key_prion]['sols']['nb_try']
                new_sol_list = []
                nb_try = -1

                if nb_already_found < nb_solution_wanted:

                    # number opf edge to have 0 degree of freedom
                    nedge_in = 2*(n_poly+nb_node_in) - 3 - n_poly
                    #print(nedge_in,n_poly)
                    #assert nedge_in >= n_poly, "not enough node"

                    # load value
                    prion_pos = db.db_pairs[key_healthy][key_prion]['struct']
                    print(prion_pos)
                    poly_prion = Polygon(prion_pos)

                    pos_prion_no_loop = prion_pos[:-1]
                    translation_test = db.db_pairs[key_healthy][key_prion]['translation']

                    # translate the healthy one
                    pos_healthy_test_trans = [0]*(len(pos_healty_no_loop))
                    for hp in translation_test.keys():
                        pos_healthy_test_trans[hp] = pos_healty_no_loop[translation_test[hp]]

                    Xs0 = [[p[0] for p in pos_healthy_test_trans], [p[1]
                                                                    for p in pos_healthy_test_trans]]
                    XsF = [[p[0] for p in pos_prion_no_loop], [p[1]
                                                               for p in pos_prion_no_loop]]

                    # Find position of center node
                    nb_try, nb_found = 0, nb_already_found
                    new_sol_list = []
                    
                    conn = None

                    # we stop either when enough solution is found
                    while nb_try < max_nb_try and nb_found < nb_solution_wanted:
                        

                        if conn_func is None:
                            conn_unformated = get_random_conn(
                                n_poly, nedge_in, nb_node_in)
                            conn = []

                            for in_node_id, connected_to in enumerate(conn_unformated):
                                real_id = n_poly + in_node_id + 1
                                for node in connected_to:
                                    conn.append([node+1, real_id])
                        elif conn is None:
                            conn = conn_func((key_healthy, key_prion))
                    

                        okay = True
                        x0 = [[random.random() for i in range(nb_node_in)], [
                            random.random() for i in range(nb_node_in)]]
                        
                
                        #transform in matlab
                        Xs0_m = matlab.double(Xs0)
                        XsF_m = matlab.double(XsF)
                        x0_m = matlab.double(x0)
                        conn_m = matlab.double(conn)

                        Xu_m, fu_m = eng.construct_network_c(
                            Xs0_m, XsF_m, x0_m, conn_m, 1, nargout=2)

                        fu = np.array(fu_m._data)
                        Xu = np.array(Xu_m._data).reshape(Xu_m.size[::-1]).T
                        avg_dist = np.average(fu)

                        if avg_dist > 0.001:
                            okay = False
                            print("did not converge")

                        Xu_init_x, Xu_init_y = Xu[0], Xu[1]
                        Xu_final_x, Xu_final_y = Xu[2], Xu[3]

                        # check if center node all inside the polygon
                        for x, y in zip(Xu_init_x, Xu_init_y):
                            if not poly_healthy.contains(Point((x, y))):
                                okay = False
                        for x, y in zip(Xu_final_x, Xu_final_y):
                            if not poly_prion.contains(Point((x, y))):
                                okay = False
                        if nb_try%100==0:
                            print(nb_try)
                        nb_try += 1
                        if okay:
                            print("sol found: " + str(nb_try))
                            # put all know solution in a list
                            all_sol = copy.copy(already_found_sol)
                            all_sol.extend(copy.copy(new_sol_list))
                            # check if the solution is different enough
                            if not check_if_solution_already_found(all_sol, Xu, conn, min_avg_distance_sol):
                                nb_found += 1
                                new_sol_list.append((Xu, conn))
                                print(x0)
                                print("sol unique")
                            else:
                                print("sol already exist")
                else:
                    print("enough solution already found ")

                if len(new_sol_list) > 0:  # add new solution if we found some
                    db.add_solution(key_healthy, key_prion, new_sol_list,nb_try )
                    print("sol added in db for " + str(i_h)+","+str(i_p))
                else:
                    db.add_solution(key_healthy, key_prion, [],nb_try )
                    print("no solution to add")
                print(len(new_sol_list), i_h, i_p, key_healthy, key_prion)
            if lock:
                db.unlock(int(key_healthy))  # unlock the healthy


if __name__ == "__main__":

    #keys to find internal nodes. 
    #None will do them all
    list_key_dict={1737806584325749351: [ -8771402139612191995]}

    nb_solution_wanted  = 1000 #number of target solution for each pair prion-healthy of external nodes
    path = "/Users/mathieuouellet/Desktop/n=5_5deg_2_fig4b" #path of the database
    max_nb_try = 10000 #number of try to find a solution 
    min_avg_distance_sol = 0.01 #accepted error 
    db = database.Database(path)
    nb_node_in = 2 #number of internal nodes

    #function to get the connection pattern 
    def get_conn(pair_2_key):
        trans_pair_1 = db.db_pairs[pair_2_key[0]][pair_2_key[1]]['translation']
        conn_p_nat_pair_1 = db.db_pairs[pair_2_key[0]][pair_2_key[1]]['sols'][pair_2_key[2]]['conn']
        conn_h_nat_pair_1 = [[[i for i in trans_pair_1 if trans_pair_1[i]==n1-1][0]+1,n2] for n1,n2 in conn_p_nat_pair_1]
        return conn_h_nat_pair_1



    eng = start_matlab_engine() #matlab engine
    s = eng.genpath('AMP/src/polygon/matlab/') #path to matlab script for the posititioning of internal nodes
    eng.addpath(s, nargout=0)
    find_solution(db, eng, 
        list_key_dict=list_key_dict, 
        nb_solution_wanted=nb_solution_wanted,
        max_nb_try = max_nb_try,
        min_avg_distance_sol = min_avg_distance_sol,
        nb_node_in=nb_node_in,
        conn_func = get_conn
        )
