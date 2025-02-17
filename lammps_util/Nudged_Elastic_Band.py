
import itertools
from collections import Counter
import os
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import lammps_util.protein_template as protein_template
import polygon_util as pu
import numpy as np
import math
import copy
import subprocess


ETOL_MINIMIZE = 1.0e-10

# method that compute the


def create_from_protein_template(protein_template):
    print("print not implemented yet")

# Create the final atom position file needed to do a NEB analysis


def create_final_file(path, atom_list):
    try:
        os.remove(path)
    except:
        print("new file")
    sep = "   "
    with open(path, "a") as f:
        f.write(str(len(atom_list)) + "\n")
        for i in range(len(atom_list)):
            f.write(str(
                i+1) + sep + str(atom_list[i, 1]) + sep + str(atom_list[i, 2]) + sep + "0.0 \n")


def count_binding(binding_list):
    count_1 = Counter([p[0] for p in binding_list])
    count_2 = Counter([p[1] for p in binding_list])
    return max([max(list(count_1.values())), max(list(count_2.values()))])


def avg_dist(list_pos, np_array_pos):

    return sum(
        math.sqrt(
            (list_pos[i][0] - np_array_pos[i, 0]) ** 2
            + (list_pos[i][1] - np_array_pos[i, 1]) ** 2
        )
        for i in range(len(list_pos))
    )


def get_distance(pos1, pos2):
    return math.sqrt((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2)


def add_center_node(pos_h, pos_p, inside_node):
    pos_p_2 = (copy.copy(pos_p))
    pos_p_2.extend(list(zip(inside_node[2], inside_node[3])))
    all_p = np.array(pos_p_2)
    # get position with inside node for healthy
    pos_h_2 = (copy.copy(pos_h))
    pos_h_2.extend(list(zip(inside_node[0], inside_node[1])))
    all_h = np.array(pos_h_2)

    return all_h, all_p


def add_edge(dict_edge, dict_node, key1, key2, node_pos, binded, n_out):

    if key1[0] == n_out:
        key1 = (0, key1[1])
    if key2[0] == n_out:
        key2 = (0, key2[1])

    # if none of them were deleted
    if key1 in dict_node and key2 in dict_node:
        dict_edge[(key1, key2)] = get_distance(
            node_pos[key1[0]], node_pos[key2[0]])
        return (key1, key2)
    # if one of them was deleted because of the binding we still want to add the edge
    else:
        if key1 not in dict_node and key2 in dict_node:
            dict_edge[((binded[key1[0] % n_out], 'p'), key2)] = get_distance(
                node_pos[key1[0] % n_out], node_pos[key2[0]])
            return ((binded[key1[0] % n_out], 'p'), key2)
        if key1 in dict_node and key2 not in dict_node:
            dict_edge[(key1, (binded[key2[0] % n_out], 'p'))] = get_distance(
                node_pos[key1[0]], node_pos[key2[0] % n_out])
            return (key1, (binded[key2[0] % n_out], 'p'))
            


def translate_position(non_translated_pos, translation):

    translated_pos = [0]*len(non_translated_pos)
    for i in range(len(non_translated_pos)):
        if i in translation:
            translated_pos[i] = (
                non_translated_pos[translation[i], 0], non_translated_pos[translation[i], 1])
        else:
            translated_pos[i] = (non_translated_pos[i, 0],
                                 non_translated_pos[i, 1])

    return translated_pos


def transform_all_point_back(binded_from_sol, binded_to_move, solution_poly, full_point_to_transform, inverse):
    # transform the inside point by refitting the subset of bonded point for HP
    trans1_hp, data1_hp, trans2_hp, data2_hp = pu.find_good_initial_guess(
        set_pt1=binded_from_sol, set_pt2=binded_to_move, data_to_shift=full_point_to_transform)
    dist1_hp, dist2_hp = avg_dist(list(solution_poly.exterior.coords), data1_hp), avg_dist(
        list(solution_poly.exterior.coords), data2_hp)
    if not inverse:
        full_point_transformed, tr = data1_hp, trans1_hp
        if dist2_hp < dist1_hp:
            full_point_transformed, tr = data2_hp, trans2_hp
    else:
        full_point_transformed, tr = data2_hp, trans1_hp
        if dist2_hp < dist1_hp:
            full_point_transformed, tr = data1_hp, trans2_hp

    return full_point_transformed, tr


def get_dict_node(full_bottom_node, full_top_node, binded_pairs_list, n_out, translation):
    # set un the dict for the nodes
    dict_node_temp, dict_node = {}, {}
    for i in range(len(full_bottom_node)):
        dict_node_temp[(i, 'h')] = full_bottom_node[i]
    for i in range(len(full_top_node)):
        dict_node_temp[(i, 'p')] = full_top_node[i]
    # remove the loop
    del dict_node_temp[(n_out, 'h')]
    del dict_node_temp[(n_out, 'p')]
    # remvove healty binding
    for p in binded_pairs_list:
        del dict_node_temp[(p[0], 'h')]

    # translate the healthy and put in new dict
    for i in range(len(full_bottom_node)):
        if (i, 'h') in dict_node_temp:
            if i in translation:
                dict_node[(translation[i], 'h')] = dict_node_temp[(i, 'h')]
            else:
                dict_node[(i, 'h')] = dict_node_temp[(i, 'h')]
        if (i, 'p') in dict_node_temp:
            dict_node[(i, 'p')] = dict_node_temp[(i, 'p')]
    return dict_node


def get_dict_edge(full_bottom_node, dict_node, binded_dict, connection,  n_out, translation):

    key_in, key_out = [], []
    # get edge list
    dict_edge = {}
    # translate the healthy of hp to fit pp
    full_healthy_hp_translated = translate_position(
        full_bottom_node, translation)
    # put the outside edge for both
    for i in range(n_out+1):
        # get the position of the two consecutive edge of the outside
        key1, key2 = (i, 'h'), ((i+1) % n_out, 'h')
        key_in.append(add_edge(dict_edge, dict_node, key1, key2,
                 full_healthy_hp_translated, binded_dict, n_out))
        
        # the second protein doesn't have any edge removed therefore we can add them all
        key1p, key2p = (i, 'p'), ((i+1) % n_out, 'p')
        key_in.append(add_edge(dict_edge, dict_node, key1p, key2p,
                 full_healthy_hp_translated, binded_dict,  n_out))
        
    # edge to the inside
    for p in connection:
        f_, t_ = p[0], p[1]
        # for the botom protein some edge are removed so we need to check
        key1, key2 = (f_, 'h'), (t_, 'h')
        key_out.append(add_edge(dict_edge, dict_node, key1, key2,
                 full_healthy_hp_translated, binded_dict, n_out))
        # for the top one no edge are removed
        key1p, key2p = (f_, 'p'), (t_, 'p')
        key_out.append(add_edge(dict_edge, dict_node, key1p, key2p,
                 full_healthy_hp_translated, binded_dict, n_out))

    return dict_edge, (key_in, key_out)


def _get_position_in_protein_dictionarry(dict_hp):
    return {key: i for i, key in enumerate(list(dict_hp.keys()))}


def _initialize_protein_pair(dict_hp, dict_pp, dict_edge_pp, n_out):

    protein_pair = protein_template.Protein_Template()
    dict_hp_key = list(dict_hp.keys())

    healthy_pos_protein_format, prion_pos_protein_format = [], []
    dict_pos_in_protein = {}

    # translate_healthy_position
    for i, key in enumerate(dict_hp_key):
        typ = 1  # inside node
        if key[0] <= n_out:
            typ = 0  # outside node
        pos_hp, pos_pp = dict_hp[key], dict_pp[key]
        dict_pos_in_protein[key] = i
        healthy_pos_protein_format.append([typ, pos_hp[0], pos_hp[1]])
        prion_pos_protein_format.append([typ, pos_pp[0], pos_pp[1]])

    conn_protein_format, conn_length_list = [], []

    # define all the edges/interaction
    for key in dict_edge_pp:
        from_id = dict_pos_in_protein[key[0]]
        to_id = dict_pos_in_protein[key[1]]
        conn_protein_format.append([from_id, to_id])
        conn_length_list.append(dict_edge_pp[key])

    # format in numpy array
    healthy_pos_protein_format = np.array(healthy_pos_protein_format)
    prion_pos_protein_format = np.array(prion_pos_protein_format)
    conn_protein_format = np.array(conn_protein_format)

    # create the protein template
    protein_pair.set_healthy_structure(
        position=healthy_pos_protein_format, connection=conn_protein_format)
    protein_pair.set_prion_structure(position=prion_pos_protein_format)
    protein_pair.set_bond_length(conn_length_list)

    return protein_pair, dict_pos_in_protein


def get_protein_pairs_binded_bypass(pos_h, pos_p, translation, inside_node, conn, hp_vec, pp_vec):
    
    n_out = len(pos_h)-1
    best_sol_hp, binded_hp =  hp_vec
    best_sol_pp, binded_pp =  pp_vec
    # correct connection (fit the indices and not matlab)
    conn = [(c[0]-1, c[1]) for c in conn]

    # get back the binding PH and the binding PP
    healthy_poly, prion_poly = Polygon(pos_h), Polygon(pos_p)

    # get position with inside node for prion / healthy
    all_h, all_p = add_center_node(pos_h, pos_p, inside_node)
    full_healthy_hp, full_prion_pp_p2 = np.copy(all_h), np.copy(all_p)

    # get binded node from the prion that do not move
    healthy_binded, prion_binded = [p[0] for p in binded_hp], [p[1] for p in binded_hp]
    init_prion_3 = [list(prion_poly.exterior.coords)[i] for i in prion_binded]
    # transform the inside point by refitting the subset of bonded point for HP
    sol_hp_3 = [list(best_sol_hp.exterior.coords)[i]for i in prion_binded]
    full_prion_hp, trans_hp = transform_all_point_back( sol_hp_3, init_prion_3, best_sol_hp, all_p, False)
    # transform the inside point by refitting the subset of bonded point for PP
    sol_pp_3 = [list(best_sol_pp.exterior.coords)[i] for i in prion_binded]
    full_prion_pp_p1, trans_pp = transform_all_point_back( sol_pp_3, init_prion_3, best_sol_pp, all_p, False)

    # set un the dict for the nodes
    binded_translated = {translation[p[0]]: p[1] for p in binded_hp}
    dict_hp = get_dict_node(full_healthy_hp, full_prion_hp, binded_hp, n_out, translation)
    dict_edge_hp, key_io_hp = get_dict_edge( full_healthy_hp, dict_hp, binded_translated, conn, n_out, translation)

    # set un the dict for the nodes
    identity_translation = {i: i for i in range(n_out)}
    dict_pp = get_dict_node(full_prion_pp_p2, full_prion_pp_p1, binded_pp, n_out, identity_translation)
    dict_edge_pp, key_io_pp = get_dict_edge(full_prion_pp_p2, dict_pp, binded_translated, conn, n_out, identity_translation)

    #get dict to translate the position in the molecule from internal key ie. (3,p)
    tranlate_pos = _get_position_in_protein_dictionarry(dict_hp)
 
    in_edge = [ (tranlate_pos[p[0]],tranlate_pos[p[1]]) for p in key_io_hp[0] if p is not None]
    out_edge = [(tranlate_pos[p[0]],tranlate_pos[p[1]]) for p in key_io_hp[1] if p is not None]


    assert set(dict_edge_pp.keys()) == set(
        dict_edge_hp.keys()), "edge set should be the same"
    assert set(dict_pp.keys()) == set(
        dict_hp.keys()), "node set should be the same"

    protein_pair, dict_pos_in_protein = _initialize_protein_pair(
        dict_hp, dict_pp, dict_edge_pp, n_out)

    return protein_pair, dict_pos_in_protein, (in_edge,out_edge)

    

def get_protein_pairs_binded(pos_h, pos_p, translation, binded_hp, inside_node, conn, target_previous_pp_fitness, target_previous_hp_fitness, min_r_squared_default, n_point_max):

    n_out = len(pos_h)-1

    # correct connection (fit the indices and not matlab)
    conn = [(c[0]-1, c[1]) for c in conn]

    # get back the binding PH and the binding PP
    healthy_poly, prion_poly = Polygon(pos_h), Polygon(pos_p)


    fit_ordered_hp, sol_ordered_hp, output_index_hp, _ = pu.estimate_stacking_icp(
        healthy_poly, prion_poly, avg_dist_tol=0.3, n_point_max=n_point_max,  multiple_bind=False, output_index=True)
    best_fit_hp, best_sol_hp = fit_ordered_hp[-1], sol_ordered_hp[-1]
    hp_fit_max = pu.rescale_fit(best_fit_hp, min_r_squared_default)


    fit_ordered_pp, sol_ordered_pp, output_index_pp, _ = pu.estimate_stacking_icp(
        prion_poly, prion_poly, avg_dist_tol=0.3, n_point_max=n_point_max,  multiple_bind=False, output_index=True)
    best_fit_pp, best_sol_pp = fit_ordered_pp[-1], sol_ordered_pp[-1]
    pp_fit_max = pu.rescale_fit(best_fit_pp, min_r_squared_default)


    # todo assert if they are kinda the same
    assert (hp_fit_max-target_previous_hp_fitness) < 0.01, "there is a big difference in the fitness" + \
        str(hp_fit_max-target_previous_hp_fitness)
    assert (pp_fit_max-target_previous_pp_fitness) < 0.01, "there is a big difference in the fitness" + \
        str(pp_fit_max-target_previous_pp_fitness)

    # get position with inside node for prion / healthy
    all_h, all_p = add_center_node(pos_h, pos_p, inside_node)
    full_healthy_hp, full_prion_pp_p2 = np.copy(all_h), np.copy(all_p)

    # get binded node from the prion that do not move
    healthy_binded, prion_binded = [p[0] for p in binded_hp], [p[1] for p in binded_hp]
    init_prion_3 = [list(prion_poly.exterior.coords)[i] for i in prion_binded]
    # transform the inside point by refitting the subset of bonded point for HP
    sol_hp_3 = [list(best_sol_hp.exterior.coords)[i]for i in prion_binded]
    full_prion_hp, trans_hp = transform_all_point_back( sol_hp_3, init_prion_3, best_sol_hp, all_p, False)
    # transform the inside point by refitting the subset of bonded point for PP
    sol_pp_3 = [list(best_sol_pp.exterior.coords)[i] for i in prion_binded]
    full_prion_pp_p1, trans_pp = transform_all_point_back( sol_pp_3, init_prion_3, best_sol_pp, all_p, False)

    # set un the dict for the nodes
    binded_translated = {translation[p[0]]: p[1] for p in binded_hp}
    dict_hp = get_dict_node(full_healthy_hp, full_prion_hp, binded_hp, n_out, translation)
    dict_edge_hp, key_io_hp = get_dict_edge( full_healthy_hp, dict_hp, binded_translated, conn, n_out, translation)

    # set un the dict for the nodes
    identity_translation = {i: i for i in range(n_out)}
    binded_pp = [(translation[p[0]], p[1]) for p in binded_hp]
    dict_pp = get_dict_node(full_prion_pp_p2, full_prion_pp_p1, binded_pp, n_out, identity_translation)
    dict_edge_pp, key_io_pp = get_dict_edge(full_prion_pp_p2, dict_pp, binded_translated, conn, n_out, identity_translation)

    #get dict to translate the position in the molecule from internal key ie. (3,p)
    tranlate_pos = _get_position_in_protein_dictionarry(dict_hp)
 
    in_edge = [ (tranlate_pos[p[0]],tranlate_pos[p[1]]) for p in key_io_hp[0] if p is not None]
    out_edge = [(tranlate_pos[p[0]],tranlate_pos[p[1]]) for p in key_io_hp[1] if p is not None]


    assert set(dict_edge_pp.keys()) == set(
        dict_edge_hp.keys()), "edge set should be the same"
    assert set(dict_pp.keys()) == set(
        dict_hp.keys()), "node set should be the same"

    protein_pair, dict_pos_in_protein = _initialize_protein_pair(
        dict_hp, dict_pp, dict_edge_pp, n_out)

    return protein_pair, dict_pos_in_protein, (in_edge,out_edge)

# Create the task  file needed to do a NEB analysis


def create_neb_task_file(path, init_path, final_path, path_dir_dump, timestep, spring_inter_replica=1.0, MaxE=0.0, nbMin=200000, nbClimb=200000, maxF=0.0, prnt=-1, minimize=False):
    try:
        os.remove(path)
    except Exception:
        print("new file")

    sep = "     "
    with open(path, "a") as f:
        f.write("# LAMMPS task file for neb \n \n ")

        f.write(f'{sep}dimension   2' + "\n")
        f.write(f'{sep}atom_style bond' + "\n")
        f.write(f'{sep}boundary   f f p' + "\n")
        f.write(f'{sep}bond_style harmonic' + "\n")
        f.write(sep + "atom_modify	map array sort 0 0.0"+"\n")
        f.write("\n")

        f.write(sep + "variable  u uloop 100\n")
        f.write(f'{sep}read_data  {init_path}' + "\n")

        if minimize:
            f.write(
                f"minimize {str(ETOL_MINIMIZE)}  {str(ETOL_MINIMIZE/100)}  "
                + str(nbMin // 10)
                + "  "
                + str(nbMin)
                + " \n"
            )

        f.write(f'{sep}timestep  {str(timestep)}' + "\n")
        f.write(
            f'{sep}fix 1 all neb {str(spring_inter_replica)} parallel ideal'
            + "\n"
        )

        f.write("\n")
        f.write(f'{sep}thermo 100' + "\n")
        f.write(f'{sep}fix 2 all enforce2d' + "\n")
        f.write(f'{sep}dump 1 all atom 10 {path_dir_dump}dump.neb.$u' + "\n")
        f.write(f'{sep} min_style quickmin' + "\n")

        if MaxE == 0.0 and maxF == 0.0:
            MaxE = 0.00000000000001
        if prnt == -1:
            prnt = int(nbMin/50)

        f.write(
            f'{sep}neb {str(MaxE)} {str(maxF)} {str(nbMin)} {str(nbClimb)} '
            + str(prnt)
            + " final "
            + final_path
            + "\n"
        )


def restructure_output(neb_final_py_path, dump_dir, output_file="dump.restructured.final"):

    str_to_run = (
        f"python {str(neb_final_py_path)} -o {str(output_file)} -r {dump_dir}"
        + "dump.neb.*"
    )

    test = subprocess.Popen([str_to_run], stdout=subprocess.PIPE)
    output = test.communicate()[0]
    print(output)

    return str_to_run


def plot_position(output_file_path, index_list=None, marker=False):

    if index_list is None:
        index_list = []
    list_to_plot = get_atom_position(output_file_path)

    for i in range(len(index_list)):
        xpos = list_to_plot[index_list[i]][0]
        ypos = list_to_plot[index_list[i]][1]
        plt.scatter(xpos, ypos)

        if marker:
            for j in range(len(xpos)):
                plt.text(xpos[j], ypos[j], str(j))

    plt.show()


def get_atom_position(output_file_path):

    file1 = open(output_file_path, 'r')
    Lines = file1.readlines()

    timestep_id = [i for i, line in enumerate(Lines) if 'TIMESTEP' in line]

    atom_per_sim = int(Lines[timestep_id[0]+3])

    pos_string_list_list = [[] for _ in range(len(timestep_id))]
    atom_pos = [[[], []] for _ in range(len(timestep_id))]

    for i, j in itertools.product(range(len(timestep_id)), range(atom_per_sim)):

        string = Lines[timestep_id[i]+j+9]

        splitted = string.split()
        map_object = map(float, splitted)
        list_of_float = list(map_object)

        pos_string_list_list[i].append((list_of_float[2], list_of_float[3]))
        atom_pos[i][0].append(list_of_float[2])
        atom_pos[i][1].append(list_of_float[3])

    return atom_pos


def get_position_healthy_prion(all_position, dict_pos_in_protein, binded_pp, n_node, n_out):
    pos_bottom, pos_top = [], []

    for i in range(n_node+1):
        if i != n_out:  # remove the node that was previously used to close the outside polygon

            if (i, 'h') in dict_pos_in_protein:

                pos = all_position[dict_pos_in_protein[(i, 'h')]]
            else:

                pos = all_position[dict_pos_in_protein[(binded_pp[i], 'p')]]

            if len(all_position[0]) == 3:
                pos_bottom.append([pos[1], pos[2]])
            else:
                pos_bottom.append([pos[0], pos[1]])
            if (i, 'p') in dict_pos_in_protein:
                pos = all_position[dict_pos_in_protein[(i, 'p')]]
                if len(all_position[0]) == 3:
                    pos_top.append([pos[1], pos[2]])
                else:
                    pos_top.append([pos[0], pos[1]])

    return pos_bottom, pos_top


def plot_energy(log_file_path):

    file1 = open(log_file_path, 'r')
    Lines = file1.readlines()

    timestep_id = []
    last_line = Lines[-1]
    splitted = last_line.split()
    map_object = map(float, splitted)
    list_of_float = list(map_object)

    list_of_float = list_of_float[9:]
    n_replica = len(list_of_float) // 2

    rc_list = []
    energy_list = []

    for i in range(len(list_of_float)):
        if i % 2 == 0:
            rc_list.append(list_of_float[i])
        else:
            energy_list.append(list_of_float[i])

    plt.plot(rc_list, energy_list)
    plt.show()


def get_energy(log_file_path):

    file1 = open(log_file_path, 'r')
    Lines = file1.readlines()

    timestep_id = []
    last_line = Lines[-1]
    splitted = last_line.split()
    map_object = map(float, splitted)
    list_of_float = list(map_object)

    list_of_float = list_of_float[9:]
    n_replica = len(list_of_float) // 2

    rc_list = []
    energy_list = []

    for i in range(len(list_of_float)):
        if i % 2 == 0:
            rc_list.append(list_of_float[i])
        else:
            energy_list.append(list_of_float[i])

    return rc_list, energy_list





# Create the task  file needed to do a NEB analysis
