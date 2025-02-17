
import sqlitedict as sql
import lammps_util.protein_template as protein_template
import numpy as np


class Database():

    min_r_squared_default = 0.001
    relaxed_min_r = 0.01
    param = ["n_node", "length_list", "min_angle", "is_simple", "max_fitness_healthy",
             "min_fitness_prion", "distance_binded", "hp_pp_dif", "n_point_max"]

    def __init__(self, filename):
        self.db_parameter = sql.SqliteDict(
            f'{str(filename)}.sqlite', tablename="parameter"
        )

        self.db_polygon = sql.SqliteDict(
            f'{str(filename)}.sqlite', tablename="polygon"
        )

        self.db_pairs = sql.SqliteDict(
            f'{str(filename)}.sqlite', tablename="pairs")
        self.db_lock = sql.SqliteDict(
            f'{str(filename)}.sqlite', tablename="lock")

    ##########################
    #         Polygon        #
    ##########################

    def set_up_polygon(self, n_node, length_list,  min_angle, is_simple):
        self.db_parameter["length_list"] = length_list
        self.db_parameter["n_node"] = n_node
        self.db_parameter["min_angle"] = min_angle
        self.db_parameter["is_simple"] = is_simple
        self.db_parameter.commit()

    def add_polygon(self, node_pos, dict_data = None, commit=True):

        if dict_data is None:
            dict_data = {}
        assert node_pos[0] == node_pos[-1], "should be in polygon format and loop back"
        #assert not key in self.db.keys(), "this polygon is already in the database"
        dict_data["struct"] = node_pos

        self.db_polygon[hash(tuple(node_pos))] = dict_data
        if commit:
            self.db_polygon.commit()

    def get_list_pos(self, hash):
        return self.db_polygon[hash]['struct']

    def get_list_pos_fit(self, hash):
        dict_data = self.db_polygon[hash]
        return dict_data['struct'], dict_data['fitness']

    def get_hash_prion_healthy(self):
        # COND1: healthy bind less then max_healthy
        # COND2: prion bind more then min_prion

        # TODO kinda dumb should just create a conditional iterator
        assert "max_fitness_healthy" in self.db_parameter, "max_fitness_healthy not set"
        assert "min_fitness_prion" in self.db_parameter, "min_fitness_prion not set"
        max_fitness_healthy, min_fitness_prion = self.db_parameter[
            "max_fitness_healthy"], self.db_parameter["min_fitness_prion"]
        healthy_hash, prion_hash = [], []

        for key, poly_dic in self.db_polygon.items():
            fit = poly_dic['fitness']
            # COND1
            if fit < max_fitness_healthy:
                healthy_hash.append(key)
            # COND2
            if fit > min_fitness_prion:
                prion_hash.append(key)

        return healthy_hash, prion_hash

    def add_polygon_comparions_list(self, hash_healthy, list_hash_prion: list[int], commit=True):

        assert all(
            type(h) == int for h in list_hash_prion
        ), "list_hash_prion should be a list of int"

        dict_data_actual = self.db_polygon[hash_healthy]
        # if the list doesn't exsist
        if 'compared_to' not in dict_data_actual:
            dict_data_actual['compared_to'] = set(list_hash_prion)
        else:
            already_checked: set = dict_data_actual['compared_to']
            for prion_hash in list_hash_prion:
                if prion_hash in already_checked:
                    print("prion was already in the checked list")
                already_checked.add(prion_hash)
            dict_data_actual['compared_to'] = already_checked
        self.db_polygon[hash_healthy] = dict_data_actual
        if commit:
            self.db_polygon.commit()

    def commit_polygon(self):
        list_keys = list(self.db_parameter.keys())
        assert "n_node" in list_keys, "n_node not in parameters"
        assert "min_angle" in list_keys, "min_angle not in parameters"
        assert "is_simple" in list_keys, "is_simple not in parameters"
        assert "length_list" in list_keys, "length_list not in parameters"

        self.db_polygon.commit()

    def get_polygon(self, node_pos):
        return self.db_polygon[hash(tuple(node_pos))]

    ##########################
    #          Pairs         #
    ##########################

    def set_up_pairs(self, max_fitness_healthy, min_fitness_prion, distance_binded, hp_pp_dif, n_point_max):
        self.db_parameter["max_fitness_healthy"] = max_fitness_healthy
        self.db_parameter["min_fitness_prion"] = min_fitness_prion
        self.db_parameter["distance_binded"] = distance_binded
        self.db_parameter["hp_pp_dif"] = hp_pp_dif
        self.db_parameter["n_point_max"] = n_point_max
        self.db_parameter.commit()

    def check_if_pair_already_tested(self, hash_healthy, hash_prion):
        # if nothing was compared to it then return false
        if 'compared_to' not in self.db_polygon[hash_healthy]:
            return False
        # if there was some prion comparison done then check
        return int(hash_prion) in self.db_polygon[hash_healthy]['compared_to']

    def add_pair(self, healthy_pos, list_prion_pos, list_data_prion, commit=True):

        assert healthy_pos[0] == healthy_pos[-1], "should be in polygon format and loop back"
        # if not already there then we can just add the new pair
        if hash(tuple(healthy_pos)) not in self.db_pairs:
            dic_all_prions = {"struct_healthy": healthy_pos}
            for prion, data in zip(list_prion_pos, list_data_prion):
                data['struct'] = prion
                dic_all_prions[hash(tuple(prion))] = data
        else:
            dic_all_prions = self.db_pairs[hash(tuple(healthy_pos))]
            for prion, data in zip(list_prion_pos, list_data_prion):
                if hash(tuple(prion)) in dic_all_prions:
                    print("prion was already in the pair, it will be overwritten")
                data['struct'] = prion
                dic_all_prions[hash(tuple(prion))] = data
        self.db_pairs[hash(tuple(healthy_pos))] = dic_all_prions
        if commit:
            self.db_pairs.commit()

    def add_data_to_pairs(self, healthy_pos, prion_pos, data,  commit=True):
        assert healthy_pos[0] == healthy_pos[-1], "should be in polygon format and loop back"
        assert prion_pos[0] == prion_pos[-1], "should be in polygon format and loop back"

        dict_data_pair_healthy = self.db_pairs[hash(tuple(healthy_pos))]
        dict_prion = dict_data_pair_healthy[hash(tuple(prion_pos))]

        for key, dat in data.items():
            dict_prion[key] = dat

        dict_data_pair_healthy[hash(tuple(prion_pos))] = dict_prion
        self.db_pairs[hash(tuple(healthy_pos))] = dict_data_pair_healthy

        if commit:
            self.db_pairs.commit()

    def get_hash_paired_healthy(self):
        return list(self.db_pairs.keys())

    def get_healthy_from_pair(self, key_healthy):
        return self.db_pairs[key_healthy]['struct_healthy']

    def get_prion_list_from_pair(self, key_healthy):
        list_prion_key = []
        list_prion_pos = []
        for key in self.db_pairs[key_healthy].keys():
            if key != 'struct_healthy':
                list_prion_key.append(key)
                list_prion_pos.append(
                    self.db_pairs[key_healthy][key]['struct'])

        return list_prion_key, list_prion_pos

    ##########################
    #    Pairs  solution     #
    ##########################

    def set_up_internal_node(self, nb_node_in, max_nb_try, nb_solution_wanted, min_avg_distance_sol):

        self.db_parameter['max_nb_try'] = max_nb_try
        self.db_parameter['nb_node_in'] = nb_node_in
        self.db_parameter['nb_solution_wanted'] = nb_solution_wanted
        self.db_parameter['min_avg_distance_sol'] = min_avg_distance_sol
        self.db_parameter.commit()

    # sol list [(Xu, conn)]
    def add_solution(self, key_healthy, key_prion, sol_list, nb_tested):

        # if solution are already in the dict
        if 'sols' in self.db_pairs[key_healthy][key_prion]:
            dict_sol = self.db_pairs[key_healthy][key_prion]['sols']
        else:
            dict_sol = {}
        for inside_node, connection in sol_list:
            hash_key = hash(str(inside_node))
            dict_sol[hash_key] = {
                'inside_node': inside_node, 'conn': connection}

        if 'nb_try' in dict_sol:
            dict_sol['nb_try'] = dict_sol['nb_try']+nb_tested
        else:
            dict_sol['nb_try'] = nb_tested

        # just the last update would probably be needed
        pair_dict = self.db_pairs[key_healthy]
        pair_data = pair_dict[key_prion]
        pair_data['sols'] = dict_sol

        # this one should be the only one needed to update the dict
        self.db_pairs[key_healthy] = pair_dict
        self.db_pairs.commit()

    def get_all_sol(self, key_healthy, key_prion):
        list_key = []
        if "sols" in self.db_pairs[key_healthy][key_prion]:
            list_sol_key = list(
                self.db_pairs[key_healthy][key_prion]['sols'].keys())

            if 'nb_try' in list_sol_key:
                list_sol_key.remove('nb_try')
            if list_sol_key:
                list_key.extend(
                    key
                    for key in list_sol_key
                    if len(
                        self.db_pairs[key_healthy][key_prion]["sols"][key][
                            'inside_node'
                        ]
                    )
                    > 0
                )

        return list_key, [self.db_pairs[key_healthy][key_prion]["sols"][key] for key in list_key]

    def load_protein_pair(self, key_healthy, key_prion, sol_key, translate = True):

        pair_data = self.db_pairs[key_healthy][key_prion]
        prion_pos = pair_data['struct']
        translation = pair_data['translation']
        healthy_pos = self.db_pairs[key_healthy]['struct_healthy']
        inside_node = pair_data['sols'][sol_key]['inside_node']
        conn = pair_data['sols'][sol_key]['conn']
        protein_pair = protein_template.Protein_Template()

        # translate_healthy_position
        if translate:
            pos_healthy_trans = [0]*(len(healthy_pos)-1)
            for hp in translation.keys():
                # useless [:-1] I know but I like it
                pos_healthy_trans[hp] = healthy_pos[:-1][translation[hp]]
        else:
            pos_healthy_trans = healthy_pos[:-1]



        healthy_pos_protein_format = []
        prion_pos_protein_format = []
        # add outside node
        for i in range(len(pos_healthy_trans)):
            healthy_pos_protein_format.append(
                [0, pos_healthy_trans[i][0], pos_healthy_trans[i][1]])
            prion_pos_protein_format.append(
                [0, prion_pos[i][0], prion_pos[i][1]])
  
        in_healthy_x, in_healthy_y = inside_node[0], inside_node[1]
        in_prion_x, in_prion_y = inside_node[2], inside_node[3]
        # add inside node
        for i in range(len(in_healthy_x)):
            healthy_pos_protein_format.append(
                [1, in_healthy_x[i], in_healthy_y[i]])
            prion_pos_protein_format.append([1, in_prion_x[i], in_prion_y[i]])
        conn_protein_format = [[pair[0]-1, pair[1]-1] for pair in conn]
        conn_protein_format.extend(
            [i, (i + 1) % (len(healthy_pos) - 1)]
            for i in range(len(healthy_pos) - 1)
        )

        # format in numpy array
        healthy_pos_protein_format = np.array(healthy_pos_protein_format)
        prion_pos_protein_format = np.array(prion_pos_protein_format)
        conn_protein_format = np.array(conn_protein_format)

        # create the protein template
        protein_pair.set_healthy_structure(
            position=healthy_pos_protein_format, connection=conn_protein_format)
        protein_pair.set_prion_structure(position=prion_pos_protein_format)

        return protein_pair

    ##########################
    #     Pairs stability    #
    ##########################

    def set_up_stability(self, n_test_temperature, start_temp, max_temp, maxtime):

        self.db_parameter['n_test_temperature'] = n_test_temperature
        self.db_parameter['start_temp'] = start_temp
        self.db_parameter['max_temp'] = max_temp
        self.db_parameter['maxtime'] = maxtime
        self.db_parameter.commit()

    def set_up_interaction(self, mass, K):
        self.db_parameter['mass'] = mass
        self.db_parameter['K'] = K
        self.db_parameter.commit()

    def add_stability(self, key_healthy, key_prion, sol_key, mbr_triplet_healthy, mbr_triplet_prion):
        # no need to open it slowly like that
        pair_dict = self.db_pairs[key_healthy]
        pair_prion = pair_dict[key_prion]
        pair_solution = pair_prion['sols']
        solution = pair_solution[sol_key]
        # store the triplet
        solution['stability_prion'] = mbr_triplet_prion
        solution['stability_healthy'] = mbr_triplet_healthy
        # waste of time only the last one is probably needed
        pair_solution[sol_key] = solution
        pair_prion['sols'] = pair_solution
        pair_dict[key_prion] = pair_prion
        self.db_pairs[key_healthy] = pair_dict
        # commit
        self.db_pairs.commit()

    def add_RMS_stability(self, key_healthy, key_prion, sol_key, temperature_list,  healthy_angle_list, prion_angle_list, RMS_hp):
        # no need to open it slowly like that
        pair_dict = self.db_pairs[key_healthy]
        pair_prion = pair_dict[key_prion]
        pair_solution = pair_prion['sols']
        solution = pair_solution[sol_key]
        # store the triplet
        solution['RMS_stability_temp'] = temperature_list
        solution['RMS_stability_healthy'] = healthy_angle_list
        solution['RMS_stability_prion'] = prion_angle_list
        solution['RMS_hp'] = RMS_hp
        # waste of time only the last one is probably needed
        pair_solution[sol_key] = solution
        pair_prion['sols'] = pair_solution
        pair_dict[key_prion] = pair_prion
        self.db_pairs[key_healthy] = pair_dict
        # commit
        self.db_pairs.commit()

    ##########################
    #     NEB computation    #
    ##########################

    def _add_NEB(self, key_healthy, key_prion, sol_key, reaction_coordinate, energy_list, coord_key, energy_key, atom_list_list = None):
        # no need to open it slowly like that
        pair_dict = self.db_pairs[key_healthy]
        pair_prion = pair_dict[key_prion]
        pair_solution = pair_prion['sols']
        solution = pair_solution[sol_key]
        # store the triplet
        solution[coord_key] = reaction_coordinate
        solution[energy_key] = energy_list
        # waste of time only the last one is probably needed
        pair_solution[sol_key] = solution
        pair_prion['sols'] = pair_solution
        pair_dict[key_prion] = pair_prion
        self.db_pairs[key_healthy] = pair_dict
        # commit
        self.db_pairs.commit()

    def add_NEB_HP(self, key_healthy, key_prion, sol_key, reaction_coordinate, energy_list, atom_list_list = None):
        self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, energy_list, 'NEB_HP_reaction_coordinate', 'NEB_HP_energy')

        if not atom_list_list is None:
            self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, atom_list_list, 'NEB_HP_reaction_coordinate', 'NEB_HP_atom_list_list')

        
    def delete_NEB_HP(self, key_healthy, key_prion, sol_key):
        # no need to open it slowly like that
        pair_dict = self.db_pairs[key_healthy]
        pair_prion = pair_dict[key_prion]
        pair_solution = pair_prion['sols']
        solution = pair_solution[sol_key]
        # store the triplet
        del solution['NEB_HP_reaction_coordinate']
        del solution['NEB_HP_energy']
        # waste of time only the last one is probably needed
        pair_solution[sol_key] = solution
        pair_prion['sols'] = pair_solution
        pair_dict[key_prion] = pair_prion
        self.db_pairs[key_healthy] = pair_dict
        # commit
        self.db_pairs.commit()

    def add_NEB_HP_PP(self, key_healthy, key_prion, sol_key, reaction_coordinate, energy_list, pos_list, angle_list, atom_list_list = None ):
        self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, energy_list, 'NEB_HP_PP_reaction_coordinate', 'NEB_HP_PP_energy')
        self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, pos_list, 'NEB_HP_PP_reaction_coordinate', 'NEB_HP_PP_pos_list')
        self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, angle_list, 'NEB_HP_PP_reaction_coordinate', 'NEB_HP_PP_angle_list')

        if not atom_list_list is None:
            self._add_NEB(key_healthy, key_prion, sol_key, reaction_coordinate, atom_list_list, 'NEB_HP_PP_reaction_coordinate', 'NEB_HP_PP_atom_list_list')


    ##########################
    #          lock          #
    ##########################

    def set_up_lock(self):
        self.db_lock["locked"] = set()
        self.db_lock.commit()

    def clear_all_lock(self):
        print(str(len(self.db_lock["locked"])) + " locks cleared")
        self.db_lock["locked"] = set()
        self.db_lock.commit()

    def lock(self, hash):
        assert type(hash) == int, "please pass a int"
        if not self.is_locked(hash):
            active_lock: set = self.db_lock["locked"]
            active_lock.add(hash)
            self.db_lock["locked"] = active_lock
            self.db_lock.commit()
            return True
        return False

    def unlock(self, hash):
        assert type(hash) == int, "please pass a int"
        if self.is_locked(hash):
            active_lock: set = self.db_lock["locked"]
            active_lock.remove(hash)
            self.db_lock["locked"] = active_lock
            self.db_lock.commit()
            print(f"hash: {str(hash)} lock cleared")
        return True

    def is_locked(self, hash):
        assert type(hash) == int, "please pass a int"
        print(list(self.db_lock.keys()))
        if hash in self.db_lock["locked"]:
            return True
        return False

    ##########################
    #          Utils         #
    ##########################

    def close(self):
        self.db_parameter.close()
        self.db_polygon.close()
        self.db_pairs.close()
