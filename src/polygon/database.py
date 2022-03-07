import os
import json
import sqlitedict  as sql


class Database():

    
    min_r_squared_default = 0.001
    relaxed_min_r =  0.01
    param = ["n_node","length_list","min_angle","is_simple", "max_fitness_healthy","min_fitness_prion","distance_binded","hp_pp_dif","n_point_max"]
    


    def __init__(self, filename):
        self.db_parameter= sql.SqliteDict(str(filename)+".sqlite",tablename="parameter")
        self.db_polygon = sql.SqliteDict(str(filename)+".sqlite",tablename="polygon")
        self.db_pairs = sql.SqliteDict(str(filename)+".sqlite",tablename="pairs")


    ##########################
    #         Polygon        #
    ##########################

    def set_up_polygon(self, n_node, length_list,  min_angle, is_simple):
        self.db_parameter["is_polygon_locked"] = False
        self.db_parameter["length_list"] = length_list
        self.db_parameter["n_node"] = n_node
        self.db_parameter["min_angle"] = min_angle
        self.db_parameter["is_simple"] = is_simple
        self.db_parameter.commit()


    def lock_polygon(self):
        self.db_parameter["is_polygon_locked"] = True

    def add_polygon(self, node_pos, dict_data={}, commit = True):
        
        assert node_pos[0]==node_pos[-1], "should be in polygon format and loop back"
        #assert not key in self.db.keys(), "this polygon is already in the database"
        dict_data["struct"] = node_pos

        self.db_polygon[hash(tuple(node_pos))]=dict_data
        if commit:
            self.db_polygon.commit()

    def commit_polygon(self):
        list_keys = list(self.db_parameter.keys())
        assert "is_polygon_locked" in list_keys, "is_polygon_locked not in parameters"
        assert "n_node" in list_keys, "n_node not in parameters"
        assert "min_angle" in list_keys, "min_angle not in parameters"
        assert "is_simple" in list_keys, "is_simple not in parameters"
        assert "length_list" in list_keys, "length_list not in parameters"
        assert not self.db_parameter["is_polygon_locked"] in list_keys, "The polygon list is locked"
        
        self.db_polygon.commit()
        
        
    def get_polygon(self, node_pos):
        return self.db_polygon[hash(tuple(node_pos))]

    ##########################
    #          Pairs         #
    ##########################


    def set_up_pairs(self, max_fitness_healthy, min_fitness_prion, distance_binded, hp_pp_dif, n_point_max ):
        self.db_parameter["max_fitness_healthy"] = max_fitness_healthy
        self.db_parameter["min_fitness_prion"] = min_fitness_prion
        self.db_parameter["distance_binded"] = distance_binded
        self.db_parameter["hp_pp_dif"] = hp_pp_dif
        self.db_parameter["n_point_max"] = n_point_max
        self.db_parameter.commit()



    def add_pair(self, healthy_pos, list_prion_pos, list_data_prion, commit = True):

        assert healthy_pos[0]==healthy_pos[-1], "should be in polygon format and loop back"

        dic_all_prions = {}
        dic_all_prions["struct_healthy"] = healthy_pos
        for prion, data in zip(list_prion_pos, list_data_prion):
            data['struct'] = prion
            dic_all_prions[hash(tuple(prion))] = data
        self.db_pairs[hash(tuple(healthy_pos))] = dic_all_prions

        if commit:
                self.db_pairs.commit()

    def add_data_to_pairs(self, healthy_pos, prion_pos, data,  commit = True):
        assert healthy_pos[0]==healthy_pos[-1], "should be in polygon format and loop back"
        assert prion_pos[0]==prion_pos[-1], "should be in polygon format and loop back"

        dict_data_pair_healthy  = self.db_pairs[hash(tuple(healthy_pos))]   
        dict_prion = dict_data_pair_healthy[hash(tuple(prion_pos))]
        
        for key, dat in data.items():
            dict_prion[key]=dat
            
        dict_data_pair_healthy[hash(tuple(prion_pos))] = dict_prion
        self.db_pairs[hash(tuple(healthy_pos))] = dict_data_pair_healthy
        
        if commit:
            self.db_pairs.commit() 


    ##########################
    #          Utils         #
    ##########################


    def close(self):
        self.db_parameter.close()
        self.db_polygon.close()
        self.db_pairs.close()



  