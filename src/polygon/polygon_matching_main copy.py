import polygon_util as pu
from shapely.geometry import Polygon
import database  


#COND1: healthy bind less then max_healthy 
#COND2: prion bind more then min_prion
#COND3: check if the fitness of hp is smaller of at least hp_pp_dif than for pp
#COND4 check if the rounded fitness of hp is bigger or equal than the rounded fitness of  hh
#COND5: check if h and p bind to the same in p
#COND6: check if we can morph the healty into the prion while binded 
def select_pairs(db : database.Database):

    length_list =  db.db_parameter["length_list"]
    distance_binded = db.db_parameter["distance_binded"]
    hp_pp_dif = db.db_parameter["hp_pp_dif"]
    n_point_max = db.db_parameter["n_point_max"]
    min_r_squared_default = db.min_r_squared_default
    relaxed_min_r = db.relaxed_min_r

    healthy_hash, prion_hash = db.get_hash_prion_healthy()

    print("potential prion count:" + str(len(prion_hash)))
    print("potential healthy count:" + str(len(healthy_hash)))
    
    list_pp_fit_max, list_pp_sol_max = [], []
    
    for key_p in prion_hash:
        if False:
            prion_pos_list = db.get_list_pos(key_p)
            prion_poly = Polygon(prion_pos_list)

            #get binding of prion prion and more importantly the binded pairs 
            fit_ordered_pp, sol_ordered_pp, output_index_pp, _= pu.estimate_stacking_icp(prion_poly, prion_poly, n_point_max = n_point_max,  multiple_bind = False , output_index = True )
            pp_fit_max, pp_sol_max = pu.rescale_fit(fit_ordered_pp[-1], min_r_squared_default), sol_ordered_pp[-1]
            list_pp_fit_max.append(pp_fit_max)
            list_pp_sol_max.append(pp_sol_max)

    print("ready")

    for  key_h in db.db_pairs.keys():
        for  key_p in db.db_pairs[key_h].keys():
            
            if key_p !='struct_healthy':
                
            #check if the healthy is locked

                healthy_pos, hh_fit_max = db.get_list_pos_fit(key_h)
                healthy_poly = Polygon(healthy_pos)
                prion_pos_list, data_list = [], []
                hash_tested = []
                old_hp = db.db_pairs[key_h][key_p]['hp_fitness']
                old_pp = db.db_pairs[key_h][key_p]['pp_fitness']

                #check if healthy,prion pair was already checked 
                #add the prion to the list of new tested pair
                hash_tested.append(int(key_p))
                prion_pos = db.get_list_pos(key_p)
                prion_poly = Polygon(prion_pos)
                #get binding of healthy with prion and get the binded pairs 
                fit_ordered_hp, sol_ordered_hp, output_index_hp, _= pu.estimate_stacking_icp(healthy_poly, prion_poly, avg_dist_tol=0.3, n_point_max = n_point_max,  multiple_bind = False , output_index = True )
            
                #define distance accepted for binding
                time_min_binded = distance_binded/min_r_squared_default

                hp_fit_max, hp_sol_max = pu.rescale_fit(fit_ordered_hp[-1], min_r_squared_default), sol_ordered_hp[-1]
                #pp_fit_max, pp_sol_max = list_pp_fit_max[p_i], list_pp_sol_max[p_i]

                _, binded_hp, _ = pu.get_fitness(healthy_poly, hp_sol_max, min_r_squared = min_r_squared_default, time_min_binded=time_min_binded, output_index = True, compute_intersection=False)
                #_, binded_pp, _ = pu.get_fitness(prion_poly, pp_sol_max, min_r_squared = min_r_squared_default, time_min_binded=time_min_binded, output_index = True, compute_intersection=False)

                hp_fit_relaxed, _, _ = pu.get_fitness(healthy_poly, hp_sol_max, min_r_squared = relaxed_min_r, time_min_binded=time_min_binded, output_index = True, compute_intersection=False)
                hp_fit_relaxed = pu.rescale_fit(hp_fit_relaxed, relaxed_min_r)
                
                print(hp_fit_max, old_hp)


if __name__ == "__main__":
    db = database.Database("src/polygon/data/n=5_5deg_2")
    select_pairs(db)