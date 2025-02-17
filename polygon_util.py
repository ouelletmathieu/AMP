import itertools
from multiprocessing.dummy import Array
from typing import Tuple
from matplotlib import pyplot
from pandas import DataFrame, array
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import random
import shapely
import math
import numpy as np
from sympy import Integer
from collections import Counter
from typing import NewType
import re
import icp as icp


def euclid_squared_minimum(points: list[list], min_d_squared: float) -> bool:
    """Check if the minimum of the set of points is at least min_d_squared distance.
    Return a boolean.

    Args:
        points (list(list)): ex. [[x1,y1,...,],[x2,y2,...,],...,[xn,yn,...,]]
        min_d_squared (float): min distance required

    Returns:
        bool: if all points are at least min_d_squared distance from each other
    """
    min_okay = True
    for i in range(len(points)-1):
        for j in range(len(points)-1):

            if i<j:
                tot = sum((points[i][p]-points[j][p])**2 for p in range(len(points[i])))
                if tot < min_d_squared:
                    min_okay = False
                    break
        if not min_okay:
            break

    return min_okay

def find_random_poly(n: int, min_d_squared = 0.0) -> Polygon:
    """Create a random polygon with n points and with min_d_squared within all points 
        The random polygon is in the space [0,1]^2. There is no fixed distance requirement 

    Args:
        n (int): number of vertices
        min_d_squared (float, optional): Minimum distance between points. Defaults to 0.

    Returns:
        Polygon: a valid 2d polygon that can be non simple 
    """
    found = False
    poly = []

    while not found: 
        polygon1 = []
        for _ in range(n):
            pt = (random.random(),random.random())
            polygon1.append(pt)
        poly = Polygon(polygon1)

        if poly.is_simple:
            if min_d_squared == 0:
                found = True

            elif euclid_squared_minimum(list(poly.exterior.coords),min_d_squared):
                found = True
    return poly

def move_poly(poly1: Polygon, rx: float, ry: float, radian: float) -> Polygon:
    """move the polygon (create a new one) the old one is unchanged

    Args:
        poly1 (Polygon): polygon before movement. This polygon wont be changed 
        rx (float): distance to move in the x direction
        ry (float): distance to move in the y direction 
        radian (float): radian for rotation 

    Returns:
        Polygon: new polygon (slow) moved, 
    """

    #those two should commute 
    geom1 = shapely.affinity.translate(poly1, xoff=rx, yoff= ry,zoff =0 )
    return shapely.affinity.rotate(
        geom1, radian, origin='centroid', use_radians=True
    )

def get_poly_angle(poly1: Polygon) -> float:
    """get the angle defined by the line between the centroid and the first point

    Args:
        poly1 (Polygon): Polygon

    Returns:
        float: angle defined by the line between the centroid and the first point
    """
    centroid = list(poly1.centroid.coords)[0]
    pt1 = list(poly1.exterior.coords)[0]
    new_pt = (pt1[0] - centroid[0],pt1[1] - centroid[1])
    return math.atan2(new_pt[1],new_pt[0])

def move_to(poly1: Polygon, newx: float, newy: float, newtheta: float) -> Polygon:
    """move a point to the specified x,y and set the angle (see get poly angle)

    Args:
        poly1 (Polygon): polygon to move
        newx (float): x position to move to
        newy (float): y position to move to
        newtheta (float): theta of the new position (center to first point)

    Returns:
        Polygon: return new polygon at the new x,y coordinate with the 
        angle newtheta between the center and the first vertex
    """
    x, y =  list(poly1.centroid.coords)[0]
    theta = get_poly_angle(poly1)
    return move_poly(poly1, newx-x, newy-y, newtheta-theta)

def get_scaled_points(poly: Polygon, sx: float, sy: float) -> list[tuple]:
    """return a list of points scaled by (sx,sy) 

    Args:
        poly (Polygon): polygon to scale
        sx (float): scale in the x direction 
        sy (float): scale in the y direction 

    Returns:
        list[list]: list of scaled point  [(x1,y1),(x2,y2),...,(xn,yn)]
    """
    return [(sx*pt[0], sy*pt[1]) for pt in list(poly.exterior.coords)]

def get_fitness(poly1: Polygon, poly2: Polygon, min_r_squared = 0.001, time_min_binded = 1, output_index = False, compute_intersection = True)-> tuple[float,list,float]:
    """return the fitness 1/r^2 and min_r_squared is the cap for small r

    Args:
        poly1 (Polygon): Polygon1 
        poly2 (Polygon): Polygon2
        min_r_squared (float, optional): minimum distance before capping the 1/r2 potential. Defaults to 0.001.
        time_min_binded (float, optional): ratio of min_r_squared to consider being binded . Defaults to 1.5.
        output_index (bool, optional): If true return pair that are closer than time_min_binded*min_r_squared together as being binded. Defaults to False.

    Returns:
        If output_index is false then only return fitness else it return:
        Tuple[float,list,float]: fitness, list of binded pairs, overlaping factor
        
        fit: fitness
        set_pair_binded: list of pair ((poly1_i,poly2_i),...) where  poly1_i is binded with poly2_i
        inter: overlaping factor
    """
    set_pair_binded = []        #pairs of binded point (r < min_r_squared) (poly1, poly2)
    pt1, pt2 = list(poly1.exterior.coords), list(poly2.exterior.coords)
    fit = 0
    inter, inter2 = 1, 1
    valid = 1

    #check if they are one on top of eachother. If they are return 0
    full_dist = 0
    for i in range(len(pt1)):
        full_dist += (pt1[i][0]-pt2[i][0])**2 + (pt1[i][1]-pt2[i][1])**2
    if full_dist < min_r_squared:
        valid = 0

    #check the amount of intersection. Can return false empty intersection 
    #therefore we move the polylgon and check it a second time
    if poly1.intersects(poly2):
        inter -= 2*poly1.intersection(poly2).area/(poly1.area + poly2.area)
        if inter<0:
            inter=0
    #move poly a little bit since if htey share the same points we get false negative
    if poly1.intersects(move_poly(poly2,0.0001,0.0001,0.001)):
        inter2 -= 2*poly1.intersection(poly2).area/(poly1.area + poly2.area)
        if inter2<0:
            inter2=0
    if inter2<inter:
        inter=inter2

    for i in  range(len(pt1)-1):
        for j in  range(len(pt2)-1):
            d = ((pt1[i][0]-pt2[j][0])**2 + (pt1[i][1]-pt2[j][1])**2)
            if d<=min_r_squared:
                d = min_r_squared
            if d<=time_min_binded*min_r_squared:
                set_pair_binded.append((i,j))
            fit+= 1.0/(d**2)

    if fit<0.00000000001:
        fit = 0.00000000001
    if compute_intersection:
        full_fit = inter*math.log(fit)
    else:
        full_fit = math.log(fit)

    if output_index:
        return full_fit*valid, set_pair_binded, inter

    return full_fit*valid

def center_poly(poly1: Polygon) -> Polygon:
    """center the polygon at (0,0)

    Args:
        poly1 (Polygon): Polygon to move

    Returns:
        Polygon: new Polygon
    """

    return move_to(poly1, 0, 0, 0)

def farthest_from_centroid(poly1: Polygon) -> float:
    """return the length farthest point from centroid 

    Args:
        poly1 (Polygon): Polygon

    Returns:
        float: the distance of the point the farthest to the centroid
    """
    pt_list= list(poly1.exterior.coords)
    centroid = list(poly1.centroid.coords)[0]
    max_dist_2 = 0
    for pt in pt_list:
        d = ((pt[0]-centroid[0])**2 + (pt[1]-centroid[1])**2)
        if  d > max_dist_2:
            max_dist_2 = d

    return math.sqrt(max_dist_2)

########################################################
#      Potential based binding  (behave poorly)        #
########################################################

def find_a_good_start(poly1: Polygon, poly2: Polygon, min_r_squared = 0.001, r_list = None, n_center = 40, n_second = 40) -> tuple[Polygon, list[float], list[Polygon], list[list[tuple]]]:
    """try to place polygon 2 at n_center angle arround polygon 1 at all distance in r_list with the polygon rotated 
    Args:
        poly1 (Polygon): Polygon1
        poly2 (Polygon): Polygon2
        min_r_squared (float, optional): distance where the potential 1/r2 is capped. Defaults to 0.001.
        r_list (list, optional): list of distance centroid centroid to try from polygon1 to polygon2 . Defaults to [0.7,0.85,1,1.15,1.3,1.5].

    Returns:
        Tuple[Polygon, list[float], list[Polygon], list[list[tuple]]]: 
            first polygon, list of fitness for all starting pos, polygon for all starting pos, list of pos
    """
    if r_list is None:
        r_list = [0.7,0.85,1,1.15,1.3,1.5]
    poly1 = center_poly(poly1)
    main_angle_center = [i*2*math.pi/n_center for i in range(n_center)]
    angle_second = [i*2*math.pi/n_second for i in range(n_second)]

    print(main_angle_center)

    d1 = farthest_from_centroid(poly1)
    d2 = farthest_from_centroid(poly2)
    distance = d1/2 + d2/2

    pos_list = []
    fit_list = []
    poly_list = []

    for r in r_list:
        for a1 in main_angle_center:
            for a2 in angle_second:
                poly_new = move_to(poly2, r*distance*math.cos(a1), r*distance*math.sin(a1), a2)
                poly_list.append(poly_new)
                fit_list.append(get_fitness(poly1, poly_new, min_r_squared))
                pos_list.append((r*distance*math.cos(a1), r*distance*math.sin(a1), a2))
    return poly1, fit_list, poly_list, pos_list

def estimate_derivative(poly1: Polygon, poly2: Polygon, h_xy: float, h_theta: float, min_r_squared = 0.001) -> tuple[float,float,float]:
    """estimate the derivative of the fitness for the capped 1/r2 potential given a spatial delta of h_xy and a rotation delta of h_theta. 

    Args:
        poly1 (Polygon): Polygon1
        poly2 (Polygon): Polygon2
        h_xy (float):  spatial delta 
        h_theta (float): angle theta (in rad)
        min_r_squared (float, optional):  distance where the potential 1/r2 is capped. Defaults to 0.001.

    Returns:
        Tuple[float,float,float]: partial_f_x, partial_f_y, partial_f_t
    """

    fxph = get_fitness(poly1, move_poly(poly2, h_xy, 0, 0), min_r_squared)
    fxmh = get_fitness(poly1, move_poly(poly2, -h_xy, 0, 0), min_r_squared)
    fyph = get_fitness(poly1, move_poly(poly2, 0, h_xy,  0), min_r_squared)
    fymh = get_fitness(poly1, move_poly(poly2, 0, -h_xy, 0), min_r_squared)
    ftph = get_fitness(poly1, move_poly(poly2, 0, 0,  h_theta), min_r_squared)
    ftmh = get_fitness(poly1, move_poly(poly2, 0, 0, -h_theta), min_r_squared)

    partial_f_x  = (fxph - fxmh)/(2*h_xy)
    partial_f_y  = (fyph - fymh)/(2*h_xy)
    partial_f_t  = (ftph - ftmh)/(2*h_xy)

    return partial_f_x, partial_f_y, partial_f_t

def move_in_gradient(poly1: Polygon, poly2: Polygon, rate: float, h_xy: float, h_theta: float, min_r_squared = 0.001) -> tuple[Polygon, Polygon]:
    """move polygon 2 in gradient by rate by estimating the derivative (See estimate_derivative) with delta h_xy and h_theta  
        return the position of the found local minima as a new polygon2. Polygon 1 is returned too
    Args:
        poly1 (Polygon): Polygon1
        poly2 (Polygon): Polygon2
        rate (float): the ratio of the length of the gradient that will be used to move
        h_xy (float): spatial delta
        h_theta (float): angle theta (in rad)
        min_r_squared (float, optional): distance where the potential 1/r2 is capped. Defaults to 0.001.

    Returns:
        Tuple[Polygon, Polygon]: polygon1 inputed is returned and a new polygon 2 at the local minima
    """
    partial_f_x, partial_f_y, partial_f_t = estimate_derivative(poly1, poly2, h_xy, h_theta, min_r_squared)
    poly2_new = move_poly(poly2, rate*partial_f_x, rate*partial_f_y, rate*partial_f_t)

    return poly1, poly2_new

def estimate_stacking(poly1: Polygon, poly2: Polygon, n_keep = 12, rate = 0.001, h_xy = 0.001, h_theta = 0.001,  min_r_squared = 0.001 ) -> tuple[float, float, Polygon, Polygon]:
    """estimate the best way to stackl two polygon using a gradient and 1/r2 potential approach. It does not work really well. 

    Args:
        poly1 (Polygon): Polygon1
        poly2 (Polygon): Polygon2
        n_keep (int, optional): _number of initial position kept for gradient descent. Defaults to 12.
        rate (float, optional): the ratio of the length of the gradient that will be used to move. Defaults to 0.001.
        h_xy (float, optional): spatial delta. Defaults to 0.001.
        h_theta (float, optional):angle theta (in rad). Defaults to 0.001.
        min_r_squared (float, optional): distance where the potential 1/r2 is capped. Defaults to 0.001.

    Returns:
        Tuple[float, float, Polygon, Polygon]: best initial position fitness, best final position fitness,  
            best initial position polygon, best final position polygon (all for polygon 2)
    """
    poly1, fit_list, poly_list, _ = find_a_good_start(poly1, poly2,  min_r_squared)
    ind_and_fit = [(i,fit_list[i]) for i in range(len(fit_list))]
    ind_and_fit.sort(key = lambda x: x[1])
    inds = ind_and_fit[-n_keep:]

    best_before = inds[-1][1]
    best_found, best_poly = 0, 0

    for ind, fit_b in inds:
        best_for_ind, best_poly_for_ind = 0,0
        poly2 = poly_list[ind]
        for _ in range(40):
            fit = get_fitness(poly1, poly2, min_r_squared)
            if fit >best_for_ind: 
                best_for_ind = fit
                best_poly_for_ind = poly2
            poly1, poly2 = move_in_gradient(poly1, poly2, rate, h_xy, h_theta, min_r_squared)

        if best_for_ind > best_found:
            best_found = best_for_ind
            best_poly = best_poly_for_ind

    return best_before, best_found, poly1, best_poly

########################################################
#                   Point matching                     #
########################################################

Sol_type = NewType('Sol_type', tuple[list[float], list[tuple[float,float]]])

def move_to_center(set_pt1):
    cg_1 = np.sum(set_pt1,0)/len(set_pt1)
    return set_pt1 - cg_1 

def find_good_initial_guess(set_pt1: np.ndarray, set_pt2: np.ndarray, data_to_shift: np.ndarray = None ) -> tuple[np.ndarray,np.ndarray]:
    """Find a good initial guess by fitting the centroid of set_pt1 and set_pt2 and keepping their main angle the same where
        by main angle we mean the best line that pass through the cloud. Two solution for the two direction of the line 
        are returned. If  data_to_shift is not None the array data_to_shift is transformed according to the two direction  

    Args:
        set_pt1 (np.ndarray): (n x dim) position for subset of points of polygon 1 
        set_pt2 (np.ndarray): (n x dim)  position for subset of points of polygon 2
        data_to_shift (np.ndarray, optional): (n x dim) array of thing that will be shifted in 
            accordance with how the second set is moved is moved.

    Returns:
        Tuple[np.ndarray,np.ndarray]: return best fit in the two direction 
        Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]: return best fit in the two direction and data_to_shift also moved in the two directions.
    """
    cg_1 = np.sum(set_pt1,0)/len(set_pt1)
    cg_2 = np.sum(set_pt2,0)/len(set_pt2)
    trans = set_pt2 - cg_2 
    transc = np.copy(trans)

    m1, b1 = np.polyfit([l[0] for l in set_pt1],  [l[1] for l in set_pt1], 1)
    m2, b2 = np.polyfit([l[0] for l in set_pt2],  [l[1] for l in set_pt2], 1)

    T = icp.rotation_matrix_2d( math.atan2(m1,1) - math.atan2(m2,1) )
    trans1 = np.dot(T, trans.T).T
    trans1 = trans1 +  cg_1
    #dist1, _ = icp.nearest_neighbor(set_pt1, trans1)

    T2 = icp.rotation_matrix_2d( math.atan2(m1,1) - math.atan2(m2,1) + math.pi ) 
    T2[:,0] = T2[:,0]
    trans2 = np.dot(T2, transc.T).T
    trans2 = trans2 +  cg_1
    #dist2, _ = icp.nearest_neighbor(set_pt1, trans2)

    if not isinstance(data_to_shift, type(None)):
        data = np.copy(data_to_shift) - cg_2 
        data1 =  np.dot(T, data.T).T +  cg_1
        data = np.copy(data_to_shift) - cg_2 
        data2 =  np.dot(T2, data.T).T +  cg_1
        return trans1, data1, trans2, data2
        
    else:
        return trans1, trans2

def transform_homogeneous_2d(array: np.ndarray, transform: np.ndarray) ->  np.ndarray:
    """transform a (nx2) np.ndarray given a transformation matrix transform

    Args:
        array (np.ndarray): _description_
        transform (np.ndarray): _description_

    Returns:
        np.ndarray: transformed np.array (nx2) of point 
    """
    C1 = np.ones((len(array), 3))
    C1[:,0:2] = array
    C1 = np.dot(transform, C1.T).T
    return C1[:,0:2]

def translate_np_array(src_array: np.ndarray, id_src: int, target_array: np.ndarray, id_target: int)-> np.ndarray:
    """translate src_array based on the position of its vertex id_src to the position of id_target vertex of  target_array

    Args:
        src_array (np.ndarray): (n x dim) array to translate
        id_src (int): id of vertex of src_array to use for position
        target_array (np.ndarray): (n x dim) array to act as target
        id_target (int): id of the target vertex of target_array

    Returns:
        np.ndarray: translated array 
    """
    pts_mapped = np.array( (target_array[id_target,0], target_array[id_target,1]) )
    pts_poly = np.array( ( src_array[id_src,0], src_array[id_src,1]) )

    return  src_array + pts_mapped-pts_poly

def check_self_matching(n_point: int, poly1: Polygon,  avg_dist_tol=0.1, ratio_max_overlap_area = 0.02) -> Polygon:
    """ (Prefer using the check_matching method with two time the same polygon). Find the best matching by trying to fit n_point points using the ICP method ans initial 
        guess using linear line fit. This method consider a subset of point smalller than n_point. 

    Args:
        n_point (int): max number of point (excluded) to consider in the matching
        poly1 (Polygon): Polygon 
        avg_dist_tol (float): average distance for the subset of n_points in the ICP algotrithm 
        ratio_max_overlap_area (float): ratio of the area of poly1 that can overlap when they are binded  

    Returns:
        Polygon: solution for polygon 
    """
    area_poly1 = poly1.area
    max_inter = ratio_max_overlap_area*area_poly1
    real_point = list(poly1.exterior.coords)[:-1]

    list_of_index = []
    list_of_subset = []
    list_found_poly = []

    for i in range(len(real_point)):
        lst_sel_pt = []
        lst_sel_ind = []
        for j in range(n_point):
            ind = (i + j)%len(real_point)
            lst_sel_pt.append(real_point[ind])
            lst_sel_ind.append(ind)
        list_of_subset.append(lst_sel_pt)
        list_of_index.append(lst_sel_ind)

    for i in range(len(list_of_subset)):
        for j in range(i+1,len(list_of_subset)):

            first_index = list_of_index[j][0]
            A = np.array(list_of_subset[i])
            B = np.array(list_of_subset[j])
            all_pts_np = np.array(real_point)

            B1, all_pts_new1, B2, all_pts_new2 = find_good_initial_guess(A, B, data_to_shift = all_pts_np)

            T1, dist1, nb_step1 = icp.icp(B1, A)
            T2, dist2, nb_step2 = icp.icp(B2, A)

            B1 = transform_homogeneous_2d(B1, T1)
            B2 = transform_homogeneous_2d(B2, T2)

            #if found mathching is good enough
            if np.average(dist1) < avg_dist_tol :

                poly2_pts_1 = transform_homogeneous_2d(all_pts_new1, T1)
                poly2_pts_1 = translate_np_array(poly2_pts_1, first_index, B1, 0)
                poly2_1 = Polygon(tuple(map(tuple, poly2_pts_1)) )
                inter_1 = poly1.intersection(poly2_1).area

                if  inter_1<max_inter:
                    list_found_poly.append(poly2_1)

                """               
                plt.plot([l[0] for l in A], [l[1] for l in A])
                plt.plot([l[0] for l in B1], [l[1] for l in B1])
                plt.plot([l[0] for l in poly2_pts_1], [l[1] for l in poly2_pts_1])
                plt.show()
                """ 

            if np.average(dist2) < avg_dist_tol :
                poly2_pts_2 = transform_homogeneous_2d(all_pts_new2, T2)
                poly2_pts_2 = translate_np_array(poly2_pts_2, first_index, B2, 0)
                poly2_2 = Polygon(tuple(map(tuple, poly2_pts_2)) )
                inter_2 = poly1.intersection(poly2_2).area

                #check if they overlap      
                if  inter_2<max_inter:
                    list_found_poly.append(poly2_2)

                """
                plt.plot([l[0] for l in A], [l[1] for l in A])
                plt.plot([l[0] for l in B2], [l[1] for l in B2])
                plt.plot([l[0] for l in poly2_pts_2], [l[1] for l in poly2_pts_2])
                plt.show()
                """

    return list_found_poly

def check_matching(n_point: int, poly1: Polygon, poly2: Polygon, no_icp_threshold=0.01, avg_dist_tol=0.1, ratio_max_overlap_area = 0.02):
    """Find the best matching  between poly1 and poly2 by trying to fit n_point points using the ICP method ans initial 
        guess using linear line fit. This method consider a subset of point smalller than n_point. 

    Args:
        n_point (int): max number of point to try to fit (excluded)
        poly1 (Polygon): polygon 1
        poly2 (Polygon): polygon 2
        avg_dist_tol (float): average distance for the subset of n_points in the ICP algotrithm 
        ratio_max_overlap_area (float): ratio of the area of poly1 that can overlap when they are binded  


    Returns:
        Polygon: return the best polygon fit
    """
    #compute max overlap allowed
    area_poly = (poly1.area + poly2.area)/2
    max_inter = ratio_max_overlap_area*area_poly
    #get stationarry point of poly1 as numpy
    real_point_1, real_point_2 = (
        list(poly1.exterior.coords)[:-1],
        list(poly2.exterior.coords)[:-1],
    )


    #only works if it is the same size
    if len(real_point_1)!=len(real_point_2):
        print("polygon are not of the same size")
        return []

    list_of_index = []                              #index of point for each subset for both poly to test the mating 
    list_of_subset_1,list_of_subset_2 = [], []      #index of point for each subset for both poly
    list_found_poly = []                            #list of solution
    list_of_index = []                              #pairs ((i,j,k), (l,m,n)) of matching index

    for i in range(len(real_point_1)):                                      #
        lst_sel_pt_1,lst_sel_pt_2  = [], []
        lst_sel_ind = []
        for j in range(n_point):
            lst_sel_pt_1.append(real_point_1[(i + j)%len(real_point_1)])
            lst_sel_pt_2.append(real_point_2[(i + j)%len(real_point_2)])
            lst_sel_ind.append((i + j)%len(real_point_1))

        list_of_subset_1.append(lst_sel_pt_1)
        list_of_subset_2.append(lst_sel_pt_2)
        list_of_index.append(lst_sel_ind)

    for i, j in itertools.product(range(len(list_of_subset_1)), range(len(list_of_subset_2))):
        all_pts_2_np =  np.array(real_point_2)
        first_index_2  = list_of_index[j][0]            #index of the fist point for each subset of prion
        index_set_1, index_set_2 = list_of_index[i], list_of_index[j]
        A, B  = np.array(list_of_subset_1[i]), np.array(list_of_subset_2[j])


        B1, all_pts_2_np_new1, B2, all_pts_2_np_new2 = find_good_initial_guess(A, B, data_to_shift = all_pts_2_np)

        if n_point!=2:

            dist0 = icp.nearest_neighbor(B1, A)[0]

            if np.mean(dist0)>no_icp_threshold:
                T1, dist1, nb_step1 = icp.icp(B1, A)
                T2, dist2, nb_step2 = icp.icp(B2, A)
                B1, B2 = transform_homogeneous_2d(B1, T1), transform_homogeneous_2d(B2, T2)  
            else:
                dist1,dist2 = dist0, dist0
                T1,T2 = np.identity(3),np.identity(3)
        else:
            dist1,dist2 = 0, 0
            T1,T2 = np.identity(3),np.identity(3)



        #if found mathching is good enough for transformation 1 
        if np.average(dist1) < avg_dist_tol : 
            poly2_pts_1 = transform_homogeneous_2d(all_pts_2_np_new1, T1)
            poly2_pts_1 = translate_np_array(poly2_pts_1, first_index_2, B1, 0)
            poly2_1 = Polygon(tuple(map(tuple, poly2_pts_1)) )
            inter_1 = poly1.intersection(poly2_1).area
            #check if they overlap
            if  inter_1<max_inter and not poly1.within(poly2_1) and not poly2_1.within(poly1):
                list_found_poly.append(poly2_1)
        #if found mathching is good enough for transformation 2 
        if np.average(dist2) < avg_dist_tol :
            poly2_pts_2 = transform_homogeneous_2d(all_pts_2_np_new2, T2)
            poly2_pts_2 = translate_np_array(poly2_pts_2, first_index_2, B2, 0)
            poly2_2 = Polygon(tuple(map(tuple, poly2_pts_2)) )
            inter_2 = poly1.intersection(poly2_2).area
            #check if they overlap      
            if  inter_2<max_inter and not poly1.within(poly2_2) and not poly2_2.within(poly1) :
                list_found_poly.append(poly2_2)

    return list_found_poly

def is_one_to_one_binding(binding_pairs: tuple[list[Integer],list[Integer]]) -> bool:
    # sourcery skip: simplify-len-comparison
    """check if the binding between two polygon is one to one

    Args:
        binding_pairs (list[list[int],list[int]]): ex. ([1,6],[7,6],[8,9]) where here 1 and 7 
        of polygon 1 is binded with 6 and 8 with 9

    Returns:
        bool: return false if more than one vertice is binded to the same vertice or the opposite
    """
    if len(binding_pairs)==0:
        return False

    bound1, bound2 = Counter([x[0] for x in binding_pairs]),Counter([x[1] for x in binding_pairs])

    if bound1.most_common()[0][1] > 1 : #check if the most common binded atom has more than one binded to it:
        return False
    if bound2.most_common()[0][1] > 1 : #check if the most common binded atom has more than one binded to it:
        return False

    return True

def estimate_stacking_icp(poly1: Polygon, poly2: Polygon, n_point_max = 4, avg_dist_tol=0.1, ratio_max_overlap_area = 0.05, \
     output_index = False, multiple_bind = True ) -> tuple[list[float], list[Polygon], list[list], list[float]]:
    """ gives a list of potential staking of the two polygon poly1 and poly2 considering a maximum of n_point_max matching (included)
    use the method check_matching for each n and check for multiple binding if asked and order the solution.  

    Args:
        poly1 (Polygon): polygon 1
        poly2 (Polygon): polygon 2
        n_point_max (int, optional): maximum number of point to try to bind. Defaults to 4.
        avg_dist_tol (float, optional): average distance for the subset of n_points in the ICP algotrithm 
        ratio_max_overlap_area (float): ratio of the area of poly1 that can overlap when they are binded  
        multiple_bind (bool, optional): If false no multiple bind is accepted. Defaults to True.

    Returns:
    if output index:
        tuple[list[float], list[Polygon], list[list], list[float]]: 
        fit_ordered: list of fitness for all solution in order
        sol_ordered: list of polygon in order of fitness
        index_ordered: lisf of (list of pair ((poly1_i,poly2_i),...) where  poly1_i is binded with poly2_i)
        intersect_ordered: list of the penality for each 
    """
    sol_ordered, fit_ordered, index_ordered, intersect_list = [],[], [], []

    uncheck_sol = []
    for n_point in range(2, n_point_max+1):
        uncheck_sol.extend(check_matching(n_point, poly1, poly2,  avg_dist_tol, ratio_max_overlap_area))


    for sol in uncheck_sol:
        fit, index_pair_fit, intersect = get_fitness(poly1, sol, output_index = True)
        #if not multiple bind need to check if there is more than one atom binded on the same site. 
        okay = True
        if not multiple_bind and not is_one_to_one_binding(index_pair_fit):
            okay = False
        if okay:
            sol_ordered.append(sol)
            fit_ordered.append(fit)
            index_ordered.append(index_pair_fit)
            intersect_list.append(intersect)

    ordered_all = list(
        sorted(
            zip(fit_ordered, sol_ordered, index_ordered, intersect_list),
            key=lambda pair: pair[0],
        )
    )


    fit_ordered, sol_ordered, index_ordered, intersect_ordered  = [x[0] for x in ordered_all],  [x[1] for x in ordered_all], [x[2] for x in ordered_all], [x[3] for x in ordered_all]

    if output_index:
        return fit_ordered, sol_ordered, index_ordered, intersect_ordered

    return fit_ordered, sol_ordered 

########################################################
#        molecular distance geometry problem           #
########################################################

def update(list_length: list[float], sol: Sol_type, theta: float) -> Sol_type:
    """utility to find random polygons with same edges length. It increase the solution vector adding a new edge with
    the good edge length and with and angle theta with the previous edge
    (Sol_type = Tuple[list[float], list[Tuple[float,float]]])

    Args:
        list_length (list[float]): example [1,1,1,1,1] generates pentagon with edge length of 1 
        sol (Tuple[list[float], list[Tuple[float,float]]]): sol has the shape ((angles), (current pos))
            angles: is the list of angles from between links grow each time a link is added
            current pos: (x,y) position of the last point of the last added link. A valid polygon should go back to (0,0)
            list_length: length of edges of the polygon. They are in order of appearance from (0,0)
        theta (float): new angle to add

    Returns:
       Sol_type : return and updated solution (like input sol)
    """
    worker_id = len(sol[0])
    if worker_id == len(list_length):
        return True, sol
    new_pos = (sol[1][0] + list_length[worker_id]*math.cos(theta), sol[1][1] + list_length[worker_id]*math.sin(theta))
    return False, (sol[0] + (theta,), new_pos)

def get_pos(list_length: list[float], sol: Sol_type) -> list[tuple[float]]:
    """get the position list of each point to use in shapely (warning read full description) or do computation. Last position should be close to [0,0] but
    cannot be used directly in shapely. The last position need to be replaced by (0,0) to fully close the polygon. 


    Args:
        list_length (_type_): _ example [1,1,1,1,1] generates pentagon with edge length of 1 
        sol (Tuple[list[float], list[Tuple[float,float]]]): sol has the shape ((angles), (current pos))
            angles: is the list of angles from between links grow each time a link is added
            current pos: (x,y) position of the last point of the last added link. A valid polygon should go back to (0,0)
            list_length: length of edges of the polygon. They are in order of appearance from (0,0)

    Returns:
        list[Tuple(float)]: list of position [(x1,y1), (x2,y2), ..., (xn,yn)]
    """
    pos = [(0,0),]
    for i in range(len(list_length)):
        new_pos = (pos[-1][0] + list_length[i]*math.cos(sol[0][i]), pos[-1][1] + list_length[i]*math.sin(sol[0][i]))
        pos.append(new_pos)

    return pos

def random_sol(list_length: list[float], min_angle: float) -> Sol_type:
    """create a random solution for a given list_length (may not be good) return the sol and the distance of the last point squared. 
        Should be 0 if it was a perfect solution

    Args:
        list_length (list[float]): 
        min_angle (float): minimum angle in the solution. If 0 degenerated polygon are accepted

    Returns:
        Sol_type: return solution, shape ((angles), (current pos))
            angles: is the list of angles from between links grow each time a link is added
            current pos: (x,y) position of the last point of the last added link. A valid polygon should go back to (0,0)
            list_length: length of edges of the polygon. They are in order of appearance from (0,0)
    """

    sol = (tuple(),(0,0))
    for _ in range(len(list_length)):
        found = False
        if len(sol[0]) > 0:
            last_angle = sol[0][-1]
            while not found:
                new_theta = random.random()*2*math.pi
                delta = abs(((math.pi - last_angle)+new_theta)%(2*math.pi))
                if delta > min_angle and delta < 2*math.pi-min_angle: 
                    found = True
        else:
            new_theta = 0 #first angle always 0

        _, sol = update(list_length, sol, new_theta)
    return sol, sol[1][0]**2 + sol[1][1]**2

def find_set_of_conformer_polygon(list_length:list[float], max_error_squared=0.001, max_try = 100000, is_simple = True, min_angle= 0.174533) -> list[Sol_type]:
    """create a list of conformer polygon with all the same list_length e.x. [1,1,1,1,1] for pentagon with same edge lenght 
        and a maximal squared distance from the origin of the last vertex given by max_error_squared (perfect solution should be 0).
        The last edge is always a little bit shorter or longer.

    Args:
        list_length (list[float]):  length of edges of the polygon. They are in order of appearance from (0,0)
        max_try (int):
        max_error_squared (float): max squared distance the solution can be of respecting the length in list_length. 
        is_simple (bool): if True solution will not intersect 
        min_angle (float): minimum angle in the solution. If 0 degenerated polygon are accepted
        
    Returns:
        Sol_type: return solution, shape ((angles), (current pos))
            angles: is the list of angles from between links grow each time a link is added
            current pos: (x,y) position of the last point of the last added link. A valid polygon should go back to (0,0)
            list_length: length of edges of the polygon. They are in order of appearance from (0,0)
    """    
    list_sol = []
    for _ in range(max_try):
        sol, fit = random_sol(list_length,min_angle)
        delta = abs(((math.pi -sol[0][-1])+sol[0][0])%(2*math.pi))
        if fit < max_error_squared and delta > min_angle and delta < 2*math.pi-min_angle:
            position = get_pos(list_length, sol)[:-1] 
            position.append((0,0))

            if is_simple:
                poly = Polygon(position)
                if poly.is_simple:
                    list_sol.append(position)
            else:
                list_sol.append(position)

    return list_sol


########################################################
#                        IO                            #
########################################################

def sol_string_to_poly(sol_str: str, n: int) -> Polygon:
    """convert a solution string to a polygon

    Args:
        sol_str (str): string representation of a solution vector
        n (int): number of vertex

    Returns:
        Polygon: the corresponding polygon
    """

    list_part  = re.split('\(|\)',sol_str)
    angles = tuple(float(a) for a in list_part[1].split(','))
    last_pos = tuple(float(a) for a in list_part[3].split(','))

    pts = get_pos([1]*n, (angles,last_pos))
    return Polygon(pts)

def rescale_fit(fit: float, min_r_squared = 0.001)-> float:
    """return the rescaled fitness where 3 mean that 3 point are at the min_r_squared distance

    Args:
        fit (float): fitness in the 1/r2 form 
        min_r_squared (float, optional): the value used to cap 1/r2 in the fitness . Defaults to 0.001.

    Returns:
        float: the new fitness in the number of binded point form
    """
    return (min_r_squared**2)*math.exp(fit)

def read_txt_file_fit_pos(df: DataFrame) -> tuple[list[float], list[list[tuple[float]]]]:
    """given a dataframe with one column formatted like:
        fitness, sol 

    Args:
        df (DataFrame): dataframe see description 

    Returns:
        fit_list, pos_list : tuple(list(float), list(list(tuple(float))))
            fit_list: list of fitness for each solution
            pos_list: list of position 
    """
    column = df['fit'] if 'fit' in df else df.iloc[:, 0]
    pos_list = []
    fit_list = []
    for nb in range(len(column)):
        list_part  = re.split('\[|\]',column[nb])
        fit = float(list_part[0])
        list_pos_xy = re.split('\(|\)',list_part[1])
        position = [
            tuple(
                float(val.replace('\'', ''))
                for val in list_pos_xy[2 * i + 1].split(',')
            )
            for i in range(len(list_pos_xy) // 2)
        ]

        pos_list.append(position)
        fit_list.append(rescale_fit(fit))

    return fit_list, pos_list

def get_angle_from_pos(pos_list:list[tuple[float]]) -> list[float]:
    """transform the list of position to an angle list 

    Args:
        pos_list (list[tuple(float)]): position list ex: [(x1,y1),(x2,y2),...,(x1,y1)]

    Returns:
        angle list list(float): ex. [0.4, 3.2, 6.0,...]
    """
    angles = []
    for i in range(len(pos_list)-1):
        ind1, ind2 = i, i+1

        dx = pos_list[ind2][0]-pos_list[ind1][0]
        dy = pos_list[ind2][1]-pos_list[ind1][1]
        angles.append(math.atan2(dy,dx))
    return angles

#translation healthy index -> prion index
def find_missing_translation(length_list, partial_translation):
    #ex trans = 0:3, 3:5, 6:6
    #prion = [0,1,2,3,4,5,6]
    #healthy = [3,?,?,5,?,6]
    #index not none = 5

    index = list(partial_translation.keys())[0]
    val = partial_translation[index]

    #generate two possible polygon full translation 
    clockwise, anti_clockwise = {},{}

    for i in range(len(length_list)):
        clockwise[(index+i)%len(length_list)] = (val+i)%len(length_list)
        anti_clockwise[(index+i)%len(length_list)] = (val-i)%len(length_list)
    clockwise_okay, anti_clockwise_okay = True, True
    for key in partial_translation.keys():
        if clockwise[key]!= partial_translation[key]:
            clockwise_okay=False
        if anti_clockwise[key]!= partial_translation[key]:
            anti_clockwise_okay=False
        if not clockwise_okay and not anti_clockwise_okay:
            break

    if clockwise_okay and check_if_length_match(length_list, clockwise):
        return clockwise
    elif anti_clockwise_okay and check_if_length_match(length_list, anti_clockwise):
        return anti_clockwise
    else:
        return []

#translation healthy index -> prion index    
def check_if_length_match(length_list, complete_translation):

    for i in range(len(length_list)):
        p_ind_start = complete_translation[i]
        if length_list[i]!=length_list[p_ind_start]:
            return False
    
    return True

def to_tuple(list_pos):
    return tuple(list(list_pos))

########################################################
#                        Utils                         #
########################################################





########################################################
#                        PLOTS                         #
########################################################

# TODO
def test_plot_fit(poly1, poly2):
    poly1n = center_poly(poly1)
    poly2n = center_poly(poly2)
    fit_list = []
    poly2n = move_poly(poly2n, 2, 0, 0.4)
    for _ in range(200):
        poly2n = move_poly(poly2n, -0.01, 0, 0.4)
        fit_list.append(get_fitness(poly2n, poly1n))

    return fit_list
    

    

