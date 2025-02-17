import subprocess
import polygon_util as pu
import math
import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import icp
from functools import cmp_to_key
from sklearn import manifold
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
from matplotlib import offsetbox
import io
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
from matplotlib import offsetbox


# plot correlation between angles 
def plot_correlation_angles(fitness_list, postition_list, shifted_angle = False):

    n_edge = len(postition_list[0])-1

    a_list = [list() for _ in range(n_edge)]
    for pos in postition_list:
        if  shifted_angle:
            angles = get_smaller_angle_description(pos)
        else:
            angles = pu.get_angle_from_pos(pos)
        for i in range(len(angles)):
            a_list[i].append(angles[i])

    fig, axs = plt.subplots(len(a_list), len(a_list), sharex=True, sharey=True, figsize=(4*len(a_list),4*len(a_list)))

    for i in range(len(a_list)):
        for j in range(i,len(a_list)):
            axs[i,j].scatter(a_list[i],a_list[j], c = fitness_list, s=10)

    plt.show()

    return fig, axs

# learn a manifold to use to map other points return the manifold and a panda dataframe with column name column_i
def learn_manifold(postition_list, dim =3, n_neighbors=10):


    # fill array with the positions 
    flatten_pos = []
    for p_vec in postition_list:
        new_p_vec = []
        for pos in p_vec:
            new_p_vec.extend((pos[0], pos[1]))
        flatten_pos.append(new_p_vec)
    flatten_pos = np.array(flatten_pos)


    # create the iso_map for the norestriction 
    iso_manifold  = manifold.Isomap( n_components=dim, n_neighbors=n_neighbors)
    full_manifold = iso_manifold.fit_transform(flatten_pos)
    col = [f'Component {str(i+1)}' for i in range(dim)]
    full_manifold_pd = pd.DataFrame(full_manifold, columns=col)

    return iso_manifold, full_manifold_pd

##########################################
#                 plot                   #
##########################################

def plot_manifold_3d(manifold_df, fig= None, active_axis=None, fitness_list = None, angle_1 = 30, angle_2 = 310, dim= 3, plot_axis = False, alpha = 0.3):
    
    if active_axis is None:
        fig = plt.figure(figsize=(24, 24), dpi=440)
        ax = fig.add_subplot(111,  projection=Axes3D.name)
    else:
        ax = active_axis

    ax.view_init(angle_1, angle_2)
    xs, ys, zs = manifold_df['Component 1'], manifold_df['Component 2'], manifold_df['Component 3']

    if not plot_axis:
        ax.axis('off')
        fig.axes[0].get_xaxis().set_visible(False)
        fig.axes[0].get_yaxis().set_visible(False)

    if dim == 3:
        if fitness_list is not None:
            ax.scatter(xs, ys, zs, c=fitness_list, marker="o", alpha=alpha)
        else:
            ax.scatter(xs, ys, zs, marker="o", c="sandybrown",alpha=alpha)
    if dim == 4:
        ts = manifold_df['Component 4']
        ax.scatter(xs, ys, zs, c=ts, marker="o",alpha=alpha)

    return fig, ax

def get_position_on_manifold(postition_list, manifold,  dim=3):
    flatten_fitted = []
    for p_vec in postition_list:
        new_p_vec = []
        for pos in p_vec:
            new_p_vec.extend((pos[0], pos[1]))
        flatten_fitted.append(new_p_vec)


    flatten_fitted = np.array(flatten_fitted)
    fitted_manifold = manifold.transform(flatten_fitted)
    col = [f'Component {str(i+1)}' for i in range(dim)]
    fitted_manifold_pd = pd.DataFrame(fitted_manifold, columns=col)
    toreturn = [0]*dim
    for i, col_name in enumerate(col):
        toreturn[i] = fitted_manifold_pd[col_name]

    return toreturn


def plot_on_manifold_3d(postition_list, manifold, fig= None, active_axis= None, fitness_list= None, angle_1= 30, angle_2= 310, dim= 3, plot_axis= False, alpha= 0.3):
    
    if active_axis is None:
        fig =  plt.figure(figsize=(24, 24), dpi=440)
        ax = fig.add_subplot(111, projection=Axes3D.name)
    else:
        ax = active_axis

    ax.view_init(angle_1, angle_2)

    xyz = get_position_on_manifold(postition_list, manifold, dim)

    if dim == 3:
        if fitness_list is not None:
            ax.scatter(xyz[0], xyz[1], xyz[2], c=fitness_list, marker="o", alpha=alpha)
        else:
            ax.scatter(xyz[0], xyz[1], xyz[2], marker="o", c="sandybrown",alpha=alpha)
    if dim == 4:
        ax.scatter(xyz[0], xyz[1], xyz[2], c=xyz[3], marker="o",alpha=alpha)


    return fig, ax


def plot_polygon_on_manifold_3d(postition_list, manifold, fig, active_axis, d_min=0.004, dim=3 ):

    # Create a dummy axes to place annotations to
    ax2 = fig.add_subplot(111,frame_on=False) 
    ax2.axis("off")
    ax2.axis([0,1,0,1])

    pos_used = []

    xyz = get_position_on_manifold(postition_list, manifold, dim)

    for i in range(len(xyz[0])):
        xp,yp = [pos[0] for pos in postition_list[i]], [pos[1] for pos in postition_list[i]]
        s = (xyz[k][i] for k in range(dim))
        x,y =  proj(s, active_axis, ax2)
        valid = True
        for x_old, y_old in pos_used:
            d = (x_old-x)**2 + (y_old-y)**2
            if d < d_min:
                valid = False
                break
        if valid:
            image(ax2,c(xp,yp, 0.015),[x,y])
            pos_used.append((x,y))
    
    return fig, active_axis

##########################################
#                 util                   #
##########################################

# method to always sort the angle of the polygon starting with the smallest one. The sorting is based on the
# dictionary sort so if the first angle is the same the second angle is compared. 

def is_a_smaller(a, b, tol):
    is_smaller = True 
    i = 0 
    while is_smaller and i<len(a): 
        if a[i]-b[i] > tol:
            is_smaller =  False
        elif b[i]-a[i] > tol:
            break
        i+=1

    return is_smaller

def my_compare(item1, item2):
    TOL = 1e-09

    if is_a_smaller(item1, item2, TOL):
        return -1
    elif is_a_smaller(item2, item1, TOL):
        return 1
    else:
        return 0

def make_angle_greater_than_0(angles):
    new_angles = []
    for ang in angles:
        if ang<0:
            new_angles.append(2*math.pi+ang)
        else:
            new_angles.append(ang)
    return new_angles


def get_smaller_angle_description(pos, do_print = False):
    angle_list = []
    
    for i in range(len(pos)-1):
        edge_zeroed = i
        init, final = pos[edge_zeroed], pos[(edge_zeroed+1)]
        pos_np = np.array(pos)
        theta_rot = math.atan2(final[1]-init[1], final[0]-init[0])
        rot_mat = icp.rotation_matrix_2d(-theta_rot)
        new_pos_np = np.dot(rot_mat, pos_np.T).T
        x0,y0 = new_pos_np[edge_zeroed,0], new_pos_np[edge_zeroed,1]
        new_pos_np = new_pos_np-np.array((x0,y0))

        angles = pu.get_angle_from_pos([tuple(val) for val in new_pos_np])
        angles = [angles[(k + i)%(len(angles))] for k in range(len(angles))] 
        angles = make_angle_greater_than_0(angles)
        angle_list.append(angles)

    if do_print:
        for ang in angle_list:
            print(ang)

    sorted_angle_list = sorted(angle_list, key=cmp_to_key(my_compare))[0]
    
    if sorted_angle_list[1]>math.pi:
        pos = [(val[0],-val[1]) for val in pos]
        return get_smaller_angle_description(pos, do_print)
    return sorted_angle_list



def proj(X, ax1, ax2):
    """ From a 3D point in axes ax1, 
        calculate position in 2D in ax2 """
    x,y,z = X
    x2, y2, _ = proj3d.proj_transform(x,y,z, ax1.get_proj())
    return ax2.transData.inverted().transform(ax1.transData.transform((x2, y2)))

def image(ax,arr,xy):
    """ Place an image (arr) as annotation at position xy """
    im = offsetbox.OffsetImage(arr, zoom=0.1)
    im.image.axes = ax
    ab = offsetbox.AnnotationBbox(im, xy, xybox=(-30., 30.),
                        xycoords='data', boxcoords="offset points",
                        pad=0.3, arrowprops=dict(arrowstyle="->"))
    ax.add_artist(ab)


def get_image_poly(xp,yp,ddd ):
    img_buf = io.BytesIO()
    dpi = 300
    fig2 = plt.figure( figsize=(1, 1), dpi=dpi )
    
    #fig2.set_size_inches(0.2, 0.2)
    
    ax = fig2.add_subplot(111)
    ax.set_aspect('equal')
    ax.axis('off')
    fig2.axes[0].get_xaxis().set_visible(False)
    fig2.axes[0].get_yaxis().set_visible(False)
    ax.plot(xp,yp, linewidth=3)
    fig2.savefig(img_buf, transparent=True, bbox_inches=0,  dpi=dpi)
    im = Image.open(img_buf)
    plt.close(fig2)
    return im

def measure_distance(pos1, list_pos2):
    with open('tmp.dat', 'w') as f:

        for pos in pos1[:-1]:
            
            f.write(f'{str(pos[0])} {str(pos[1])}' + " \n")
        f.write( " \n" )

        for pos_list in list_pos2:
            for pos in pos_list[:-1]:
                
                f.write(f'{str(pos[0])} {str(pos[1])}' + " \n")
            f.write( " \n" )
    p = subprocess.Popen(['./sim_bin < tmp.dat > out.dat'],  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    with open('out.dat') as f:
        return f.readlines()