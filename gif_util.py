import imageio
import os
import matplotlib.pyplot as plt



def create_gif_multiple_molecule(path, bound_to_plot,  count_to_plot, conn_to_plot, n_molecule,  box_side):

    filenames = []
    min_x, max_x, min_y, max_y = -1*box_side, box_side, -1*box_side, box_side
    in_conn, out_conn = conn_to_plot

    # for each frame
    for frame in range(len(in_conn)):

        fig, (ax4, ax3, ax1, ax2) = plt.subplots(4, 1, figsize=(8, 12),
                                                 dpi=240, gridspec_kw={'height_ratios': [1, 1, 1, 3]})

        ax2.set_xlim(min_x-0.05, max_x+0.05)
        ax2.set_ylim(min_y-0.05, max_y+0.05)
        ax2.set_aspect('equal', adjustable='box')

        # for each molecule
        for mol_id in range(n_molecule):

            mol_frame_out_conn = out_conn[frame][mol_id]
            mol_frame_in_conn = in_conn[frame][mol_id]

            for j in range(len(mol_frame_out_conn)):
                xl, yl = mol_frame_out_conn[j]
                ax2.plot(xl, yl, c='b')

            for j in range(len(mol_frame_in_conn)):
                xl, yl = mol_frame_in_conn[j]
                ax2.plot(xl, yl, c='r')

        self_bonded, outer_bounded, specific_binded_number = bound_to_plot

        xl = [k/len(in_conn) for k in range(frame+1)]
        yl_s = self_bonded[:frame+1]
        yl_o = outer_bounded[:frame+1]
        yl_1 = [specific_binded_number[k][1] for k in range(frame+1)]
        yl_2 = [specific_binded_number[k][2] for k in range(frame+1)]
        yl_3 = [specific_binded_number[k][3] for k in range(frame+1)]

        ax1.plot(xl, yl_s, linewidth=3, label='self binding total')
        ax1.plot(xl, yl_o, linewidth=3, label='outer binding total')
        ax1.set_xlim((0, 1))
        ax1.set_ylim((0, 1.2*max([max(self_bonded), max(outer_bounded)])))
        ax1.legend()

        ax3.plot(xl, yl_1, linewidth=3, label='binding 1')
        ax3.plot(xl, yl_2, linewidth=3, label='binding 2')
        ax3.plot(xl, yl_3, linewidth=3, label='binding 3')
        ax3.set_xlim((0, 1))
        ax3.set_ylim((0, 1.2 * max(max(p)
                                   for p in specific_binded_number)))
        ax3.legend()

        count_healthy, count_prion, count_other = count_to_plot

        conf_h = [count_healthy[k] for k in range(frame+1)]
        conf_p = [count_prion[k] for k in range(frame+1)]
        conf_o = [count_other[k] for k in range(frame+1)]

        ax4.plot(xl, conf_h, linewidth=3, c='g', label='healthy count')
        ax4.plot(xl, conf_p, linewidth=3, c='r', label='prion count')
        ax4.plot(xl, conf_o, linewidth=3, c='y', label='other count')
        ax4.set_xlim((0, 1))
        ax4.legend()

        filename = f'{frame}.png'
        filenames.append(filename)

        # save frame
        plt.savefig(filename)
        plt.close()

    # build gif
    with imageio.get_writer(path, mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)

