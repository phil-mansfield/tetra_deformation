import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc

files = ["../data/grids/grid_L500_1.00000_%d.dat" % n for n in range(128)]

def main():
    palette.configure(False)

    skip = 1
    L = 500.0
    
    for fname in files:
        in_halo, rho_rho_m, v_v_init, l1, l2, l3 = np.loadtxt(fname).T

        layer_name = fname[19: -4]
        
        w = 128 // skip
        l = (L*128)/1024
        
        assert w*w == len(in_halo)
    
        in_halo = np.reshape(in_halo, (w, w))
        v_v_init = np.reshape(v_v_init, (w, w))
        rho_rho_m = np.reshape(rho_rho_m, (w, w))
        l1 = np.reshape(l1, (w, w))
        l2 = np.reshape(l2, (w, w))
        l3 = np.reshape(l3, (w, w))
    
        fig, ax = plt.subplots(1, 3, figsize=(18, 6))
        ax[0].set_title(r"${\rm Within\,}R_{\rm vir}$")
        ax[0].set_ylabel(r"$u_y$")
        ax[0].set_xlabel(r"$u_x$")
        ax[0].imshow(in_halo, origin="lower", cmap="binary",
                     extent=[0, l, 0, l])
    
        ax[1].set_xlabel(r"$u_x$")
        ax[1].set_title(r"$\log_{10}(V/V_{\rm init})$")
        ax[1].imshow(v_v_init, origin="lower", cmap="bwr",
                     extent=[0, l, 0, l], vmin=-4, vmax=4)

        ax[2].set_xlabel(r"$u_x$")
        ax[2].set_title(r"$\rho/\langle\rho_m\rangle$")
        im = ax[2].imshow(rho_rho_m, origin="lower", cmap="bwr",
                          extent=[0, l, 0, l], vmin=-4, vmax=4)

        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        
        ax[1].get_yaxis().set_visible(False)

        plt.savefig("../plots/density_%s.png" % layer_name)
        
        cmap, vmin, vmax = "RdGy_r", -5, 5
        fig, ax = plt.subplots(1, 3, figsize=(18, 6))

        ax[0].set_title(r"$\lambda_1$")
        ax[0].set_ylabel(r"$u_y$")
        ax[0].set_xlabel(r"$u_x$")
        ax[0].imshow(l1, origin="lower", cmap=cmap,
                     extent=[0, l, 0, l], vmin=vmin, vmax=vmax)
        
        ax[1].set_title(r"$\lambda_2$")
        ax[1].set_xlabel(r"$u_x$")
        ax[1].imshow(l2, origin="lower", cmap=cmap,
                    extent=[0, l, 0, l], vmin=vmin, vmax=vmax)
        
        ax[2].set_title(r"$\lambda_3$")
        ax[2].set_xlabel(r"$u_x$")
        im = ax[2].imshow(l3, origin="lower", cmap=cmap,
                          extent=[0, l, 0, l], vmin=vmin, vmax=vmax)

        ax[1].get_yaxis().set_visible(False)
        ax[2].get_yaxis().set_visible(False)

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        
        plt.savefig("../plots/lambda_%s.png" % layer_name)
    
if __name__ == "__main__": main()
