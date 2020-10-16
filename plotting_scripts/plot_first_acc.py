import numpy as np
import matplotlib.pyplot as plt
import array
import palette
from palette import pc
import scipy.stats as stats
import numpy.random as random

L = 62.5
l = L / 8
data_dir = "/home/phil/code/src/github.com/phil-mansfield/tetra_deformation/data/L63/Lagrangian_info"

def main():
    palette.configure(True)
    
    host_id = array.array("l")
    first_acc = array.array("f")
    density = array.array("f")
    lambda1 = array.array("f")
    lambda2 = array.array("f")
    lambda3 = array.array("f")

    first_acc_file = "%s/first_acc000.dat" % data_dir
    density_file = "%s/L_density000.dat" % data_dir

    w = 128
    n = w*w*w
    
    with open(first_acc_file, "rb") as fp:
        host_id.fromfile(fp, n)
        first_acc.fromfile(fp, n)

    with open(density_file, "rb") as fp:
        density.fromfile(fp, n)
        lambda1.fromfile(fp, n)
        lambda2.fromfile(fp, n)
        lambda3.fromfile(fp, n)

    first_acc = np.reshape(first_acc, (w, w, w))
    density = np.log10(np.reshape(density, (w, w, w)))
    lambda1 = np.reshape(lambda1, (w, w, w))
    lambda2 = np.reshape(lambda2, (w, w, w))
    lambda3 = np.reshape(lambda3, (w, w, w))

    choice = random.choice(np.arange(w*w*w), 10000)
    d_flat = density[:,:,:].flatten()[choice]
    l1_flat = lambda1[:,:,:].flatten()[choice]
    fa_flat = first_acc[:,:,:].flatten()[choice]
    ok = fa_flat < 2

    #fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    fig, ax = plt.subplots(1, 1) 
    ax.scatter(d_flat[ok], np.log10(fa_flat[ok]),
               s=15, c=l1_flat[ok])
    ax.set_xlabel(r"$\log_{10}(\rho_{\rm IC})$")
    ax.set_ylabel(r"$\log_{10}(a_{\rm acc})$")

    z = 100
    
    plt.figure()
    plt.imshow(np.log10(first_acc[:, :, z]), origin="lower",
               extent=[0, l, 0, l], cmap="afmhot")
    plt.title(r"$\log_{10}(a_{\rm acc})$")
    plt.xlabel(r"$u_x\,({\rm Mpc}/h)$")
    plt.ylabel(r"$u_y\,({\rm Mpc}/h)$")
    plt.colorbar()

    plt.figure()
    plt.imshow(density[:, :, z], origin="lower", extent=[0, l, 0, l],
               cmap="bwr", vmin=-0.35, vmax=0.35)
    plt.title(r"$\log_{10}(\rho_{\rm IC})$")
    plt.xlabel(r"$u_x\,({\rm Mpc}/h)$")
    plt.ylabel(r"$u_y\,({\rm Mpc}/h)$")
    plt.colorbar()

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    ax[0].set_title(r"$\lambda_{1,{\rm IC}}$")
    ax[0].imshow(lambda1[:, :, z], origin="lower", extent=[0, l, 0, l],
                 cmap="RdGy_r", vmin=-0.4, vmax=0.4)
    ax[0].set_xlabel(r"$u_x\,({\rm Mpc}/h)$")
    ax[0].set_ylabel(r"$u_y\,({\rm Mpc}/h)$")
    
    ax[1].set_title(r"$\lambda_{2,{\rm IC}}$")
    ax[1].imshow(lambda2[:, :, z], origin="lower", extent=[0, l, 0, l],
                 cmap="RdGy_r", vmin=-0.4, vmax=0.4)
    ax[1].set_xlabel(r"$u_x\,({\rm Mpc}/h)$")
    ax[2].set_title(r"$\lambda_{3,{\rm IC}}$")
    im = ax[2].imshow(lambda3[:, :, z], origin="lower", extent=[0, l, 0, l],
                      cmap="RdGy_r", vmin=-0.4, vmax=0.4)
    ax[2].set_xlabel(r"$u_x\,({\rm Mpc}/h)$")
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    ax[1].get_yaxis().set_visible(False)
    ax[2].get_yaxis().set_visible(False)    

    plt.show()
    
if __name__ == "__main__": main()
