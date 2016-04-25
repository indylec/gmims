import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
#from setup_matplotlib import *
import sys

map_in=sys.argv[1]
vmin=float(sys.argv[2])
vmax=float(sys.argv[3])
use_planck_cmap=int(sys.argv[4])

map=hp.read_map(map_in,nest=True)


#out_filename = input_filename.split(".")[0]
if use_planck_cmap:
    ############### CMB colormap
    from matplotlib.colors import ListedColormap
    import numpy as np
    colombi1_cmap = ListedColormap(np.loadtxt("/home/leclercq/repos/gmims/Planck_Parchment_RGB.txt")/255.)
    colombi1_cmap.set_bad("gray") # color of missing pixels
    colombi1_cmap.set_under("white") # color of background, necessary if you want to use
    # this colormap directly with hp.mollview(m, cmap=colombi1_cmap)
    cmap = colombi1_cmap
    #out_filename += "_planck_cmap"
    hp.mollview(map,nest=True,min=vmin,max=vmax,cmap=cmap)

else:
    hp.mollview(map,nest=True)

plt.show()
