import numpy as np
import healpy as hp
import sys


map_in=sys.argv[1]

map=hp.read_map(map_in,nest=True)

map[np.where(map!=hp.UNSEEN)]*=1.55

map_out=map_in.split(".")[0]+"_1.55.fits"

hp.write_map(map_out,map,nest=True)
