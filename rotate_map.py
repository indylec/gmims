import numpy as np
import healpy as hp
import sys
import matplotlib.pyplot as plt

infile=sys.argv[1]
outfile=sys.argv[2]

inmap=hp.read_map(infile,nest=True)

size=inmap.shape[0]

pix_array=np.arange(size)

invecs_x, invecs_y, invecs_z =hp.pix2vec(512,pix_array,nest=True)

r = hp.rotator.Rotator(coord=['G','C'])

outvecs_x, outvecs_y, outvecs_z = r(invecs_x, invecs_y, invecs_z)

new_pix_array=hp.vec2pix(512,outvecs_x,outvecs_y,outvecs_z,nest=True)

outmap = inmap[new_pix_array]

hp.write_map(outfile,outmap,nest=True)


#fig1=plt.figure(figsize=(8,5),dpi=150)
#hp.mollview(outmap, nest=True,fig=fig1.number, coord='G',min=0.,max=1.)
#hp.graticule(dpar=15., dmer=40., coord=['C'])
#plt.savefig(outfile,dpi=150, bbox_inches="tight")
#plt.show()
