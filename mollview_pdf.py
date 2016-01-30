import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys

map_in=sys.argv[1]
mymin=float(sys.argv[2])
mymax=float(sys.argv[3])

map=hp.read_map(map_in,nest=True)

map_pdf=map_in.rsplit('.',1)[0]+'.pdf'

fig1=plt.figure(figsize=(8,5),dpi=150)
hp.mollview(map,nest=True,fig=fig1.number,coord='G', xsize=8*150,min=mymin,max=mymax,title=map_in)
hp.graticule(dpar=15., dmer=40., coord=['C'],local=True)
plt.show()
#plt.savefig(map_pdf,dpi=150, bbox_inches="tight")
