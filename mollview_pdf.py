import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from setup_matplotlib import *
import sys

map_in=sys.argv[1]
vmin=float(sys.argv[2])
vmax=float(sys.argv[3])


map=hp.read_map(map_in,nest=True)

map=hp.reorder(map,n2r=True)

map[np.where(map==hp.UNSEEN)]=np.nan
map[np.where(map==0.0)]=np.nan

mask=np.ones(hp.nside2npix(512))
mask[np.where(~np.isnan(map))]=0

map_pdf=map_in.rsplit('.',1)[0]+'.pdf'

nside=512

xsize = 2000
ysize = xsize/2.

unit = r"$\mathrm{K}$"

theta = np.linspace(np.pi, 0, ysize)
phi   = np.linspace(-np.pi, np.pi, xsize)

#r = hp.rotator.Rotator(coord=['G','C'])

longitude = np.radians(np.linspace(-180, 180, xsize))
latitude = np.radians(np.linspace(-90, 90, ysize))

# project the map to a rectangular matrix xsize x ysize
PHI, THETA = np.meshgrid(phi, theta)

#GALTHETA,GALPHI=r(THETA,PHI)

#GALPHI, GALTHETA=  np.meshgrid(galphi, galtheta)

grid_pix = hp.ang2pix(nside, THETA, PHI)

#grid_pix = hp.ang2pix(nside, GALTHETA, GALPHI)

grid_mask=mask[grid_pix]
grid_map=np.ma.MaskedArray(map[grid_pix], grid_mask)
#grid_map = map[grid_pix]

from matplotlib.projections.geo import GeoAxes

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)

width=20.

fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)/(3./2.)))
    # matplotlib is doing the mollveide projection
ax = fig.add_subplot(111,projection='mollweide')

from matplotlib import cm
myafm=cm.get_cmap()
#myafm=cm.afmhot
myafm.set_bad(color='grey')
#myafm.set_under("white")
#myafm.set_over("white")


# rasterized makes the map bitmap while the labels remain vectorial
# flip longitude to the astro convention
image = plt.pcolormesh(longitude[::-1], latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=True,cmap=myafm)

# graticule
ax.set_longitude_grid(60)
ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
if width < 10:
    ax.set_latitude_grid(45)
    ax.set_longitude_grid_ends(90)


# colorbar
cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05, ticks=[vmin, vmax])
cb.ax.xaxis.set_label_text(unit)
cb.ax.xaxis.labelpad = -8
# workaround for issue with viewers, see colorbar docstring
cb.solids.set_edgecolor("face")

ax.tick_params(axis='x', labelsize=10,color='w',labelcolor='w')
ax.tick_params(axis='y', labelsize=10)

ax.grid(color='w')

# remove tick labels
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticklabels([])
# remove grid
#ax.xaxis.set_ticks([])
#ax.yaxis.set_ticks([])

# remove white space around figure
#spacing = 0.05
#plt.subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)

plt.grid(True)

#fig1=plt.figure(figsize=(8,5),dpi=150)
#hp.mollview(map,nest=True,fig=fig1.number,coord='C' ,xsize=8*150,min=mymin,max=mymax,title='',unit='Brightness (K)')
#hp.graticule(dpar=15., dmer=10., coord=['G'])

#plt.show()
plt.savefig(map_pdf,bbox_inches="tight")#,dpi=150)
