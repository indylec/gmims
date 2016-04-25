import healpy as hp 
import numpy as np
import sys

inmap=sys.argv[1]
out=sys.argv[2]
res_out=int(sys.argv[3])

map_in=hp.read_map(inmap, nest=True)

map_out=hp.ud_grade(map_in,res_out,order_in='NESTED',order_out='NESTED',pess=True)

hp.write_map(out,map_out,nest=True)
