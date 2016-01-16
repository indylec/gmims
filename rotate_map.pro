pro rotate_map,inmap,outmap,nside,coords,SET_COORDS=set_coords

sizein=size(inmap)

pix_array=lindgen(sizein[1])

pix2vec_nest,nside,pix_array,invecs

if (keyword_set(set_coords) && (n_elements(coords) NE 0)) then begin outvecs=rotate_coord(invecs,inco=coords[1],outco=coords[0])
endif else begin outvecs=rotate_coord(invecs,inco='G',outco='C')
endelse

invecs=0

vec2pix_nest,nside,outvecs,new_pix_array

outvecs=0

outmap = inmap[new_pix_array]



end
   
