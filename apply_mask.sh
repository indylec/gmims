#!/bin/bash

python ~/repos/gmims/apply_mask.py TP_gal_bin0_nest_sm_512.fits 0 gal_kq85.fits

python ~/repos/gmims/apply_mask.py TP_gal_bin1_nest_sm_512.fits 1 gal_kq85.fits

python ~/repos/gmims/apply_mask.py TP_gal_bin4_nest_sm_512.fits 4 gal_kq85.fits

python ~/repos/gmims/apply_mask.py TP_gal_bin7_nest_sm_512.fits 7 gal_kq85.fits

python ~/repos/gmims/apply_mask.py TP_gal_bin8_nest_sm_512.fits 8 gal_kq85.fits