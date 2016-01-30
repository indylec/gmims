#!/bin/bash

python ~/repos/gmims/coeff_all.py haslam_adjusted.fits g0g.fits haslam_adjusted g0g
python ~/repos/gmims/coeff_all.py haslam_adjusted.fits g1g.fits haslam_adjusted g1g
python ~/repos/gmims/coeff_all.py haslam_adjusted.fits g7g.fits haslam_adjusted g7g
python ~/repos/gmims/coeff_all.py haslam_adjusted.fits g8g.fits haslam_adjusted g8g
python ~/repos/gmims/coeff_all.py haslam_adjusted.fits g0.fits haslam_adjusted g0

python ~/repos/gmims/coeff_all.py haslam.fits g0g.fits haslam g0g
python ~/repos/gmims/coeff_all.py haslam.fits g0.fits haslam g0

python ~/repos/gmims/coeff_all.py g0g.fits haslam_adjusted.fits g0g haslam_adjusted
python ~/repos/gmims/coeff_all.py g0.fits haslam_adjusted.fits g0 haslam_adjusted


python ~/repos/gmims/coeff_all.py stockert_adjusted.fits g0g.fits stockert_adjusted g0g
python ~/repos/gmims/coeff_all.py stockert_adjusted.fits g1g.fits stockert_adjusted g1g
python ~/repos/gmims/coeff_all.py stockert_adjusted.fits g7g.fits stockert_adjusted g7g
python ~/repos/gmims/coeff_all.py stockert_adjusted.fits g8g.fits stockert_adjusted g8g
python ~/repos/gmims/coeff_all.py stockert_adjusted.fits g0.fits stockert_adjusted g0

python ~/repos/gmims/coeff_all.py stockert.fits g0g.fits stockert g0g
python ~/repos/gmims/coeff_all.py stockert.fits g0.fits stockert g0

python ~/repos/gmims/coeff_all.py g0g.fits stockert_adjusted.fits g0g stockert_adjusted
python ~/repos/gmims/coeff_all.py g0.fits stockert_adjusted.fits g0 stockert_adjusted

python ~/repos/gmims/coeff_all.py g0g.fits g8g.fits g0g g8g
python ~/repos/gmims/coeff_all.py g8g.fits g0g.fits g8g g0g
python ~/repos/gmims/coeff_all.py g1g.fits g8g.fits g1g g8g
python ~/repos/gmims/coeff_all.py g1g.fits g7g.fits g1g g7g
python ~/repos/gmims/coeff_all.py g0g.fits g7g.fits g0g g7g

python ~/repos/gmims/coeff_all.py g0.fits g8.fits g0 g8
python ~/repos/gmims/coeff_all.py g8.fits g0.fits g8 g0