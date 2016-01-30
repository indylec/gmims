#!/bin/bash

python ~/repos/gmims/apply_mask.py g1.fits eq2gal_kq85.fits
python ~/repos/gmims/apply_mask.py g7.fits eq2gal_kq85.fits
python ~/repos/gmims/apply_mask.py g8.fits eq2gal_kq85.fits

python ~/repos/gmims/apply_mask.py g0g.fits gal_kq85.fits
python ~/repos/gmims/apply_mask.py g1g.fits gal_kq85.fits
python ~/repos/gmims/apply_mask.py g7g.fits gal_kq85.fits
python ~/repos/gmims/apply_mask.py g8g.fits gal_kq85.fits