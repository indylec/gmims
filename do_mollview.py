#!/bin/bash

cd ~/gmims/ttplots

for f in *spec_index*; do
    echo "$f";
    python ~/repos/gmims/mollview_pdf.py "$f" -6 -1
done
