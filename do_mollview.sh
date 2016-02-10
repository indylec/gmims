#!/bin/bash

cd ~/gmims/ttplots

for f in g*g*spec_index*.fits; do
    echo "$f";
    python ~/repos/gmims/mollview_pdf.py "$f" -4 -1
done
