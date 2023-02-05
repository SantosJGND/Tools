#!/bin/bash

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate ./env

DIR=/home/bioinf/Desktop/OTHER/KAUST/INVERSIONS/

for gff in $DIR/*.bed; do
    echo $gff
    python -u Distribution_simple.py --input $gff
    
done
