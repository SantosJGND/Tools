#!/bin/bash

### variables 
# MAX_SIZE: get files under that size
# seqtk script: script that takes  fastq filepath, proportion, outdir. Samples from file to gzipped subset in subdir.
# FILE_DIR: directory with fastq files.

MAX_SIZE=400000
SEQTK_SCRIPT="seqtk_sample_prop.sh"
FILE_DIR=`pwd`"/"

### output directory
OUTDIR=$FILE_DIR"compressed/"

mkdir -p $OUTDIR

######
######

for file_path in $FILE_DIR*gz; do

echo $file_path

file_size=$(du -s $file_path | awk '{print $1}')
prop=`bc -l <<< $MAX_SIZE' / '$file_size`

echo $prop $file_size

sh $SEQTK_SCRIPT $file_path $prop $OUTDIR

done
