#!/bin/bash


SEQTK_DIR="/mnt/sdb/televir_work/televir_envs/remap/remap/bin/"

file=$1
prop_sample=$2
OUTDIR=$3

name=${file%".fastq.gz"}
name=`basename $name`
new_file=$name"_subset.fastq.gz"

echo $file $new_file $prop_sample >> $OUTDIR"sample_report.txt"

$SEQTK_DIR"seqtk" sample $file $prop_sample | bgzip -c > $OUTDIR$new_file
