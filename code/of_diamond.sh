#!/bin/bash

# Load modules
module load OrthoFinder diamond mcl fastme

# Define paths
workdir='/home/faalm/projects/cogeqc_benchmark'
outdir='/home/faalm/projects/cogeqc_benchmark/products/result_files'

# Run OrthoFinder
orthofinder -f "$workdir/data" -S diamond -I 1.5 -o "$outdir/default_1_5" -og
orthofinder -f "$workdir/data" -S diamond -I 2 -o "$outdir/default_2" -og
orthofinder -f "$workdir/data" -S diamond -I 3 -o "$outdir/default_3" -og
orthofinder -f "$workdir/data" -S diamond -I 1 -o "$outdir/default_1" -og
