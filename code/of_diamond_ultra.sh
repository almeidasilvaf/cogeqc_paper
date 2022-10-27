#!/bin/bash

# Load modules
module load OrthoFinder diamond mcl fastme

# Define paths
workdir='/home/faalm/projects/cogeqc_benchmark'
outdir='/home/faalm/projects/cogeqc_benchmark/products/result_files'

# Run OrthoFinder
orthofinder -f "$workdir/data" -S diamond_ultra_sens -I 1.5 -t 8 -o "$outdir/ultra_1_5" -og
orthofinder -f "$workdir/data" -S diamond_ultra_sens -I 2 -t 8 -o "$outdir/ultra_2" -og
orthofinder -f "$workdir/data" -S diamond_ultra_sens -I 3 -t 8 -o "$outdir/ultra_3" -og
orthofinder -f "$workdir/data" -S diamond_ultra_sens -I 1 -t 8 -o "$outdir/ultra_1" -og
