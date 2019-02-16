#!/bin/bash

name=$(basename $1 .fasta)
./split_genome.sh $1 $name

for f in ${name}/*.fasta; do
  echo sbatch --qos=normal --mem-per-cpu=4096 --output="slurm-${name}-$(basename $f .fasta)-%j.out" --wrap="ribosor.py $f"
done
