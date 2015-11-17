#!/bin/bash

# this runs snakemake on cluster, requesting 2GB ram 
# per job and running a maximum of 50 jobs at a time
# needs to be run on head node and should probably be 
# run in screen session so that can be left to
# run in background

snakemake --rerun-incomplete --cluster "qsub -l h_vmem=16g -V" --jobs 30
# snakemake --rerun-incomplete --drmaa " -l h_vmem=8g -V" --jobs 30
`
