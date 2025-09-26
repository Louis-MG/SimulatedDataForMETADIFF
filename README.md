# Simulated Data For META-DIFF

Repository to replicate our simulation used in the paper METADIFF

# Structure

This repository is has the following structure: 
 - functional analysis (10 paired-end samples per condition, 3 conditions, 5 replicates) 
    - configs files 
    - modified reference genomes 
    - id to genome file with modified reference genomes 
 - taxonomic analysis (10 paired-end samples per condition, 12 conditions, 5 replicates) 
    - configs files 
    - reference genomes 
    - id to genome file 
    - abundance profiles
 - figures:
    - tables for taxonomy and functional
    - script to generate the figures from the tables
    - scripts (R'n'bash) to generate intermediate files and tables
 - global files: 
    - metadata 
    - scripts to replicate abundance profiles (`traitement_acne_2022.R`)
    - env yml for replication

# Tools and scripts

The abundances were generated using the Rscript `traitement_acne_2022.R`
. The abundances and config files (with seeds) were used with [CAMISIM
V1.3](https://github.com/CAMI-challenge/CAMISIM) using the following
command:

``` bash
#create a mamba env using the camisim.yml file
#activate env
#start the simulation
python3 metagenomesimulation.py {config_file}
```
NOTE: the functional simulation uses abundance profiles from the
taxonomic simulation.

To get insight into the -R and bash- scripts, use a `--help` argument.
