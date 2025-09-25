# Simulated Data For META-DIFF
Repository to replicate our simulation used in the paper METADIFF

# Structure

This repository is has the following structure:
 - functionnal analysis (10 paired-end samples per condition, 3 conditions, 5 replicates)
 - taxonomic analysis (10 paired-end samples per condition, 12 conditions, 5 replicates)
 - global files:
     - genomes used
     - global config

# Tools and scripts

The abundances were generated using the Rscript `traitement_acne_2022.R` .
The abundances and config files (with seeds) were used with [CAMISIM V1.3](https://github.com/CAMI-challenge/CAMISIM) using the following command:

```bash
python3 metagenomesimulation.py {config_file}
```
The yaml to reproduce the camisim environment is `camisim.yml`. The github version of CAMISIM was used within the environment.

NOTE: the functional simulation uses abundances from the taxonomic simulation.
