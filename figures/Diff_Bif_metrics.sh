#!/bin/bash

set -e

threads=1
skip1=false
correction=benjamini
export TMPDIR=$HOME/tmp

Help(){
	echo "
This script runs the pipeline to get kmers differentially enriched in the case samples. The different steps are: preparing the fof (kmdiff_fof_prep.sh), finding the kmers of interest (Kmdiff.sh), finding their belonging (FindingKmero.sh) and getting a few metrics forthe species of interest (metrics.R). Its argumens are:
	-c --cases <PATH> directory containing cases (FC) sample files.
	-C --controls <PATH> directory containing controls (noise) sample files.
	-o --output <PATH> path to the output directory. It will be created if it exists.
	-g --genomes <PATH> path to the directory containg the genomes used for camisim read generation (1 per file)
	-s --species <PATH> path to a directory with the genomes of the species of interest.
	-k --correction <STRING> correction method to use with kmdiff. {benjamini}
	-t --threads <INT> number of threads to use. {1}
	--skip1) skips the kmdiff step if it is already done. Expects files to be in the output argument.
	-h --help displays this help message and exits.
	"
}

if [[ $# -eq 0 ]]
then
	Help
	exit
fi

while [[ $# -gt 0 ]]
do
	case $1 in
	-c | --cases) cases="$2"
	shift 2;;
	-C | --controls) controls="$2"
	shift 2;;
	-o | --output) output="$2"
	shift 2;;
	-g | --genomes) genomes="$2"
	shift 2;;
	-s | --species) species="$2"
	shift 2;;
	-t | --threads) threads="$2"
	shift 2;;
	-k | --correction) correction="$2"
	shift 2;;
	--skip1) skip1=true
	shift;;
	-h | --help) Help; exit;;
	-* | --*) unknown="$1"; echo -e "Unkown option: ${unknown}"; Help; exit;;
	*) shift;;
	esac
done


if [ -d "$output" ]
then
	if [ ! "$skip1" = true ]
	then
		echo -e "WARNING: ${output} directory already exists."
	fi
else
	mkdir "$output"
fi

if [ ! -d "$cases" ]
then
	echo -e "ERROR: ${cases} directory does not exist. Exiting."
fi

if [ ! -d "$controls" ]
then
	echo -e "ERROR: ${controls} directory does not exist. Exiting."
fi

if [ ! "$skip1" = true ]
then
	/mnt/projects_tn01/KMER_analysis/gitlab/kmdiff_fof_prep.sh --cases "$cases" --controls "$controls" --output "$output"/fof.txt
	/mnt/projects_tn01/KMER_analysis/gitlab/Kmdiff.sh --file "$output"/fof.txt --output "$output"/kmdiff --threads "$threads" --correction "$correction"
fi

#stops the pipeline if no kmer is found to be enriched in case
diff_kmers=$(cat "$output"/kmdiff/case_kmers.fasta | wc -l)
if [[ "$diff_kmers" -eq 0 ]]
then
	echo -e "STOP: kmdiff didn't find any kmer differentially enriched between case and control"
	exit 1
fi

/mnt/projects_tn01/KMER_analysis/gitlab/FindingKmero.sh --output "$output" --genomes "$genomes" --species "$species" --kmdiff "$output"/kmdiff/case_kmers.fasta --threads "$threads"

Rscript /mnt/projects_tn01/KMER_analysis/gitlab/metrics.R "$output" > "$output"/results_model.txt
mv "$output"/kmdiff/case_kmers.fasta "$output"/; mv "$output"/kmdiff/control_kmers.fasta "$output"/
