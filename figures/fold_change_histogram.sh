#!/bin/bash

set -e

Help() {
echo "
This script samples kmers from a result and calculates fold changes for each. Then displays the values in a histogram to help the analysis. Its arguments are :
	-i --input <PATH> path to the input case_kmers.fasta (output of kmdiff).
	-o --output <PATH> path to the output directory.
	-c --case_or_control self-explanatory, needed for the figure [case,control].
	-h --help displays this help message and exits.
"
}

if [ $# -eq 0 ]
then
	Help
	exit 0
fi

while [ $# -ge 1 ]
do
	case $1 in
	-i | --input) input="$2"
	shift 2;;
	-o | --output) output="$2"
	shift 2;;
	-c | --case_or_control) case_control="$2"
	shift 2;;
	-h | --help) Help; exit;;
	-* | --*) unknown="$1"; echo -e "ERROR: unkown option $unknown"; Help; exit;;
	*) shift;;
	esac
done

grep -F '>' "$input" > "$output"/kmers_FC.tsv
#garde uniquement les chiffres
sed -i 's/.*control=//g;s/_case=/\t/g' "$output"/kmers_FC.tsv

Rscript /mnt/projects_tn01/KMER_analysis/gitlab/fold_figure.R "$output" "$case_control"
#rm "$output"/kmers_FC.tsv
