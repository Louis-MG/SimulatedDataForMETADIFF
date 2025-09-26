#!/bin/bash

set -e

Help() {
	echo "
This script uses Bifrost to get the presence/absence patterns of Kmdiff's kmers in the genomes used for reads generation. It first finds kmers unique to each genomne of interest relatively to the other genomes. It then indexes the genomes and unique kmers, and finally finds the patterns. Its arguments are :
	-g --genomes <PATH> path to the directory containing the genomes used for read generation (1 per file).
	-s --species <PATH> path to a directory with the genomes of the species of interest.
	-k --kmdiff <PATH> path to the kmdiff output.
	-o --output <PATH> path to the output directory. Will be created if it does not exist.
	-t --threads <INT> number of threads to use. [1]
	-h --help displays this help message and exits.

The main output is:
	- a query_results.tsv file, tab-separated, with kmers as lines and genomes as columns. Each kmer is either present (1) or absent (0) for each genome.
	"
}

threads=1

if [ $# -eq 0 ]
then
	Help
	exit
fi

while [[ $# -ge 1 ]]
do
	case $1 in
	-g | --genomes) genomes="$2"
	shift 2;;
	-s | --species) species="$2"
	shift 2;;
	-k | --kmdiff) kmdiff="$2"
	shift 2;;
	-o | --output) output="$2"
	shift 2;;
	-t | --threads) threads="$2"
	shift 2;;
	-h | --help) Help; exit;;
	-* | --*) unknown="$1"; echo -e "Unkown option: ${unknown}"; Help; exit;;
	*) shift;;
	esac
done


if [ ! -f "$kmdiff" ]
then
	echo -e "ERROR: $kmdiff does not exist. Exiting."
	exit 1
fi

if [ ! -d "$genomes" ]
then
	echo -e "ERROR: $genomes does not exist. Exiting."
	exit 1
fi

if [ ! -d "$species" ]
then
        echo -e "ERROR: $species does not exist. Exiting."
	exit 1
fi

#######################################
#
#	ETAPE 1 : KMER UNIQUES
#
#######################################

echo -e "Concatening genome ..."
mkdir $output/tmp_genomes
echo $genomes/* | xargs cat > $output/tmp_genomes/all.fasta

echo -ne "\033[2K \rObtaining unique kmers ..."
#on obtient les kmers uniques a chaque genome
parallel_uniques() {
	local sp_target=$1
	local file_sp_target=$( cut -d "/" -f6 <<< $1)
	local name_sp_target=${file_sp_target//.fna/}
	local output=$3
	cat $output/tmp_genomes/all.fasta > $output/tmp_genomes/all_except_$name_sp_target.fasta
	cat $(find $2/* ! -name $file_sp_target) >> $output/tmp_genomes/all_except_$name_sp_target.fasta
	/mnt/projects_tn01/KMER_analysis/tools/UniqueKMER/uniquekmer -f $sp_target -k 31 -r $output/tmp_genomes/all_except_${name_sp_target}.fasta -e 0  -o $output/tmp_genomes/unique_kmers_${name_sp_target}
	mv $output/tmp_genomes/unique_kmers_${name_sp_target}/kmercollection.fasta  $output/unique_kmers_${name_sp_target}.fasta
	awk -i inplace '{if (/>/) {line=$0; sum=0} else {sum+=1; KMER=$0; print line ":kmer:" sum "\n" KMER} }' $output/unique_kmers_${name_sp_target}.fasta
}
export -f parallel_uniques
parallel --tmpdir ~/tmp_parallel --jobs $threads parallel_uniques ::: $(readlink -f $species/*) ::: $species ::: $output


#######################################
#
#       ETAPE 2 : BIFROST
#
#######################################

echo -ne "\033[2K \rPreparing bifrost job..."
readlink -f $genomes/* > $output/tmp_genomes/fof_bifrost_build.txt
readlink -f $species/* >> $output/tmp_genomes/fof_bifrost_build.txt
readlink -f $output/unique_kmers_* >> $output/tmp_genomes/fof_bifrost_build.txt

echo -ne "\033[2K \rBifrost running for kmers enriched in cases, on reference genomes ....\n"
#on construit un BDG colore avec les kmers uniques a chaque genome et on query les kmers de kmdiff dessus
/mnt/projects_tn01/KMER_analysis/tools/bifrost/build/src/Bifrost build -t $threads --colors --input-ref-file $output/tmp_genomes/fof_bifrost_build.txt -o $output/tmp_genomes/bifrost_colored_dbg
/mnt/projects_tn01/KMER_analysis/tools/bifrost/build/src/Bifrost query -t $threads -e 1 --input-graph-file $output/tmp_genomes/bifrost_colored_dbg.gfa.gz --input-query-file $kmdiff --input-color-file $output/tmp_genomes/bifrost_colored_dbg.color.bfg -o $output/result_genomes_bifrost_query

sed -i 's/\/mnt\/scratch\/LM\/genomes_minus_4\///g' $output/result_genomes_bifrost_query.tsv
sed -i 's/\/mnt\/scratch\/LM\/bifrost\///g' $output/result_genomes_bifrost_query.tsv
sed -i 's/\/mnt\/scratch\/LM\/4genomes\///g' $output/result_genomes_bifrost_query.tsv
sed -i 's/.fna//g' $output/result_genomes_bifrost_query.tsv
sed -i 's/.fasta//g' $output/result_genomes_bifrost_query.tsv
string_output=${output////\\/}
sed -i -e "s/$string_output//g; s/\///g" $output/result_genomes_bifrost_query.tsv

echo -ne "\033[2K \rBifrost running enriched in controls, on reference genomes ....\n"

kmdiff=${kmdiff//case/control}
/mnt/projects_tn01/KMER_analysis/tools/bifrost/build/src/Bifrost query -t $threads -e 1 --input-graph-file $output/tmp_genomes/bifrost_colored_dbg.gfa.gz --input-query-file $kmdiff --input-color-file $output/tmp_genomes/bifrost_colored_dbg.color.bfg -o $output/controls_result_genomes_bifrost_query

rm -dr $output/tmp_genomes
