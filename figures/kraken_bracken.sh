#!/bin/bash

mkdir temp
/mnt/software/bioinfo/Kraken2-daydream/bin/kraken2 --db /mnt/software/dataset/db/kraken2/2022/ --threads 80 --output output_kraken --report temp/test_kraken_report.tsv  $(ls /mnt/loreal_tn04/projects/LOREAL_TN04/de_scratch_is6/AD_LO_shotgun_acne/preprocessing/good/SRR806970*-28046.good.fq)

exit
/mnt/software/bioinfo/Kraken2-2.1.2/bin/kraken2 --db /mnt/software/dataset/db/kraken2/2022/ --threads 80 --output temp --report temp/paired_kraken_report.tsv  /mnt/loreal_tn04/projects/LOREAL_TN04/de_scratch_is6/AD_LO_shotgun_acne/preprocessing/good/*_R1* /mnt/loreal_tn04/projects/LOREAL_TN04/de_scratch_is6/AD_LO_shotgun_acne/preprocessing/good/*_R2*


/mnt/software/bioinfo/Bracken-2.6/Bracken/bracken -r 150 -d /mnt/software/dataset/db/kraken2/2022/ \
	-i ./temp/kraken_report.tsv \
	-o bracken_output \
	-l S \
	-w bracken_S_report.tsv \
	-t 10 \

/mnt/software/bioinfo/Bracken-2.6/Bracken/bracken -r 150 -d /mnt/software/dataset/db/kraken2/2022/ \
        -i ./temp/paired_kraken_report.tsv \
        -o paired_bracken_output \
        -l S \
        -w paired_bracken_S_report.tsv \
        -t 10 \

ls *bracken_S_report.tsv > input_script_alban.txt
python3 ./resume_taxonomy_bracken_v2.py \
        -l input_script_alban.txt \
        -o temp/output.otu \
        -S temp/output.S \
        -G temp/output.G \
        -C temp/output.C \
        -F temp/output.F \
        -O temp/output.O \
        -P temp/output.P \
        -K temp/output.K \
        -na /mnt/software/dataset/db/kraken2/2022/taxonomy/names.dmp \
        -no /mnt/software/dataset/db/kraken2/2022/taxonomy/nodes.dmp \
        -me /mnt/software/dataset/db/kraken2/2022/taxonomy/merged.dmp \

mv temp/output.S ./acne_kraken_2022.tsv
