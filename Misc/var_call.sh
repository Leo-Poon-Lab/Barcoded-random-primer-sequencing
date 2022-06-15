#!/bin/sh

# conda activate ngs
cd ../data

for dir in `ls -d */`
do
    cd $dir
	for fwdread in `ls | grep '_1.fastq.gz'`
	do
		sample=$(echo $fwdread | cut -d"_" -f 1)
		rwsread=$sample"_2.fastq.gz"

		# begin varcall
		## freebayes
		~/softwares/freebayes/scripts/freebayes-parallel <(~/softwares/freebayes/scripts/fasta_generate_regions.py reference.fasta.fai 650) 47 -f reference.fasta -F 0.01 $bamfile > ../../results/"freebayes_"$sample".vcf"				
		## VarDict
		bedtools makewindows -w 650 -b reference.bed > ref_slide_window.bed ##TODO
		cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 47 -k "~/softwares/VarDictJava/build/install/VarDict/bin/VarDict -G reference.fasta -f 0.01 -N samplename -th 48 -b $bamfile -R {} " | ~/softwares/VarDictJava/build/install/VarDict/bin/teststrandbias.R | ~/softwares/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N samplename -E -f 0.01 > ../../results/"vardict_"$sample".vcf"
		## lofreq
		/home/hggu/softwares/lofreq/src/lofreq/lofreq call-parallel --pp-threads 47 -f reference.fasta -o ../../results/"lofreq_"$sample".vcf" $bamfile

	done

done


