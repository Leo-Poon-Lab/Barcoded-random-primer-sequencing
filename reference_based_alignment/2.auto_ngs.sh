#!/bin/sh

#conda activate
#cp -a ../scripts/ORF_WT.csv ../data/
cp -a ../scripts/reference.fasta ../data/

cd ../data
mkdir ../results
mkdir ../FastQC
mkdir ../results/MEGAHIT
mkdir ../results/ivar_consensus

for dir in `ls -d */`
do
    cd $dir
	for fwdread in `ls | grep '_1.fastq.gz'`
	do
		cp ../reference.fasta ./
		~/softwares/samtools-1.12/bin/samtools faidx reference.fasta
		~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index reference.fasta
		#~/softwares/bwa/bwa index reference.fasta
		java -jar ~/softwares/picard.jar CreateSequenceDictionary -R reference.fasta -O reference.dict
		~/softwares/bioawk/bioawk -c fastx '{print $name"\t1\t"length($seq)"\t"$name}' reference.fasta > reference.bed

		sample=$(echo $fwdread | cut -d"_" -f 1)
		rwsread=$sample"_2.fastq.gz"

		# # de novo assembly	
		# ~/softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -t 46 -1 $fwdread -2 $rwsread -o ../../results/MEGAHIT/$sample

		# resequencing
		~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 46 reference.fasta $fwdread $rwsread | ~/softwares/samtools-1.12/bin/samtools view -b --threads 46 - | ~/softwares/samtools-1.12/bin/samtools sort -o aln.sorted.bam --threads 46 -m 4G - 
		mv aln.sorted.bam $sample"_sorted.bam"
		~/softwares/samtools-1.12/bin/samtools index $sample"_sorted.bam"

		~/softwares/samtools-1.12/bin/samtools mpileup -aa -d 1000000000 -f reference.fasta $sample"_sorted.bam" > $sample"_mpileup.txt"			
		cat $sample"_mpileup.txt" | ivar consensus -p $sample"_consensus" -n N -m 100 
		
		mv $sample"_consensus.fa" ../../results/ivar_consensus/$sample"_consensus.fa"
		
		## readcount
		cat $sample"_mpileup.txt" | ~/softwares/mpileup2readcounts/build/mpileup2readcounts > ../../FastQC/$sample"_bam.readcount.mpileup.txt"	
		
		## FastQC
		# fastqc $sample"_sorted.bam" -t 48 -o ../../FastQC/
		rm $sample"_mpileup.txt"
		rm $sample"_mpileup_ori.txt"
		rm *.bam
		rm *.bai
		rm $sample"_consensus"*			
	done
	cd ..
done


