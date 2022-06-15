#!/bin/sh

#de novo
cd ../data

for dir in `ls -d */`
do
    cd $dir
	for fwdread in `ls | grep '_1.fastq'`
	do
		sample=$(echo $fwdread | cut -d"_" -f 1)
		rwsread=$sample"_2.fastq.gz"
		# java -jar /home/hggu/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 48 -phred33 $fwdread $rwsread output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
		/home/alison/softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 $fwdread -2 $rwsread -t 10 -o "megahit_"$sample
	done
	cd ..
done
