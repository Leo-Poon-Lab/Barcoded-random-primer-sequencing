#!/bin/sh

cp -a barcodes.tsv ../data/
cd ../data

for fwdread in `ls | grep '_1.fastq.gz'`
do
	sample=$(echo $fwdread | cut -d"_" -f 1)
	rwsread=$sample"_2.fastq.gz"
	mkdir $sample
	mv $fwdread $sample
	mv $rwsread $sample
done

for dir in `ls -d */`
do
    cd $dir
	fwdread=`ls | grep '_1.fastq.gz'`
	sample=$(echo $fwdread | cut -d"_" -f 1)
	Rscript ../../scripts/0.barcode_split.R $sample
	cd ..
done