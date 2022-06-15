#!/bin/sh

cd ../data


for dir in `ls -d */`
do
    cd $dir
	for fwdread in `ls | grep '_1.fastq.gz'`
	do
		echo $fwdread
		zcat $fwdread | paste - - - - | wc -l 
	done
	cd ..
done

