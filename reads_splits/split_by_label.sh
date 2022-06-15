#!/bin/sh

awk -v FS="\t" -v OFS="\n" '
{
	label = $9; 
	{print $1,$2,$3,$4 | "gzip > "label"_1.fastq.gz"; print $5,$6,$7,$8 | "gzip > "label"_2.fastq.gz"}
}
' - 