#!/bin/sh

awk -v FS="\t" -v OFS="\n" '{head_1 = substr($2,1,5); head_2 = substr($6,1,5); {print head_1>"tmp_1.seq"; print head_2>"tmp_2.seq"}}' - 