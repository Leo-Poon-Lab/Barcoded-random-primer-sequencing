#!/usr/bin/Rscript
library(Biostrings)
library(parallel)
library(tidyverse)
# library(future.apply)
# plan(multiprocess, gc = TRUE)
# options(future.globals.maxSize = 2 * 1000 * 1024^2) # 2GB
# ulimit::memory_limit(2000)

# sample = "1-210-2233"
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

wd_ori <- getwd()
dir = paste0("../data/", sample, "/")
print(dir)
# dir = dirs[1]
setwd("../../scripts")
barcodes <- readr::read_tsv("./barcodes.tsv", col_names = F)
system(paste("/bin/bash -c", shQuote(paste0("cp -a ./to_barcode.sh ", dir))))
system(paste("/bin/bash -c", shQuote(paste0("cp -a ./split_by_label.sh ", dir))))
setwd(dir)

# step 1
print("running step 1: extracting barcodes")

read_1_path <- list.files(".", "_1.fastq.gz", full.names = T)
read_2_path <- gsub("_1.fastq.gz", "_2.fastq.gz", read_1_path)
read_1_path <- gsub("//", "/", read_1_path)
read_2_path <- gsub("//", "/", read_2_path)
sample_name <- strsplit(read_1_path, "/", fixed=T)[[1]]
sample_name <- sample_name[length(sample_name)]
sample_name <- strsplit(sample_name, "_", fixed=T)[[1]][1]
print(sample_name)

command <- paste0("paste <(zcat ", read_1_path, "|paste - - - -) <(zcat ", read_2_path, "| paste - - - - )", " | ", "./to_barcode.sh")
system(paste("/bin/bash -c", shQuote(command)))

# step 2
print("running step 2: labeling barcode")
read_1 <- readLines("./tmp_1.seq")
read_2 <- readLines("./tmp_2.seq")

read_1 <- DNAStringSet(read_1)	
read_2 <- DNAStringSet(read_2)	
read_1_com <- complement(read_1)
read_2_com <- complement(read_2)

read_1 <- as.character(read_1)
read_2 <- as.character(read_2)
read_1_com <- as.character(read_1_com)
read_2_com <- as.character(read_2_com)

total_i <- length(read_1)
#format(object.size(read_1), units = "Mb")

fastq_label <- mclapply(seq_along(read_1), function(i){
	# print(round(i/total_i*100, 2))
	rst_all <- c(read_1[i], read_2[i], read_1_com[i], read_2_com[i])
	check <- rst_all %in% barcodes$X2
	if(any(check)){
		b_tmp <- barcodes$X1[barcodes$X2 == rst_all[check][1]]
		return(b_tmp)
	} else {
		return("NA")
	}
}, mc.cores = 24)

#gc()
fastq_label <- unlist(fastq_label)
writeLines(fastq_label, "./fastq_label")

count_t <- table(fastq_label)
hous1 <- head(sort(table(read_1[fastq_label=="NA"]), decreasing = T), 20)
hous2 <- head(sort(table(read_1[fastq_label=="NA"]), decreasing = T), 20)
write_csv(tibble(seq = names(c(hous1, hous2)), num = c(hous1, hous2)), paste0("./head_of_uassigned_seqs_", sample, ".csv"))
write_csv(tibble(seq = names(count_t), num = count_t), paste0("./barcode_stat_", sample, ".csv"))


# step 3
print("running step 3: spliting samples")
command <- paste0("paste <(zcat ", read_1_path, "|paste - - - -) <(zcat ", read_2_path, "| paste - - - - ) <(cat fastq_label) ", "| ./split_by_label.sh")
system(paste("/bin/bash -c", shQuote(command)))

system("mkdir archived")
system(paste0("mv -f ", sample_name, "* ",  "./archived"))
file_tmp <- list.files(pattern = "fastq.gz")
file.rename(file_tmp, paste0(sample_name, "-", file_tmp))

system(paste("/bin/bash -c", shQuote("gzip -f ./fastq_label ./fastq_label")))
file.rename("./fastq_label.gz", paste0(sample_name, "-", "fastq_label.gz"))
system(paste("/bin/bash -c", shQuote("gzip -f ./tmp_1.seq ./tmp_1.seq")))
file.rename("./tmp_1.seq.gz", paste0(sample_name, "-", "tmp_1.seq.gz"))
system(paste("/bin/bash -c", shQuote("gzip -f ./tmp_2.seq ./tmp_2.seq")))
file.rename("./tmp_2.seq.gz", paste0(sample_name, "-", "tmp_2.seq.gz"))
setwd(wd_ori)	

