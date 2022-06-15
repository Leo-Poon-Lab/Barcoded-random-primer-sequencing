library(Biostrings)

dir.create("../results")
dir.create("../results/MEGAHIT/")
sample_name <- list.files("../data/", "final.contigs.fa$", recursive = T, full.names = T)
sample_name <- sample_name[!grepl("intermediate", sample_name)]

for(i in seq_along(sample_name)){
    sample_tmp <- sample_name[i]
	sample_i <- strsplit(sample_tmp, "_")[[1]]
    sample_i <- sample_i[length(sample_i)]
	sample_i <- strsplit(sample_i, "/", fixed = T)[[1]][1]
	print(sample_i)
	
	data <- readDNAStringSet(sample_tmp)
    multi <- sapply(names(data), function(x){
		tmp = strsplit(x, " len=")[[1]][1]
		tmp = strsplit(tmp, " multi=")[[1]][2]
		return(as.numeric(tmp))
	})
	if(length(data)> 100){
		data_out <- data[order(multi, decreasing = T)][1:100]
	} else {
		data_out <- data[order(multi, decreasing = T)]
	}

    writeXStringSet(data_out, paste0("../results/MEGAHIT/",sample_i, "_multi100.fasta"))
    
}



