# Barcoded-random-primer-sequencing

---

- `./reads_splits/`: FASTQ paired-end raw reads were divided into groups according to their barcode sequence if available. The reads division/assignment was done with in-house scripts. Essentially, we extracted the heading and tailing five nucleotide bases (barcode region) from each read, and afterwards split original sample-specific fastq reads files to different sample-barcode-specific fastq reads files.

- `./denovo_assembly/`: For each sample-barcode-specific fastq reads files, reference genome-free de novo assembly was carried out using MEGAHIT (v1.2.9) to the raw reads to generate consensus. The most relevant reference genome was determined by comparing the assembled consensus against global database using NCBI BLASTN programme. 

- `./reference_based_alignment/`: Basing on the identified reference genome, reference-based alignment for raw reads were performed using BWA-MEM2 (v2.0pre2). Consensus sequence was obtained using Samtools mpileup (v1.12) with iVar consensus (v1.3.1), and the depths of each nucleotide position was estimated using mpileup2readcounts (https://github.com/gatoravi/mpileup2readcounts). 

- `./Misc/`: Coverage plots were generated according to the depths by ggplot2 in R. 

---