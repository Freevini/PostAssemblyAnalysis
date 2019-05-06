#!/usr/bin/env Rscript

##Usage example
#./PostAssemblyAnalysis.R -a data/testGenome.fasta -F data/testReads_R1.fastq -R data/testReads_R2.fastq -o testOut/
#install.packages("optparse")
library("optparse")

option_list = list(
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Input Assembly", metavar="fasta"),
  
  make_option(c("-i", "--index"), type="character", default=NULL, 
              help="Is an Index needed Y/N [default Y]", metavar="Y"),
  
  make_option(c("-F", "--forward"), type="character", default="NA", 
              help="Illumina Forward raw read in fastq.gz", metavar="fastq"),
  
  make_option(c("-R", "--reverse"), type="character", default="NA", 
              help="Illumina Reverse raw read in fastq.gz", metavar="fastq"),
  
  make_option(c("-l", "--LongRead"), type="character", default="NA", 
              help="Long reads raw in fastq.gz", metavar="fastq"),
  
  make_option(c("-o", "--output"), type="character", default="Stdout", 
              help="Output directory for all files created", metavar="PATH"),
  
  make_option(c("-t", "--threads"), type="numeric", default="5", 
              help="number of threads", metavar="threads"),
  
  make_option(c("-db", "--BlastDB"), type="character", default="Stdout", 
              help="Blast database used for blast search", metavar="blast")
  
  
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$assembly)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat(paste(c("Genome Assembly:",opt$assembly), collapse='\t'), '\n')
cat(paste(c("Index: ",opt$index), collapse='\t'), '\n')
cat(paste(c("forward Read:",opt$forward), collapse='\t'), '\n')
cat(paste(c("reverse Read:",opt$reverse), collapse='\t'), '\n')
cat(paste(c("Long Reads:",opt$LongRead), collapse='\t'), '\n')
cat(paste(c("output directory:",opt$output), collapse='\t'), '\n')
cat(paste(c("number of threads:",opt$threads), collapse='\t'), '\n')
# 
# if (opt$threads=="NA"){opt$threads <- 5}
# cat(opt$threads)
cat("##--------------------------------",'\n')   
cat("##01 create ouput directories",'\n')
cat("##--------------------------------",'\n') 

cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')

system(paste0("mkdir -p ",opt$output,"/RawReadsMerged"))

cat(paste(c("End:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')

cat("##--------------------------------",'\n')   
cat("##01 rawRead mapping minimap2 to assembly",'\n')
cat("##--------------------------------",'\n') 
if (opt$LongRead!="NA") {
  system(paste0("mkdir -p ",opt$output,"/rawReads2assembly_minimap2/bamqc"))
  
  cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
  
  #system(paste0("cat ",opt$forward, " ",opt$forward," > ", opt$output,"/RawReadsMerged/rawReadsMerged.fastq"))
  
 
  #system(paste0("minimap2 -t ",opt$threads," -ax map-ont ",opt$assembly," ", opt$LongRead, " > ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam"))
  #system(paste0("samtools view -bS " ,opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam | samtools sort - ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted "))
  

 system(paste0("minimap2 -t ",opt$threads," -ax map-ont ",opt$assembly," ", opt$LongRead, " | samtools sort -@8 -O BAM -o ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted.bam - "))



  ##Unmapped
  system(paste0("samtools view -b -f 4 ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted.bam > ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted_unmapped.bam" ))
  system(paste0("bedtools bamtofastq -i ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted_unmapped.bam -fq ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted_unmapped.fastq"))
  
  ##Qualimap
  
  system(paste0("qualimap bamqc -bam ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted.bam -outdir ",opt$output,"/rawReads2assembly_minimap2/bamqc  --java-mem-size=2G"))
  
  system(paste0("samtools index ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted.bam "))
  
  # system(paste0("awk -f /home/bioinf/extradrv/fastq2fasta-master/fastq2fasta.awk \
  #  $Overall_output_directory/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted_unmapped.fastq > \
  #  $Overall_output_directory/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted_unmapped.fasta
  # 

  cat(paste(c("End:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
  
  
}else {
  cat("No long reads supplied: skipping Minimap2", '\n')
  
}


cat("##--------------------------------",'\n')   
cat("##002 BWA Index",'\n')
cat("##--------------------------------",'\n') 
if (opt$index=="Y") {
  system(paste0("bwa index ",opt$assembly))
  
} else if (opt$index=="N") {
  cat("No BWA index needed: skipping indexing", '\n')
  
} else {
  cat("Index parameter not recognised", '\n')
  
}

cat("##--------------------------------",'\n')   
cat("##02 AllIllumina Reads 2 Assembly",'\n')
cat("##--------------------------------",'\n') 


if (opt$forward!="NA") {


cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
system(paste0("mkdir -p ",opt$output,"/Allillumina2assembly/bamqc"))

#system(paste0("bwa index ",opt$assembly))

#system(paste0("bwa mem -t ",opt$threads," ",opt$assembly ," \'<zcat ",opt$forward, " ",opt$forward  , "\' > ",opt$output, "/Allillumina2assembly/rawIllumina_bwamem_Scaffolds.sam"))

system(paste0("bwa mem -t ",opt$threads," ",opt$assembly ," \'<zcat ",opt$forward, " ",opt$reverse  , " \' | samtools sort -@8 -O BAM -o ",opt$output, "/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted.bam - "))
# | samtools sort -@8 -O BAM -o $Overall_output_directory/0${num}_${run}Freebayes/rawIlluminaMapped/rawIllumina_${name}_bwamem_Assembly_sorted.bam - 

#system(paste0("samtools view -bS " ,opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds.sam | samtools sort - ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted "))



##Unmapped
system(paste0("samtools view -b -f 4 ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted.bam > ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted_unmapped.bam" ))
system(paste0("bedtools bamtofastq -i ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted_unmapped.bam -fq ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted_unmapped.fastq"))

##Qualimap

system(paste0("qualimap bamqc -bam ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted.bam -outdir ",opt$output,"/Allillumina2assembly/bamqc  --java-mem-size=2G"))
system(paste0("samtools index ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted.bam "))



} else {
  cat("No Illumina data given", '\n')

}



if (opt$BlastDB!="NA") {

  

cat("##--------------------------------",'\n')   
cat("##03 Blast",'\n')
cat("##--------------------------------",'\n') 
cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
system(paste0("mkdir -p ",opt$output,"/Blast/"))


paste0("blastn -num_threads " ,opt$threads," -max_hsps 1 -max_target_seqs 1 -task megablast -show_gis -query ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted.bam " \
" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen\" " , \
" -db ", opt$db ," -out " ,opt$output,"/Blast/rawIllumina_bwamem_blast.txt "," -evalue 0.01 -word_size 12")

}else {
  cat("No BLASTdb: skipping blastn", '\n')
  
}


#print()

# print(paste("Genome Assembly:",opt$assembly))

# 
# 
# print(opt$file)
# #dat <- read.csv(file = "/home/vincent/Desktop/Projects/oldRMKs_20181031/02_script/names.txt", header = FALSE)
# # 
# main <- function() {
#   # args <- commandArgs(trailingOnly = TRUE)
#   filename <- opt$file
#   #top <- args[2]
#   
#   
#   dat <- read.csv(file = filename, header = FALSE)
#   
#   
#   list <- as.vector(dat)$V1
#   
#   for (id in list) {
#     # print("===================================================")
#     # print(top)
#     # print("===================================================")
#     # 
#     
#     print("===================================================")
#     print(id)
#     print("===================================================")
#     
#     
#     system(paste0("grep 's__' /home/vincent/Desktop/Projects/oldRMKs_20181031/03_metaphlan2/",id,"_raw/metaphlan_profile/profiled_metagenome_",id,"_raw.txt | grep -v 't__'| cut -d '|' -f 7
#                   "))
#     
#   }
#   
#   
# }
# 
# main()


# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$file)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# 
#id <- "hallo"
# list <- as.vector(args[1])

# ####============quick overview================
# for (id in list) {
#   
#   
#   print("===================================================") 
#       print(id)
#   print("===================================================") 
#   d
#   
#   system(paste0("grep 's__' /home/vincent/Desktop/Projects/oldRMKs_20181031/03_metaphlan2/",id,"_raw/metaphlan_profile/profiled_metagenome_",id,"_raw.txt | grep -v 't__'| cut -d '|' -f 7
# "))
#   
# }



