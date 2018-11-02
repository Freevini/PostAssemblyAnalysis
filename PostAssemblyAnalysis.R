#!/usr/bin/env Rscript
#./PostAssemblyAnalysis.R -a data/testGenome.fasta -F data/testReads_R1.fastq -R data/testReads_R2.fastq -o testOut/

library("optparse")

option_list = list(
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Input Assembly", metavar="fasta"),
  
  make_option(c("-F", "--forward"), type="character", default="NA", 
              help="Illumina Forward raw read in fastq.gz", metavar="fastq"),
  
  make_option(c("-R", "--reverse"), type="character", default="NA", 
              help="Illumina Reverse raw read in fastq.gz", metavar="fastq"),
  
  make_option(c("-l", "--LongRead"), type="character", default="NA", 
              help="Long reads raw in fastq.gz", metavar="fastq"),
  
  make_option(c("-o", "--output"), type="character", default="Stdout", 
              help="Output directory for all files created", metavar="PATH"),
  
  make_option(c("-t", "--threads"), type="numeric", default="5", 
              help="number of threads", metavar="threads")
  
  
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$assembly)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat(paste(c("Genome Assembly:",opt$assembly), collapse='\t'), '\n')
cat(paste(c("forward Read:",opt$forward), collapse='\t'), '\n')
cat(paste(c("reverse Read:",opt$reverse), collapse='\t'), '\n')
cat(paste(c("Long Reads:",opt$LongRead), collapse='\t'), '\n')
cat(paste(c("output directory:",opt$output), collapse='\t'), '\n')
cat(paste(c("number of threads:",opt$threads), collapse='\t'), '\n')




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
  
  system(paste0("cat ",opt$forward, " ",opt$forward," > ", opt$output,"/RawReadsMerged/rawReadsMerged.fastq"))
  
  system(paste0("minimap2 -t ",opt$threads," -ax map-pb ",opt$assembly," ", opt$output, "/RawReadsMerged/rawReadsMerged.fastq > ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam"))
  system(paste0("samtools view -bS " ,opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam | samtools sort - ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted "))
  
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
cat("##02 AllIllumina Reads 2 Assembly",'\n')
cat("##--------------------------------",'\n') 
cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
system(paste0("mkdir -p ",opt$output,"/Allillumina2assembly/bamqc"))

system(paste0("bwa index ",opt$assembly))

system(paste0("bwa mem $Assembly_Input \'<zcat ",opt$forward, opt$reverse ,"\' -t " ,opt$threads, "> ",opt$output, "/Allillumina2assembly/rawIllumina_bwamem_Scaffolds.sam"))

#|samtools sort -@8 -O BAM -o $Overall_output_directory/Allillumina2assembly/rawIllumina_${species}_${name}_bwamem_Scaffolds_sorted.bam - | tee $Overall_output_directory/Allillumina2assembly/rawIllumina_${species}_${name}_bwamem_Scaffolds_sorted.log

  ###--------------
  ##03_AllIllumina2Assembly
  ###--------------


  # 
  # $QUAL_ROOT/qualimap bamqc -bam $Overall_output_directory/Allillumina2assembly/rawIllumina_${species}_${name}_bwamem_Scaffolds_sorted.bam -outdir $Overall_output_directory/Allillumina2assembly/bamqc --java-mem-size=2G
  # 
  # samtools index $Overall_output_directory/Allillumina2assembly/rawIllumina_${species}_${name}_bwamem_Scaffolds_sorted.bam


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
#   
#   
#   system(paste0("grep 's__' /home/vincent/Desktop/Projects/oldRMKs_20181031/03_metaphlan2/",id,"_raw/metaphlan_profile/profiled_metagenome_",id,"_raw.txt | grep -v 't__'| cut -d '|' -f 7
# "))
#   
# }


