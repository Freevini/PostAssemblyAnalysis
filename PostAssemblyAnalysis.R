#!/usr/bin/env Rscript

##Usage example
#./PostAssemblyAnalysis.R -a data/testGenome.fasta -F data/testReads_R1.fastq -R data/testReads_R2.fastq -o testOut/
#install.packages("optparse")
library("optparse")

option_list = list(
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Input Assembly", metavar="fasta"),
  
  make_option(c("-x", "--name"), type="character", default=NULL, 
              help="sample Name", metavar="NA"),
  
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
  
  make_option(c("-s", "--sniffles"), type="character", default="NULL", 
              help="Do you want to sniffle for rearrangments Y/N [default N]", metavar="N"),
  
  make_option(c("-n", "--nucmer"), type="character", default="NULL", 
              help="Do you want to run nucmer Y/N [default N]", metavar="N"),
  
  make_option(c("-db", "--BlastDB"), type="character", default="Stdout", 
              help="Blast database used for blast search", metavar="blast")
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$assembly)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat(paste(c("Sample name:",opt$name), collapse='\t'), '\n')
cat(paste(c("Genome Assembly:",opt$assembly), collapse='\t'), '\n')
cat(paste(c("Index: ",opt$index), collapse='\t'), '\n')
cat(paste(c("forward Read:",opt$forward), collapse='\t'), '\n')
cat(paste(c("reverse Read:",opt$reverse), collapse='\t'), '\n')
cat(paste(c("Long Reads:",opt$LongRead), collapse='\t'), '\n')
cat(paste(c("output directory:",opt$output), collapse='\t'), '\n')
cat(paste(c("number of threads:",opt$threads), collapse='\t'), '\n')
cat(paste(c("Sniffles:",opt$sniffles), collapse='\t'), '\n')
cat(paste(c("Nucmer:",opt$nucmer), collapse='\t'), '\n')
cat(paste(c("Blast DB:",opt$BlastDB), collapse='\t'), '\n')


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

  cat(paste(c("Start and go:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')

  #system(paste0("cat ",opt$forward, " ",opt$forward," > ", opt$output,"/RawReadsMerged/rawReadsMerged.fastq"))


  #system(paste0("minimap2 -t ",opt$threads," -ax map-ont ",opt$assembly," ", opt$LongRead, " > ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam"))
  #system(paste0("samtools view -bS " ,opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds.sam | samtools sort - ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted "))

#cat("goGogo")
 system(paste0("minimap2 -a -t ",opt$threads," -ax map-ont ",opt$assembly," ", opt$LongRead, " | samtools sort -@",opt$threads," -O BAM -o ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted.bam - "))
 #system(paste0("minimap2 -t ",opt$threads," -ax map-ont ",opt$assembly," ", opt$LongRead, " | samtools sort -@8 -O BAM -o ",opt$output,"/rawReads2assembly_minimap2/rawReads_minimap2Mapped_Scaffolds_sorted.bam - "))



  ##Unmapped
  system(paste0("samtools view -b -f 4 ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted.bam > ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted_unmapped.bam" ))
  system(paste0("bedtools bamtofastq -i ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted_unmapped.bam -fq ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted_unmapped.fastq"))

  ##Qualimap

  system(paste0("qualimap bamqc -bam ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted.bam -outdir ",opt$output,"/rawReads2assembly_minimap2/bamqc  --java-mem-size=2G"))

  system(paste0("samtools index ",opt$output,"/rawReads2assembly_minimap2/",opt$name,"_rawReads_minimap2Mapped_sorted.bam "))

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

system(paste0("bwa mem -t ",opt$threads," ",opt$assembly ," \'<zcat ",opt$forward, " ",opt$reverse  , " \' | samtools sort -@",opt$threads," -O BAM -o ",opt$output, "/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_Scaffolds_sorted.bam - "))
# | samtools sort -@8 -O BAM -o $Overall_output_directory/0${num}_${run}Freebayes/rawIlluminaMapped/rawIllumina_${name}_bwamem_Assembly_sorted.bam -

#system(paste0("samtools view -bS " ,opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds.sam | samtools sort - ",opt$output,"/Allillumina2assembly/rawIllumina_bwamem_Scaffolds_sorted "))



##Unmapped
system(paste0("samtools view -b -f 4 ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted.bam > ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted_unmapped.bam" ))
system(paste0("bedtools bamtofastq -i ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted_unmapped.bam -fq ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted_unmapped.fastq"))

##Qualimap

system(paste0("qualimap bamqc -bam ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted.bam -outdir ",opt$output,"/Allillumina2assembly/bamqc  --java-mem-size=2G"))
system(paste0("samtools index ",opt$output,"/Allillumina2assembly/",opt$name,"_rawIllumina_bwamem_sorted.bam "))



} else {
  cat("No Illumina data given", '\n')

}



if (opt$BlastDB!="NA") {

  

cat("##--------------------------------",'\n')   
cat("##03 Blast",'\n')
cat("##--------------------------------",'\n') 
cat(paste(c("Start:",cat(as.character(Sys.time()[1]))), collapse='\t'), '\n')
system(paste0("mkdir -p ",opt$output,"/Blast/"))
#cat(paste(c("hello","world_01")))


system(paste0("blastn -num_threads " ,opt$threads," -max_hsps 5 -max_target_seqs 5 -task megablast -show_gis -query ",opt$assembly, " -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen \" -db ", opt$BlastDB ," -out " ,opt$output,"/Blast/",opt$name,"_assembly_blast.txt -evalue 0.01 -word_size 12"))

}else {
  cat("No BLASTdb: skipping blastn", '\n')
  
}


##-------   
##NGMLR 
##-------   

if (opt$sniffles=="Y") {
  
  system(paste0("mkdir -p ",opt$output,"/rawReads2assembly_Sniffles/{mapping,sniffles}" ))
  
  system(paste0("ngmlr -t ",opt$threads,
                " -r ", opt$assembly ,
                " -q ", opt$LongRead , " -x ont |samtools sort -@",opt$threads," -O BAM -o ", opt$output ,"/rawReads2assembly_Sniffles/mapping/",opt$name,"_rawReads_nglmrMapped_sorted.bam - " ))
  
  system(paste0("samtools index " ,opt$output ,"/rawReads2assembly_Sniffles/mapping/",opt$name,"_rawReads_nglmrMapped_sorted.bam"))
  
  
  system(paste0("sniffles -t ", opt$threads , " -m " , opt$output ,"/rawReads2assembly_Sniffles/mapping/",opt$name,"_rawReads_nglmrMapped_sorted.bam -v " ,opt$output ,"/rawReads2assembly_Sniffles/sniffles/",opt$name,"_sniffles2assembly.vcf"))
  

} else {
  cat("Sniffles step skipped", '\n')
  
}


###==========================
##nucmer
###==========================



if (opt$nucmer=="Y") {
  
  sytem(paste0("mkdir -p ",opt$output,"/nucmer/" ))
  
  system(paste0("cd " ,opt$output ,"/nucmer/"))
   
  system(paste0("nucmer -maxmatch -nosimplify ",opt$assembly , " " ,opt$assembly))

  system(paste0("delta-filter -i 95 -l 50000 ",opt$output ,"/nucmer/out.delta > " ,opt$output ,"/nucmer/filtered_095_5000bp.txt "))
  
  system(paste0("grep \"[[:space:]]\" ",opt$output ,"/nucmer/filtered_095_5000bp.txt >  | grep -v "/" > " ,opt$output ,"/nucmer/filtered_095_5000bp_CLEANED.txt "))
  
  
} else {
  cat("nucmer step skipped", '\n')
  
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



