# ======== Packages required =========
# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
# Rcran
packages<-c('Xmisc','gsubfn')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio.packages<-c('Shortread')
for (bio.package in bio.packages){
  if(bio.package %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(bio.package)}
}
# === setting environment ===
parser <- ArgumentParser$new()
parser$add_usage('merge_fastqgz_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line. Merge multiple lanes fastq.gz data into one')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
parser$helpme()

# === variables ====

dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
setwd(dirPath)
if (dir.exists(dirPath)){
  cat(paste0("Setting ",dirPath," as the working directory\n"))
} else {
  cat("Directory is not existing!\n")
  quit()
}
# =============================
R1.gz.files<-list.files(dirPath, 'L00[1-9]_R1_001.fastq.gz$')
R2.gz.files<-list.files(dirPath, 'L00[1-9]_R2_001.fastq.gz$')
# merge R1 files
require(ShortRead)
require(gsubfn)
R1.sample_names<-NULL
for(gz.file in R1.gz.files){
  R1.sample_name<-regmatches(gz.file, gregexpr(".*?_L", gz.file))[[1]]
  R1.sample_names<-c(R1.sample_names,R1.sample_name)
}
R1.group_names<-row.names(table(R1.sample_names))

for (group_name in R1.group_names){
  fls<-NULL
  for (gz.file in R1.gz.files){
    if (regmatches(gz.file, gregexpr(".*?_L", gz.file))[[1]]==group_name){
      fls<-c(fls,gz.file)
    }
  }
  fout = paste0(gsub(x = group_name,pattern = '_L',replacement = '_R1'),".fastq.gz")
  for (fl in fls) {
    fq = readFastq(fl)
    writeFastq(fq, fout, mode="a", compress=TRUE)
  }
}
# merge R2 files if paired-end
if(!rapportools::is.empty(R2.gz.files)){
  R2.sample_names<-NULL
  for(gz.file in R2.gz.files){
    R2.sample_name<-regmatches(gz.file, gregexpr(".*?_L", gz.file))[[1]]
    R2.sample_names<-c(R2.sample_names,R2.sample_name)
  }
  R2.group_names<-row.names(table(R2.sample_names))
  for (group_name in R2.group_names){
    fls<-NULL
    for (gz.file in R2.gz.files){
      if (regmatches(gz.file, gregexpr(".*?_L", gz.file))[[1]]==group_name){
        fls<-c(fls,gz.file)
      }
    }
    fout = paste0(gsub(x = group_name,pattern = '_L',replacement = '_R2'),".fastq.gz")
    for (fl in fls) {
      fq = readFastq(fl)
      writeFastq(fq, fout, mode="a", compress=TRUE)
    }
  }
}
