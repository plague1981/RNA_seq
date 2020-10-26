# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
if("Xmisc" %in% rownames(installed.packages()) == FALSE) {
  install.packages('Xmisc')}
library(Xmisc)
# === setting environment ===
parser <- ArgumentParser$new()
parser$add_usage('edgeR_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line.')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
parser$helpme()

# === variables ====

dirPath <- dir
setwd(dirPath)
dirPath <-gsub ('\\\\','/',dirPath)
if (dir.exists(dirPath)){
  cat(paste0("Setting ",dirPath," as the working directory\n"))
} else {
  cat("Directory is not existing!\n")
  quit()
}
# === packages requirement ====
if (!check.packages('ShortRead')){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install('ShortRead')} 
if (!check.packages('gsubfn')){install.packages('gsubfn')}
library('ShortRead')
library('gsubfn')
# =============================
all.gz.files<-list.files(dirPath, 'L00[1-9]_R[1,2]_001.fastq.gz$')

sample_names<-NULL
for(gz.file in all.gz.files){
  sample_name<-regmatches(gz.file, gregexpr(".*?_L", gz.file))[[1]]
  sample_names<-c(sample_names,sample_name)
}
group_names<-row.names(table(sample_names))


for (group_name in group_names){
  fls<-NULL
  for (gz.file in all.gz.files){
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

