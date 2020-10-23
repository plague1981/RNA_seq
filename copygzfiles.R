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
parser$add_argument('--fdir', type = 'character', default = getwd(), help = '"directory",Enter your source directory')
parser$add_argument('--tdir', type = 'character', default = getwd(), help = '"directory",Enter your destination directory')
parser$helpme()

# Copy *.fastq.gz files to destination fold
fdirPath <-gsub ('\\\\','/',fdir)
tdirPath <-gsub ('\\\\','/',tdir)
if (dir.exists(fdirPath) && dir.exists(tdirPath)){
  cat(paste0("Copying *.gz files from ",fdirPath,'to',tdirPath,"\n"))
} else {
  cat("At least one of Directories is not existing!\n")
  quit()
}
dirs<-list.dirs(fdirPath)

gz.files<-NULL
for (dir in dirs){
  gz.file<-list.files(dir,'gz$')
  gz.file<-file.path(dir,gz.file)
  gz.files<-c(gz.files,gz.file)
}

file.copy(gz.files, tdirPath)
