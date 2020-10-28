# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
# Rcran
packages<-c('Xmisc')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio.packages<-c("Rhisat2")
for (bio.package in bio.packages){
  if(bio.package %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(bio.package)}
}

# === setting environment ===
parser <- Xmisc::ArgumentParser$new()
parser$add_usage('Rsubread_featureCount_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line.')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
parser$add_argument('--ref', type = 'character', default = 'genome.fa', help = '"reference",Enter your reference filename (eg. *.fa), not compressed!')
parser$add_argument('--i', type = 'character', default = 'hg', help = '"index",Enter your index prefix')
parser$add_argument('--od', type = 'character', default = getwd(), help = '"outdir",Enter your output directory')
parser$helpme()

# === variables ====
dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
od <-gsub ('\\\\','/',od)
if (dir.exists(dirPath) && dir.exists(od)){
  setwd(dirPath)
  cat(paste0("Setting ",dirPath," as the working directory\n"))
} else if (!dir.exists(dirPath) && !dir.exists(od)) {
  cat("Both directories are not existing!\n")
  quit()
} else if (!dir.exists(dirPath)){
  cat("Reference directory is not existing!\n")
  quit()
} else if (!dir.exists(od)){
  cat("Output directory is not existing!\n")
  quit()
}
# ======

Rhisat2::hisat2_build(references = ref, outdir = od, prefix = i, force = TRUE, strict = TRUE, execute = TRUE)
