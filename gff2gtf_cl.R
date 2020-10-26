# ======== Packages required =========
if("Xmisc" %in% rownames(installed.packages()) == FALSE) {
  install.packages(Xmisc)}
# === setting environment ===
parser <- Xmisc::ArgumentParser$new()
Xmisc::parser$add_usage('Rsubread_index_cl.R [options]')
Xmisc::parser$add_description('An executable R script parsing arguments from Unix-like command line.')
Xmisc::parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
Xmisc::parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
Xmisc::parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
Xmisc::parser$add_argument('--f', type = 'character', help = '"filename",Enter your *.gff filename')

# === variables ====
dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
if (dir.exists(dirPath)){
  print(paste0("Setting ",dirPath," as the working directory"))
} else {
  print("Directory is not existing!")
  quit()
}

# ======== Packages required (cont.) =========
if (!Xmisc::check.packages(rtracklayer)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rtracklayer")} 
if (!Xmisc::check.packages(tkWidgets)){
  BiocManager::install("tkWidgets")} 

# === execute align command ===
# import gff data
gff_data <- rtracklayer::import(f)

# export gtf data
chSgz<-tkWidgets::hasChar(".gz", what = "suffix")
if (chSgz(f)){
  gtf_filename<-gsub(f,pattern = '.gff.gz',replacement = '.gtf')
  rtracklayer::export(gff_data,gtf_filename,"gtf")
} else {
  gtf_filename<-gsub(f,pattern = '.gff',replacement = '.gtf')
  rtracklayer::export(gff_data,gtf_filename,"gtf")
}
