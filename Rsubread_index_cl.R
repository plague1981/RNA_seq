# ======== Packages required =========
if("Xmisc" %in% rownames(installed.packages()) == FALSE) {
  install.packages(Xmisc)}
library(Xmisc)

# === setting environment ===
parser <- ArgumentParser$new()
parser$add_usage('Rsubread_index_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line.')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
parser$add_argument('--i', type = 'character', default = 'hg', help = '"index",Name your index')
parser$add_argument('--db', type = 'character', default = 'genome.fa', help = '"database",basename of created index files')
parser$add_argument('--gI', type = 'logical', default = FALSE, help = '"gappedIndex", logical  indicating  if  a  gapped  index  or  a  full  index  will  be  built.A gappedindex contains 16mers (subreads) that are extracted every three bases from areference genome, whereas a full index contains subreads extracted from everychromosomal location of a genome.   The index contains a hash table,  whichincludes sequences of subreads and their corresponding chromosomal locations.Default value of this argument isFALSE(ie. a full index is built).')
parser$add_argument('--iS', type = 'logical', default = FALSE, help = '"indexSplit", logical indicating if an index can be split into multiple blocks. The block size isdetermined by value of parametermemory.FALSEby default (ie. a single-blockindex is generated).')
parser$add_argument('--m', type = 'numeric', default = 8000, help = '"memory", a numeric value specifying the amount of memory (in megabytes) used for storing the index during read mapping. 8000 MB by default. Note that this optionis ignored when "iS" is FALSE')
parser$add_argument('--TH', type = 'numeric', default = 100, help = '"threshold",a numeric value specifying the threshold for removing highly repetitive subreads(16bp mers).  100 by default.  Subreads will be excluded from the index if theyoccur more than threshold number of times in the genome.')
parser$add_argument('--cs', type = 'logical', default = FALSE, help = '"colorspace",logical specifying the mode of the index. If TRUE, a color space index will bebuilt. Otherwise, a base space index will be built.')
parser$helpme()
# === variables ====
dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
  if (dir.exists(dirPath)){
    print(paste0("Setting ",dirPath," as the working directory"))
  } else {
    print("Directory is not existing!")
    quit()
  }
basename <- i
database <- db

# ======== Packages required (cont.) =========
if (!check.packages(Rsubread)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("Rsubread")}  

library(Rsubread)

# ======== build index (biuld index in server if possible) ==========
if (!file.exists(paste0(basename,".00.b.array"))){
    if (file.exists(database)){
      ref<-file.path(dirPath, database)
      print(paste0("Setting ",database," as the index file for Rsubread"))
    }else {
      print(paste0(database," does not exist, please check the filename and try again!"))
      quit()
    }
buildindex(basename = basename, reference = ref, gappedIndex = gI, indexSplit = iS, memory = m, TH_subread = TH, colorspace = cs)
} else {
  print(paste(basename,"index exists. Please run Rsubread_align.R"))
}
