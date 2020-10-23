# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
if("Xmisc" %in% rownames(installed.packages()) == FALSE) {
  install.packages(Xmisc)}
library(Xmisc)
if (!check.packages(BiocManager)) {install.packages("BiocManager")}
if (!check.packages(edgeR)) {BiocManager::install("edgeR")}
if (!check.packages(statmod)) {install.packages("statmod")}
# === setting environment ===
parser <- ArgumentParser$new()
parser$add_usage('Rsubread_featureCount_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line.')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
# annotation
parser$add_argument('--ai', type = 'character', default = "hg38", help = 'a character string specifying an in-built annotation used for read summarization. It has four possible values including"mm10","mm9","hg38"and"hg19"')
parser$add_argument('--ae', type = 'character', default = 'NULL', help = 'A character string giving name of a user-provided annotation file or a data frame including user-provided annotation data. If the annotation is in GTF format, it can only be provided as a file')
parser$add_argument('--iG',type='logical', default = FALSE, help='"isGTFAnnotationFile", if the annotation provided via the annot.extargument is in GTF format or  not. This option is only applicable when --ae is not NULL')
parser$add_argument('--Gf', type = 'character', default = "exon" , help = 'a character string denoting the type of features that will be extracted from a GTF annotation. This argument is only applicable when --iG==TRUE')
parser$add_argument('--GT', type = 'character', default = "gene_id", help = 'a character string denoting the type of attributes in a GTF annotation that will be used to group features.This argument is only applicable when --iG==TRUE')
parser$add_argument('--GTe', type = 'character', default = 'NULL', help = 'a character vector specifying extra GTF attribute types that will also be includedin the counting output. These attribute types will not be used to group features.')
parser$add_argument('--cA', type = 'character', default = 'NULL', help = 'the name of a comma-delimited text file that includes aliases of chromosome names')
# level of summarization
parser$add_argument('--uMF',type='logical', default = TRUE, help='"useMetaFeatures", if the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg. genes)')
# overlap between reads and features
parser$add_argument('--aMO',type='logical', default = FALSE, help='"allowMultiOverlap", a read is allowed to be assigned to more than one feature(or meta-feature), when it is found to overlap with more than one feature (or meta-feature).')
parser$add_argument('--mO', type = 'numeric', default = 1, help = '"minOverlap", a minimum number of overlapped bases required for assigninga read to a feature (or a meta-feature)')
parser$add_argument('--fO', type = 'numeric', default = 0, help = '"fracOverlap", a minimum fraction of overlapping bases in a read that is requiredfor read assignment. Value should be within range [0,1].')
parser$add_argument('--fOF', type = 'numeric', default = 0, help = '"fracOverlapFeature",a minimum fraction of bases included in a feature that is required. for overlapping with a read or a read pair. Value should be within range [0,1].')
parser$add_argument('--lO',type='logical', default = FALSE, help='"largestOverlap",a read (or read pair) will be assigned to the feature (or meta-feature) thathas the largest number of overlapping bases, if the read (or read pair) overlapswith multiple features (or meta-features).')
parser$add_argument('--nO', type = 'numeric', default = NULL, help = '"nonOverlap",the maximum number of bases in a read (or a read pair) that are allowed not to overlap the assigned feature. "NULL" by default (ie. no limit is set)')
parser$add_argument('--nOF', type = 'numeric', default = NULL, help = '"nonOverlapFeature",the maximum number of non-overlapping bases allowed in afeature during read assignment. "NULL" by default (ie. no limit is set)')
# read shift, extension and reduction
parser$add_argument('--rST', type = 'character', default = 'upstream', help = '"readShiftType",a character string indicating how reads should be shifted.  It has four possiblevalues  includingupstream,downstream,leftandright.')
parser$add_argument('--rSS', type = 'numeric', default = 0, help = '"readShiftSize", the number of bases the reads will be shifted by. Negative value is not allowed.')
parser$add_argument('--rE5', type = 'numeric', default = 0, help = '"readExtension5", a number of bases extended upstream from 5 end of each read.Negative value is not allowed.')
parser$add_argument('--rE3', type = 'numeric', default = 0, help = '"readExtension3", a number of bases extended upstream from 3 end of each read.Negative value is not allowed.')
parser$add_argument('--r2p', type = 'numeric', default = NULL, help = '"read2pos", each read should be reduced to its 5 most base or 3 mostbase. It has three possible values:NULL,5(denoting 5 most base) and3(denoting 3 most base). Default value is NULL, ie. no read reduction will be performed.')
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
if (ae=='NULL'){ae <-NULL}
if (GTe=='NULL'){GTe <-NULL}
if (cA=='NULL'){cA <-NULL}

bam_files <- list.files(dirPath,pattern = ".bam$")

# ======== Packages required (cont.) =========
if (!check.packages(Rsubread)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Rsubread")}  
library(Rsubread)

# ========= featurCount ===========
Rsub_fc<-function(bam_files){
  featureCounts(bam_files,
                # annotation
                annot.inbuilt = ai, # c("mm10","mm9","hg38","hg19")
                annot.ext = ae,
                isGTFAnnotationFile = iG,
                GTF.featureType = Gf,
                GTF.attrType = GT,
                GTF.attrType.extra = GTe,
                chrAliases = cA,
                # level of summarization
                useMetaFeatures = uMF,
                # overlap between reads and features
                allowMultiOverlap = aMO,
                minOverlap = mO,
                fracOverlap = fO,
                fracOverlapFeature = fOF,
                largestOverlap = lO,
                nonOverlap = nO,
                nonOverlapFeature = nOF,
                # Read shift, extension and reduction
                readShiftType = rST,
                readShiftSize = rSS,
                readExtension5 = rE5,
                readExtension3 = rE3,
                read2pos = r2p,
                # multi-mapping reads
                countMultiMappingReads = TRUE,
                # fractional counting
                fraction = FALSE,
                # long reads
                isLongRead = FALSE,
                # read filtering
                minMQS = 0,
                splitOnly = FALSE,
                nonSplitOnly = FALSE,
                primaryOnly = FALSE,
                ignoreDup = FALSE,
                # strandness
                strandSpecific = 0,
                # exon-exon junctions
                juncCounts = FALSE,
                genome = NULL,
                # parameters specific to paired end reads
                isPairedEnd = TRUE,  # change to "FLASE" if not.
                requireBothEndsMapped = TRUE,
                checkFragLength = TRUE,
                minFragLength = 50,
                maxFragLength = 600,
                countChimericFragments = TRUE,
                autosort = TRUE,
                # number of CPU threads
                nthreads = 1,
                # read group
                byReadGroup = FALSE,
                # report assignment result for each read
                reportReads = "BAM", #c(CORE, SAM and BAM)
                reportReadsPath = getwd(),
                # miscellaneous
                sampleSheet = NULL,
                cellBarcodeList = NULL,
                maxMOp = 10,
                tmpDir = ".",
                verbose = FALSE)
}

# ===== Import group names from groups.txt or groups.xlsx
print("Importing groups.xlsx")
if (file.exists("groups.xlsx")){
  require(readxl)
  print("groups.xlsx was found")
  groups <- factor(unlist(read_excel("groups.xlsx",range=cell_cols("B"))))
} else if (file.exists("groups.txt")){
  print("groups.xlsx was not found, loading data from groups.txt")
  groups<- read.delim("groups.txt")
  groups<- unlist(groups[,"Groups"])
} else {
  cat("Neither groups.xlsx nor groups.txt was found")
  quit()
}

# ===== Options before analysis =========================
# Processing data if you don't have counts.txt and genes.txt
repeat{
  cat("Analyze data form *.BAM? Yes or No:\n")
  data_source_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", data_source_inp, "\n"))
  data_source_answer<-tolower(data_source_inp)
  if (data_source_answer=="y"|data_source_answer=="yes"|data_source_answer=="n"|data_source_answer=="no"){
    break()
  } else {cat("Please enter yes or no [y/n]!")}
}
# Statistics analysis
repeat{
  cat("Single facotr for statistics? Yes or No:\n")
  factor_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", factor_inp, 'for single facotr statistics',"\n"))
  factor_answer<-tolower(factor_inp)
  if (factor_answer=="y"|factor_answer=="yes"|factor_answer=="n"|factor_answer=="no"){
    break()  
  } else {cat("Please enter yes or no [y/n]!")}
}
# Overall statistics
repeat{
  cat("Wanna Overall Differential expression analysis? Yes or No:\n")
  overall_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", overall_inp, 'for overall Differential expression analysis',"\n"))
  overall_answer<-tolower(overall_inp)
  if (overall_answer=="y"|overall_answer=="yes"|overall_answer=="n"|overall_answer=="no"){
    break() 
  } else {cat("Please enter yes or no [y/n]!")}
}
# Comparison statistics
repeat{
  cat("Wanna Compare Two groups Differential expression analysis? Yes or No:\n")
  compare_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", compare_inp, 'for Differential expression analysis',"\n"))
  compare_answer<-tolower(compare_inp)
  if (compare_answer=="y"|compare_answer=="yes"){
    con<-vector()
    for (n in 1:length(levels(groups))){
      con<-c(con,0)
    }
    groups_list<-matrix(levels(groups))
    colnames(groups_list)<-"group"
    print(groups_list)
    repeat{
      n<-as.integer(readline(prompt = paste0("Please enter the 1st group number [1-",length(groups_list),"]:")))
      m<-as.integer(readline(prompt = paste0("Please enter the 2nd group number [1-",length(groups_list),"]:")))
      if (!n %in% 1:length(groups_list) && !m %in% 1:length(groups_list)){
        print("The given groups are not the group list. Please enter correct groups.")
      } else if (!n %in% 1:length(groups_list)) {
        print("The 1nd group is not the group list. Please re-enter correct groups.")
      } else if (!m %in% 1:length(groups_list)){
        print("The 2nd group is not the group list. Please re-enter correct groups.")
      } else if (n==m){
        print("The given groups are the same one! Please re-enter two different groups.")
      } else {
        break()
      }
    }
    group_1 <- groups_list[n,]
    con[n]<- -1
    group_2<- groups_list[m,]
    con[m]<- 1
    break()  
  } else if (compare_answer=="n"|compare_answer=="no"){
    break()
  } else {
    cat("Please enter yes or no [y/n]!")
  }
}
# ==========================================================
# Processing data if you don't have counts.txt and genes.txt 

if (character_to_logical(data_source_answer)){
  # analyzing data
  fc<-Rsub_fc(bam_files)
  # Create a DGEList object from fc
  y<-DGEList(counts = fc$counts[,1:length(bam_files)],genes = fc$annotation[,c("GeneID","Length")], group = groups)
  # export data
  write.table(y$counts, file="counts.txt", row.names = TRUE, col.names = TRUE, sep = "\t" )
  write.table(y$genes, file="genes.txt", row.names = TRUE, col.names = TRUE, sep = "\t" )
  setwd('C:\\Users\\Changyi.Lin\\Documents\\DNA sequencing template')
} else if (!character_to_logical(data_source_answer)) {
  print("Analyzing data from counts.txt and genes.txt")
  # import data
  counts.txt <- "counts.txt"
  raw.data.counts <- read.delim(counts.txt)
  genes.txt <- "genes.txt"
  raw.data.genes <- read.delim(genes.txt)
  # Create a DGEList object from counts.txt and genes.txt
  require(edgeR)
  y <- DGEList(counts = raw.data.counts, genes = data.frame(raw.data.genes), group = groups)
} 

# calculate normal factors
y<-calcNormFactors(y)

# edit design
design <- model.matrix(~ 0 + groups, data = y$samples)
colnames(design) <- levels(groups)

# ====== Statistics (edgeR) ==================
n_occur <- data.frame(table(groups))
if (1 %in% n_occur[,"Freq"]){
  print("One or more groups have No Duplicates")
  y<-estimateGLMCommonDisp(y, method = "deviance", robust=TRUE, subset = NULL)
  # ====== There are three options if no replicates ============
  # Option 1: asign a dispersion
  # well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
  # Assign Biological coefficient of variation, which is the square root of the dispersion parameter under the negative binomial mode
  #bcv <- 0.4
  #et <- exactTest(y, dispersion=bcv^2, pair = c("C","A"))
  # Option 2: use estimateGLMCommonDisp()
  #y<-estimateGLMCommonDisp(y, method = "deviance", robust=TRUE, subset = NULL)
  # Option 3: use housekeeping gene(s)
  # Assign which gene is houskeeping gene which is not affected by treatment
  #housekeeping<-"ACTB"
  # create a copy of the data object with only one treatment group
  #y1<-y
  # Assign all samples as a group in new dataset
  #y1$samples$group <- 1
  # Get the Dispersion from new dataset
  #y0 <- estimateDisp(y1[housekeeping,], trend="none", tagwise=FALSE)
  # Assign the dispersion to the original dataset
  #y$common.dispersion <- y0$common.dispersion
  
} else {
  if (character_to_logical(factor_answer,ignore.case = TRUE)){
    # ====== Pairwise comparisons between two or more groups (classic)
    # = quantile-adjusted conditional maximum likelihood (qCML)
    # = qCML method is only applicable on datasets with "a single factor" design
    # = This method proves to be accurate and nearly unbiased even for small counts and small numbers of replicates
    require(statmod)
    y <- estimateDisp(y, robust=TRUE)
    # ===== Alternatives below to get y$pseudo.counts,
    #y <- estimateCommonDisp(y)
    #y <- estimateTagwiseDisp(y)
    #row.names(y$pseudo.counts)<-gene_names
  } else if (!character_to_logical(factor_answer,ignore.case = TRUE)) {
    # ====== For More complex experiments (glm functionality)
    y <- estimateDisp(y, design = design, robust=TRUE)
    # ===== Alternatives below to get y$pseudo.counts,
    #y <- estimateGLMCommonDisp(y, design)
    #y <- estimateGLMTrendedDisp(y, design)
    #y <- estimateGLMTagwiseDisp(y, design)
  } 
}
library(DESeq2)
x<-DESeqDataSetFromMatrix(countData = raw.data.counts,colData = y$samples, design = design)
x
