# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
packages<-c('Xmisc','readxl','xlsx')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio.packages<-c('Rsubread','tkWidgets','edgeR')
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
parser$add_argument('--gf', type = 'character', default = 'groups.xlsx', help = '"filename",Enter the filename of the group file')
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
parser$add_argument('--nO', type = 'integer', default = 0, help = '"nonOverlap",the maximum number of bases in a read (or a read pair) that are allowed not to overlap the assigned feature. "NULL" by default (ie. no limit is set)')
parser$add_argument('--nOF', type = 'integer', default = 0, help = '"nonOverlapFeature",the maximum number of non-overlapping bases allowed in afeature during read assignment. "NULL" by default (ie. no limit is set)')
# read shift, extension and reduction
parser$add_argument('--rST', type = 'character', default = 'upstream', help = '"readShiftType",a character string indicating how reads should be shifted.  It has four possiblevalues  includingupstream,downstream,leftandright.')
parser$add_argument('--rSS', type = 'numeric', default = 0, help = '"readShiftSize", the number of bases the reads will be shifted by. Negative value is not allowed.')
parser$add_argument('--rE5', type = 'numeric', default = 0, help = '"readExtension5", a number of bases extended upstream from 5 end of each read.Negative value is not allowed.')
parser$add_argument('--rE3', type = 'numeric', default = 0, help = '"readExtension3", a number of bases extended upstream from 3 end of each read.Negative value is not allowed.')
parser$add_argument('--r2p', type = 'integer', default = 0, help = '"read2pos", each read should be reduced to its 5 most base or 3 mostbase. It has three possible values:NULL,5(denoting 5 most base) and3(denoting 3 most base). Default value is NULL, ie. no read reduction will be performed.')
# multi-mapping reads
parser$add_argument('--MMR',type='logical', default = TRUE, help='"countMultiMappingReads",if multi-mapping reads/fragments should be counted. "NH" tag is used to located multi-mapping reads in the input BAM/SAMfiles.')
# fractional counting
parser$add_argument('--fr',type='logical', default = FALSE, help='"fraction", if fractional counts are produced for multi-mapping readsand/or multi-overlapping reads')
# long reads
parser$add_argument('--iL',type='logical', default = FALSE, help='"isLongRead", if input data contain long reads.  This option should be set toTRUEif counting Nanopore or PacBio long reads.')
# read filtering
parser$add_argument('--mM',type ='integer', default = 0, help = '"minMQS", the minimum mapping quality score a read must satisfy in orderto be counted. For paired-end reads, at least one end should satisfy this criteria.')
parser$add_argument('--sO',type='logical', default = FALSE, help='"splitOnly", if  only split alignments (their CIGAR strings containletter "N") should be included for summarization.')
parser$add_argument('--nS',type='logical', default = FALSE, help='"nonSplitOnly", if only non-split alignments (their CIGAR strings donot contain letter "N") should be included for summarization')
parser$add_argument('--pO',type='logical', default = FALSE, help='"primaryOnly", if only primary alignments should be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files')
parser$add_argument('--iD',type='logical', default = FALSE, help='"ignoreDup", if reads marked as duplicates should be ignored.')
# strandness
parser$add_argument('--sS',type ='integer', default = 0, help = '"strandSpecific", an integer vector indicating if strand-specific read counting should be performed.Length of the vector should be either1(meaning that the value is applied to allinput files), or equal to the total number of input files provided. Each vector ele-ment should have one of the following three values:0(unstranded),1(stranded)and2(reversely stranded). Default value of this parameter is0(ie. unstrandedread counting is performed for all input files)')
# exon-exon junctions
parser$add_argument('--jC',type='logical', default = FALSE, help='"juncCounts", if number of reads supporting each exon-exon junction will bereported. Junctions are identified from those exon-spanning reads in input data.')
parser$add_argument('--g', type = 'character', default = 'NULL', help = '"genome",the name of a FASTA-format file that includes the ref-erence sequences used in read mapping that produced the provided SAM/BAMfiles.')
# parameters specific to paired end reads
parser$add_argument('--iPE',type='logical', default = FALSE, help='"isPairedEnd", if counting should be performed on read pairs or reads.')
parser$add_argument('--BEM',type='logical', default = FALSE, help='"requireBothEndsMapped", if both ends from the same fragment are required to be suc-cessfully  aligned  before  the  fragment can be assigned to a feature or meta-feature. This parameter is only appliable when isPairedEnd is TRUE')
parser$add_argument('--cFL',type='logical', default = FALSE, help='"checkFragLength", if the two ends from the same fragment are required to satisify the fragment length criteria before the fragment can be assigned to a feature or meta-feature.')
parser$add_argument('--miF',type ='integer', default = 50, help = '"minFragLength", the minimum fragment length for paired-end reads.')
parser$add_argument('--maF',type ='integer', default = 600, help = '"maxFragLength", the maximum fragment length for paired-end reads.')
parser$add_argument('--cCF',type='logical', default = TRUE, help='"countChimericFragments", if a chimeric fragment, which has its two reads mappedto different chromosomes, should be counted or not.')
parser$add_argument('--a',type='logical', default = TRUE, help='"autosort", if the automatic read sorting is enabled.')
# number of CPU threads
parser$add_argument('--nt', type = 'integer', default = 1, help = '"nthreads",integer giving  the  number  of  threads  used  for  running  this  function.')
# read group
parser$add_argument('--bRG',type='logical', default = FALSE, help='"byReadGroup", read counting will be performed for each individual read group.')
# report assignment result for each read
parser$add_argument('--rR', type = 'character', default = 'NULL', help = '"reportReads",output detailed read assignment results for each read (or fragment if paired end).')
parser$add_argument('--rRP', type = 'character', default = 'NULL', help = '"reportReadsPath", the directory where files including detailed assign-ment results are saved.  If NULL, the results will be saved to the current working directory.')
# miscellaneous
parser$add_argument('--s', type = 'character', default = 'NULL', help = '"sampleSheet", a  character  string  specifying  the  single-cell  RNA  sample  sheet  file.  If NULL, featureCounts runs on the bulk RNAseq mode.')
parser$add_argument('--cBL', type = 'character', default = 'NULL', help = '"cellBarcodeList", the  file name containing  the  list  of  cell  barcodes  forscRNA sample preparation.')
parser$add_argument('--mM', type = 'number', default = 10, help = '"maxMOp",the maximum number of "M" operations (matches or mis-matches)allowed in a CIGAR string.')
parser$add_argument('--tD', type = 'character', default = '.', help = '"tmpDir", the  directory  under  which  intermediate  files  aresaved (later removed). By default, current working directory is used.')
parser$add_argument('--v',type='logical', default = FALSE, help='"verbose", if verbose information for debugging will be generated.')

parser$helpme()
# === variables ====
dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
if (dir.exists(dirPath)){
  setwd(dirPath)
  print(paste0("Setting ",dirPath," as the working directory"))
} else {
  print("Directory is not existing!")
  quit()
}
if (ae=='NULL'){ae <-NULL}
if (GTe=='NULL'){GTe <-NULL}
if (cA=='NULL'){cA <-NULL}
if (g=='NULL'){g <-NULL}
if (rR=='NULL'){rR <-NULL}
if (rRP=='NULL'){rRP <-NULL}
if (s=='NULL'){s <-NULL}
if (cBL=='NULL'){cBL <-NULL}
# ===== Import group names from groups.txt or groups.xlsx

chSxlsx<-tkWidgets::hasChar(".xlsx", what = "suffix")
chStxt<-tkWidgets::hasChar(".txt", what = "suffix")
cat(paste0("Importing ",gf,'\n'))
if (file.exists(gf) && chSxlsx(gf)){
  cat(paste0(gf," was found\n"))
  groups <- factor(unlist(readxl::read_excel("groups.xlsx",range=readxl::cell_cols("B"))))
  groups_name <- unlist(readxl::read_excel("groups.xlsx",range=readxl::cell_cols("A")))
} else if (file.exists(gf) && chStxt(gf)){
  groups<- read.delim(gf)
  groups<- factor(unlist(groups[,"Groups"]))
  groups_name<- unlist(groups[,"Sample_names"])
} else {
  cat(paste0("Cannot find ",gf,'\n'))
  quit()
}
names(groups)<-groups_name
# ==========================================================
# Get Count table
# Processing for *.BAM

BAM.files<-list.files(dirPath,pattern = ".BAM$",ignore.case = TRUE)
path.BAM.files <- file.path(dirPath, BAM.files)
# re-order groups based on BAM.files
group_order<-NULL
for (n in 1:length(BAM.files)){
  if (!is.element(BAM.files[n], names(groups))){
    cat(paste0(BAM.files[n], ' dose not exist in groups\n'))
    cat(paste0('Please check sample names and ',gf, ' \n'))
    quit()
  }
  for (m in 1:length(names(groups))){
    if (BAM.files[n]==names(groups)[m]){
      group_order<-c(group_order,m)
      break()
    } 
  }
}
groups<-groups[group_order]
# analyzing data
fc<-Rsubread::featureCounts(path.BAM.files
                            # annotation
                            ,annot.inbuilt = ai # c("mm10","mm9","hg38","hg19")
                            ,annot.ext = ae
                            ,isGTFAnnotationFile = iG
                            ,GTF.featureType = Gf
                            ,GTF.attrType = GT
                            ,GTF.attrType.extra = GTe
                            ,chrAliases = cA
                            # level of summarization
                            ,useMetaFeatures = uMF
                            # overlap between reads and features
                            ,allowMultiOverlap = aMO
                            ,minOverlap = mO
                            ,fracOverlap = fO
                            ,fracOverlapFeature = fOF
                            ,largestOverlap = lO
                            #,nonOverlap = nO
                            #,nonOverlapFeature = nOF
                            # Read shift, extension and reduction
                            ,readShiftType = rST
                            ,readShiftSize = rSS
                            ,readExtension5 = rE5
                            ,readExtension3 = rE3
                            #,read2pos = r2p
                            # multi-mapping reads
                            ,countMultiMappingReads = MMR
                            # fractional counting
                            ,fraction = fr
                            # long reads
                            ,isLongRead = iL
                            # read filtering
                            ,minMQS = mM
                            ,splitOnly = sO
                            ,nonSplitOnly = nS
                            ,primaryOnly = pO
                            ,ignoreDup = iD
                            # strandness
                            ,strandSpecific = sS
                            # exon-exon junctions
                            ,juncCounts = jC
                            ,genome = g
                            # parameters specific to paired end reads
                            ,isPairedEnd = iPE
                            ,requireBothEndsMapped = BEM
                            ,checkFragLength = cFL
                            ,minFragLength = miF
                            ,maxFragLength = maF
                            ,countChimericFragments = cCF
                            ,autosort = a
                            # number of CPU threads
                            ,nthreads = nt
                            # read group
                            ,byReadGroup = bRG
                            # report assignment result for each read
                            ,reportReads = rR #c(CORE, SAM and BAM)
                            ,reportReadsPath = rRP
                            # miscellaneous
                            ,sampleSheet = s
                            ,cellBarcodeList = cBL
                            ,maxMOp = mM
                            ,tmpDir = tD
                            ,verbose = v
)
# Create a DGEList object from fc
y<-edgeR::DGEList(counts = fc$counts[,1:length(path.BAM.files)],genes = fc$annotation[,c("GeneID","Length")], group = groups)
# export data
require(xlsx)
cat('Creating "counts.txt"...')
write.table(y$counts, file="counts.txt", row.names = TRUE, col.names = TRUE, sep = "\t" )
cat('Done!\n')
cat('Creating "genes.txt"...')
write.table(y$genes, file="genes.txt", row.names = TRUE, col.names = TRUE, sep = "\t" )
cat('Done!\n')
