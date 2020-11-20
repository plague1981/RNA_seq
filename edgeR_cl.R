# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
# Rcran
packages<-c('Xmisc','readxl','xlsx','statmod','stringr','rJava')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio.packages<-c('Rsubread','edgeR',"topGO")
for (bio.package in bio.packages){
  if(bio.package %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(bio.package)}
}
# === setting environment ===
parser <- Xmisc::ArgumentParser$new()
parser$add_usage('edgeR_cl.R [options]')
parser$add_description('An executable R script parsing arguments from Unix-like command line.')
parser$add_argument('--h',type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--dir', type = 'character', default = getwd(), help = '"directory",Enter your working directory')
parser$helpme()
# === variables ====
#dir<-'C:/Users/Changyi.Lin/Desktop/Vincent/merged'
dirPath <- dir
dirPath <-gsub ('\\\\','/',dirPath)
if (dir.exists(dirPath)){
  setwd(dirPath)
  cat(paste0("Setting ",dirPath," as the working directory\n"))
} else {
  cat("Directory is not existing!\n")
  quit()
}

# ===== Import group names from groups.txt or groups.xlsx
cat("Importing files needed (eg. groups.xlsx/groups.txt, counts.txt and genes.txt)\n")

if ((file.exists("groups.xlsx")|file.exists("groups.txt"))&& file.exists('counts.txt') && file.exists('genes.txt')){
  suppressPackageStartupMessages(require(readxl))
  cat("All files were found!\n")
  groups <- factor(unlist(read_excel("groups.xlsx",range=cell_cols("B"))))
  groups_name <- unlist(read_excel("groups.xlsx",range=cell_cols("A")))
  groups_name<-gsub('-','.',groups_name)
  names(groups)<-groups_name
  raw.data.counts <- read.delim("counts.txt")
  raw.data.genes <- read.delim("genes.txt")
  colnames(raw.data.counts)
  group_order<-NULL
  for (n in 1:length(colnames(raw.data.counts))){
    if (!is.element(colnames(raw.data.counts)[n], names(groups))){
      cat(paste0(colnames(raw.data.counts)[n], ' dose not exist in groups\n'))
      cat(paste0('Please check sample names and ',gf, ' \n'))
      quit()
    }
    for (m in 1:length(names(groups))){
      if (colnames(raw.data.counts)[n]==names(groups)[m]){
        group_order<-c(group_order,m)
        break()
      } 
    }
  }
  groups<-groups[group_order]
} else {
  cat("At least one file is missing\n")
  quit()
}

# ===== Options before analysis =========================
# Need of converting gene id to symbol
repeat{
  cat("Do you need to convert gene_id to symbol?Yes or No:\n")
  convert_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", convert_inp, 'for converting gene_id to symbol\n'))
  convert_inp<-tolower(convert_inp)
  if (convert_inp=="y"|convert_inp=="yes"|convert_inp=="n"|convert_inp=="no"){
    break() 
  } else {cat("Please enter yes or no [y/n]!")}
}
if (convert_inp=="y"|convert_inp=="yes"){
  repeat{
    cat("Which database are you using?\n")
    cat('1. NCBI\n')
    cat('2. Ensembl\n')
    db_inp <- readLines("stdin", n = 1L)
    if (db_inp==1){
      db_inp<-'NCBI'
      db_keytype<-"ENTREZID"
      break()} 
    else if (db_inp==2){
      db_inp<-'Ensembl'
      db_keytype<-"ENSEMBL"
      break()} else {cat("Please enter the number![1-2]\n")}
  }
  cat(paste("you entered", db_inp, 'as the database\n'))
}
# Sample info
Symbol.packages<-c("org.Hs.eg.db","org.Mm.eg.db",'org.Rn.eg.db',"org.EcK12.eg.db")
s_species<-c('Human','Mouse','Rat','E. coli strain K12')
repeat{
  cat("What kind species are your samples:\n")
  cat('1. Human\n')
  cat('2. Mouse\n')
  cat('3. Rat\n')
  cat('4. E. coli strain K12\n')
  species_inp <- readLines("stdin", n = 1L)
  if(as.numeric(species_inp) %in% c(1:4)){
    suppressPackageStartupMessages(require(stringr))
    pattern <- "org\\.\\s*(.*?)\\s*\\.eg\\.db"
    sample_species<- regmatches(Symbol.packages[1], regexec(pattern, Symbol.packages[1]))[[1]][2]
    break()
  } else {cat("Please enter the number![1-4]\n")}
}

if(Symbol.packages[as.numeric(species_inp)] %in% rownames(installed.packages()) == FALSE) {
  BiocManager::install(Symbol.packages[as.numeric(species_inp)])}
if (species_inp==1){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  database<-org.Hs.eg.db
} else if (species_inp==2){
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  database<-org.Mm.eg.db
} else if (species_inp==3){
  suppressPackageStartupMessages(require(org.Rn.eg.db))
  database<-org.Rn.eg.db
} else if (species_inp==4){
  suppressPackageStartupMessages(require(org.EcK12.eg.db))
  database<-org.EcK12.eg.db  
}
cat(paste("you entered", s_species[species_inp], 'as the species\n'))

# Statistics analysis
repeat{
  cat("Single factor for statistics? Yes or No:\n")
  factor_inp <- readLines("stdin", n = 1L)
  cat(paste("you entered", factor_inp, 'for single facotr statistics',"\n"))
  factor_answer<-tolower(factor_inp)
  if (factor_answer=="y"|factor_answer=="yes"|factor_answer=="n"|factor_answer=="no"){
    break()  
  } else {cat("Please enter yes or no [y/n]!\n")}
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
      cat(paste0("Please enter the 1st group number [1-",length(groups_list),"]:"))
      n<-readLines("stdin", n = 1L)
      cat(paste0("Please enter the 2nd group number [1-",length(groups_list),"]:"))
      m<-readLines("stdin", n = 1L)
      n<-as.integer(n)
      m<-as.integer(m)
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
if ((overall_answer=="n"|overall_answer=="no")&&(compare_inp=="n"|compare_inp=="no")){
  cat('No files were created!...')
  quit()
}
# ==========================================================
# Processing data if you don't have counts.txt and genes.txt 

cat("Analyzing data from counts.txt and genes.txt\n")
# Create a DGEList object from counts.txt and genes.txt
suppressPackageStartupMessages(require(edgeR))
y <- DGEList(counts = raw.data.counts, genes = data.frame(raw.data.genes), group = groups)
table(rowSums(y$counts==0)==12)
# remove worthless genes
cat('Removing worthless genes by "filterByExpr"...')
keep.exprs <- filterByExpr(y, group=groups)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
cat('Done!\n')
# calculate normal factors
cat('Calculating normal factors!...')
y<-calcNormFactors(y)
cat('Done!\n')
# edit design
design <- model.matrix(~ 0 + groups, data = y$samples)
colnames(design) <- levels(groups)

# ====== Statistics (edgeR) ==================
n_occur <- data.frame(table(groups))
if (1 %in% n_occur[,"Freq"]){
  cat("One or more groups have No Duplicates\n")
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
  if (Xmisc::character_to_logical(factor_answer,ignore.case = TRUE)){
    # ====== Pairwise comparisons between two or more groups (classic)
    # = quantile-adjusted conditional maximum likelihood (qCML)
    # = qCML method is only applicable on datasets with "a single factor" design
    # = This method proves to be accurate and nearly unbiased even for small counts and small numbers of replicates
    suppressPackageStartupMessages(require(statmod))
    y <- estimateDisp(y, robust=TRUE)
    # ===== Alternatives below to get y$pseudo.counts,
    #y <- estimateCommonDisp(y)
    #y <- estimateTagwiseDisp(y)
    #row.names(y$pseudo.counts)<-gene_names
  } else if (!Xmisc::character_to_logical(factor_answer,ignore.case = TRUE)) {
    # ====== For More complex experiments (glm functionality)
    y <- estimateDisp(y, design = design, robust=TRUE)
    # ===== Alternatives below to get y$pseudo.counts,
    #y <- estimateGLMCommonDisp(y, design)
    #y <- estimateGLMTrendedDisp(y, design)
    #y <- estimateGLMTagwiseDisp(y, design)
  } 
}
# Add gene symbols and length if needed
# ADD Symbols

if (convert_inp=="y"|convert_inp=="yes"){
    # NCBI reference
    if (suppressWarnings(!is.na(as.numeric(row.names(y)[1])))){
      Symbol<- mapIds(database,keys=rownames(y), keytype=db_keytype, column="SYMBOL")
      gene_id_names <- data.frame(Symbol=Symbol)
    # Convert GeneID to Symbols
      gene_names<-NULL
      for (gene_id in as.character(row.names(y))) {
        if (is.na(gene_id_names[gene_id,"Symbol"])) {
         gene_names<-c(gene_names, gene_id)
       } else {
         gene_name<-as.character(gene_id_names[gene_id,"Symbol"])
         gene_names<-c(gene_names, gene_name)
       }
      }
    row.names(y)<-gene_names
    }
}

# ======= Overall Differential expression analysis ======
# Genewise Negative Binomial Generalized Linear Model with Quasi-likelihood

if (overall_answer=="y"|overall_answer=="yes"){
  fit <- glmQLFit(y, design = design, robust = TRUE)
  anova.like.result <- glmQLFTest(fit, coef=2:length(levels(groups)))
  out <- topTags(anova.like.result, n = "Inf")$table
  total.cpm.table<-cpm(y)
  total.rpkm.table<-rpkm(y)
  for (n in 3:ncol(out)){
    total.rpkm.table<-cbind(total.rpkm.table, out[,n][match(rownames(total.rpkm.table), rownames(out))])
    total.cpm.table<-cbind(total.cpm.table, out[,n][match(rownames(total.cpm.table), rownames(out))])
  }
  
  sample.names.rpkm<-colnames(rpkm(y))
  sample.names.cpm<-colnames(cpm(y))
  colnames(total.rpkm.table)<-c(sample.names.rpkm,colnames(out[,3:ncol(out)]))
  colnames(total.cpm.table)<-c(sample.names.cpm,colnames(out[,3:ncol(out)]))
  suppressPackageStartupMessages(require(xlsx))
  cat('Creating "total.genes.rpkm.table.xlsx"...')
  write.xlsx(total.rpkm.table, "total.genes.rpkm.table.xlsx")
  cat("Done!\n")
  cat('Creating "total.genes.cpm.table.xlsx"...')
  write.xlsx(total.cpm.table, "total.genes.cpm.table.xlsx")
  cat('Done!\n')
  
} else if (overall_answer=="n"|overall_answer=="no"){
  cat("No Overall file was created!\n")
} 

# ======== Two groups Differential expression =========
# Select two groups for comparison

if (compare_inp=="y"|compare_inp=="yes"){
  if (factor_answer=="y"|factor_answer=="yes"){
  # the exact test is only applicable to experiments with a single factor
    test.result <- exactTest(y, pair = c(group_1,group_2) )
  } else {
    # Genewise Negative Binomial Generalized Linear Model with Quasi-likelihood
    fit <- glmQLFit(y, design = design) #"DGEGLM"
    #con<-c(-1,0,0,1)
    test.result <- glmQLFTest(fit, contrast = con) # "DGELRT"
    ##Alternative method, Genewise Negative Binomial Generalized Linear Model
    ##Fit a negative binomial generalized log-linear model to the read counts for each gene
    #fit <- glmFit(y,design = design) #"DGEGLM"
    #result <- glmLRT(fit,contrast = con) # "DGELRT" 
  }
  # export the statistic results
  out <- topTags(test.result, n = "Inf")$table
  
  # Get the sample names in two groups
  sample.names.group_1<-row.names(y$samples[y$samples[,"group"]==group_1,])
  sample.names.group_2<-row.names(y$samples[y$samples[,"group"]==group_2,])
  sample.names<-c(sample.names.group_1,sample.names.group_2)

  # Get the rpkm and cpm of samples in two groups
  rpkm.table<-rpkm(y)[,sample.names]
  cpm.table<-cpm(y)[,sample.names]
  # Combine rpkm table and statistic results
  out.col<-c("GeneID","Length","logFC","logCPM","PValue","FDR")
  for (coln in out.col){
    rpkm.table<-cbind(rpkm.table, out[,coln][match(rownames(rpkm.table), rownames(out))])
    cpm.table<-cbind(cpm.table, out[,coln][match(rownames(cpm.table), rownames(out))])
  }
  colnames(rpkm.table)<-c(sample.names,out.col)
  colnames(cpm.table)<-c(sample.names,out.col)
  rpkm.table<-rpkm.table[,c("GeneID","Length",sample.names,"logFC","logCPM","PValue","FDR")]
  cpm.table<-cpm.table[,c("GeneID","Length",sample.names,"logFC","logCPM","PValue","FDR")]

  # Export the results as an xlsx file.
  suppressPackageStartupMessages(require(xlsx))
  cat(paste0('Creating ',group_1,group_2,"_genes.rpkm.table.xlsx..."))
  write.xlsx(rpkm.table, paste0(group_1,group_2,"_genes.rpkm.table.xlsx"))
  cat('Done!\n')
  cat(paste0('Creating ',group_1,group_2,"_genes.cpm.table.xlsx..."))
  write.xlsx(cpm.table, paste0(group_1,group_2,"_genes.cpkm.table.xlsx"))
  cat('Done!')
  # Intepret the differential expression results
  pdf(file="LogFCvsPvalue.pdf")
  plot(out[,"logFC"],-log2(out[,"PValue"]))
  dev.off()
}
