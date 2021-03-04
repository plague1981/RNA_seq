# Retrive gene data
gene_data<-function(mygene, tar_gene,s1,s2){
  gene.fpkm <- cummeRbund::fpkm(mygene)
  s1.fpkm.tar_gene<-gene.fpkm[gene.fpkm$sample_name==s1,]
  row.names(s1.fpkm.tar_gene)<-tar_gene
  s2.fpkm.tar_gene<-gene.fpkm[gene.fpkm$sample_name==s2,]
  df.tar_gene<-rbind(s1.fpkm.tar_gene,s2.fpkm.tar_gene)
  return(df.tar_gene)
}
# Retrive isoform data
isoform_data<-function(mygene,s1,s2){
  isoform.fpkm <- cummeRbund::fpkm(cummeRbund::isoforms(mygene))
  s1.fpkm.tar_gene.isoform<-isoform.fpkm[isoform.fpkm$sample_name==s1,]
  s2.fpkm.tar_gene.isoform<-isoform.fpkm[isoform.fpkm$sample_name==s2,]
  df.tar_gene.isoform.table<-rbind(s1.fpkm.tar_gene.isoform,s2.fpkm.tar_gene.isoform)
  df.tar_gene.isoform<-split(df.tar_gene.isoform.table, list(df.tar_gene.isoform.table$isoform_id))
  return(df.tar_gene.isoform)
}
# Draw isoform plot
draw_isoform_plot<-function(isoform.data){
  isoform.plots<-	
    ggplot(isoform.data,aes(x=sample_name,y=fpkm))+
    geom_bar(stat ="identity", color="black", width=0.5, position=position_dodge())+
    geom_errorbar(aes(ymin=conf_lo, ymax=conf_hi), width=.2,position=position_dodge(0.5))+
    labs(title=paste0(names(isoform.data)," differetial expression"), x="Group name", y = "FPKM")+
    theme_classic() + scale_fill_manual(values=c('#123456','#345678'))
}

# DESeq2: Instead of genes, we could also extract exons, coding sequences, or transcripts with the same function.
generateCountTable <- function(files, transcripts="TxDb.Hsapiens.UCSC.hg38.knownGene",overlapto="gene") {
  require(transcripts, character.only=TRUE)
  require(GenomicRanges)
  require(Rsamtools)
  require(GenomicAlignments)
  txdb<-transcriptsBy(get(transcripts,
                          envir=.GlobalEnv),
                      overlapto)
  l<- vector("list", length(files))
  for(i in 1:length(files)) {
    alns <- readGappedReads(files[i])
    strand(alns) <- "*"
    hits <- countOverlaps(alns,txdb)
    l[[i]] <- countOverlaps(txdb, alns[hits==1])
    names(l) <- gsub("\\.bam", "", files)
  }
  ct<-as.data.frame(l)
  ct
}
# Conver Gene ID to SYMBOL
getMatrixWithSymbols <- function(df){
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  
  geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(df), column="SYMBOL", keytype="ENTREZID", multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}
# draw a heat map
drawheatmap<-function(DESfile){
  require("RColorBrewer")
  require("gplots")
  DESfile_eset<-ExpressionSet(counts(DESfile, normalized = TRUE))
  sel<-order(rowMeans(counts(DESfile, normalized = TRUE)),decreasing = TRUE)[1:input$n_heatmap]
  hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
  heatmap.2(exprs(DESfile_eset)[sel,],col=hmcol,trace="none",margin=c(10,6))
}
# draw a volcano plot
volcano_plot<-function(res_cases){
  df<-data.frame(res_cases$log2FoldChange,-log(res_cases$pvalue), res_cases$padj, row.names = row.names(res_cases))
  
  # Make a basic volcano plot
  plot(df$res_cases.log2FoldChange, df$X.log.res_cases.pvalue., pch=20, main="Volcano plot", xlab="log2_FoldChange", ylab="-log10_pvalue")
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  df_padj<-subset(df, res_cases.padj<.05)
  points(df_padj$res_cases.log2FoldChange, df_padj$X.log.res_cases.pvalue., pch=20, col="red")
  
  df_log2fc<-subset(df, abs(res_cases.log2FoldChange)>1)
  points(df_log2fc$res_cases.log2FoldChange, df_log2fc$X.log.res_cases.pvalue., pch=20, col="orange")
  
  df_Sig<-subset(df, res_cases.padj<.05 & abs(res_cases.log2FoldChange)>1)
  points(df_Sig$res_cases.log2FoldChange, df_Sig$X.log.res_cases.pvalue., pch=20, col="green")
  
  # Label points with the textxy function from the calibrate plot
  require("calibrate")
  textxy(df_Sig$res_cases.log2FoldChange, df_Sig$X.log.res_cases.pvalue., labs=row.names(df_Sig), cex=.8)
}
# COunt plot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
counts_dotplotly<-function(df,genename){
  p<-ggplot(df, aes(x=condition.., y=count)) + 
    ggtitle(label = genename) +
    geom_count(stat="identity") +
    stat_summary(fun.data=data_summary, geom='errorbar', color="red",position = position_dodge(.9)) 
  return(ggplotly(p))
}
sum_table<-function(df){
  df_mean <- ddply(.data = df, .variables = 'condition..',summarise, length = mean(count))
  df_sd <- ddply(.data = df, .variables = 'condition..', summarise, length = sd(count))
  df_sum <- data.frame(df_mean, df_sd$length)
  df_sum <- rename(df_sum, c('condition'='Group','length'='Mean',"df_sd.length" = "SD"))
  return(df_sum)
}
# get groups info
read_group<-function(input){
  require(readxl)
  file2table<-read_excel(input, col_types="text", range=cell_cols(c("A","B")), col_names = TRUE, sheet = 1)
  return(file2table)
}

# statistics
y_estimate<-function(y){
  if (1 %in% n_occur()[,"Freq"]){
    require(edgeR)
    y_estimate<-estimateGLMCommonDisp(y, method = "deviance", robust=TRUE, subset = NULL)
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
    if (length(levels(group_factors()))==2){
      # ====== Pairwise comparisons between two or more groups (classic)
      # = quantile-adjusted conditional maximum likelihood (qCML)
      # = qCML method is only applicable on datasets with "a single factor" design
      # = This method proves to be accurate and nearly unbiased even for small counts and small numbers of replicates
      require(statmod)
      y_estimate <- estimateDisp(y, robust=TRUE)
      # ===== Alternatives below to get y$pseudo.counts,
      #y <- estimateCommonDisp(y)
      #y <- estimateTagwiseDisp(y)
      #row.names(y$pseudo.counts)<-gene_names
    } else {
      # ====== For More complex experiments (glm functionality)
      y_estimate <- estimateDisp(y, design = design(), robust=TRUE)
      # ===== Alternatives below to get y$pseudo.counts,
      #y <- estimateGLMCommonDisp(y, design)
      #y <- estimateGLMTrendedDisp(y, design)
      #y <- estimateGLMTagwiseDisp(y, design)
    } 
  }
  return(y_estimate)
}
con<-function(ref, exp){
  con<-vector()
  for (n in 1:length(levels(group_factors()))){
    con<-c(con,0)
  }
  con[which(levels(group_factors())==ref)]<- -1
  con[which(levels(group_factors())==exp)]<- 1
  return(con)
}
volcano_plot<-function(x){
  v<-volcanoly(x, col = c("#252525"), point_size = 5, 
            effect_size_line = c(input$logFC_left, input$logFC_right), effect_size_line_color = "grey", effect_size_line_width = 0.5, effect_size_line_type = 2, 
            genomewideline = input$log10P, genomewideline_color = "grey", genomewideline_width = 0.5,
            genomewideline_type = 2, highlight = NULL, highlight_color = "red",
            xlab = NULL, ylab = "-log10(p)", title = paste(input$ref_edgeR, 'vs',input$contrast_edgeR, "Volcano Plot"))
  return(v)
}
volcano_plot_FDR<-function(x){
  v<-volcanoly(x, col = c("#252525"), point_size = 5, 
               effect_size_line = c(input$logFC_left, input$logFC_right), effect_size_line_color = "grey", effect_size_line_width = 0.5, effect_size_line_type = 2, 
               genomewideline = input$FDR, genomewideline_color = "grey", genomewideline_width = 0.5,
               genomewideline_type = 2, highlight = NULL, highlight_color = "red",
               xlab = NULL, ylab = "-log10(FDR)", title = paste(input$ref_edgeR, 'vs',input$contrast_edgeR, "Volcano Plot"))
  return(v)
}
heatmap_plot<-function(r,c){
  m<-data.matrix(cpm.table()[r,c])
  fig <- plot_ly(
    x = c, y = r,
    z = m, type = "heatmap"
  )
  return(fig)
}
