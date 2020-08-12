
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
counts_dotplot<-function(df){
  df$condition..<-as.factor(df$condition..)
  p<-ggplot(df, aes(x=condition.., y=count)) + 
    geom_dotplot(binaxis='y', stackdir='center')
  return(p)
}
