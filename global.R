
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
