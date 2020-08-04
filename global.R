# Retrive gene data
gene_data<-function(mygene, tar_gene,s1,s2){
  gene.fpkm <- fpkm(mygene)
  s1.fpkm.tar_gene<-gene.fpkm[gene.fpkm$sample_name==s1,]
  row.names(s1.fpkm.tar_gene)<-tar_gene
  s2.fpkm.tar_gene<-gene.fpkm[gene.fpkm$sample_name==s2,]
  df.tar_gene<-rbind(s1.fpkm.tar_gene,s2.fpkm.tar_gene)
  return(df.tar_gene)
}
# Retrive isoform data
isoform_data<-function(mygene,s1,s2){
  isoform.fpkm <- fpkm(isoforms(mygene))
  s1.fpkm.tar_gene.isoform<-isoform.fpkm[isoform.fpkm$sample_name==s1,]
  s2.fpkm.tar_gene.isoform<-isoform.fpkm[isoform.fpkm$sample_name==s2,]
  df.tar_gene.isoform.table<-rbind(s1.fpkm.tar_gene.isoform,s2.fpkm.tar_gene.isoform)
  df.tar_gene.isoform<-split(df.tar_gene.isoform.table, list(df.tar_gene.isoform.table$isoform_id))
  isoform.table.column.name<-c("isoform_id","sample_name","fpkm","conf_hi","conf_lo","quant_status")
  
  for (n in 1:length(df.tar_gene.isoform)){
    assign(paste0("isoform.table_",n),df.tar_gene.isoform[n])
    write(paste0('isoform.table_',n,'<-data.frame(isoform.table_',n,')'),"tmp.R")
    write(paste0('colnames(isoform.table_',n,')<-isoform.table.column.name'),"tmp.R", append = TRUE)
    source("tmp.R", local = TRUE)
    file.remove("tmp.R")
  }
  
}
