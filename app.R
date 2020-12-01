# $sudo apt-get install r-cran-rjava (for rJava package)
# ======== Packages required =========
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16384m"))
# Rcran
packages<-c('shiny','shinythemes','shinydashboard','shinycssloaders','plyr','dplyr','gridExtra','plotly','rapportools',"calibrate",'gplots','rJava','xlsx','readxl')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio.packages<-c('Rsubread','edgeR','DESeq2','cummeRbund')
for (bio.package in bio.packages){
  if(bio.package %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(bio.package)}
}
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(plyr)
library(dplyr)
library(cummeRbund)
library(gridExtra)
library(plotly)
library(DESeq2)
library(rapportools)

gc()
ui<- dashboardPage(
  dashboardHeader(title = 'NGS data analysis'),
  dashboardSidebar(
    sidebarMenu(
      menuItem('CummeRbund',tabName = 'CummeRbund', icon = icon('chart-line')),
      menuItem('DESeq2',tabName = 'DESeq2', icon = icon('chart-line')),
      menuItem('edgeR',tabName = 'edgeR', icon = icon('chart-line'))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = 'CummeRbund', 
              navbarPage(title = 'Analyze data with CummeRbund',
                         tabPanel('Groups Summary', icon = icon('calendar-plus'),
                                  textInput('directory',label = 'Please copy and paste the directory:',value = getwd()),
                                  actionButton('get_summary',label = 'Summarize'),
                                  shiny::tags$hr(),
                                  shiny::tags$p('Info summary'),
                                  tableOutput('summary') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('Dispersion plot'),
                                  plotOutput('dispersion') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csScatterMatrix plot'),
                                  plotOutput('csScatterMatrix') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csVolcanoMatrix plot'),
                                  plotOutput('csVolcanoMatrix') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csDensity plot'),                                      
                                  radioButtons(inputId = 'select_csDensity', label = 'Please select output csDensity plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                  plotOutput('csDensity') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csBoxplot plot'), 
                                  radioButtons(inputId = 'select_csBoxplot', label = 'Please select output csBoxplot plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                  plotOutput('csBoxplot') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csDendro plot'),                                       
                                  radioButtons(inputId = 'select_csDendro', label = 'Please select output csDendro plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                  plotOutput('csDendro') %>% withSpinner(color="#0dc5c1"),
                         ),
                         tabPanel('Groups Expression', icon = icon('calendar-plus'),
                                  textInput(inputId = 'gene_name',label = 'Please enter a gene:',value = 'example:GAPDH'),
                                  actionButton('get_gene',label = 'Analyze'),
                                  shiny::tags$hr(),
                                  shiny::tags$p('Gene Expression in groups'),
                                  plotOutput('genePlots') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  tableOutput('groups_genes_rpkm'),
                                  tableOutput('groups_isoforms_rpkm')
                         ),
                         tabPanel('Pairwise Summary', icon = icon('calendar-plus'),
                                  uiOutput('group1'),
                                  uiOutput('group2'),
                                  actionButton('get_plots',label = 'Get Plots'),
                                  shiny::tags$hr(),
                                  shiny::tags$p('MA_plot'),
                                  plotOutput('MA_plot') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csScatter plot'),
                                  plotOutput('csScatter') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  shiny::tags$p('csVolcano plot'),
                                  plotOutput('csVolcano') %>% withSpinner(color="#0dc5c1")
                         ),
                         tabPanel('Pairwise Expression', icon = icon('calendar-plus'),
                                  textInput('2group_gene_name', label = 'Please enter a gene', value = 'example:GAPDH'),
                                  actionButton('get_2group_gene',label = 'Analyze'),
                                  shiny::tags$hr(),
                                  plotOutput('2group_gene_plot') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  tableOutput('2group_gene_rpkm') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  uiOutput('2group_isoforms_plot') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$hr(),
                                  uiOutput('2group_isoforms_rpkm') %>% withSpinner(color="#0dc5c1")
                         )
              )
      ),
      tabItem(tabName = 'DESeq2',
              navbarPage('Analyze data with DESeq2',
                         tabPanel('bamfiles', icon = icon('calendar-plus'), 
                                  textInput('directory_DESeq2',label = 'Please copy and paste the directory of *.bam files:',value = getwd()),
                                  tableOutput('fileslist')
                         ),
                         tabPanel('Count Table',icon = icon('calendar-plus'),
                                  fileInput(inputId = 'file',label = 'Select raw data file:',multiple = FALSE),
                                  fileInput(inputId = 'condition',label = 'Select condition file:',multiple = FALSE),
                                  actionButton(inputId = 'get_counttable',label = 'Analyze'),
                                  tableOutput('counttable') %>% withSpinner(color="#0dc5c1"),
                                  
                                  actionButton(inputId = 'get_dds',label = 'DESeqDataSet'),
                                  tableOutput('dds') %>% withSpinner(color="#0dc5c1")
                         ),
                         tabPanel('Sum-Plots', icon = icon('map'),
                                  shiny::tags$p('PCA'),
                                  plotlyOutput('plot_PCA') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('DispEsts'),
                                  plotOutput('plot_DispEsts') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('MA'),
                                  plotOutput('plot_MA') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('Sparsity'),
                                  plotOutput('plot_Sparsity') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('Heat map'),
                                  sliderInput(inputId = 'n_heatmap',label = 'How many genes in the heat map?',min = 2,max = 100,value = 10),
                                  plotOutput('heatmap') %>% withSpinner(color="#0dc5c1")
                                  
                         ),
                         tabPanel('Differetial expression', icon = icon('calendar-plus'),
                                  textInput(inputId = 'gene',label = 'Please enter a gene:',value = 'GAPDH'),
                                  actionButton(inputId = 'gene_analyze',label = 'Get data'),
                                  plotlyOutput('DE_plot') %>% withSpinner(color="#0dc5c1"),
                                  tableOutput('DE_table') %>% withSpinner(color="#0dc5c1"),
                                  tableOutput('DE_sum_table') %>% withSpinner(color="#0dc5c1")
                         ),
                         tabPanel('Pairwise',icon = icon('calendar-plus'),
                                  uiOutput('ref_group'),
                                  uiOutput('contrast_groups'),
                                  actionButton(inputId = 'pairwise',label = 'Get data'),
                                  tableOutput('pw_table') %>% withSpinner(color="#0dc5c1")
                         )
              ) #navbarPage: DESeq2
      ), # tabItem:DESeq2
      tabItem(tabName = 'edgeR',
              navbarPage(title = 'Analyze data with edgeR',
                         tabPanel('Raw Data',icon = icon('calendar-plus'),
                                  fileInput(inputId = 'counts_file',label = 'Please select your counts file(eg. counts.txt)',multiple = FALSE),
                                  fileInput(inputId = 'genes_file',label = 'Please select your genes file(eg. genes.txt)',multiple = FALSE),
                                  fileInput(inputId = 'groups_file',label = 'Please select your groups file(eg. groups.xlsx or .txt)',multiple = FALSE),
                                  actionButton(inputId = 'upload_files',label = 'Upload'),
                                  #shiny::tags$p('Groups info'),
                                  #tableOutput('groups_table') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('Counts table preview'),
                                  tableOutput('counts_table') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('Genes table preview'),
                                  tableOutput('genes_table') %>% withSpinner(color="#0dc5c1"),
                                  shiny::tags$p('lib size calibration'),
                                  tableOutput('y_norm_factor') %>% withSpinner(color="#0dc5c1"),
                                  downloadLink("sample_info.download", "Download")
                         ),
                         tabPanel(title = 'Overall analysis', icon = icon('calendar-plus'),
                                  actionButton(inputId = 'get_cpm_table',label = 'cpm table'),
                                  tableOutput('cpm.table') %>% withSpinner(color="#0dc5c1"),
                                  actionButton(inputId = 'get_rpkm_table',label = 'rpkm table'),
                                  tableOutput('rpkm.table') %>% withSpinner(color="#0dc5c1"),
                                  actionButton(inputId = 'get_out_table',label = 'statistic table'),
                                  tableOutput('out_table') %>% withSpinner(color="#0dc5c1"),
                                  downloadLink("overall.download", "Download")
                         ),
                         tabPanel(title = 'Two groups analysis', icon = icon('calendar-plus'),
                                  uiOutput('ref_group_edgeR'),
                                  uiOutput('contrast_group_edgeR'),
                                  actionButton(inputId = 'get_out42groups_table',label = 'statistic table'),
                                  tableOutput('out42groups_table') %>% withSpinner(color="#0dc5c1"),
                                  downloadLink("twogroup.download", "Download")
                         ),
                         tabPanel(title = 'Visualization', icon = icon('map'),
                                  sliderInput(inputId = 'logFC_left',label = 'Meaningful logFC below:',min = -5,max = 5,value = -1, step = 0.25),
                                  sliderInput(inputId = 'logFC_right',label = 'Meaningful logFC above:',min = -5,max = 5,value = 1, step = 0.25),
                                  sliderInput(inputId = 'log10P',label = 'Meaningful -log10PPvalue above:',min = 1,max = 5,value = 2, step = 0.25),
                                  plotlyOutput('volcano_plots') %>% withSpinner(color="#0dc5c1")
                                  
                         )
              ) #navbarPage: edgeR
      ) # tabItem:edgeR
    ) # tabItems
  ) # dashboardBody
) #dashboardPage

server <- function(input, output, session){
  source('global.R', local = TRUE)
  #   
  # variables
  cuff<-eventReactive(input$get_summary,{
    cummeRbund::readCufflinks(input$directory)
  })
  mygene<-eventReactive(input$get_gene,{
    cummeRbund::getGene(cuff(),input$gene_name)
  })
  
  groups<-reactive({
    if (is.null(cuff())){
      return(NULL)
    } else
      Info<-cummeRbund::replicates(cuff())
    return(levels(factor(Info$sample_name)))
  })
  
  output$summary<-renderTable({
    cuff()
    Info<-cummeRbund::replicates(cuff())
    return(Info)
    #grid.table(Info, theme=ttheme_minimal(base_size = 16))
  })
  output$dispersion<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else 
      cummeRbund::dispersionPlot(genes(cuff()))
  })
  output$csScatterMatrix<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else 
      cummeRbund::csScatterMatrix(genes(cuff()))
  })
  output$csVolcanoMatrix<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else 
      cummeRbund::csVolcanoMatrix(genes(cuff()))
  })
  output$csDensity<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else if (input$select_csDensity=='Group'){
      cummeRbund::csDensity(genes(cuff()))
    } else if (input$select_csDensity=='Individual'){
      cummeRbund::csDensity(genes(cuff()),replicates=T)
    }
  })
  output$csBoxplot<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else if (input$select_csBoxplot=='Group'){
      cummeRbund::csBoxplot(genes(cuff()))
    } else if (input$select_csBoxplot=='Individual'){
      cummeRbund::csBoxplot(genes(cuff()),replicates=T)
    }
  })
  output$csDendro<-renderPlot({
    if (is.null(cuff())){
      return(NULL)
    } else if (input$select_csDendro=='Group'){
      cummeRbund::csDendro(genes(cuff()))
    } else if (input$select_csDendro=='Individual'){
      cummeRbund::csDendro(genes(cuff()),replicates=T)
    }
  })
  # tabPanel groups gene and isoforms analysis
  output$groups_genes_rpkm<-renderTable({
    if (is.null(mygene())){
      return(NULL)
    } else
      cummeRbund::fpkm(mygene())
  })
  output$groups_isoforms_rpkm<-renderTable({
    if (is.null(mygene())){
      return(NULL)
    } else
      cummeRbund::fpkm(isoforms(mygene()))
  })
  output$genePlots<-renderPlot({
    if (is.null(mygene())){
      return(NULL)
    } else
      g1 <- expressionBarplot(mygene())
    g2 <- expressionBarplot(cummeRbund::isoforms(mygene()))
    grid.arrange(g1, g2, nrow=1)
  })
  # tabPanel Pairwise summary
  output$group1<-renderUI({
    group_name_list<-groups()
    selectInput(inputId = 'group1', label = 'Please select first group:', choices = group_name_list ,selected = group_name_list[1])
  })
  output$group2<-renderUI({
    group_name_list<-groups()
    selectInput(inputId = 'group2', label = 'Please select second group:', choices = group_name_list ,selected = group_name_list[2])
  })
  
  MA_plot<-eventReactive(input$get_plots,{
    cummeRbund::MAplot(cummeRbund::genes(cuff()), input$group1, input$group2)
  })
  csScatter_plot<-eventReactive(input$get_plots,{
    cummeRbund::csScatter(cummeRbund::genes(cuff()), input$group1, input$group2)
  })
  csVolcano_plot<-eventReactive(input$get_plots,{
    cummeRbund::csVolcano(cummeRbund::genes(cuff()), input$group1, input$group2)
  })
  output$MA_plot<-renderPlot({
    if (is.null(MA_plot())){
      return(NULL)
    } else
      MA_plot()
  })
  output$csScatter<-renderPlot({
    if (is.null(csScatter_plot())){
      return(NULL)
    } else
      csScatter_plot()
  })
  output$csVolcano<-renderPlot({
    if (is.null(csVolcano_plot())){
      return(NULL)
    } else
      csVolcano_plot()
  })
  # tabPanel Pairwise gene and isoforms analysis
  mygene_2group<-eventReactive(input$get_2group_gene,{
    cummeRbund::getGene(cuff(),input$`2group_gene_name`)
  })
  df.tar_gene<-eventReactive(input$get_2group_gene,{
    gene_data(mygene_2group(), input$`2group_gene_name`,input$group1, input$group2)
  })
  df.tar_gene.isoform<-eventReactive(input$get_2group_gene,{
    isoform_data(mygene_2group(),input$group1, input$group2)
  })
  output$`2group_gene_plot`<-renderPlot({
    if (is.null(mygene_2group())){
      return(NULL)
    } else
      g <- ggplot(df.tar_gene(), aes(x=sample_name,y=fpkm))+
        geom_bar(stat = "identity", color="black", width=0.5, position=position_dodge())+
        geom_errorbar(aes(ymin=conf_lo, ymax=conf_hi), width=.2,position=position_dodge(0.5))+
        labs(title=paste(input$`2group_gene_name`, "differetial expression"), x="Group name", y = "FPKM")+theme_classic() + 
        scale_fill_manual(values=c('#123456','#345678'))
    return(g)
  })
  output$`2group_gene_rpkm`<-renderTable({
    if (is.null(mygene_2group())){
      return(NULL)
    } else
      df.tar_gene()
  })
  output$`2group_isoforms_plot`<-renderUI({
    if (is.null(df.tar_gene.isoform())){
      return(NULL)
    } else
      plot_output_list <- lapply(1:length(df.tar_gene.isoform()), function(i) {
        plotname <- paste0("plot", i)
        plot_output_object <- plotlyOutput(plotname)
        plot_output_object <- renderPlotly({
          isoform.table.column.name<-c("isoform_id","sample_name","fpkm","conf_hi","conf_lo","quant_status")
          isoform.data<-data.frame(df.tar_gene.isoform()[i])
          colnames(isoform.data)<-isoform.table.column.name
          draw_isoform_plot(isoform.data)
        })
      })
    return(plot_output_list)
  })
  output$`2group_isoforms_rpkm`<-renderUI({
    if (is.null(df.tar_gene.isoform())){
      return(NULL)
    } else
      table_output_list <- lapply(1:length(df.tar_gene.isoform()), function(i) {
        tablename <- paste0("table", i)
        table_output_object <- plotlyOutput(tablename)
        table_output_object <- renderTable({
          isoform.table.column.name<-c("isoform_id","sample_name","fpkm","conf_hi","conf_lo","quant_status")
          isoform.data<-data.frame(df.tar_gene.isoform()[i])
          colnames(isoform.data)<-isoform.table.column.name
          isoform.data
        })
      })
    return(table_output_list)
  })
  # DESeq2
  bamfiles<-reactive({
    bamfiles<-list.files(input$directory_DESeq2, pattern = ".bam$")
    bamfiles<-file.path(input$directory_DESeq2, bamfiles )
  })
  
  counttable<-eventReactive(input$get_counttable,{
    if (is.null(bamfiles()) & is.null(input$file)){
      return(NULL)
    } else if (!is.empty(bamfiles())){
      return(generateCountTable(bamfiles()))
    } else if (!is.null(input$file)){
      counttable<-read.csv(file = input$file$datapath ,sep = '\t')
      return(counttable)
    }
  })
  condition<-reactive({
    if (is.null(input$condition)){
      return(NULL)
    } else
      grouptable<-read.csv(file = input$condition$datapath ,sep = '\t')
    return(unlist(matrix(grouptable[1,])))
  })
  
  des<-reactive({
    if (is.null(counttable())|is.null(condition())){
      return(NULL)
    } else
      des<-DESeqDataSetFromMatrix(counttable(), DataFrame(condition()), ~ condition())
    des<-estimateSizeFactors(des)
    return(des)
  })
  M1symb<-reactive({
    if (is.null(des())){
      return(NULL)
    } else
      return(getMatrixWithSymbols(des()))
  })
  # DESeqResults extration
  output$ref_group<-renderUI({
    if (is.null(condition())){
      return(NULL) 
    } else
      group_choices<-row.names(table(condition()))
    #group_choices<-row.names(table(M1symb()$groups))
    selectInput(inputId = 'ref_DESeq2', label = 'Please select the reference group:', choices = group_choices,selected = group_choices[1])
  })
  output$contrast_groups<-renderUI({
    if (is.null(condition())){
      return(NULL) 
    } else
      group_choices<-row.names(table(condition()))
    selectInput(inputId = 'contrast_DESeq2', label = 'Please select the contrast group:',choices = group_choices,selected = group_choices[2])
  })
  dds<-eventReactive(input$get_dds,{
    if (is.null(M1symb())){
      return(NULL)
    } else
      dds<-DESeq(M1symb()) # Convert to DESeqDataSet format
    #dds$condtion.. <- factor(condition())
    dds$condition..<-as.factor(dds$condition..)
    #dds$condition..<-relevel(dds$condition.., ref = input$ref_DESeq2)
    design(dds) <- ~ condition.. + 0 # they still have the _vs_ in the names
    dds<-DESeq(dds, betaPrior=FALSE)
    dds<-nbinomWaldTest(dds) # test for significance of change in deviance. c(nbinomLRT(),nbinomWaldTest())
    return(dds)
  })
  
  output$dds<-renderTable(rownames = TRUE,{
    if (is.null(dds())){
      return(NULL)
    } else 
      return(data.frame(resultsNames(dds())))
  })
  
  # Getting Significant data (padj<0.1)
  resSig<-reactive({
    if (is.null(res_cases())){
      return(NULL)
    } else
      return(res_cases[res_cases$padj<0.1,])
  })
  
  output$fileslist<-renderTable({
    if (is.null(bamfiles())){
      return(NULL)
    } else
      return(row.names(table(bamfiles())))
  })
  output$counttable<-renderTable(rownames = TRUE,spacing = 'xs',{
    if (is.null(counttable())){
      return(NULL)
    } else
      return(head(counts(M1symb(), normalized = TRUE)))
  })
  output$plot_PCA<-renderPlotly({
    if (is.null(dds())){
      return(NULL)
    } else
      vstcounts <- varianceStabilizingTransformation(dds(), blind=TRUE)
    g<-plotPCA(vstcounts, intgroup=head(colnames(colData(dds())),-1))
    return(ggplotly(g))
  })
  output$plot_DispEsts<-renderPlot({
    if (is.null(dds())){
      return(NULL)
    } else
      return(plotDispEsts(dds(), genecol="black", fitcol="red",finalcol="dodgerblue",cex=.2))
  })
  
  output$plot_MA<-renderPlot({
    if (is.null(dds())){
      return(NULL)
    } else
      return(plotMA(dds(),alpha=0.05, xlab="mean of normalized count"))
  })
  
  output$plot_Sparsity<-renderPlot({
    if (is.null(dds())){
      return(NULL)
    } else
      return(plotSparsity(dds(),normalized = TRUE))
  })
  output$heatmap<-renderPlot({
    if (is.null(dds())){
      return(NULL)
    } else
      return(drawheatmap(dds()))
  })
  # Differential expression section
  DE_data<-eventReactive(input$gene_analyze,{
    return(plotCounts(dds(), gene = input$gene, intgroup = 'condition..',returnData = TRUE,normalized = TRUE))
  })
  output$DE_plot<-renderPlotly({
    if (is.null(DE_data())){
      return(NULL)
    } else
      counts_dotplotly(DE_data(), input$gene)
  })
  output$DE_table<-renderTable(colnames = TRUE,rownames = TRUE, {
    if (is.null(DE_data())){
      return(NULL)
    } else
      DE_data()
  })
  output$DE_sum_table<-renderTable(colnames = TRUE,rownames = TRUE, {
    if (is.null(DE_data())){
      return(NULL)
    } else
      return(sum_table(DE_data()))
  })
  # Pairwise section
  
  # Remove NA containing rows
  res_cases<-eventReactive(input$pairwise,{
    contrast_list<-NULL
    for (n in 1:length(row.names(table(condition)))){
      if (row.names(table(condition))[n]==input$ref_DESeq2){
        contrast_list<- c(contrast_list,1)
      } else if (row.names(table(condition))[n]==input$contrast_DESeq2){
        contrast_list<- c(contrast_list,-1)
      } else
        contrast_list<-c(contrast_list,0)
    }
    res<-results(dds(), contrast = contrast_list)
    return(res[complete.cases(res),])
  })
  output$pw_table<-renderTable(rownames = TRUE,{
    #return(resultsNames(dds()))
    if (is.null(res_cases())){
      return(NULL)
    } else
      return(head(data.frame(res_cases())))
  })
  # ===== edgeR
  # variants
  readgroup<- eventReactive(input$upload_files,{
    read_group(input$groups_file$datapath)
  })
  group_nrow<- reactive({
    nrow(readgroup())
  })
  group_factors<-reactive({
    factor(unlist(readgroup()$Groups))
  })
  readcounts<- eventReactive(input$upload_files,{
    read.table(input$counts_file$datapath,row.names = 1)
  })
  readgenes<- eventReactive(input$upload_files,{
    read.table(input$genes_file$datapath,row.names = 1)
  })
  y<- eventReactive(input$upload_files,{
    require(edgeR)
    y<-DGEList(counts = data.frame(readcounts()), genes = data.frame(readgenes()), group = group_factors())
    keep.exprs <- filterByExpr(y, group=group_factors())
    y <- y[keep.exprs,, keep.lib.sizes=FALSE]
    y<-calcNormFactors(y)
    return(y)
  })
  design <- eventReactive(input$upload_files,{
    design <-model.matrix(~ 0 + group_factors(), data = y()$samples)
    colnames(design) <- levels(group_factors())
    return(design)
  })
  n_occur <- reactive({
    data.frame(table(group_factors()))
  })
  y_estimate_result<-reactive({
    return(y_estimate(y()))
  })
  cpm.table<-eventReactive(input$get_cpm_table,{
    return(cpm(y_estimate_result()))
  })
  rpkm.table<-eventReactive(input$get_rpkm_table,{
    return(rpkm(y_estimate_result()))
  })
  out<-eventReactive(input$get_out_table,{
    fit <- glmQLFit(y_estimate_result(), design = design(), robust = TRUE)
    anova.like.result <- glmQLFTest(fit, coef=2:length(levels(group_factors())))
    out <- topTags(anova.like.result, n = "Inf")$table
    return(out)
  })
  con_seq<-reactive({
    return(con(input$ref_edgeR, input$contrast_edgeR))
  })
  out42groups<-eventReactive(input$get_out42groups_table,{
    if (length(levels(group_factors()))==2){
      # the exact test is only applicable to experiments with a single factor
      test.result <- exactTest(y_estimate_result(), pair = c(input$ref_edgeR,input$contrast_edgeR) )
    } else {
      # Genewise Negative Binomial Generalized Linear Model with Quasi-likelihood
      fit <- glmQLFit(y_estimate_result(), design = design()) #"DGEGLM"
      #con<-c(-1,0,0,1)
      test.result <- glmQLFTest(fit, contrast = con_seq()) # "DGELRT"
      ##Alternative method, Genewise Negative Binomial Generalized Linear Model
      ##Fit a negative binomial generalized log-linear model to the read counts for each gene
      #fit <- glmFit(y,design = design) #"DGEGLM"
      #result <- glmLRT(fit,contrast = con) # "DGELRT" 
    }
    # export the statistic results
    out <- topTags(test.result, n = "Inf")$table
    return(out)
  })
  
  # Output
  #output$groups_table<-renderTable({
  #  readgroup()
  #})
  output$counts_table<-renderTable(rownames = TRUE,{
    if (is.null(readcounts())){
      return(NULL)
    } else
      utils::head(readcounts())
  })
  output$genes_table<-renderTable(rownames = TRUE,{
    if (is.null(readgenes())){
      return(NULL)
    } else
      utils::head(readgenes())
  })
  output$y_norm_factor<-renderTable(rownames = TRUE,{
    if (is.null(y())){
      return(NULL)
    } else
      y()$sample
  })
  output$sample_info.download <- downloadHandler(
    filename =  "norm_factor.xlsx"
    ,
    content = function(file) {
      require(xlsx)
      write.xlsx(y()$sample, sheetName = 'norm_factor',file)
    })
  output$cpm.table<-renderTable(rownames = TRUE,digits = 6, spacing = 'xs',{
    if (is.null(cpm.table())){
      return(NULL)
    } else
      utils::head(cpm.table())
  })
  output$rpkm.table<-renderTable(rownames = TRUE,digits = 6, spacing = 'xs',{
    if (is.null(rpkm.table())){
      return(NULL)
    } else
      utils::head(rpkm.table())
  })
  output$out_table<-renderTable(rownames = TRUE,digits = 6, spacing = 'xs',{
    if (is.null(out())){
      return(NULL)
    } else
      utils::head(data.frame(out()))
  })
  output$overall.download <- downloadHandler(
    filename =  "result_tables.xlsx"
    ,
    content = function(file) {
      require(xlsx)
      write.xlsx(cpm.table(), sheetName = 'cpm',file)
      write.xlsx(rpkm.table(), sheetName = 'rpkm',file, append = TRUE)
      write.xlsx(out(), sheetName = 'out',file, append = TRUE)
    })
  # two groups comparison
  output$ref_group_edgeR<-renderUI({
    if (is.null(group_factors())){
      return(NULL) 
    } else
      group_choices<-levels(group_factors())
    selectInput(inputId = 'ref_edgeR', label = 'Please select the reference group:', choices = group_choices,selected = group_choices[1])
  })
  output$contrast_group_edgeR<-renderUI({
    if (is.null(group_factors())){
      return(NULL) 
    } else
      group_choices<-levels(group_factors())
    selectInput(inputId = 'contrast_edgeR', label = 'Please select the contrast group:',choices = group_choices,selected = group_choices[2])
  })
  output$out42groups_table<-renderTable(rownames = TRUE,digits = 6, spacing = 'xs',{
    if (is.null(out())){
      return(NULL)
    } else
      utils::head(data.frame(out42groups()))
  })
  output$twogroup.download <- downloadHandler(
    filename =  "2groups_result_table.xlsx"
    ,
    content = function(file) {
      require(xlsx)
      write.xlsx(data.frame(out42groups()), sheetName = 'out',file)
  })
  output$volcano_plots<-renderPlotly({
    require(manhattanly)
    v_data<-data.frame(out42groups()[,c('GeneID','logFC','PValue')])
    v_data[,'logFC']<-as.numeric(v_data[,'logFC'])
    v_data[,'PValue']<-as.numeric(v_data[,'PValue'])
    x<-volcanor(v_data,p = 'PValue',effect_size = 'logFC',snp = 'GeneID')
    x$effectName<-'EFFECTSIZE'
    return(volcano_plot(x))
  })
  
}
shinyApp(ui, server)
