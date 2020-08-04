library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(cummeRbund)
library(gridExtra)

ui<- dashboardPage(
      dashboardHeader(title = 'NGS data analysis'),
      dashboardSidebar(
        sidebarMenu(
          menuItem('CummeRbund',tabName = 'CummeRbund', icon = icon('chart-line'))
        )
      ),
      dashboardBody(
        tabItems(
          tabItem(tabName = 'CummeRbund', 
                  navbarPage(title = 'Analyze data with CummeRbund',
                             tabPanel('Groups Summary', icon = icon('calendar-plus'),
                                      textInput('directory',label = 'Please copy and paste the directory:',value = getwd()),
                                      uiOutput('groups_choices'),
                                      actionButton('get_summary',label = 'Summarize'),
                                      shiny::tags$hr(),
                                      shiny::tags$p('Info summary'),
                                      plotOutput('summary') %>% withSpinner(color="#0dc5c1"),
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
                             tabPanel('2Group Summary', icon = icon('calendar-plus'),
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
                             tabPanel('2Group Expression', icon = icon('calendar-plus'),
                                      textInput('2group_gene_name', label = 'Please enter a gene', value = 'example:GAPDH'),
                                      actionButton('get_2group_gene',label = 'Analyze'),
                                      shiny::tags$hr(),
                                      plotOutput('2group_gene_plot'),
                                      shiny::tags$hr(),
                                      tableOutput('2group_gene_rpkm'),
                                      uiOutput('2group_isoforms_rpkm')
                             )
                            )
                  
                  )
        )
      )
)

server <- function(input, output, session){
source('global.R', local = TRUE)
  
# variables
cuff<-eventReactive(input$get_summary,{
  readCufflinks(input$directory)
})
mygene<-eventReactive(input$get_gene,{
  getGene(cuff(),input$gene_name)
})

groups<-eventReactive(input$get_summary,{
  Info<-replicates(cuff())
  return(levels(factor(Info$sample_name)))
})

output$groups_choices<- renderUI({
  checkboxGroupInput(inputId = "groups_name", label = "Please select groups you want to compare with:", choices = groups(), selected = groups())
})

output$summary<-renderPlot({
  cuff()
  Info<-replicates(cuff())
  grid.table(Info, theme=ttheme_minimal(base_size = 16))
})
output$dispersion<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else 
      dispersionPlot(genes(cuff()))
})
output$csScatterMatrix<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else 
    csScatterMatrix(genes(cuff()))
})
output$csVolcanoMatrix<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else 
    csVolcanoMatrix(genes(cuff()))
})
output$csDensity<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else if (input$select_csDensity=='Group'){
    csDensity(genes(cuff()))
  } else if (input$select_csDensity=='Individual'){
    csDensity(genes(cuff()),replicates=T)
  }
})
output$csBoxplot<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else if (input$select_csBoxplot=='Group'){
    csBoxplot(genes(cuff()))
  } else if (input$select_csBoxplot=='Individual'){
    csBoxplot(genes(cuff()),replicates=T)
  }
})
output$csDendro<-renderPlot({
  if (is.null(cuff())){
    return(NULL)
  } else if (input$select_csDendro=='Group'){
    csDendro(genes(cuff()))
  } else if (input$select_csDendro=='Individual'){
    csDendro(genes(cuff()),replicates=T)
  }
})

output$groups_genes_rpkm<-renderTable({
  if (is.null(mygene())){
    return(NULL)
  } else
  fpkm(mygene())
})
output$groups_isoforms_rpkm<-renderTable({
  if (is.null(mygene())){
    return(NULL)
  } else
  fpkm(isoforms(mygene()))
})
output$genePlots<-renderPlot({
  if (is.null(mygene())){
    return(NULL)
  } else
  g1 <- expressionBarplot(mygene())
  g2 <- expressionBarplot(isoforms(mygene()))
  grid.arrange(g1, g2, nrow=1)
})
# tabPanel 2Group summary
output$group1<-renderUI({
  group_name_list<-groups()
  selectInput(inputId = 'group1', label = 'Please select first group:', choices = group_name_list ,selected = group_name_list[1])
})
output$group2<-renderUI({
  group_name_list<-groups()
  selectInput(inputId = 'group2', label = 'Please select second group:', choices = group_name_list ,selected = group_name_list[2])
})

MA_plot<-eventReactive(input$get_plots,{
  MAplot(genes(cuff()), input$group1, input$group2)
})
csScatter_plot<-eventReactive(input$get_plots,{
  csScatter(genes(cuff()), input$group1, input$group2)
})
csVolcano_plot<-eventReactive(input$get_plots,{
  csVolcano(genes(cuff()), input$group1, input$group2)
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
# tabPanel 2Group gene and isoforms analysis
mygene_2group<-eventReactive(input$get_2group_gene,{
  getGene(cuff(),input$`2group_gene_name`)
})
df.tar_gene<-eventReactive(input$get_2group_gene,{
  gene_data(mygene_2group(), input$`2group_gene_name`,input$group1, input$group2)
})


output$`2group_gene_plot`<-renderPlot({
  if (is.null(mygene_2group())){
    return(NULL)
  } else
  g1 <- tableGrob(df.tar_gene(),theme=ttheme_minimal(base_size = 5))
  g2 <- ggplot(df.tar_gene(), aes(x=sample_name,y=fpkm))+
    geom_bar(stat = "identity", color="black", width=0.5, position=position_dodge())+
    geom_errorbar(aes(ymin=conf_lo, ymax=conf_hi), width=.2,position=position_dodge(0.5))+
    labs(title=paste(input$`2group_gene_name`, "differetial expression"), x="Group name", y = "FPKM")+theme_classic() + 
    scale_fill_manual(values=c('#123456','#345678'))
})
output$`2group_gene_rpkm`<-renderTable({
  if (is.null(mygene_2group())){
    return(NULL)
  } else
  df.tar_gene()
})
output$`2group_isoforms_rpkm`<-renderUI({
  table_output_list <- lapply(1:length(input$show_results), function(i) {
    tablename <- paste0("table", i)
    table_output_object <- plotlyOutput(tablename)
    table_output_object <- renderTable({
      source('global.R', local = TRUE)
      gene_expression_table(input$show_results[i]) 
    })
  })
  return(table_output_list)
})

}


shinyApp(ui, server)

