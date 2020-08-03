library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(GenomicAlignments)
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
                                      textInput('directory',label = 'Please copy and paste the directory:'),
                                      actionButton('get_summary',label = 'Summarize'),
                                      plotOutput('summary'),
                                      plotOutput('dispersion') %>% withSpinner(color="#0dc5c1"),
                                      plotOutput('csScatterMatrix') %>% withSpinner(color="#0dc5c1"),
                                      plotOutput('csVolcanoMatrix') %>% withSpinner(color="#0dc5c1"),
                                      radioButtons(inputId = 'select_csDensity', label = 'Please select output csDensity plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                      plotOutput('csDensity') %>% withSpinner(color="#0dc5c1"),
                                      radioButtons(inputId = 'select_csBoxplot', label = 'Please select output csBoxplot plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                      plotOutput('csBoxplot') %>% withSpinner(color="#0dc5c1"),
                                      radioButtons(inputId = 'select_csDendro', label = 'Please select output csDendro plot type:',choices = c('Group','Individual'), selected = 'Group'),
                                      plotOutput('csDendro') %>% withSpinner(color="#0dc5c1"),
                             ),
                             tabPanel('Groups Expression', icon = icon('calendar-plus'),
                                      textInput(inputId = 'gene_name',label = 'Please enter a gene:',value = 'example:GAPDH'),
                                      actionButton('get_gene',label = 'Analyze'),
                                      plotOutput('genePlots') %>% withSpinner(color="#0dc5c1"),
                                      tableOutput('groups_genes_rpkm'),
                                      tableOutput('groups_isoforms_rpkm')
                             ),
                             tabPanel('2Group Summary', icon = icon('calendar-plus'),
                                      uiOutput('group1'),
                                      uiOutput('group2'),
                                      actionButton('get_plots',label = 'Get Plots'),
                                      #tags$div('MA_plot'),
                                      plotOutput('MA_plot') %>% withSpinner(color="#0dc5c1"),
                                      #tags$div('csScatter plot'),
                                      plotOutput('csScatter') %>% withSpinner(color="#0dc5c1"),
                                      #tags$div('csVolcano plot'),
                                      plotOutput('csVolcano') %>% withSpinner(color="#0dc5c1")
                             ),
                             tabPanel('2Group Expression', icon = icon('calendar-plus')
                                      
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

output$group1<-renderUI({
  group_name_list<-groups()
  selectInput(inputId = 'group1', label = 'Please select first group:', choices = group_name_list ,selected = group_name_list[1])
})
output$group2<-renderUI({
  group_name_list<-groups()
  selectInput(inputId = 'group2', label = 'Please select second group:', choices = group_name_list ,selected = group_name_list[2])
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
}
shinyApp(ui, server)

