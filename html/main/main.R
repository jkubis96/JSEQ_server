library(shiny)
library(shinydashboard)
library(shinyvalidate)


options(shiny.maxRequestSize = 300000*1024^2)

t <- getwd()

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "JSEQ_scRNAseq"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Description", tabName = "description", icon = icon("book-open")),
      menuItem("Analysis", tabName = "analysis", icon = icon("calculator")),
      menuItem("Validation", tabName = "validate", icon = icon("folder-open")),
      menuItem("Data exploration", tabName = "exp", icon = icon("list")),
      menuItem("Contact", tabName = "contact", icon = icon("address-book"))
      
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "home",
              img(src = "logo_jbs.PNG", height = 100, width = 100, align = "right"),
              #img(src = "ichb.jpg", height = 100, width = 400, align = "right"),
              h1('JSEQ_scRNAseq - tool'),
              br(),
              h4("Welcome to the JSEQ Single-Cell RNA Seq tool", align = "left"),
              h4("This approach and software was tested on AAV vectors.", align = "left"),
              br(),
              br(),
              h3('Autors: Jakub Kubis, Maciej Figiel', align = "center"),
              h3('Institute of Bioorganic Chemistry', br(),'Polish Academy of Sciences', br() , 'Department of Molecular Neurobiology', br(), align = "center"),
              
              
      ),
      
      tabItem(tabName = "analysis",
              h1("Single cell analysis by JSEQ", align = "left"),
              fluidRow(
                
                column(3,
                       h4("Account & Project", align = "left"),
                       textInput("user", "Username", 
                                 value = ""),
                       textInput("email", "E-mail", 
                                 value = ""),
                       textInput("project_name", "Project name", 
                                 value = ""),
                       textInput("description", "Short project description", 
                                 value = ""),
                       
                       
                ),
                
                column(3,
                       h4("Data", align = "left"),
                       textInput("source", "Data source", 
                                 value = "https://"),
                       selectInput("species", "Species", 
                                   choices = list('human' = 'human', 'mouse' = 'mouse')),
                       selectInput("tissue", "Tissue type", 
                                   choices = list('bladder' = 'bladder', 'blood' = 'blood', 'blood vessels' = 'blood vessels',
                                                  'bone marrow' = 'bone marrow', 'brain' = 'brain', 'eye' = 'eye',
                                                  'gut' = 'gut', 'heart' = 'heart', 'kidney' = 'kidney', 'liver' = 'liver',
                                                  'lung' = 'lung', 'pancreas' = 'pancreas', 'skin' = 'skin', 'stomach' = 'stomach',
                                                  'organoids' = 'organoids', 'other' = 'other')),
                       textInput("affiliation", "Tissue affiliation", 
                                 value = "", placeholder = 'cortex, beta-cells, etc.'),
                       
                       selectInput("status", "Cells status", 
                                   choices = list('healthy' = 'healthy', 'disease' = 'disease', 'unknow' = 'unknow')),
                       
                       selectInput("development", "Development status", 
                                   choices = list('mature' = 'mature', 'development' = 'development', 'cell culture' = 'cell_culture')),
                       
                ),
                column(3,
                       h4("Analysis settings", align = "left"),
                       
                       selectInput("input", "Input type", 
                                   choices = list('fastq' = 'fastq', 'sparse_matrix' = 'sparse', 'data_frame [.csv, .tsv, .txt]' = 'matrix')),
                       
                       selectInput("library", "Library type", 
                                   choices = list('dropseq' = 'dropseq', '10x_V5' = '10xv5', '10x_V3' = '10xv3')),
                       
                       
                       sliderInput("cell_n", "Estimated number of cells",
                                   min = 100, max = 200000, value = 1000),
                       
                       selectInput("reads", "Read length [+-/ 20 nt]", 
                                   choices = list('75' = 75, '100' = 100, '150' = 150, '200' = 200, '250' = 250, '300' = 300), selected = '100'),
                       
                       selectInput("naming", "Naming type", 
                                   choices = list('canonical' = 'canonical', 'non-canonical' = 'non_canonical')),
                       fileInput("files", "File input", multiple = TRUE, buttonLabel = 'Files'),
                       
                ),
                
                column(3,
                       h4("Action", align = "left"),
                       actionButton("check", 'Check', icon = icon('check')),
                       h4(""),
                       actionButton("start", "Run", icon = icon('play'))
                ))
              
              
      ),
      tabItem(tabName = "description",
              h2("Widgets tab content", align = "justify")
              
              
      ),
      tabItem(tabName = "validate",
              h2("Widgets tab content", align = "justify")
              
              
      ),
      tabItem(tabName = "exp",
              h2("Widgets tab content")
              
      ),
      tabItem(tabName = "contact",
              h2("Widgets tab content")
              
      )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$start, {
    
      
    CHECK <- TRUE
    
    while(CHECK) {
      samp<-c(2:9,letters,LETTERS,"!", "#", "%")
      code <- paste(sample(samp,10),collapse="")
      
      if (TRUE %in% grepl(pattern = code, x = list.files(file.path('../../projects')))) {
        CHECK = TRUE
      } else {
        CHECK = FALSE
      }
      
    }
    
    directory = paste0('project_name-', as.character(input$project_name), '__',
                       'user-', as.character(input$user) , '__', 
                       'id-', as.character(code) , '__',
                       'species-', as.character(input$species), '__',
                       'tissue-', as.character(input$tissue), '__',
                       'tissue_affiliation-', as.character(input$affiliation), '__',
                       'cell_disease_status-', as.character(input$status), '__',
                       'cell_development_status-', as.character(input$development))
    
    
    dir.create(file.path( paste0('../../projects/',directory )))
    dir.create(file.path( paste0('../../projects/',directory, '/sc_data'  )))
    dir.create(file.path(paste0('../../projects/',directory, '/results'  )))
    dir.create(file.path(paste0('../../projects/',directory, '/tmp'  )))
    
    
    
    if (input$input == 'fastq') {
      dir.create(file.path(paste0('../../projects/',directory, '/fast_data' )))
      
      destDir <- file.path(paste0('../../projects/',directory, '/fast_data' ))
      
      inFile <- input$files
      file.copy( inFile$datapath, file.path(destDir, inFile$name) )
      
    } else {
      destDir <- file.path(paste0('../../projects/',directory, '/sc_data' ))
      
      inFile <- input$files
      file.copy(inFile$datapath, file.path(destDir, inFile$name) )
    }
    
    conf = paste0(
      'user=', as.character(input$user), '\n',
      'email=', as.character(input$email), '\n',
      'project_name=', as.character(input$project_name), '\n',
      'source=', as.character(input$source), '\n',
      'species=', as.character(input$species), '\n',
      'tissue=', as.character(input$tissue), '\n',
      'tissue_affiliation=', as.character(input$affiliation), '\n',
      'cell_disease_status=', as.character(input$status), '\n',
      'cell_development_status=', as.character(input$development) ,' \n',
      'data_format=', as.character(input$input) ,' \n',
      'library=', as.character(input$library) ,' \n',
      'marker_type=', as.character(input$naming) ,' \n',
      'READS_LENGHT=', as.character(input$reads) ,' \n',
      'id=', as.character(code) ,'\n',
      'directory_name=',directory, '\n',
      'cell=',input$cell_n
      
    )
    
    conf <- gsub(pattern = ' ', replacement = '', conf)
    
    writeLines(conf, file.path(paste0('../../projects/',directory, "/config")))
    
    
    c = paste0('echo "' ,directory,'" >> ../../tasks')
    
    system(c)

    
  })
  
  observeEvent(input$check, {
    
    iv <- InputValidator$new()
    iv$add_rule("user", sv_required())
    iv$add_rule("project_name", sv_required())
    iv$add_rule("description", sv_required())
    iv$add_rule("affiliation", sv_required())
    iv$add_rule("email", sv_required())
    iv$add_rule("email", sv_email())
    iv$enable()
    
  })
  
  
  
  
}


shinyApp(ui = ui, server = server)