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
      id = "tabs",
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
                       textInput("user", "Author", 
                                 value = "",  placeholder = 'name & affiliation'),
                       textInput("email", "E-mail", 
                                 value = "", placeholder = 'user@mail.com'),
                       textInput("project_name", "Project name", 
                                 value = "", placeholder = 'name'),
                       textInput("description", "Short project description", 
                                 value = "", placeholder = 'information about data content'),
                       textInput("source", "Data source", value = '', placeholder = 'DOI, link, GEO, unpublished, etc.')
                       
                      
                       
                       
                ),
                
                column(3,
                       h4("Data", align = "left"),
                       
                       selectInput("species", "Species", 
                                   choices = list('human' = 'human', 'mouse' = 'mouse')),
                       
                       
                       selectInput("sex", "Sex", 
                                   choices = list('male' = 'male', 'female' = 'female', 'other' = 'other')),
                       
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

                       
                       numericInput("cell_n", "Estimated number of cells:", value = 1000 , min = 100, max = 200000),
                       
                       selectInput("reads", "Read length [+-/ 20 nt]", 
                                   choices = list('75' = 75, '100' = 100, '150' = 150, '200' = 200, '250' = 250, '300' = 300), selected = '100'),
                       
                       selectInput("naming", "Naming type", 
                                   choices = list('canonical' = 'canonical', 'non-canonical' = 'non_canonical')),
                       fileInput("files", "Files input", multiple = TRUE, accept = c('.csv', '.tsv', '.txt', '.fastq', '.mtx'), buttonLabel = 'Files'),
                       
                ),
                
                column(3,
                       h4("Action", align = "left"),
                       actionButton("check", 'Check', icon = icon('check')),
                       h4(""),
                       uiOutput("run_buttom"),
                       uiOutput("accepted")
                       
                ))
              
              
      ),
      tabItem(tabName = "description",
              h2("Widgets tab content", align = "justify")
              
              
      ),
      tabItem(tabName = "validate",
              h2("Data validation", align = "justify"),
              
              fluidRow(
                
                column(3, 
                       h4("Project", align = "left"),
                       textInput("id_code", "Project ID", 
                                 value = "", placeholder = '10-character code'),
                       
                       actionButton("search", 'Check', icon = icon('print')),
                       selectInput("val_choose", 
                                          h3("Data validation:"), 
                                          choices = list("Yes" = 1, 
                                                         "No" = 2
                                          ), selected = 1),
                       
                       
                      
                        actionButton("val", "Validate", icon = icon('edit')),
                  
                  
                        HTML(paste('','',
                             'Data validation.',
                             'If the output of the analysis is appropriate you can share data for other users. If you choose "Yes" the data will be available in "Data exploration" tab where you can obtain more data [analysis output files in hdf5 formats]. The data can be always removed after contact with us. If you choose "No" data will be removed. If you will not make choice the data will remove in a few days.', sep="<br/>"))
                  
                       
                      
                       
                       
                ),
                
                column(9,style = "background-color:#FFFFFF",
                       h4("Results", align = "left"),
                       uiOutput("project_validation", inline = TRUE)
                    
                       
                ))
              
              
      ),
      tabItem(tabName = "exp",
              h2("Widgets tab content"),
              
              
              column(2, 
                     h4("Project", align = "left"),
                     textInput("id_code", "Project ID", 
                               value = "", placeholder = '10-character code'),
                     
                     actionButton("search", 'Check', icon = icon('print')),
                     selectInput("val_choose", 
                                 h3("Data validation:"), 
                                 choices = list("Yes" = 1, 
                                                "No" = 2
                                 ), selected = 1),
                     
                     
                     
                     actionButton("val", "Validate", icon = icon('edit')),
                     
                     
                     HTML(paste('','',
                                'Data validation.',
                                'If the output of the analysis is appropriate you can share data for other users. If you choose "Yes" the data will be available in "Data exploration" tab where you can obtain more data [analysis output files in hdf5 formats]. The data can be always removed after contact with us. If you choose "No" data will be removed. If you will not make choice the data will remove in a few days.', sep="<br/>"))
                     
                     
                     
                     
                     
              ),
              
              column(10,style = "background-color:#FFFFFF",
                     h4("Results", align = "left"),
                     DT::DTOutput("my_datatable")
                     
                     
              )
              
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
      samp<-c(2:9,letters,LETTERS,"!", "#")
      code <- paste(sample(samp,10),collapse="")
      
      if (TRUE %in% grepl(pattern = code, x = list.files(file.path('../../projects'))) && TRUE %in% grepl(pattern = code, x = list.files(file.path('../../validated')))) {
        CHECK = TRUE
      } else {
        CHECK = FALSE
      }
      
    }
    
    directory = paste0('project_name-', as.character(gsub(pattern = ' ', '_', x =  input$project_name)), '__',
                       'id-', as.character(code) , '__',
                       'species-', as.character(input$species), '__',
                       'sex-', as.character(input$sex), '__',
                       'tissue-', as.character(input$tissue), '__',
                       'tissue_affiliation-', as.character(gsub(pattern = ' ', '_', x =  input$affiliation)), '__',
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
      'author=', as.character(input$user), '\n',
      'email=', as.character(input$email), '\n',
      'project_name=', as.character(gsub(pattern = ' ', '_', x =  input$project_name)), '\n',
      'source=', as.character(input$source), '\n',
      'species=', as.character(input$species), '\n',
      'sex=', as.character(input$sex), '\n',
      'tissue=', as.character(input$tissue), '\n',
      'tissue_affiliation=', as.character(gsub(pattern = ' ', '_', x =  input$affiliation)), '\n',
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
    
    output$accepted <- renderUI({ 
      HTML(paste('','',
                 'Your analysis has started...',
                 'Analysis can last several hours.',
                 'If the analysis is finished, you will receive ',
                 'an email with information about the finish', 
                 'and access code to the analysis.',
                 'Then the results you can check in the Validation tab', sep="<br/>"))
    
       })

    
  })
  
  observeEvent(input$check, {
    
    iv <- InputValidator$new()
    iv$add_rule("user", sv_required())
    iv$add_rule("project_name", sv_required())
    iv$add_rule("description", sv_required())
    iv$add_rule("affiliation", sv_required())
    iv$add_rule("email", sv_required())
    iv$add_rule("email", sv_email())
    iv$add_rule('source', sv_required())
    iv$add_rule('files', sv_required())
    iv$add_rule('cell_n', sv_required())
    iv$enable()
    
    if (length(input$cell_n) > 0 && length(input$files) > 0 && length(input$user) > 0 && length(input$project_name) > 0 && length(input$description) > 0 && length(input$affiliation) > 0 && length(input$email) > 0 && grepl('@', input$email) && length(input$source) > 0) {
        
      
      output$run_buttom <- renderUI({
        actionButton("start", "Run", icon = icon('play')) })
      output$accepted <- renderUI({ 
        'All data correct. Run can run analysis.'
        
      })
    } else {
   
      
      output$accepted <- renderUI({ 
        'Provide all required data ...'
        
      })
      
    }
    
  
    
    
  })
  
  
  observeEvent(input$search, {

    if (nchar(input$id_code) == 10 && TRUE %in% grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))) {
      output$project_validation <- renderUI({
        includeHTML(file.path('../../projects/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))], 'Report.html'))
        
        })
      


    
    } else {
      
      
      output$project_validation <- renderUI({ 
        HTML(paste('','',
                   'No projects found.',
                   'Possible scenarios:',
                   '-wrong project ID',
                   '-analysis failed due to poor data quality',
                   '- the input data was in the wrong format', 
                   '-selected analysis options were wrongly matched to the dataset',
                   'You can try to re-enable analytics or contact us for guidance on a failed analytics. Contact details can be found in the Contact tab.', sep="<br/>"))
        


        
      })
      
    }
    
    
  })
  
  
  observeEvent(input$val, {
    

    
      if (nchar(input$id_code) < 10) {
      
      
      iv <- InputValidator$new()
      iv$add_rule("id_code", sv_required())
      iv$enable()
      
      output$project_validation <- renderUI({ 
        HTML(paste('','',
                   'Too short ID', sep="<br/>"))
        
        
        
        
      })
        

      
    } else if (input$val_choose != 0 && nchar(input$id_code) == 10 && TRUE %in% grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))) {
      
      if(input$val_choose == 1) {
        
        
        mv = paste0('cd ../../; mv ' ,paste0('$(pwd)/projects/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))]),' ', paste0('$(pwd)/validated/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))]))
        
        system(mv)
        
        
        
        
      } else {
        
        rm = paste0('cd ../../; rm -rf ' ,paste0('$(pwd)/projects/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))]))
        
        system(rm)        
      }
      
      
    } else if (nchar(input$id_code) == 10 && FALSE  == unique(grepl(pattern = input$id_code, x = list.files(file.path('../../projects'))))) {
      
      
      output$project_validation <- renderUI({ 
        HTML(paste('','',
                   'No projects found.',
                   'Wrong ID', sep="<br/>"))
        
        
        
        
      })
      
      
    }
    
  })
  
  

  
   values <- reactiveValues()
    
    
    
    
    list <- list()
    
    n = 0
    
    for (file in list.files('../../projects/')) {
      n = n + 1
      list$project_name[n] <- sub('__.*', '', gsub('.*project_name-', '', file))
      list$project_id[n] <- sub('__.*', '', gsub('.*id-', '', file))
      list$species[n] <- sub('__.*', '', gsub('.*species-', '', file))
      list$tissue[n] <- sub('__.*', '', gsub('.*tissue-', '', file))
      list$tissue_affiliation[n] <- sub('__.*', '', gsub('.*tissue_affiliation-', '', file))
      list$disease_status[n] <- sub('__.*', '', gsub('.*disease_status-', '', file))
      list$sex[n] <- sub('__.*', '', gsub('.*sex-', '', file))
      list$cell_development_status[n] <- sub('__.*', '', gsub('.*cell_development_status-', '', file))
      
      
      
      
    }
    
    values$df <- as.data.frame(list)
    output$my_datatable <- DT::renderDataTable(as.data.frame(list),
                        selection = 'single')

    # output$my_datatable = DT::renderDataTable(
    #   values$df, escape = FALSE, selection = 'none', server = FALSE, 
    #   options = list(dom = 't', paging = FALSE, ordering = FALSE),
    #   editable = TRUE,
    #   callback = JS("table.rows().every(function(i, tab, row) {
    #     var $this = $(this.node());
    #     $this.attr('id', this.data()[0]);
    #     $this.addClass('shiny-input-container');
    #   });
    #   Shiny.unbindAll(table.table().node());
    #   Shiny.bindAll(table.table().node());")
    # )
    
    
  
  
  
}


shinyApp(ui = ui, server = server)