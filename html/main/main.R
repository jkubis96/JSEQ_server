library(shiny)
library(shinydashboard)
library(shinyvalidate)
library(DT)
library(rhdf5)
library(plotly)


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
              tabsetPanel(type = "tabs",
                          tabPanel("Data", column(12,style = "background-color:#FFFFFF",
                                                  
                                                  uiOutput('species'),
                                                  uiOutput('sex'),
                                                  uiOutput('tissue'),
                                                  uiOutput('affiliation'),
                                                  uiOutput('disease'),
                                                  uiOutput('status'),
                                                  DT::DTOutput("my_datatable"),
                                                  
                                                 
                                                  
                                                  
                                                  
                          )
                          
                          
                          
                          ),
                          tabPanel("Info", 
                                   column(12,style = "background-color:#FFFFFF",
                                          
                                  
                                    
                                          uiOutput("explor1"))
                                   
                                   
                                   ),
                         
                          tabPanel("Visualization",
                                  
                                   column(8,style = "background-color:#FFFFFF",

                                          plotlyOutput('explor4', width =  '1200px', height = '650px')

                                        ), column(2,offset = 2,
                                                  
                                          uiOutput('explor3')
                                                  
                                        ),
                                   


                                   
                          ),
                          tabPanel("Markers", 
                                   
                              
                                    column(10,style = "background-color:#FFFFFF",
                                           
                                           DTOutput('explor6', width =  '1100px')
                                           
                                    ), column(2,offset = 0,
                                              
                                              uiOutput('explor5')
                                              
                                    ),
                                   
                                   
                          ),
                          tabPanel("Gene_viewer", 
                                   column(10,style = "background-color:#FFFFFF",
                                          
                                          plotlyOutput('explor9', width =  '1100px', height = '650px')
                                          
                                   ), column(2,offset = 0,
                                             
                                             uiOutput('explor7'),
                                             uiOutput('explor8'),
                                             uiOutput('explor7.1')
                                             
                                   ),
                                   
                                   
                          ),
                          tabPanel("Interactions", 
                                   column(2,style = "background-color:#FFFFFF",
                                          
                                          
                                          
                                         )
                                   
                                   
                          ), tabPanel("Report", 
                                      column(12,style = "background-color:#FFFFFF",
                                             
                                             # htmlOutput("explor2"),
                                             
                                             
                                      )
                                      
                                      
                          ),
                                   
                               
                         
              ),
              
              
              
      ),
      
      tabItem(tabName = "contact",
            
              
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
        
        
        destDir <- file.path('../../validated')
        
        inFile <- file.path('../../projects', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))])
        

        file.copy(inFile, destDir, recursive=TRUE)
        
        unlink(inFile, recursive=TRUE)
        
        output$project_validation <- renderUI({ 
          HTML(paste('','',
                     'Project accepted for publication.', sep="<br/>"))
        })

        
        
       
        
        
        
      } else {
        
        inFile <- file.path('../../projects', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))])
        
        unlink(inFile, recursive=TRUE)  
        
        output$project_validation <- renderUI({ 
          HTML(paste('','',
                     'Project removed succedfully', sep="<br/>"))
        })
      }
      
    } else if (nchar(input$id_code) == 10 && FALSE  == unique(grepl(pattern = input$id_code, x = list.files(file.path('../../projects'))))) {
      
      
      output$project_validation <- renderUI({ 
        HTML(paste('','',
                   'No projects found.',
                   'Wrong ID', sep="<br/>"))
        
        
        
        
      })
      
      
    }
    
  })
  
    
    

    list <- list()
    
    
    n = 0
    
    for (file in list.files('../../validated/')) {
      n = n + 1
      list$id[n] <- sub('__.*', '', gsub('.*id-', '', file))
      list$species[n] <- sub('__.*', '', gsub('.*species-', '', file))
      list$sex[n] <- sub('__.*', '', gsub('.*sex-', '', file))
      list$tissue[n] <- sub('__.*', '', gsub('.*tissue-', '', file))
      list$affiliation[n] <- sub('__.*', '', gsub('.*tissue_affiliation-', '', file))
      list$disease[n] <- sub('__.*', '', gsub('.*disease_status-', '', file))
      list$sex[n] <- sub('__.*', '', gsub('.*sex-', '', file))
      list$status[n] <- sub('__.*', '', gsub('.*cell_development_status-', '', file))
      tmp <-  read.csv(file.path('../../validated',file, 'config'), sep = '=', header = FALSE)
      list$cells[n] <- gsub('/', '', tmp$V2[tmp$V1 == 'cells_after_qc'])
      list$subclasses[n] <- gsub('/', '', tmp$V2[tmp$V1 == 'subclasses_number'])
      list$subtypes[n] <- gsub('/', '', tmp$V2[tmp$V1 == 'subtypes_number'])
      list$action[n] <- sub('__.*', '', gsub('.*id-', '', file))

      
      
      
      
    }
    
    
    
 
    values <- reactiveValues()
    df <- as.data.frame(list)
    
    button <- function(tbl){
      function(i){
        sprintf(
          '<button id="button_%s_%s" type="button" onclick="%s">Explore</button>', 
          tbl, i, "Shiny.setInputValue('button', this.id);")
        
      }
    }
    
   
    df <- cbind(df, 
                  button = sapply(df$action, button("df")), 
                  stringsAsFactors = FALSE)
    
    df <- df[colnames(df)]
    

    species <- c()
    sex <- c()
    affiliation <- c()
    tissue <- c()
    disease <- c()
    status <- c()
    
    for (file in list.files('../../validated/')) {

      
      species <- c(species, sub('__.*', '', gsub('.*species-', '', file)))
      sex <- c(sex, sub('__.*', '', gsub('.*sex-', '', file)))
      tissue <- c(tissue, sub('__.*', '', gsub('.*tissue-', '', file)))
      affiliation <- c(affiliation, sub('__.*', '', gsub('.*tissue_affiliation-', '', file)))
      disease <- c(disease, sub('__.*', '', gsub('.*disease_status-', '', file)))
      status <- c(status, sub('__.*', '', gsub('.*cell_development_status-', '', file)))
      
    }
    
    species <- c('all', species)
    sex <- c('all', sex)
    tissue <- c('all', tissue)
    affiliation <- c('all', affiliation)
    disease <- c('all', disease)
    status <- c('all',status)
    
    species <- unique(species)
    sex <- unique(sex)
    tissue <- unique(tissue)
    affiliation <- unique(affiliation)
    disease <- unique(disease)
    status <- unique(status)
    
    names(species) <- species
    names(sex) <- sex
    names(tissue) <- tissue
    names(affiliation) <- affiliation
    names(disease) <- disease
    names(status) <- status
    
    output$species <- renderUI({ 
     
      
    column(2, selectInput("species_filter", 
                          h4("Species"), 
                          choices = species, selected = 'all')) })
    
    output$sex <- renderUI({
    column(2, selectInput("sex_filter", 
                          h4("Sex"), 
                          choices = sex, selected = 'all')) })
    
    output$tissue <- renderUI({
    column(2, selectInput("tissue_filter", 
                          h4("Tissue"), 
                          choices = tissue, selected = 'all')) })
    
    output$affiliation <- renderUI({
    column(2, selectInput("affiliation_filter", 
                          h4("Affiliation"), 
                          choices = affiliation, selected = 'all')) })
    
    output$disease <- renderUI({
    column(2, selectInput("disease_filter", 
                          h4("Disease"), 
                          choices = disease, selected = 'all')) })
    
    output$status <- renderUI({
    column(2, selectInput("status_filter", 
                          h4("Status"), 
                          choices = status, selected = 'all')) })

    
    observeEvent(input$species_filter, {
    
    if (input$species_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$species %in% input$species_filter,]
  
      
    } else {
      values$df <- df
      
     
      
      }
    })
    
    
    observeEvent(input$sex_filter, {
    if (input$sex_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$sex %in% input$sex_filter,]
    } else {
      values$df <- df
      }
    })
    
    observeEvent(input$tissue_filter, {
    if (input$tissue_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$tissue %in% input$tissue_filter,]
    } else {
      values$df <- df
      }
    })
    
    observeEvent(input$affiliation_filter, {
    if (input$affiliation_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$affiliation %in% input$affiliation_filter,]
    } else {
      values$df <- df
      }
    })
    
    observeEvent(input$disease_filter, {
    if (input$disease_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$disease %in% input$disease_filter,]
    } else { 
      values$df <- df}
    })
    
    observeEvent(input$status_filter, {
    if (input$status_filter != 'all') {
      values$df <- df
      values$df <- values$df[values$df$status %in% input$status_filter,]
    } else {
      values$df <- df
      }
    })
    
    
    
   
       output$my_datatable <- DT::renderDataTable(values$df,
                                               server = FALSE, escape = FALSE, selection = 'none')



       
       observeEvent(input[["button"]], {
         splitID <- strsplit(input[["button"]], "_")[[1]]
         tbl <- splitID[2]
         row <- splitID[3]
        
       
         txt <- read.csv(file.path('../../validated',list.files('../../validated/')[grep(row, list.files('../../validated/'))], 'config'), sep = '=', header = FALSE)
         txt <- txt[txt$V1 %in% c('author', 'project_name','email', 'source', 'sex', 'species', 'tissue', 'tissue_affiliation', 'cell_disease_status', 'cell_development_status', 'library', 'READS_LENGTH', 'cells_after_qc', 'subclasses_number', 'subtypes_number'),]
         text = ''
         for (d in 1:nrow(txt)) {
           text <- paste0(text, as.character(txt$V1[d]), ':', as.character(txt$V2[d]), '<br/>')
         }
         
         # output$info_out 
         
         output$explor1 <- renderText({HTML(text, sep="<br/>") })
         output$explor2 <- renderUI({includeHTML(file.path('../../validated/', list.files('../../validated/')[grep(row, list.files('../../validated/'))], 'Report.html'))})

         
         path = file.path('../../validated/', list.files('../../validated/')[grep(row, list.files('../../validated/'))], 'data.h5')
         cells_meta <- h5read(path, "metadata/cells_meta")
         
         output$explor3 <- renderUI({
           selectInput('plot_var', 'Choose data presentation:', choices = list('UMAP_subclasses' = 1, 'UMAP_subtypes' = 2, 'HDMAP_subtypes' = 3), selected = 1, multiple = FALSE, selectize = TRUE)
           
         })
         
         
         observeEvent(input$plot_var, {
           if (input$plot_var == 1) {
             
             
             plot <- ggplot(data = cells_meta, mapping = aes(UMAP_1, UMAP_2, color = subclass)) + geom_point()  + theme_classic() 
             
             output$explor4 <- renderPlotly({
               ggplotly(plot, originalData = FALSE, dynamicTicks = FALSE)
               
             })
             
             
             
           } else if (input$plot_var == 2) {
             
             
             plot <- ggplot(data = cells_meta, mapping = aes(UMAP_1, UMAP_2, color = subtypes)) + geom_point()  + theme_classic() 
             
             output$explor4 <- renderPlotly({
               ggplotly(plot, originalData = FALSE, dynamicTicks = FALSE)
               
             })
             
             
             
           } else if (input$plot_var == 3) {
             
             
             plot <- ggplot(data = cells_meta, mapping = aes(HDMAP_1, HDMAP_2, color = subtypes)) + geom_point()  + theme_classic() 
             
             output$explor4 <- renderPlotly({
               ggplotly(plot, originalData = FALSE, dynamicTicks = FALSE)
               
             })
             
             
             
           }
           
         
        
                  
         })
         
         
         
         output$explor5 <- renderUI({
           selectInput('markers_var', 'Choose markers:', choices = list('Subclasses' = 1, 'Subtypes' = 2, 'CSSG_markers' = 3), selected = 1, multiple = FALSE, selectize = TRUE)

         })


         observeEvent(input$markers_var, {
           if (input$markers_var == 1) {



             output$explor6 <- DT::renderDataTable({
                    datatable(h5read(path, "markers/subclass_markers"), options = list(pageLength = 15),
                              rownames= FALSE)

             })



           } else if (input$markers_var == 2) {


             output$explor6 <- DT::renderDataTable({
               datatable(h5read(path, "markers/subtypes_markers"), options = list(pageLength = 15),
                         rownames= FALSE)

             })



           } else if (input$markers_var == 3) {


            output$explor6 <- DT::renderDataTable({
               datatable(h5read(path, "markers/CSSG"), options = list(pageLength = 15),
                         rownames= FALSE)

             })



           }




         })
         
         
         output$explor7 <- renderUI({
           selectInput('genes_var', 'Choose data:', choices = list('Subclasses' = 1, 'Subtypes' = 2), selected = 2, multiple = FALSE, selectize = TRUE)
           
         })
         
         
         output$explor7.1 <- renderUI({
           actionButton("update_plot", 'Upadte', icon = icon('print'))
           
         })
           
         output$explor8 <- renderUI({
           textInput('genes_input', 'Enter genes:', value = "", width = NULL, placeholder = 'KIT, ednrb, Pax3, ...')


         })
         
         
         observeEvent(input$genes_var, {
           if (input$genes_var == 1) {
             
             tmp <- h5read(path, "frames/subclass_avg_norm_expression")
             genes <- h5read(path, "frames/subclass_avg_norm_expression_rows")
             
             genes_cssg <- gsub(pattern = '*.- ', '', colnames(tmp))
             rownames(tmp) <- genes
             
             pheat <- pheatmap::pheatmap(tmp, 
                                         clustering_method = 'ward.D',
                                         angle_col = 270, fontsize_row = 20, fontsize_col = 20)
             
             
             
             
             
      
             
             
             
           } else if (input$genes_var == 2) {
             
             
             tmp <- h5read(path, "frames/subtypes_avg_norm_expression")
             genes <- h5read(path, "frames/subtypes_avg_norm_expression_rows")
             
             genes_select <- gsub(pattern = '.*- ', '', colnames(tmp))
             
             rownames(tmp) <- genes
             

               
            
            
             
             plot <- t(as.matrix(tmp[toupper(rownames(tmp)) %in% toupper(genes_select),]))
             
             output$explor9 <- renderPlotly({plot_ly(
               x = colnames(plot), y = rownames(plot),
               z = plot, type = "heatmap",
             ) 
               
             })
             
             
             observeEvent(input$update_plot, {
               genes_select <- c(genes_select, c(strsplit(input$genes_input,split=", ",fixed=TRUE)[[1]]))
               
               plot <- t(as.matrix(tmp[toupper(rownames(tmp)) %in% toupper(genes_select),]))
               
               output$explor9 <- renderPlotly({plot_ly(
                 x = colnames(plot), y = rownames(plot),
                 z = plot, type = "heatmap",
               )
               })
               
             })
             
             
             
             
           } 
           
           
           
           
         })
         
          
        
  
         
         showModal(modalDialog(
           title = paste0("Dataset ", row, ' was chosen. Got to Explore tab'),
           size = 's',
           easyClose = TRUE,
           footer = NULL
         ))
       })

  
  
}


shinyApp(ui = ui, server = server)