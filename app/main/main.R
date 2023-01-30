library(shiny)
library(shinydashboard)
library(shinyvalidate)
library(DT)
library(rhdf5)
library(plotly)
library(rsconnect)
library(rlang)
library(shinyauthr)
library(shinyWidgets)
library("RPostgreSQL")
library("DBI")
library(sodium)
library(shinydashboardPlus)
library(leaflet)
library(maps)



options(shiny.maxRequestSize = 300000*1024^2)



sys <- read.csv('../../requirements_file/system_files', sep = '=', header = F)





ui <- dashboardPage(
  
  skin = "purple",
  dashboardHeader(title = ("JBSDA"), tags$li(class = 'dropdown' ,uiOutput('acc'))),
  footer = dashboardFooter(
    left = "By Jakub Kubis - JBioSystem",
    right = paste0("Institute of Bioorganic Chemistry, PAS, Poznan, Poland,", format(Sys.Date(), format="%Y"))
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Description", tabName = "description", icon = icon("book-open")),
      menuItem("Analysis", tabName = "analysis", icon = icon("calculator"),
               menuSubItem('scRNAseq', tabName='scRNAseq'),
               menuSubItem('RNAseq', tabName='RNAseq'),
               menuSubItem('Methylation', tabName='Methylation')),
      menuItem("Data exploration", tabName = "exp", icon = icon("list"),
               menuSubItem('scRNAseq', tabName='scRNAseq_explor'),
               menuSubItem('RNAseq', tabName='RNAseq_explor'),
               menuSubItem('Methylation', tabName='Methylation_explor')),
      menuItem("Validation", tabName = "validate", icon = icon("folder-open")),
      menuItem("Contact", tabName = "contact", icon = icon("address-book"))
      
    )
  ),
  dashboardBody(
          
      tabItems(
      # First tab content
      tabItem(tabName = "home",
              img(src = "logo_jbs.PNG", height = 100, width = 100, align = "right"),
              h1('JBDB - database'),
              br(),
              h4("", align = "left"),
              h4("", align = "left"),
              br(),
              br(),
              h3('', align = "center"),
              h3('Institute of Bioorganic Chemistry', br(),'Polish Academy of Sciences', br() , 'Department of Molecular Neurobiology', br(), align = "center"),
              
              
      ),
      
      tabItem(tabName = "scRNAseq",
             fluidRow(
                
                column(3,
                       h4("Account & Project", align = "left"),
                       br(),
                       uiOutput('login_element1'),
                       br(),
                       uiOutput('login_element2'),
                       br(),
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
                       
                       selectInput("development", "Development status", 
                                   choices = list('mature' = 'mature', 'development' = 'development', 'cell culture' = 'cell_culture')),
                       
                       selectInput("status", "Cells status", 
                                   choices = list('healthy' = 'healthy', 'disease' = 'disease', 'unknow' = 'unknow')),
                       
                       uiOutput('disease_name')
                       
                      
                       
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
      tabItem(tabName = "RNAseq",
              HTML({'Currently not avaiable'})
            
              
              
              
      ),
      tabItem(tabName = "Methylation",
              HTML({'Currently not avaiable'})
              
              
      ),
      tabItem(tabName = "RNAseq_explor",
              HTML({'Currently not avaiable'})
              
              
              
              
      ),
      tabItem(tabName = "Methylation_explor",
              HTML({'Currently not avaiable'})
              
              
      ),
      tabItem(tabName = "description",
              h2("Widgets tab content", align = "justify")
              
              
      ),
      tabItem(tabName = "validate",
              h2("Data validation", align = "justify"),
              
              fluidRow(
                
                column(3, 
                       br(),
                       uiOutput('login_element1_validation'),
                       br(),
                       uiOutput('login_element2_validation'),
                       br(),
                       textInput("id_code", "Project ID", 
                                 value = "", placeholder = '10-character code'),
                       
                       actionButton("search", 'Check', icon = icon('print')),
                       selectInput("val_choose", 
                                          h3("Data validation:"), 
                                          choices = list("Yes" = 1, 
                                                         "No" = 2
                                          ), selected = 1),
                       
                       
                      
                        actionButton("val", "Validate", icon = icon('edit')),
                  
                  
                   
                       
                       
                ),
                
                column(4,style = "background-color:#FFFFFF",
                      uiOutput("project_validation"),
                      
                    
                       
                ),
                column(4, 
                       uiOutput("project_validation_report", inline = TRUE),
                       br(),
                       uiOutput('disclaimer')
                       
                       
                ))
              
              
      ),
      tabItem(tabName = "scRNAseq_explor",
              tabsetPanel(type = "tabs",
                          tabPanel("Data", column(12,style = "height:600px;background-color:#FFFFFF",
                                                  
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
                                   column(4,style = "height:600px;background-color:#FFFFFF",
                                          
                                  
                                          HTML('<br/>'),
                                          HTML('<br/>'),
                                          tableOutput("explor1")),
                                   
                                   column(5,style = "height:600px;background-color:#FFFFFF",
                                         
                                          
                                          uiOutput('pd'),
                                          
                                          uiOutput('au'),
                                          
                                          uiOutput('user_check_button'),
                                          
                                         
                                          textOutput("project_description")
                                       
                                      
                                         
                                          ),
                                   
                                   column(3,style = "height:600px;background-color:#FFFFFF",
                                          uiOutput('ar'),
                                          
                                          
                                          HTML('<br/>'),
                                          uiOutput('report_b'),
                                          
                                          br(),
                                          uiOutput('login_element1_download'),
                                          br(),
                                          uiOutput('login_element2_download'),
                                          br(),
                                          
                                          uiOutput('down_butt'),
                                          uiOutput('down_info'),
                                          

                                          )
                                   
                                   
                                   ),
                         
                          tabPanel("Visualization",
                                  
                                   column(10, style = "height:600px;background-color:#FFFFFF",

                                          plotlyOutput('explor4', height =  '600px')

                                        ), column(2,offset = 0,
                                                  
                                          uiOutput('explor3')
                                                  
                                        ),
                                   


                                   
                          ),
                          tabPanel("Markers", 
                                   
                              
                                    column(10, style = "height:600px;background-color:#FFFFFF",
                                           
                                           DTOutput('explor6')
                                           
                                    ), column(2,offset = 0,
                                              
                                              uiOutput('explor5'),
                                              HTML('<br/>'),
                                              uiOutput('sub1'),
                                              uiOutput('sub2'),
                                              uiOutput('sub3')
                                              
                                              
                                              
                                              
                                              
                                    ),
                                   
                                   
                          ),
                          tabPanel("Gene_viewer", 
                                   column(10,style = "height:600px;background-color:#FFFFFF",
                                          
                                          plotlyOutput('explor9',  height =  '600px')
                                          
                                   ), column(2,offset = 0,
                                             
                                             uiOutput('explor7'),
                                             uiOutput('explor8'),
                                             uiOutput('explor7.1')
                                             
                                   ),
                                   
                                   
                          ),
                          tabPanel("Interactions", 
                                   column(2,style = "height:600px;background-color:#FFFFFF",
                                          
                                          
                                          
                                         )
                                   
                                   
                          )
                                   
                               
                         
              ),
              
              
              
      ),
      
      tabItem(tabName = "contact",
              fluidRow(
                column(4,),
                column(5,offset = 3,
                       leafletOutput("mymap")),
              
              )
              
      )
    )
  )
)

server <- function(input, output, session) {
    
  #location

    cities <- as.data.frame(world.cities)
    output$mymap <- renderLeaflet({
      leaflet(options = leafletOptions(maxZoom = 6)) %>%
        addProviderTiles(providers$Stamen.TonerLite,
                         options = providerTileOptions(noWrap = TRUE)
        ) %>%
        addMarkers(data = cities[cities$name %in% 'Poznan',])
    })
  
    account_data <- reactiveValues()
    account_data$USER = NULL
    
  
    #Account
  
    output$acc <- renderUI({dropdownButton(
    tags$h4("Account menu"),
    actionBttn('login_header', 'Login', icon = icon('arrow-pointer'), size = 'md'),
    br(),
    actionBttn('register', 'Register', icon = icon('id-card'), size = 'md'),
    br(),
    circle = TRUE, status = 'succes',
    icon = icon("user"), width = "200px", size = 'default', right = TRUE)
     
     })
    
    
    observeEvent(input$login_header, {
      
      
      showModal(modalDialog(
        
        title = 'Login account',
        br(),
        textInput("log_login", "Login", 
                  value = "", placeholder = 'User'),
        passwordInput("log_password", "Password", 
                      value = "", placeholder = 'min. 8-character password'),
        
        checkboxInput('robot', 'I am not a robot!', value = FALSE),
        br(),
        
        uiOutput('decision_log'),
        actionButton('login_button', 'Login', width = '300px'), 
        br(),
        br(),
        actionButton('remind', 'Remind password', width = '300px'), size = 'xl',
        
        
        
        
        
      ) 
      )
      
      
      
      
      
    })
    
    #login before analysis 
    
    output$login_element1 <- renderUI({
      
      HTML({'Login is required to start analysis'})
      
    })
    
    output$login_element2 <- renderUI({
     
      actionButton('log_analysis', 'Login', icon = icon('arrow-pointer'), width = '200px')
      
    })
    
    observeEvent(input$log_analysis, {
      
      
      showModal(modalDialog(
        
        title = 'Login account',
        br(),
        textInput("log_login", "Login", 
                  value = "", placeholder = 'User'),
        passwordInput("log_password", "Password", 
                      value = "", placeholder = 'min. 8-character password'),
        
        checkboxInput('robot', 'I am not a robot!', value = FALSE),
        br(),
        
        uiOutput('decision_log'),
        actionButton('login_button', 'Login', width = '300px'), 
        br(),
        br(),
        actionButton('remind', 'Remind password', width = '300px'), size = 'xl',
        
        
        
        
        
      ) 
      )
      
    })
    
    
    
    
    #login before validation 
    
    output$login_element1_validation <- renderUI({
      
      HTML({'Login is required to validate'})
      
    })
    
    output$login_element2_validation <- renderUI({
      
      actionButton('log_validation', 'Login', icon = icon('arrow-pointer'), width = '200px')
      
    })
    
    observeEvent(input$log_validation, {
      
      
      showModal(modalDialog(
        
        title = 'Login account',
        br(),
        textInput("log_login", "Login", 
                  value = "", placeholder = 'User'),
        passwordInput("log_password", "Password", 
                      value = "", placeholder = 'min. 8-character password'),
        
        checkboxInput('robot', 'I am not a robot!', value = FALSE),
        br(),
        
        uiOutput('decision_log'),
        actionButton('login_button', 'Login', width = '300px'), 
        br(),
        br(),
        actionButton('remind', 'Remind password', width = '300px'), size = 'xl',
        
        
        
        
        
      ) 
      )
      
    })
    
    
    #login before download
    
    login_values <- reactiveValues()
    login_values$row <- NULL
    login_values$login_element1_download <- renderUI({
    
      
      HTML({'Login is required to download'})
      
    })
    
    login_values$login_element2_download <- renderUI({
      
      actionButton('log_download', 'Login', icon = icon('arrow-pointer'), width = '200px')
      
    })
    
    
    
    
    
    observeEvent(input$log_download, {
      
      
      showModal(modalDialog(
        
        title = 'Login account',
        br(),
        textInput("log_login", "Login", 
                  value = "", placeholder = 'User'),
        passwordInput("log_password", "Password", 
                      value = "", placeholder = 'min. 8-character password'),
        
        checkboxInput('robot', 'I am not a robot!', value = FALSE),
        br(),
        
        uiOutput('decision_log'),
        actionButton('login_button', 'Login', width = '300px'), 
        br(),
        br(),
        actionButton('remind', 'Remind password', width = '300px'), size = 'xl',
        
        
        
        
        
      ) 
      )
      
    })
    
    ########################################################################
     
  
   
   observeEvent(input$login_button, {
     if (!toupper(input$log_login) %in% toupper(dbGetQuery(connec, 'SELECT user_login FROM users;')[[1]])) {
       output$decision_log <- renderUI({HTML('User does not exist')})
     } else if (!password_verify(dbGetQuery(connec, paste0("SELECT password FROM users WHERE user_login = '", toupper(input$log_login) ,"';"))[[1]], input$log_password)) {
       output$decision_log <- renderUI({HTML('Wrong password!')})
     } else if (input$robot == FALSE) {
       output$decision_log <- renderUI({HTML('Are you a robot?')})
     } else  {
       #LOGIN FUNCTION
       output$decision_log <- renderUI({HTML('Logged successfully. You can start your analysis with data exploration')})
       account_data$USER <- dbGetQuery(connec, paste0("SELECT user_name FROM users WHERE user_login = '", toupper(input$log_login) ,"';"))[[1]]
       
       output$login_element1 <- renderUI({HTML( paste('<strong> User:', account_data$USER, '</strong>'))}) 
       output$login_element2 <- renderUI({HTML(paste('<strong> Email:', dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong>'))}) 
       
       output$login_element1_validation <- renderUI({HTML( paste('<strong> User:', account_data$USER, '</strong>'))}) 
       output$login_element2_validation <- renderUI({HTML(paste('<strong> Email:', dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong>'))}) 
       
       login_values$login_element1_download <- renderUI({HTML( paste('<strong> User:', account_data$USER, '</strong>'))}) 
       login_values$login_element2_download <- renderUI({HTML('You are logged in. You can download the files.')}) 
       #TU
       if (is.null(login_values$row)) {
         login_values$down_butt <- renderUI({HTML(paste0('<br>', downloadButton('downloadData', label = "Download data")))})
         login_values$login_element1_download <- renderUI({HTML( paste('<strong> User:', account_data$USER, '</strong>'))}) 
         login_values$login_element2_download <- renderUI({HTML('You are logged in. You can download the files.')}) 

       } else {
         output$down_butt <- renderUI({HTML(paste0('<br>', downloadButton('downloadData', label = "Download data")))})
         output$login_element1_download <- renderUI({HTML( paste('<strong> User:', account_data$USER, '</strong>'))}) 
         output$login_element2_download <- renderUI({HTML('You are logged in. You can download the files.')}) 

         
         
       }
       output$acc <- renderUI({dropdownButton(
         tags$h4("Account menu"),
         actionBttn('user_profil', 'Profil', icon = icon('circle-user'), size = 'md'),
         br(),
         actionBttn('user_options', 'Options', icon = icon('gear'), size = 'md'),
         br(),
         actionBttn('user_logout', 'Logout', icon = icon('arrow-right'), size = 'md'),
         br(),
         circle = FALSE, status = 'succes',
         icon = icon("user"), label = paste('Hi!',account_data$USER), width = "200px", size = 'default', right = TRUE)
         
       })
       
       
     }
     
     
     iv <- InputValidator$new()
     iv$add_rule("user_login", sv_required())
     iv$add_rule("user_password", sv_required())
     iv$add_rule("log_login", sv_required())
     iv$add_rule("log_password", sv_required())
     iv$add_rule("user_email", sv_required())
     iv$add_rule("user_email", sv_email())
     iv$add_rule("rodo", sv_required())
     iv$enable()
     
     
   })
   
   
   #Logout, profil & options##########################################################
   
   
   observeEvent(input$user_profil, {
     
     showModal(modalDialog(
       
       tags$h2('User profil'),
       
       tags$h3(account_data$USER),
       br(),
       HTML({paste0('<strong> Name:  ', dbGetQuery(connec, paste0("SELECT name FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Surname:  ', dbGetQuery(connec, paste0("SELECT surname FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Contact:  ', dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Affiliation:  ', dbGetQuery(connec, paste0("SELECT affiliation FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Country:  ', dbGetQuery(connec, paste0("SELECT country FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> ORCID:  ', dbGetQuery(connec, paste0("SELECT ORCID FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Other links:  ', dbGetQuery(connec, paste0("SELECT links FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}),
       HTML({paste0('<strong> Account type:  ', dbGetQuery(connec, paste0("SELECT acc_type FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], '</strong> </br>')}), size = 'm',
       
       
  ) 
     )
     
     
   })
   
   
   observeEvent(input$user_logout, {
     
     output$acc <- renderUI({dropdownButton(
       tags$h4("Account menu"),
       actionBttn('login_header', 'Login', icon = icon('arrow-pointer'), size = 'md'),
       br(),
       actionBttn('register', 'Register', icon = icon('id-card'), size = 'md'),
       br(),
       circle = TRUE, status = 'succes',
       icon = icon("user"), width = "200px", size = 'default', right = TRUE)
       
     })
     
     output$login_element1 <- renderUI({
       
       HTML({'Login is required to start analysis'})
       
     })
     
     output$login_element2 <- renderUI({
       
       actionButton('log_analysis', 'Login', icon = icon('arrow-pointer'), width = '200px')
       
     })
     
     login_values$login_element1_download <- renderUI({
       
       HTML({'Login is required to download'})
       
     })
     
     login_values$login_element2_download <- renderUI({
       
       actionButton('log_download', 'Login', icon = icon('arrow-pointer'), width = '200px')
       
     })
     
     output$login_element1_validation <- renderUI({
       
       HTML({'Login is required to validate'})
       
     })
     
     output$login_element2_validation <- renderUI({
       
       actionButton('log_validation', 'Login', icon = icon('arrow-pointer'), width = '200px')
       
     })
     
     account_data$USER = NULL
     output$decision_log <- renderUI({HTML('Logout successful. You can login again')})
     
     output$down_butt <- renderUI({HTML('')})
     
     
     if (!is.null(login_values$row)) {
     output$login_element1_download <- renderUI({
       
       
       HTML({'Login is required to download'})
       
     })
     
     output$login_element2_download <- renderUI({
       
       actionButton('log_download', 'Login', icon = icon('arrow-pointer'), width = '200px')
       
     })
     
     }
     
     login_values$row <- NA
     
     
     
     
     
  
     })
     
     
     
   
   
   
   observeEvent(input$user_options, {
     
     showModal(modalDialog(
       
       title = 'Settings',
       br(),
    
       actionButton('change_data', 'Change your profil data', width = '300px'), 
       br(),
       br(),
       actionButton('change_password', 'Change your password', width = '300px'), size = 'xl',
       
       
       
       
       
     ) 
     )
     
     
   })
   
   
   observeEvent(input$change_data, {
     
     showModal(modalDialog(
       
       title = 'Data change',
       br(),
       
      
       textInput("update_email", "E-mail", 
                 value = dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_name", "Name", 
                 value = dbGetQuery(connec, paste0("SELECT name FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_surname", "Surname", 
                 value = dbGetQuery(connec, paste0("SELECT surname FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_affiliation", "Anffiliation", 
                 value = dbGetQuery(connec, paste0("SELECT affiliation FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_country", "Country", 
                 value = dbGetQuery(connec, paste0("SELECT country FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_ORCID", "ORCID", 
                 value = dbGetQuery(connec, paste0("SELECT ORCID FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       textInput("update_link", "Other links", 
                 value = dbGetQuery(connec, paste0("SELECT links FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]),
       
       br(),
       
       uiOutput('changes_suc'),
       actionButton('update_data', 'Update', width = '300px'), size = 'xl',
       
       
       
       
       
     ) 
     )
     
     
   })
   
   
   

   observeEvent(input$update_data, {
     
     dbGetQuery(connec, paste0("UPDATE users SET email = '",input$update_email,"', name = '",input$update_name,"', surname = '",input$update_surname,"', affiliation = '",input$update_affiliation,"', country = '",input$update_country,"', ORCID = '",input$update_ORCID,"', links = '",input$update_link,"' WHERE user_login = '", toupper(account_data$USER) ,"';"))
     output$changes_suc <- renderUI({HTML('Data has changed successfully')})
     
   
     
   })
   
   
   
   observeEvent(input$change_password, {
     
     showModal(modalDialog(
       
       title = 'Change password',
       br(),
      
       passwordInput("old_password", "Old password", 
                     value = "", placeholder = 'min. 8-character password'),
       
       passwordInput("new_password1", "New password 1st", 
                     value = "", placeholder = 'min. 8-character password'),
       
       passwordInput("new_password2", "New password 2nd", 
                     value = "", placeholder = 'min. 8-character password'),
      
       
       uiOutput('change_pass'),
       br(),
       actionButton('update_password', 'Change', width = '300px'), size = 'xl',
       
       
       
       
       
     ) 
     )
     
     
   })
   
   
   observeEvent(input$update_password, {
     
     if (!password_verify(dbGetQuery(connec, paste0("SELECT password FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]], input$old_password)) {
       output$change_pass <- renderUI({HTML('Wrong old password')})
     } else if (nchar(input$new_password1) < 8) {
       output$change_pass <- renderUI({HTML('The new password is too short. The minimum number of characters is 8')})
     } else if (input$new_password1 != input$new_password2) {
       output$change_pass <- renderUI({HTML('The new password is not the same in both repetitions')})
     }  else if (!grepl(pattern = '[[:upper:]]', input$new_password1) && !grepl(pattern = '[[:lower:]]', input$new_password1) && !grepl(pattern = '[[:digit:]]', input$new_password1) && !grepl(pattern = '[^[:alnum:]]', input$new_password1)) {
       output$change_pass <- renderUI({HTML('New password should contain lower & upper case letter, digit and special character [eg. JiX45%5G]')})
     } else  {
       output$change_pass <- renderUI({HTML('Password changed successfully')})
       dbGetQuery(connec, paste0("UPDATE users SET password = '",password_store(input$new_password1),"' WHERE user_login = '", toupper(account_data$USER) ,"';"))
       
       
     }
     
     iv <- InputValidator$new()
     iv$add_rule("old_password", sv_required())
     iv$add_rule("new_password1", sv_required())
     iv$add_rule("new_password2", sv_required())
     iv$enable()
     
     
     
   })
   
   
   
   
   
   
   ###################################################################
   
   observeEvent(input$remind, {
     
     
     showModal(modalDialog(
       
       title = 'Remind password',
       br(),
       textInput("remind_login", "Login", 
                 value = "", placeholder = 'User'),
       textInput("remind_email", "E-mail", 
                 value = "", placeholder = 'user@mail.com'),
       
       checkboxInput('remind_robot', 'I am not a robot!', value = FALSE),
       br(),
       
       uiOutput('decision_rem'),
       actionButton('remind_action', 'Remind', width = '300px'),
       
       
       uiOutput('code_box'),
       uiOutput('code_button'), size = 'xl'
       
       
       
       
       
     ) 
     )
     
     
   })
   
   
   
   observeEvent(input$remind_action, {
      if (!toupper(input$remind_login) %in% toupper(dbGetQuery(connec, 'SELECT user_login FROM users;')[[1]])) {
       output$decision_rem <- renderUI({HTML('User does not exist')})
     } else if (!toupper(input$remind_email) %in% toupper(dbGetQuery(connec, 'SELECT email FROM users;'))[[1]]) {
       output$decision_rem <- renderUI({HTML('Wrong e-mail')})
     } else if (input$remind_robot == FALSE) {
       output$decision_rem <- renderUI({HTML('Are you a robot?')})
     } else  {
       output$decision_rem <- renderUI({HTML('A new password has been sent to your e-mail')})

       samp<-c(2:9,letters,LETTERS)
       code_remind <- as.character(paste(sample(samp,8),collapse=""))
       
       dbGetQuery(connec, paste0("UPDATE users SET change_code = '",as.character(code_remind),"' WHERE user_login = '", toupper(input$remind_login) ,"';"))

       c = paste('python3 $(pwd)/emails/remind_password.py ', as.character(sys$V2[sys$V1 == 'sys1']), as.character(sys$V2[sys$V1 == 'sys2']), as.character(sys$V2[sys$V1 == 'sys3']) , as.character(input$remind_email), as.character(input$remind_login), as.character(code_remind))

       system(c)
       
       rm(code_remind)
       
       
       
       output$code_box <- renderUI({
         
         textInput("code_val", "Code:", 
                   value = "", placeholder = '8-character code')
         
       })
       
       output$code_button <- renderUI({
         
         actionButton('code_check', 'Confirm', width = '300px')
         
       })
       
     }
     
     iv <- InputValidator$new()
     iv$add_rule("remind_login", sv_required())
     iv$add_rule("remind_email", sv_required())
     iv$add_rule("remind_email", sv_email())
     iv$enable()
     
   })
   
  
   
   
   observeEvent(input$code_check, {
     
     if (input$code_val == dbGetQuery(connec, paste0("SELECT change_code FROM users WHERE user_login = '", toupper(input$remind_login) ,"' ;"))[[1]]) {
     
     showModal(modalDialog(
       
       title = 'Change password',
       br(),
       
       passwordInput("new_password3", "New password 1st", 
                     value = "", placeholder = 'min. 8-character password'),
       
       passwordInput("new_password4", "New password 2nd", 
                     value = "", placeholder = 'min. 8-character password'),
       
       
       uiOutput('change_pass'),
       br(),
       actionButton('change_remind', 'Change', width = '300px'), size = 'xl',
       
       
       
       
       
        ) 
      )
       
     } else {
       
       output$change_pass <- renderUI({HTML('Wrong code. Check email!')})
       
       
     }
     
   })
   
   
   observeEvent(input$change_remind, {
     
       if (nchar(input$new_password3) < 8) {
       output$change_pass <- renderUI({HTML('The new password is too short. The minimum number of characters is 8')})
     } else if (input$new_password3 != input$new_password4) {
       output$change_pass <- renderUI({HTML('The new password is not the same in both repetitions')})
     }  else if (!grepl(pattern = '[[:upper:]]', input$new_password3) && !grepl(pattern = '[[:lower:]]', input$new_password3) && !grepl(pattern = '[[:digit:]]', input$new_password3) && !grepl(pattern = '[^[:alnum:]]', input$new_password3)) {
       output$change_pass <- renderUI({HTML('New password should contain lower & upper case letter, digit and special character [eg. JiX45%5G]')})
     } else  {
       output$change_pass <- renderUI({HTML('Password changed successfully')})
       dbGetQuery(connec, paste0("UPDATE users SET password = '",password_store(input$new_password3),"' WHERE user_login = '", toupper(input$remind_login) ,"';"))
  
       
     }
     
     iv <- InputValidator$new()
     iv$add_rule("new_password3", sv_required())
     iv$add_rule("new_password4", sv_required())
     iv$enable()
     
     
     
   })
   
   
   observeEvent(input$register, {
     
     
     showModal(modalDialog(
       title = 'Register account',
       h4('Provide data:'),
       br(),
       textInput("user_login", "Login", 
                 value = "", placeholder = 'User'),
       passwordInput("user_password", "Password", 
                 value = "", placeholder = 'min. 8-character password'),
       textInput("user_email", "E-mail", 
                 value = "", placeholder = 'user@mail.com'),
       textInput("user_name", "Name", 
                 value = NULL, placeholder = 'James'),
       textInput("user_surname", "Surname", 
                 value = NULL, placeholder = 'Smith'),
       textInput("user_affiliation", "Anffiliation", 
                 value = NULL, placeholder = 'company or institute name'),
       textInput("user_country", "Country", 
                 value = NULL, placeholder = 'country name'),
       textInput("user_ORCID", "ORCID", 
                 value = NULL, placeholder = '0000-0000-0000-0000'),
       textInput("user_link", "Other links", 
                 value = NULL, placeholder = 'github, linkedin, ...'),
       checkboxInput('rodo', 'Data administrator ....', value = FALSE),
       br(),
       
      uiOutput('decision'),
      actionButton('account_create', 'Create', width = '300px'), size = 'xl',

       
       
          
       
     ) 
     )
     
     observeEvent(input$account_create, {
       if (toupper(input$user_login) %in% toupper(dbGetQuery(connec, 'SELECT user_login FROM users;')[[1]])) {
         output$decision <- renderUI({HTML('Login name in use. Change user login')})
       }  else if (nchar(input$user_login) < 4) {
         output$decision <- renderUI({HTML('Login is too short. Minimum number of characters is 4')})
       } else if (nchar(input$user_password) < 8) {
         output$decision <- renderUI({HTML(paste('Password is too short. Minimum number of characters is 8',input$user_password))})
       } else if (!grepl(pattern = '[[:upper:]]', input$user_password) || !grepl(pattern = '[[:lower:]]', input$user_password) || !grepl(pattern = '[[:digit:]]', input$user_password) || !grepl(pattern = '[^[:alnum:]]', input$user_password)) {
         output$decision <- renderUI({HTML('Password should contain lower & upper case letter, digit and special character [eg. JiX45%5G]')})
       } else if (nchar(input$user_email) == 0 || !grepl(pattern = '@', input$user_email)) {
         output$decision <- renderUI({HTML('Wrong e-mail')})
       } else if (input$rodo == FALSE) {
         output$decision <- renderUI({HTML('Read the data policy and tick the box [x] if you agree')})
       } else  {
         output$decision <- renderUI({HTML('Account created successfully. You can login & start your analysis with data exploration')})
         dbGetQuery(connec, paste0("INSERT INTO users (user_name, set_id, date) VALUES ('",account_data$USER,"','" , row, "','" , Sys.time(),"');"))
         
         sys <- read.csv('../../requirements_file/system_files', sep = '=', header = F)
         c = paste('python3 $(pwd)/emails/account_create.py ', as.character(sys$V2[sys$V1 == 'sys1']), gsub('"', '', sys$V2[sys$V1 == 'sys2']), as.character(sys$V2[sys$V1 == 'sys3']) , as.character(input$user_email), as.character(input$user_login))
         
         system(c)
         
         
       }
       
       iv <- InputValidator$new()
       iv$add_rule("user_login", sv_required())
       iv$add_rule("user_password", sv_required())
       iv$add_rule("user_email", sv_required())
       iv$add_rule("rodo", sv_required())
       iv$add_rule("user_email", sv_email())
       iv$enable()
       
       
     })
     
     
     
   })
   
   
   observeEvent(input$status, {
     if (input$status == 'disease') {
       output$disease_name <- renderUI({
         textInput("disease_name", "Disease name", 
                   value = "", placeholder = 'Alzheimer, Ovarian cancer, ...')
         
       })
     } else {
       
       output$disease_name <- renderUI({
       })
       
     }
   })
   
   
   
 
 
    ################################################
  
  {
    #SINGLE CELL
    

  
  
  observeEvent(input$start, {
    
 

    CHECK <- TRUE
    
    while(CHECK) {
      samp<-c(2:9,letters,LETTERS)
      code <- paste(sample(samp,10),collapse="")
      
      if (TRUE %in% grepl(pattern = code, x = list.files(file.path('../../projects'))) && TRUE %in% grepl(pattern = code, x = list.files(file.path('../../validated')))) {
        CHECK = TRUE
      } else {
        CHECK = FALSE
      }
      
    }
    
    directory = as.character(paste0('SC_',code))
    
    
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
    
    
    if (input$status != 'disease') {
      disease_name = ''
    } else {
      disease_name = input$disease_name
    }
  
    
    conf = paste0(
      'user_name=', as.character(account_data$USER), '\n',
      'project_name=', as.character(gsub(pattern = ' ', '_', x =  input$project_name)), '\n',
      'email=', as.character(dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(account_data$USER) ,"';"))[[1]]), '\n',
      'source=', as.character(input$source), '\n',
      'species=', as.character(input$species), '\n',
      'sex=', as.character(input$sex), '\n',
      'tissue=', as.character(input$tissue), '\n',
      'tissue_affiliation=', as.character(gsub(pattern = ' ', '_', x =  input$affiliation)), '\n',
      'cell_disease_status=', as.character(input$status), '\n',
      'disease_name=', as.character(disease_name), '\n',
      'cell_development_status=', as.character(input$development) ,' \n',
      'data_format=', as.character(input$input) ,' \n',
      'library=', as.character(input$library) ,' \n',
      'marker_type=', as.character(input$naming) ,' \n',
      'read_length=', as.character(input$reads) ,' \n',
      'id=', as.character(code) ,'\n',
      'directory_name=',directory, '\n',
      'cell=',input$cell_n, '\n',
      'description=',as.character(input$description)
      
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
    iv$add_rule('disease_name', sv_required())
    iv$enable()
    
    if (!is.null(account_data$USER) && length(input$cell_n) > 0 && length(input$files) > 0 && length(input$project_name) > 0 && length(input$description) > 0 && length(input$affiliation) > 0 && length(input$source) > 0) {
        
      
      output$run_buttom <- renderUI({
        actionButton("start", "Run", icon = icon('play')) })
      output$accepted <- renderUI({ 
        'All data correct. Run can run analysis.'
        
      })
    } else if (is.null(account_data$USER)){
      
      output$accepted <- renderUI({ 
        'Run analysis required login ...'
      })
      
    } else {
   
      
      output$accepted <- renderUI({ 
        'Provide all required data ...'
        
      })
      
    }
    
  
    
    
  })
 
  #END single cell
  
  }
   
   
   
  {
  # DATA VALIDATION
    

  
  observeEvent(input$search, {

    if (!is.null(account_data$USER) && nchar(input$id_code) == 10 && TRUE %in% grepl(pattern = input$id_code, x = list.files(file.path('../../projects'))) && TRUE %in% grepl(pattern = 'SC_', x = list.files(file.path('../../projects')))) {
      
      
      
      
      
      output$project_validation <- renderTable({
        
        tmp <-  read.csv(file.path('../../projects/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))], 'config'), sep = '=', header = FALSE)
        colnames(tmp) <- c('Variables', 'Value')
        tmp$Value <- gsub('/', '', tmp$Value)
        
        tmp
        
        })
      
      output$project_validation_report <- renderUI({
        
        actionButton("report_valid", 'View report', icon = icon('bars'))
        
        
        })
      
      
      output$disclaimer <- renderUI({
        
        
        HTML(paste('Data validation.',
                   'If the output of the analysis is appropriate you can share data for other users. If you choose "Yes" the data will be available in "Data exploration" tab where you can obtain more data [analysis output files in hdf5 formats]. The data can be always removed after contact with us. If you choose "No" data will be removed. If you will not make choice the data will remove in a few days.', sep="<br/>"))
        
        
        })
      
     

    
    } else if (is.null(account_data$USER)){
      
      output$disclaimer <- renderUI({ 
        'Data validation required login ...'
      })
      
    } else {
      
      
      output$disclaimer <- renderUI({ 
        HTML(paste('No projects found.',
                   'Possible scenarios:',
                   '-wrong project ID',
                   '-analysis failed due to poor data quality',
                   '- the input data was in the wrong format', 
                   '-selected analysis options were wrongly matched to the dataset',
                   'You can try to re-enable analytics or contact us for guidance on a failed analytics. Contact details can be found in the Contact tab.', sep="<br/>"))
        


        
      })
      
    }
    
    
  })
  
    
  observeEvent(input$report_valid, {
    
    browseURL(file.path('../../projects/', list.files(file.path('../../projects'))[grepl(pattern = input$id_code, x = list.files(file.path('../../projects')))], 'Report.html'))
    
      
      
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
        
        if (grepl(pattern = 'SC_', x = inFile)) {
     
        tmp <-  read.csv(file.path(inFile, 'config'), sep = '=', header = FALSE)
    
        
        
        dbGetQuery(connec, paste0("INSERT INTO sc_rna_seq (id, user_name, project_name, source, species, sex, tissue, tissue_affiliation, cell_disease_status, disease_name, cell_development_status, data_format, library, marker_type, read_length, cells_after_qc, subclasses_number, subtypes_number, description) VALUES 
                                  ('",gsub('/', '', tmp$V2[tmp$V1 == 'id']),"','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'user_name']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'project_name']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'source']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'species']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'sex']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'tissue']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'tissue_affiliation']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'cell_disease_status']), "','" 
                                  ,stringr::str_to_title(tmp$V2[tmp$V1 == 'disease_name']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'cell_development_status']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'data_format']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'library']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'marker_type']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'read_length']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'cells_after_qc']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'subclasses_number']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'subtypes_number']), "','" 
                                  ,gsub('/', '', tmp$V2[tmp$V1 == 'description']), "');"))
        


        file.copy(inFile, destDir, recursive=TRUE)
        
        unlink(inFile, recursive=TRUE)
        
        output$project_validation <- renderUI({ 
          HTML(paste('','',
                     'Project accepted for publication.', sep="<br/>"))
        })

        } else if (grepl(pattern = 'RS_', x = input$id_code)) {
          
          NULL
          
        } else if (grepl(pattern = 'MT_', x = input$id_code)) {
          
          NULL
          
        }
        
      
   
        
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
  # END VALIDATION
  }
  
   {
    ############################# LOADING POSTGRES DATA #######################################
    # SINGLE CELL DATA
    
    
 
    values <- reactiveValues()
    df <- as.data.frame(dbGetQuery(connec, 'SELECT * FROM sc_rna_seq;'))
    
    button <- function(tbl){
      function(i){
        sprintf(
          '<button id="button_%s_%s" type="button" onclick="%s">Explore</button>', 
          tbl, i, "Shiny.setInputValue('button', this.id);")
        
      }
    }
    
    df <- df[, c('id', 'species', 'sex', 'tissue', 'tissue_affiliation', 'cell_disease_status', 'disease_name', 'cell_development_status', 'cells_after_qc', 'subclasses_number', 'subtypes_number')]
    colnames(df) <- c('id', 'species', 'sex', 'tissue', 'affiliation', 'disease', 'disease_name', 'status', 'cells', 'subclasses', 'subtypes')
    
   
    df <- cbind(df, 
                  action = sapply(df$id, button("df")), 
                  stringsAsFactors = FALSE)
    
    

    species <- df$species
    sex <- df$sex
    affiliation <- df$affiliation
    tissue <- df$tissue
    disease <- df$disease
    status <- df$status
    

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

    
    observeEvent(c(input$species_filter, input$sex_filter, input$tissue_filter, input$affiliation_filter, input$disease_filter, input$status_filter), {
    
    values$df <- df
      
    if (input$species_filter != 'all') {
      
      values$df <- values$df[values$df$species %in% input$species_filter,]
  
      
    } 
      
    if (input$sex_filter != 'all') {
     
      values$df <- values$df[values$df$sex %in% input$sex_filter,]
    }
    
   
    if (input$tissue_filter != 'all') {
      
      values$df <- values$df[values$df$tissue %in% input$tissue_filter,]
    }
    
   
    if (input$affiliation_filter != 'all') {
     
      values$df <- values$df[values$df$affiliation %in% input$affiliation_filter,]
    } 
    
    
    if (input$disease_filter != 'all') {
      
      values$df <- values$df[values$df$disease %in% input$disease_filter,]
    } 
    
  
    if (input$status_filter != 'all') {
      
      values$df <- values$df[values$df$status %in% input$status_filter,]
    } 
    
    })
    
    
    
   
       output$my_datatable <- DT::renderDataTable(values$df,
                                               server = FALSE, escape = FALSE, selection = 'none', rownames = FALSE)


       
       
          
          
       
         observeEvent(input[["button"]], {
         splitID <- strsplit(input[["button"]], "_")[[1]]
         tbl <- splitID[2]
         row <- splitID[3]
         login_values$row <- row
         
         
         output$pd <- renderUI({HTML("<h4>Project description</h4> <br>")})
         
         
         output$ar <- renderUI({h4("Analysis report")})
         
         
         
         
       
         info <- dbGetQuery(connec, paste0("SELECT * FROM sc_rna_seq where id='",row,"';"))
         
         output$user_check_button <- renderUI({HTML(paste0('Author: ', actionLink('user_check', as.character(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]))))})
         output$project_description <- renderText(dbGetQuery(connec, paste0("SELECT description FROM sc_rna_seq where id='",row,"';"))[[1]])
         
         info <- info[, c('project_name', 'id', 'source', 'species', 'sex', 'tissue', 'tissue_affiliation', 'cell_disease_status', 'disease_name', 'cell_development_status', 'cells_after_qc', 'subclasses_number', 'subtypes_number')]
         info <- as.data.frame(t(info))
         info$n <- toupper(rownames(info))

         colnames(info) <- c('#','Data information:')
         info <- info[,c('Data information:' ,'#')]
         
         # output$info_out 
         
         output$explor1 <-  renderTable({info})
         

         
         output$login_element1_download <- login_values$login_element1_download
         
         output$login_element2_download <-login_values$login_element2_download
         
         
         output$down_butt <- login_values$down_butt
         
                                                
         #downloaded function
         output$downloadData <- downloadHandler(
           filename <- function() { paste0(list.files('../../validated/')[grep(paste0('SC_',row), list.files('../../validated/'))], '.h5') },
           
           content <- function(file) {
             dbGetQuery(connec, paste0("INSERT INTO downloads (user_name, set_id, date) VALUES ('",account_data$USER,"','" , row, "','" , Sys.time(),"');"))
             
             file.copy(from = file.path('../../validated/', list.files('../../validated/')[grep(paste0('SC_',row), list.files('../../validated/'))], 'data.h5'), to = file, recursive = TRUE)
           },
          
          )
         
         
       
         
         output$report_b <- renderUI({
           
           actionButton("report", 'View report', icon = icon('bars'), width = '200px')
           
         })
                                                
         observeEvent(input$report, {
           
         browseURL(file.path('../../validated/', list.files('../../validated/')[grep(paste0('SC_',row), list.files('../../validated/'))], 'Report.html'))})
         

         
         path = file.path('../../validated/', list.files('../../validated/')[grep(paste0('SC_',row), list.files('../../validated/'))], 'data.h5')
         cells_meta <- h5read(path, "metadata/cells_meta")
         
         output$explor3 <- renderUI({
           selectInput('plot_var', 'Choose data presentation:', choices = list('UMAP_subclasses' = 1, 'UMAP_subtypes' = 2, 'HDMAP_subtypes' = 3), selected = 1, multiple = FALSE, selectize = TRUE)
           
         })
         
         
         observeEvent(input$user_check, {
           
           
           showModal(modalDialog(
             
             tags$h2('User'),
             
             tags$h3(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]),
             br(),
             HTML({paste0('<strong> Name:  ', dbGetQuery(connec, paste0("SELECT name FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Surname:  ', dbGetQuery(connec, paste0("SELECT surname FROM users WHERE user_login = '"  ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Contact:  ', dbGetQuery(connec, paste0("SELECT email FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Affiliation:  ', dbGetQuery(connec, paste0("SELECT affiliation FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Country:  ', dbGetQuery(connec, paste0("SELECT country FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> ORCID:  ', dbGetQuery(connec, paste0("SELECT ORCID FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Other links:  ', dbGetQuery(connec, paste0("SELECT links FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}),
             HTML({paste0('<strong> Account type:  ', dbGetQuery(connec, paste0("SELECT acc_type FROM users WHERE user_login = '", toupper(dbGetQuery(connec, paste0("SELECT user_name FROM sc_rna_seq WHERE id = '", row ,"';"))[[1]]) ,"';"))[[1]], '</strong> </br>')}), size = 'm',
             
             
           ) 
           )
           
           
           
         
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
             
             subclass_markers <- h5read(path, "markers/subclass_markers")
             
             clusters_sub <- c('all', unique(subclass_markers$cluster))
             names(clusters_sub) <- clusters_sub
             subclass_sub <- c('all', unique(subclass_markers$subclass))
             names(subclass_sub) <- subclass_sub
             
             
             output$sub1 <- renderUI({selectInput('cluster_sub', 'Cluster:', choices = clusters_sub, selected = 'all', multiple = FALSE, selectize = TRUE)})
             output$sub2 <- renderUI({selectInput('subclass_sub', 'Subclass:', choices = subclass_sub, selected = 'all', multiple = FALSE, selectize = TRUE)})
             output$sub3 <- renderUI({selectInput('pval_sub', 'p_value:', choices = list('0.05' = 0.05, '0.01' = 0.01, '0.001' = 0.001, '0.0001' = 0.0001), selected = 0.05, multiple = FALSE, selectize = TRUE)})
             
         
             observeEvent(c(input$cluster_sub,input$subclass_sub, input$pval_sub), {
               
               if (input$cluster_sub != 'all' && input$subclass_sub != 'all' && input$pval_sub != 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$cluster) %in% as.character(input$cluster_sub),]
                 values$subclass_markers <- values$subclass_markers[as.numeric(values$subclass_markers$p_val) <= as.numeric(input$pval_sub),]
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$subclass) %in% as.character(input$subclass_sub),]
                 
                 
               } else  if (input$cluster_sub == 'all' && input$subclass_sub != 'all' && input$pval_sub != 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.numeric(values$subclass_markers$p_val) <= as.numeric(input$pval_sub),]
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$subclass) %in% as.character(input$subclass_sub),]
                 
                 
               } else  if (input$cluster_sub == 'all' && input$subclass_sub == 'all' && input$pval_sub != 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.numeric(values$subclass_markers$p_val) <= as.numeric(input$pval_sub),]
                 
                 
               } else  if (input$cluster_sub == 'all' && input$subclass_sub != 'all' && input$pval_sub == 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$subclass) %in% as.character(input$subclass_sub),]
                 
                 
               } else  if (input$cluster_sub != 'all' && input$subclass_sub == 'all' && input$pval_sub == 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$cluster) %in% as.character(input$cluster_sub),]

               } else  if (input$cluster_sub != 'all' && input$subclass_sub != 'all' && input$pval_sub == 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$cluster) %in% as.character(input$cluster_sub),]
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$subclass) %in% as.character(input$subclass_sub),]
                 
                 
               } else  if (input$cluster_sub != 'all' && input$subclass_sub == 'all' && input$pval_sub != 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 values$subclass_markers <- values$subclass_markers[as.character(values$subclass_markers$cluster) %in% as.character(input$cluster_sub),]
                 values$subclass_markers <- values$subclass_markers[as.numeric(values$subclass_markers$p_val) <= as.numeric(input$pval_sub),]

               }  else  if (input$cluster_sub == 'all' && input$subclass_sub == 'all' && input$pval_sub == 'all') {
                 
                 values$subclass_markers <- subclass_markers
                 
               } 
             })
             
             
             
             
             
             
             output$explor6 <- DT::renderDataTable({
                    datatable(values$subclass_markers, options = list(pageLength = 10),
                              rownames= FALSE)

             })



           } else if (input$markers_var == 2) {
             
             
             subtypes_markers <- h5read(path, "markers/subtypes_markers")
             
             clusters_typ <- c('all', unique(subtypes_markers$subtypes))
             names(clusters_typ) <- clusters_typ
             
             
             
             output$sub1 <- renderUI({selectInput('subtypes_sub', 'Subtype:', choices = clusters_typ, selected = 'all', multiple = FALSE, selectize = TRUE)})
             output$sub2 <- renderUI({selectInput('pval_typ', 'p_value:', choices = list('0.05' = 0.05, '0.01' = 0.01, '0.001' = 0.001, '0.0001' = 0.0001), selected = 0.05, multiple = FALSE, selectize = TRUE)})
             output$sub3 <- renderUI({})
             
             
             observeEvent(c(input$subtypes_sub, input$pval_typ), {
               
               if (input$subtypes_sub == 'all' && input$pval_typ == 'all') {
                 
                 values$subtypes_markers <- subtypes_markers
                
               } else  if (input$subtypes_sub == 'all' && input$pval_typ != 'all') {
                 
                 values$subtypes_markers <- subtypes_markers
                 values$subtypes_markers <- values$subtypes_markers[as.numeric(values$subtypes_markers$p_val) <= as.numeric(input$pval_typ),]
                 
               } else  if (input$subtypes_sub != 'all' && input$pval_typ == 'all') {
                 
                 values$subtypes_markers <- subtypes_markers
                 values$subtypes_markers <- values$subtypes_markers[as.character(values$subtypes_markers$subtypes) %in% as.character(input$subtypes_sub),]

               } else  if (input$subtypes_sub != 'all' && input$pval_typ != 'all') {
                 values$subtypes_markers <- subtypes_markers
                 values$subtypes_markers <- values$subtypes_markers[as.character(values$subtypes_markers$subtypes) %in% as.character(input$subtypes_sub),]
                 values$subtypes_markers <- values$subtypes_markers[as.numeric(values$subtypes_markers$p_val) <= as.numeric(input$pval_typ),]
                 
               }
             })
             
             
           
            


             output$explor6 <- DT::renderDataTable({
               datatable(values$subtypes_markers, options = list(pageLength = 10),
                         rownames= FALSE)

             })



           } else if (input$markers_var == 3) {
             

             cssg_markers <- h5read(path, "markers/CSSG")
             
             clusters_c <- c('all', unique(cssg_markers$cluster))
             names(clusters_c) <- clusters_c
             subclass_c <- c('all', unique(cssg_markers$subclass))
             names(subclass_c) <- subclass_c
             
             
             output$sub1 <- renderUI({selectInput('cluster_c', 'Cluster:', choices = clusters_c, selected = 'all', multiple = FALSE, selectize = TRUE)})
             output$sub2 <- renderUI({selectInput('subclass_c', 'Subclass:', choices = subclass_c, selected = 'all', multiple = FALSE, selectize = TRUE)})
             output$sub3 <- renderUI({selectInput('pval_c', '% unexplained cells:', choices = list('5%' = 0.05, '1%' = 0.01, '0,1%' = 0.001, '0,001%' = 0.0001), selected = 0.1, multiple = FALSE, selectize = TRUE)})
             
             
             
             observeEvent(c(input$cluster_c, input$subclass_c, input$pval_c), {
           
               
               if (input$cluster_c != 'all' && input$subclass_c != 'all' && input$pval_c != 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$cluster) %in% as.character(input$cluster_c),]
                 values$cssg_markers <- values$cssg_markers[as.numeric(values$cssg_markers$loss_pval) <= as.numeric(input$pval_c),]
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$subclass) %in% as.character(input$subclass_c),]
                 
                 
               } else  if (input$cluster_c == 'all' && input$subclass_c != 'all' && input$pval_c != 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.numeric(values$cssg_markers$loss_pval) <= as.numeric(input$pval_c),]
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$subclass) %in% as.character(input$subclass_c),]
                 
                 
               } else  if (input$cluster_c == 'all' && input$subclass_c == 'all' && input$pval_c != 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.numeric(values$cssg_markers$loss_pval) <= as.numeric(input$pval_c),]
                 
                 
               } else  if (input$cluster_c == 'all' && input$subclass_c != 'all' && input$pval_c == 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$subclass) %in% as.character(input$subclass_c),]
                 
                 
               } else  if (input$cluster_c != 'all' && input$subclass_c == 'all' && input$pval_c == 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$cluster) %in% as.character(input$cluster_c),]
                 
               } else  if (input$cluster_c != 'all' && input$subclass_c != 'all' && input$pval_c == 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$cluster) %in% as.character(input$cluster_c),]
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$subclass) %in% as.character(input$subclass_c),]
                 
                 
               } else  if (input$cluster_c != 'all' && input$subclass_c == 'all' && input$pval_c != 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 values$cssg_markers <- values$cssg_markers[as.character(values$cssg_markers$cluster) %in% as.character(input$cluster_c),]
                 values$cssg_markers <- values$cssg_markers[as.numeric(values$cssg_markers$loss_pval) <= as.numeric(input$pval_c),]
                 
               }  else  if (input$cluster_c == 'all' && input$subclass_c == 'all' && input$pval_c == 'all') {
                 
                 values$cssg_markers <- cssg_markers
                 
               } 
             })
             
            output$explor6 <- DT::renderDataTable({
               datatable(values$cssg_markers, options = list(pageLength = 5),
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
             
             genes_select <- c()
             for (col in colnames(tmp)) {
             genes_select <- c(genes_select, as.vector(stringr::str_split(col, ' ', simplify = FALSE)[[1]][2:3]))
             }
             
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
           size = 'm',
           easyClose = TRUE,
           footer = NULL,
           HTML({paste0("Dataset ", row, ' was chosen.')}),
           br(),
           HTML({paste0('Go to [Info] | [Visualization] | [Markers] | [Gene_viewe] | [Interactions] tabs and start exploring the data')})
           
         ),
        
         )
       })

    
       #END OF EXPLORATION SINGLE CELL DATA   
   }
   
  
}


shinyApp(ui = ui, server = server)