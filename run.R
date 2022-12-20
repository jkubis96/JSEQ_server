library(shiny)

args <- commandArgs()


 
path <- args[6]


runApp(file.path(path, "html/main/main.R"), launch.browser = TRUE, port = 8001)
