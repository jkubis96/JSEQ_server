 
Sys.setenv(R_INSTALL_STAGED = FALSE)
install.packages('patchwork', version = '1.1.1')
install.packages('remotes')
install.packages('viridis', version = "0.5.1")
install.packages('textclean', version = "0.9.3")
install.packages('ape', version = '5.4.1')
install.packages('umap') 
BiocManager::install("MAST")
BiocManager::install("rhdf5")
install.packages("shiny")
install.packages("shinydashboard")
install.packages('blastula')
install.packages('shinyvalidate')
