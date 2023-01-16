install.packages('rlang')
remotes::install_github("stla/NestedMenu")


setwd('/mnt/c/Users/merag/Git/JSEQ_server/app/main')
library(shiny)

runApp("main.R", launch.browser = TRUE, port = 8001, )
shiny::runApp('main.R', port = 1111, host = 'localhost')



library("DBI")
#install.packages('RPostgreSQL')
library("RPostgreSQL")
drv <- dbDriver("PostgreSQL")
connec <- dbConnect(drv, 
                    dbname = 'jbdb',
                    host = 'localhost', 
                    port = '5432',
                    user = 'admin', 
                    password = 'admin')


dbGetQuery(connec, 'CREATE TABLE sc_rna_seq (
	id VARCHAR ( 100 ) NOT NULL,
  user_name VARCHAR ( 100 ) NOT NULL,
	project_name VARCHAR ( 100 ) NOT NULL,
	source VARCHAR ( 500 ) NOT NULL,
  species VARCHAR ( 50 ) NOT NULL,
	sex VARCHAR ( 50 ) NOT NULL,
	tissue VARCHAR ( 50 ) NOT NULL,
	tissue_affiliation VARCHAR ( 200 ) NOT NULL,
  cell_disease_status VARCHAR ( 50 ) NOT NULL,
  disease_name VARCHAR ( 500 ) NOT NULL,
	cell_development_status VARCHAR ( 50 ) NOT NULL,
	data_format VARCHAR ( 50 ) NOT NULL,
	library VARCHAR ( 50 ) NOT NULL,
	marker_type VARCHAR ( 50 ) NOT NULL,
  read_length VARCHAR ( 50 ) NOT NULL,
	cells_after_qc VARCHAR ( 50 ) NOT NULL,
	subclasses_number VARCHAR ( 50 ) NOT NULL,
	subtypes_number VARCHAR ( 50 ) NOT NULL,
	description VARCHAR ( 2000 ) NOT NULL


	

);')


dbGetQuery(connec, 'CREATE TABLE users (
  user_id serial PRIMARY KEY, 
	name VARCHAR ( 500 ) NULL,
	surname VARCHAR ( 500 ) NULL,
  email VARCHAR ( 500 ) NOT NULL,
  affiliation VARCHAR ( 500 ) NULL,
	country VARCHAR ( 500 ) NULL,
	ORCID VARCHAR ( 500 ) NULL,
	links VARCHAR ( 500 ) NULL,
	user_name VARCHAR ( 500 ) NOT NULL,
	user_login VARCHAR ( 500 ) NOT NULL,
	password VARCHAR ( 500 ) NOT NULL,
  RODO VARCHAR ( 500 ) NOT NULL,
  acc_type VARCHAR ( 500 ) NOT NULL,
  change_code VARCHAR ( 500 ) NULL


);')




## FOR FUTURE

df2 <- dbGetQuery(connec, 'SELECT * FROM users;')

dbGetQuery(connec, paste0("SELECT password FROM users WHERE user_name = '",g,"';"))[[1]]


dbGetQuery(connec, paste("INSERT INTO users (name, surname, email, affiliation, country, ORCID, links, user_name, password, RODO, acc_type) VALUES ('fname', 'surnamef', 'emafffil', 'affiliation', 'counfftry', 'ORCffffID', 'links', 'user_name', 'password', 'RODO', 'dd');", sep = ","))


dbGetQuery(connec, 'DROP TABLE users;')




dbDisconnect(connec)
