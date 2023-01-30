library("RPostgreSQL")
library("DBI")

sys <- read.csv('../../requirements_file/system_files', sep = '=', header = F)
drv <- dbDriver("PostgreSQL")
connec <- dbConnect(drv, 
                    dbname = 'jbsda',
                    host = 'localhost', 
                    port = '5432',
                    user = sys$V2[sys$V1 == 'sys4'], 
                    password = sys$V2[sys$V1 == 'sys4'])

