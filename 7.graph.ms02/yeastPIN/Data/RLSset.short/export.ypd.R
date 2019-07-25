# Adopted from Qin15
# export single-gene mutant RLS at 30C and YPD for MATalpha
# based on _export_rls20140709.R
# export rls to individual files. I need them for nrmca fittings. 

#################install and declare needed packaqes for survival analysis
rm(list=ls())

library('data.table')
library('RSQLite')


#####################Declare names of files and folders. Move working directory to desired folder. Establish connection to database
File = 'rls.db'
drv <- SQLite()
con <- dbConnect(drv, dbname = File)

######################Create a complete table of conditions that have data in the RLS database. 
######################(Each row in this table will be a unique combination of genotype, mating type, media, and temperature

###get unique conditions from "set" columns (experimental treatments/mutations)
conditions1 = dbGetQuery(con, "
  SELECT DISTINCT set_genotype as genotype, set_mating_type as mat, set_media as media, set_temperature as temp
  FROM result
  WHERE pooled_by = 'file' 	
  ")
### rows in the result table of the database are not mutually exclusive. 
### In some rows, data has sometimes been pooled by genotype, background strain, etc.   

### get unique conditions from "reference" columns (the control conditions/lifespan results 
### for each row of experimental lifespan results. 
### Rows are not mutually exclusive (1 control to many experimental conditions)
conditions2 = dbGetQuery(con, "
  SELECT DISTINCT ref_genotype as genotype, ref_mating_type as mat, ref_media as media, ref_temperature as temp
  FROM result 
  WHERE pooled_by = 'file'
  ")

####combine and take unique conditions from these two
conditions = rbind(conditions1, conditions2)
conditions = unique(conditions)
conditions = conditions[complete.cases(conditions),]
row.names(conditions) = NULL
###renumber the rows. 
### important because future processes will refer to a unique condition by its row number in this table
controlConditions = conditions[conditions$genotype %in% c('BY4741', 'BY4742', 'BY4743'),]			
### create a table of conditions that have WT genotypes	

head(conditions)
table(conditions$media)

### will limit conditions to single-gene, 30C, YPD, MATalpha
conditions = conditions[!is.na(conditions$genotype), ]
conditions = conditions[ conditions$media=='YPD' & conditions$temp==30 & conditions$mat=='MATalpha', ]
conditions$single_element_flag = 0
for( i in 1:length(conditions[,1])) { 
 elements = unlist( strsplit(  conditions$genotype[i], '\\s+',))
 conditions$single_element_flag[i]= length( elements )
}
summary(conditions)
conditions= conditions[conditions$single_element_flag==1 , ]

###Add columns to the conditions data frame. These columns will be filled in by their respective variable: ie. gompertz shape/rate of the lifespan data associated with a given conditions (genotype, mating type, media, temp)
conditions$n = apply(conditions, 1, function(row) 0)
conditions$avgLS = apply(conditions, 1, function(row) 0)
conditions$StddevLS = apply(conditions, 1, function(row) 0)
conditions$medianLS = apply(conditions, 1, function(row) 0)
conditions$gompShape = apply(conditions, 1, function(row) 0)
conditions$gompRate = apply(conditions, 1, function(row) 0)
conditions$gompLogLik = apply(conditions, 1, function(row) 0)
conditions$gompAIC = apply(conditions, 1, function(row) 0)
conditions$weibShape = apply(conditions, 1, function(row) 0)
conditions$weibScale = apply(conditions, 1, function(row) 0)
conditions$weibLogLik = apply(conditions, 1, function(row) 0)
conditions$weibAIC = apply(conditions, 1, function(row) 0)

#########################################Loop through the conditions to get lifespan data 
for (r in 1:length(conditions$genotype)) {
  genotypeTemp = conditions$genotype[r]
  mediaTemp = conditions$media[r]
  temperatureTemp = conditions$temp[r]
  matTemp = conditions$mat[r]
  
  conditionName = apply(conditions[r,1:4], 1, paste, collapse=" ") 				
#### create a string to name a possible output file		
  conditionName = gsub("[[:punct:]]", "", conditionName) 						
#### remove special characters from the name (ie. quotations marks, slashes, etc.)

  genotypeTemp = gsub("'", "''", genotypeTemp)
  genotypeTemp = gsub('"', '""', genotypeTemp)
  
##### Query the database to take data (including lifespan data) from every mutually exclusive row (pooled by file, not genotype, not background, etc). 
##### There will often be multiple rows (representing different experiments) for each unique condition. This analysis pools lifespan data from multiple experiments if the conditions are all the same
  queryStatementSet = paste(
    "SELECT * ",
    "FROM result ",
    "WHERE pooled_by = 'file' AND set_genotype = '", genotypeTemp, "' AND set_mating_type = '", matTemp, "' AND set_media = '", mediaTemp, "' AND set_temperature = '", temperatureTemp, "'", 
    sep = ""
  )
  
  queryStatementRef = paste(		### both the reference and set columns will be searched for matching conditions
    "SELECT * ",
    "FROM result ",
    "WHERE pooled_by = 'file' AND ref_genotype = '", genotypeTemp, "' AND ref_mating_type = '", matTemp, "' AND ref_media = '", mediaTemp, "' AND ref_temperature = '", temperatureTemp, "'", 
    sep = ""
  )
  
  dataListSet = dbGetQuery(con, queryStatementSet)
  dataListRef = dbGetQuery(con, queryStatementRef)
  lifespansChar = unique(c(dataListSet$set_lifespans, dataListRef$ref_lifespans)) ### combine lifespan values for a given condition into a single data structure. (problems of having non-mutually exclusive rows in the ref_lifespans column are overcome by only taking unique groups of lifespans. The assumption is that no two experiments produced identical lifespan curvs)
  
  ##### Database codes the lifespan data for each experiment as a string. So, lifespanChar is a vector of strings
  ##### Convert lifespanChar into a single vector of integers lifespansTemp
  lifespansTemp = c()
  if (length(lifespansChar) > 0) {
    for (s in 1:length(lifespansChar)) {
      lifespansNum = as.numeric(unlist(strsplit(lifespansChar[s], ",")))
      lifespansNum = lifespansNum[!is.na(lifespansNum)]
      if (length(lifespansNum) > 0) {
        lifespansTemp = c(lifespansTemp, lifespansNum)
      } 
    }
  }
  
  if ( length(lifespansTemp) > 3 ) {
    require(stringr)
    conditions$media[r] = str_replace( conditions$media[r], "\\/", "")
    conditions$media[r] = str_replace( conditions$media[r], "\\/", "")
    conditions$media[r] = str_replace( conditions$media[r], "\\+", "")
    conditions$genotype[r] = str_replace( conditions$genotype[r], "\\/", "-")
    conditions$genotype[r] = str_replace( conditions$genotype[r], "\\/", "-")
    conditions$genotype[r] = str_replace( conditions$genotype[r], "\\/", "-")
    conditions$genotype[r] = str_replace( conditions$genotype[r], "\\/", "-")
    conditions$n[r] = length(lifespansTemp)  	#### record number of individuals 
    filename = paste( 'ypd30matalpha/',conditions$genotype[r], '.csv',sep='')
    out = data.frame(lifespansTemp)
    names(out) = c("rls")
      write.csv(out, filename, quote=F, row.names=F)
  }
  
}#outer loop  
    
    
