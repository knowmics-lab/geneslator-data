library(jsonlite)
library(dplyr)
library(tidyr)

source("Script/process_functions.R")
source("Script/mapping_functions.R")
source("Script/merge_functions.R")

#Read metadata
metadata <- fromJSON("Input/sources.json", simplifyVector = TRUE)
global.metadata <- metadata$global_sources
species.metadata <- metadata$species_sources

#Prepare annotation data
annotation.data.list <- list()

#Read global annotation data
global.ncbi.archive.data <- readRDS("Input/ncbi_archive.rds")
global.ncbi.discontinued.data <- readRDS("Input/ncbi_discontinued.rds")
global.ncbi.orthologs <- readRDS("Input/ncbi_orthologs.rds")
global.alliance.orthologs <- readRDS("Input/alliance_orthologs.rds")

#Load GO dictionary
go.dictionary <- readRDS("Input/go_dictionary.rds")
go.dictionary <- process.go.dictionary(go.dictionary)

#Create a taxonomy table
taxonomy.table <- species.metadata[,c("species","official_name","taxid")]

#Create species' annotation folder, if it does not exist
if(!dir.exists("Annotations")){
  mkdirs("Annotations")
}

#---------GENOME ANNOTATIONS-----------
for(species in species.metadata$species){
  
  print(paste0("Build genome annotation data for ",species,"..."))
  
  #Get species DB name and taxonomy id
  speciesdb.name <- species.metadata[species.metadata$species==species,"speciesdb_name"]
  taxid <- unlist(species.metadata[species.metadata$species==species,"taxid"])
  
  ###-----NCBI DATA------
  #Process current data
  ncbi.data <- readRDS(paste0("Input/",species,"/ncbi.rds"))
  ncbi.data <- process.ncbi.data(ncbi.data,taxid,speciesdb.name)
  #Process archive data
  ncbi.archive.data <- process.ncbi.archive.data(global.ncbi.archive.data,global.ncbi.discontinued.data,taxid,speciesdb.name)
  ncbi.replaced.data <- ncbi.archive.data[[1]]
  ncbi.discontinued.data <- ncbi.archive.data[[2]]
  #Merge current and archive data
  ncbi.data <- merge.ncbi.data(ncbi.data,ncbi.replaced.data,ncbi.discontinued.data,speciesdb.name)

  ###------SpeciesDB DATA-------
  if(!is.na(speciesdb.name)){
    speciesdb.data <- readRDS(paste0("Input/",species,"/speciesdb.rds"))
    speciesdb.data <- process.speciesdb.data(speciesdb.data,species)
  }

  ###-----Ensembl DATA----
  #Read Human HCOP data for zebrafish for symbol double check
  if(species=="Zebrafish"){
    hcop.data <- readRDS(paste0("Input/Human/speciesdb_orthologs.rds"))
  } else {
    hcop.data <- NULL
  }
  #Process current data
  ensembl.data <- readRDS(paste0("Input/",species,"/ensembl.rds"))
  ensembl.data <- process.ensembl.data(ensembl.data,speciesdb.name,ncbi.data,special.data,hcop.data,taxid,is.archive=F)
  #Merge current data with archive data
  ensembl.data[["ENSEMBLOLD ENSEMBL"]] <- NA
  ensembl.archive.data <- list()
  ensembl.archive.files <- list.files(paste0("Input/",species,"/Archives"),pattern = "ensembl_[0-9]+")
  for(file in ensembl.archive.files){
    archive.version <- strsplit(file,"_|\\.rds")[[1]][2]
    ensembl.archive.data[[archive.version]] <- readRDS(paste0("Input/",species,"/Archives/",file))
  }
  ensembl.data <- merge.with.ensembl.archive.data(ensembl.data,ensembl.archive.data,
                    species,speciesdb.name,ncbi.data,special.data,hcop.data,taxid)
  #Merge with GRCh37 data (only for Human)
  if(species=="Human"){
    ensembl.grch37.data <- readRDS(paste0("Input/",species,"/Archives/ensembl_grch37.rds"))
    ensembl.grch37.data <- process.ensembl.grch37.data(ensembl.grch37.data,speciesdb.name)
    ensembl.data <- merge.with.ensembl.grch37.data(ensembl.data,ensembl.grch37.data,ncbi.data,speciesdb.name)
  }

  ###-----Uniprot DATA----
  uniprot.data <- readRDS(paste0("Input/",species,"/uniprot.rds"))
  uniprot.data <- process.uniprot.data(uniprot.data,speciesdb.name,species)

  #Merge annotations
  annotation.data <- merge.databases(ncbi.data,ensembl.data,uniprot.data,speciesdb.data,speciesdb.name,species)
  annotation.data.list[[species]] <- annotation.data
}


#-------FUNCTIONAL ANNOTATIONS---------
for(species in species.metadata$species){
  
  print(paste0("Build functional annotation for ",species,"..."))
  
  #Get species DB name and taxonomy id
  speciesdb.name <- species.metadata[species.metadata$species==species,"speciesdb_name"]
  taxid <- unlist(species.metadata[species.metadata$species==species,"taxid"])
  
  #Orthologs
  #-------NCBI orthologs---------
  ncbi.orthologs <- process.ncbi.orthologs.data(global.ncbi.orthologs,taxid,taxonomy.table,annotation.data.list)
  
  #-----ENSEMBL orthologs-------
  ensembl.orthologs <- readRDS(paste0("Input/",species,"/ensembl_orthologs.rds"))
  if(class(ensembl.orthologs)=="list"){
    ensembl.orthologs <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  }
  ensembl.orthologs <- process.ensembl.orthologs.data(ensembl.orthologs,taxid,taxonomy.table,annotation.data.list)
  
  #-----ALLIANCE orthologs------
  alliance.orthologs <- process.alliance.orthologs.data(global.alliance.orthologs,taxid,taxonomy.table)
  
  #-----SpeciesDB orthologs--------
  if(species=="Human"){
    speciesdb.orthologs <- readRDS(paste0("Input/",species,"/speciesdb_orthologs.rds"))
    speciesdb.orthologs <- process.speciesdb.orthologs.data(speciesdb.orthologs,taxid,taxonomy.table,annotation.data.list)
  } else {
    speciesdb.orthologs <- data.frame()
  }
  
  #Merge orthologs DATA
  orthologs.data <- merge.ortho.databases(ncbi.orthologs,ensembl.orthologs,alliance.orthologs,speciesdb.orthologs,taxid,taxonomy.table)

  #----Gene Ontologies-----
  go.data <- readRDS(paste0("Input/",species,"/go.rds"))
  go.data <- process.go.data(go.data,go.dictionary,annotation.data.list,species)
  
  #----REACTOME-------
  reactome.data <- readRDS(paste0("Input/",species,"/reactome.rds"))
  reactome.data <- process.reactome.data(reactome.data,annotation.data.list,species)
  
  #---WIKIPATHWAYS----
  wikipathways.data <- readRDS(paste0("Input/",species,"/wikipathways.rds"))
  wikipathways.data <- process.wikipathways.data(wikipathways.data,annotation.data.list,species)

  #Merge annotation data with orthologs, go and pathway data
  annotation.data <- annotation.data.list[[species]]
  annotation.data <- merge(annotation.data,orthologs.data,all.x=T)
  annotation.data <- merge(annotation.data,go.data,all.x=T)
  annotation.data <- merge(annotation.data,reactome.data,all.x=T)
  annotation.data <- merge(annotation.data,wikipathways.data,all.x=T)
  
  #Write annotation table to file
  write.table(annotation.data,paste0("Annotations/",species,".txt"),sep="\t",quote=F,row.names = F)
  
}
