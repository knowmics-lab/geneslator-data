library(data.table)
library(rtracklayer)
library(jsonlite)
library(curl)
library(rvest)
library(rWikiPathways)
library(R.utils)
library(ontologyIndex)
library(rentrez)
library(xml2)
library(purrr)

source("Script/download_functions.R")

#Adding new species?
new.species <- TRUE

#Load list of URLs
list.urls <- fromJSON("Input/sources.json", simplifyVector = TRUE)
global.urls <- list.urls$global_sources
species.url <- list.urls$species_sources


##-----------GLOBAL DATABASES-----------

#NCBI archive data
ncbi.archive.data <- download.tabular.data(global.urls$ncbi_archive)
ncbi.archive.data <- ncbi.archive.data[ncbi.archive.data$`#tax_id` %in% unlist(species.url$taxid),]
saveRDS(ncbi.archive.data,"Input/ncbi_archive.rds")

#NCBI data about discontinued ids
ncbi.discontinued.data <- query.ncbi.discontinued.data(ncbi.archive.data)
saveRDS(ncbi.discontinued.data,"Input/ncbi_discontinued.rds")

#Get list of Ensembl folders
ensembl.species <- download.tabular.data(paste0(global.urls$ensembl,"current/species_EnsemblVertebrates.txt"))
ensembl.folders <- paste0("gff3/",ensembl.species$`#name`)
ensembl.species <- gsub("_"," ",gsub("_core_.*$","",ensembl.species$other_alignments))
ensembl.species <- paste0(toupper(substr(ensembl.species,1,1)),substr(ensembl.species,2,nchar(ensembl.species)))
names(ensembl.folders) <- ensembl.species

#Get list of EnsemblGenomes folders
ensembl.genome.species <- download.tabular.data(paste0(global.urls$ensemblGenomes,"current/species.txt"))
ens.divisions <- tolower(gsub("Ensembl","",ensembl.genome.species$species))
ens.collections <- gsub("_core_.*$","",ensembl.genome.species$other_alignments)
ensembl.genome.folders <- ifelse(grepl("_collection$",ens.collections),paste0(ens.divisions,"/gff3/",ens.collections,"/",ensembl.genome.species$`#name`),
       paste0(ens.divisions,"/gff3/",ensembl.genome.species$`#name`))
ensembl.genome.species <- gsub("_"," ",gsub("_gca_?[0-9]+.*$","",ensembl.genome.species$`#name`))
ensembl.genome.species <- paste0(toupper(substr(ensembl.genome.species,1,1)),substr(ensembl.genome.species,2,nchar(ensembl.genome.species)))
names(ensembl.genome.folders) <- ensembl.genome.species

#Read UNIPROT proteome species data
uniprot.species <- download.tabular.data(paste0(global.urls$uniprot,"README"),skip="Proteome_ID\t")
uniprot.species <- uniprot.species[uniprot.species$Tax_ID %in% unlist(species.url$taxid),]
uniprot.species <- uniprot.species[!duplicated(uniprot.species$Tax_ID),]

#NCBI orthologs
ncbi.orthologs <- download.tabular.data(global.urls$ncbi_orthologs)
ncbi.orthologs <- ncbi.orthologs[ncbi.orthologs$`#tax_id` %in% unlist(species.url$taxid) & 
                                 ncbi.orthologs$`Other_tax_id` %in% unlist(species.url$taxid),]
saveRDS(ncbi.orthologs,"Input/ncbi_orthologs.rds")

#Alliance orthologs
alliance.orthologs <- download.tabular.data(global.urls$alliance_genome_orthologs)
alliance.orthologs <- alliance.orthologs[alliance.orthologs$Gene1SpeciesTaxonID %in% paste0("NCBITaxon:",unlist(species.url$taxid)) & 
                                         alliance.orthologs$Gene2SpeciesTaxonID %in% paste0("NCBITaxon:",unlist(species.url$taxid)),]
saveRDS(alliance.orthologs,"Input/alliance_orthologs.rds")

#GO data
go.dictionary <- download.go.dictionary(global.urls$go_dictionary)
saveRDS(go.dictionary,"Input/go_dictionary.rds")

#Reactome data
reactome.ncbi <- download.tabular.data(paste0(global.urls$reactome,"NCBI2Reactome_All_Levels.txt"),header=F)
reactome.ensembl <- download.tabular.data(paste0(global.urls$reactome,"Ensembl2Reactome_All_Levels.txt"),header=F)
reactome.plant.ncbi <- download.tabular.data(paste0(global.urls$reactomePlant,"NCBI2PlantReactome_All_Levels.txt"),header=F)
reactome.plant.ensembl <- download.tabular.data(paste0(global.urls$reactomePlant,"Ensembl2PlantReactome_All_Levels.txt"),header=F)


##-------SPECIES SPECIFIC DATABASES---------

list.species <- species.url$species
for(species in list.species)
{
  list.species.urls <- species.url[species.url$species==species,]
  species.scientific.name <- species.url[species.url$species==species,"official_name"]
  species.taxid <- species.url[species.url$species==species,"taxid"][[1]]
  
  #Create species folder, if it does not exist
  if(!dir.exists(paste0("Input/",species))){
    mkdirs(paste0("Input/",species,"/Archives"))
  }
  
  #NCBI
  ncbi.data <- download.tabular.data(list.species.urls$ncbi_current)
  ncbi.data <- ncbi.data[ncbi.data$`#tax_id` %in% species.taxid,]
  saveRDS(ncbi.data,paste0("Input/",species,"/ncbi.rds"))
  
  #ENSEMBL
  if(species.scientific.name %in% names(ensembl.folders)){
    ensembl.url <- global.urls$ensembl
    ensembl.folder <- ensembl.folders[species.scientific.name]
  } else {
    ensembl.url <- global.urls$ensemblGenomes
    ensembl.folder <- ensembl.genome.folders[species.scientific.name]
  }
  remote.ens.version <- as.character(download.tabular.data(paste0(ensembl.url,"VERSION")))
  if(!file.exists(paste0("Input/",species,"/EnsemblVersion.txt"))){
    write.table(as.numeric(remote.ens.version)-1,paste0("Input/",species,"/EnsemblVersion.txt"),
                quote=F,row.names=F,col.names=F)
  }
  local.ens.version <- as.character(read.table(paste0("Input/",species,"/EnsemblVersion.txt")))
  if(local.ens.version!=remote.ens.version){
    list.tags.ens <- c("Name","biotype","gene_id","description")
    list.filters.ens <- c("C_gene_segment","gene","J_gene_segment","lincRNA_gene","miRNA_gene",
      "mt_gene","processed_transcript","pseudogene","RNA","rRNA_gene","snoRNA_gene","snRNA_gene",
      "V_gene_segment","VD_gene_segment","ncRNA_gene")
    ensembl.path <- paste0(ensembl.url,"current/",ensembl.folder)
    link.content <- curl_fetch_memory(paste0(ensembl.path,"/?C=S;O=D"))
    list.links <- read_html(link.content$content) %>% html_elements("a") %>% html_text(trim=T)
    ensembl.file <- list.links[grep("gff3",list.links)[1]]
    ensembl.data <- download.gff.data(paste0(ensembl.path,"/",ensembl.file),list.tags.ens,list.filters.ens)
    saveRDS(ensembl.data,paste0("Input/",species,"/ensembl.rds"))
    write.table(remote.ens.version,paste0("Input/",species,"/EnsemblVersion.txt"),quote=F,row.names=F,col.names=F)
  }
  
  #ENSEMBL ARCHIVES
  ensembl.links <- filter.remote.links(ensembl.url,"release-")
  for(archive.link in ensembl.links){
    ensembl.version <- strsplit(archive.link,"-|/")[[1]][2]
    if(ensembl.version!=remote.ens.version){
      dest.file <- paste0("Input/",species,"/Archives/ensembl_",ensembl.version,".rds")
      if(!file.exists(dest.file)){
        release.folder <- paste0(ensembl.url,archive.link,ensembl.folder)
        ensembl.archive <- download.ensembl.archive.data(release.folder, ensembl.version,
          list.tags.ens, list.filters.ens, species.scientific.name)
        if(!is.null(ensembl.archive)){
          saveRDS(ensembl.archive,paste0("Input/",species,"/Archives/ensembl_",ensembl.version,".rds"))
        }
      }
    }
  }
  
  #ENSEMBL GRCH37
  if(!is.na(list.species.urls$ensembl_grch37)){
    list.tags <- c("gene_id","gene_biotype","gene_name")
    list.filters <- "gene"
    ensembl.grch37.data <- download.gff.data(list.species.urls$ensembl_grch37,list.tags,list.filters)
    saveRDS(ensembl.grch37.data,paste0("Input/",species,"/Archives/ensembl_grch37.rds"))
  }
  
  #UNIPROT
  uniprot.class <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"SUPERREGNUM"]
  uniprot.class <- paste0(toupper(substring(uniprot.class,1,1)),substring(uniprot.class,2,nchar(uniprot.class)))
  uniprot.code <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"Proteome_ID"]
  uniprot.tax <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"Tax_ID"]
  uniprot.url <- paste0(global.urls$uniprot,uniprot.class,"/",uniprot.code,"/",uniprot.code,"_",uniprot.tax,".idmapping.gz")
  uniprot.data <- download.tabular.data(uniprot.url,header=F)
  saveRDS(uniprot.data,paste0("Input/",species,"/uniprot.rds"))
  
  #SPECIES DB
  if(!is.na(list.species.urls$speciesdb)){
    if(species=="Yeast"){
      list.tags <- c("dbxref","so_term_name","Alias","Name","gene")
      list.filters <- c("gene","ncRNA_gene","tRNA_gene","snoRNA_gene","transposable_element_gene","pseudogene",
                    "telomerase_RNA_gene","snRNA_gene","rRNA_gene")
      speciesdb.data <- download.gff.data(list.species.urls$speciesdb,list.tags,list.filters)
    } else if(species %in% c("Human","Mouse")){
      speciesdb.data <- download.tabular.data(list.species.urls$speciesdb)
    } else if(species=="Rat"){
      speciesdb.data <- download.delim.data(list.species.urls$speciesdb,header=T,comment.character="#")
    } else if(species=="Zebrafish"){
      list.tags <- c("Name","Alias","so_term_name","full_name","ID","secondaryIds")
      list.filters <- "gene"
      speciesdb.data <- download.gff.data(list.species.urls$speciesdb,list.tags,list.filters)
    } else if(species=="Fly"){
      speciesdb.file <- filter.remote.links(list.species.urls$speciesdb,"fbgn_annotation_ID")
      speciesdb.data <- download.delim.data(paste0(list.species.urls$speciesdb,speciesdb.file),
                                          header=F,comment.character="#")
    } else if(species=="Worm"){
      speciesdb.data <- fread(paste0("Input/",species,"/c_elegans.PRJNA13758.current.geneIDs.txt.gz"),data.table=F,showProgress=F)
      speciesdb.data[speciesdb.data==""] <- NA
    } else if(species=="Arabidopsis"){
      list.tags <- c("symbol","locus_type","computational_description","ID")
      list.filters <- "gene"
      speciesdb.data <- as.data.frame(readGFF(paste0("Input/",species,"/Araport11_GFF3_genes_transposons.20250813.gff.gz"),columns=character(0),tags=list.tags,
                                      filter=list(type=list.filters)))
      speciesdb.data[speciesdb.data==""] <- NA
    }
    saveRDS(speciesdb.data,paste0("Input/",species,"/speciesdb.rds"))
  }
  
  #ENSEMBL ORTHOLOGS
  if(new.species || local.ens.version!=remote.ens.version){
    ensembl.orthologs.folder <- gsub("gff3/","json/",ensembl.folder)
    ensembl.orthologs.species <- strsplit(ensembl.orthologs.folder,"/")[[1]]
    ensembl.orthologs.species <- ensembl.orthologs.species[length(ensembl.orthologs.species)]
    list.species.taxid <- unlist(species.url$taxid)
    filter.string.ortho.json <- paste0('[.genes[] | .id as $gene_id | .homologues[] | select(.taxonomy_id | IN(',
    paste0(list.species.taxid[!list.species.taxid %in% species.taxid],collapse = ","),
      ')) | {ID: $gene_id, TAXID: .taxonomy_id, ENSEMBL: .stable_id}]')
    ensembl.orthologs.data <- download.json.data(paste0(ensembl.url,"current/",ensembl.orthologs.folder,"/",ensembl.orthologs.species,".json"),
      filter.string.ortho.json)
    saveRDS(ensembl.orthologs.data,paste0("Input/",species,"/ensembl_orthologs.rds"))
  }
  
  #SPECIESDB ORTHOLOGS
  if(!is.na(list.species.urls$speciesdb_orthologs)){
    speciesdb.orthologs.data <- download.tabular.data(list.species.urls$speciesdb_orthologs)
    speciesdb.orthologs.data <- speciesdb.orthologs.data[speciesdb.orthologs.data$ortholog_species %in% species.url$taxid,]
    saveRDS(speciesdb.orthologs.data,paste0("Input/",species,"/speciesdb_orthologs.rds"))
  }
  
  #GO
  go.data <- download.delim.data(list.species.urls$go,header=F,comment.character="!")
  saveRDS(go.data,paste0("Input/",species,"/go.rds"))
  
  #REACTOME
  if(species.scientific.name %in% unique(reactome.ncbi$V6) || species.scientific.name %in% unique(reactome.ensembl$V6)){
    reactome.ncbi.data <- reactome.ncbi[reactome.ncbi$V6==species.scientific.name,]
    reactome.ensembl.data <- reactome.ensembl[reactome.ensembl$V6==species.scientific.name,]
    reactome.data <- list(ncbi=reactome.ncbi.data,ensembl=reactome.ensembl.data)
    saveRDS(reactome.data,paste0("Input/",species,"/reactome.rds"))
  } else if(species.scientific.name %in% unique(reactome.plant.ncbi$V6) || species.scientific.name %in% unique(reactome.plant.ensembl$V6)) {
    reactome.ncbi.data <- reactome.plant.ncbi[reactome.plant.ncbi$V6==species.scientific.name,]
    reactome.ensembl.data <- reactome.plant.ensembl[reactome.plant.ensembl$V6==species.scientific.name,]
    reactome.data <- list(ncbi=reactome.ncbi.data,ensembl=reactome.ensembl.data)
    saveRDS(reactome.data,paste0("Input/",species,"/reactome.rds"))
  }
  
  #WIKIPATHWAYS
  wikipathways.data <- download.wikipathways.data(global.urls$wikipathways,species.scientific.name)
  saveRDS(wikipathways.data,paste0("Input/",species,"/wikipathways.rds"))
  
}

