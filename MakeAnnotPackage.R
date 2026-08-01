library(jsonlite)
library(data.table)
library(tidyr)
library(AnnotationForge)

#Read metadata
metadata <- fromJSON("Input/sources.json", simplifyVector = TRUE)
species.metadata <- metadata$species_sources

#Create annotation DBs output folder, if it does not exist
if(!dir.exists("Output")){
  mkdirs("Output")
}

#Create annotation packages
for(org in species.metadata$species){
  
  print(paste0("Creating annotation package for ",org,"..."))
  speciesdb.name <- species.metadata[species.metadata$species==org,"speciesdb_name"]
  
  #Read annotation table
  annotation.table <- fread(paste0("Annotations/",org,".txt"),data.table = F)
  annotation.table[is.na(annotation.table$SYMBOL),"SYMBOL"] <- "NA"
  annotation.table$GID <- 1:nrow(annotation.table)
  
  #Initialize list of args for makeOrgPackage() function
  build.package.args <- list()
  
  #Gene symbol
  symbol.table <- annotation.table[,c("GID","SYMBOL")]
  symbol.table <- symbol.table[!is.na(symbol.table$SYMBOL),]
  symbol.table <- symbol.table[order(symbol.table$GID),]
  if(nrow(symbol.table)>0){
    build.package.args$symbol <- symbol.table
  }
  
  #ALIAS
  alias.table <- annotation.table[,c("GID","SYMBOL","ALIAS")]
  alias.table <- as.data.frame(alias.table %>% separate_rows(all_of(c("ALIAS")),sep="\\|"))
  alias.table <- alias.table[alias.table$SYMBOL!=alias.table$ALIAS,c("GID","ALIAS")]
  if(nrow(alias.table)>0){
    build.package.args$alias <- alias.table
  }

  #Gene locus
  locus.table <- annotation.table[,c("GID","LOCUS")]
  locus.table <- locus.table[!is.na(locus.table$LOCUS),]
  locus.table <- as.data.frame(locus.table %>% separate_rows(all_of(c("LOCUS")),sep="\\|"))
  if(nrow(locus.table)>0){
    build.package.args$locus <- locus.table
  }
  
  #Gene type
  gene.type.table <- annotation.table[,c("GID","GENETYPE")]
  gene.type.table <- gene.type.table[!is.na(gene.type.table$GENETYPE),]
  gene.type.table <- as.data.frame(gene.type.table %>% separate_rows(all_of(c("GENETYPE")),sep="\\|"))
  if(nrow(gene.type.table)>0){
    build.package.args$genetype <- gene.type.table
  }

  #Gene name
  gene.name.table <- annotation.table[,c("GID","GENENAME")]
  gene.name.table <- gene.name.table[!is.na(gene.name.table$GENENAME),]
  gene.name.table <- as.data.frame(gene.name.table %>% separate_rows(all_of(c("GENENAME")),sep="\\|"))
  if(nrow(gene.name.table)>0){
    build.package.args$genename <- gene.name.table
  }
  
  #NCBI
  ncbi.table <- annotation.table[,c("GID","ENTREZID")]
  ncbi.table <- ncbi.table[!is.na(ncbi.table$ENTREZID),]
  ncbi.table <- as.data.frame(ncbi.table %>% separate_rows(all_of(c("ENTREZID")),sep="\\|"))
  build.package.args$ncbi <- ncbi.table

  #NCBI OLD
  ncbi.old.table <- annotation.table[,c("GID","ENTREZIDOLD")]
  ncbi.old.table <- ncbi.old.table[!is.na(ncbi.old.table$ENTREZIDOLD),]
  ncbi.old.table <- as.data.frame(ncbi.old.table %>% separate_rows(all_of(c("ENTREZIDOLD")),sep="\\|"))
  build.package.args$ncbiOld <- ncbi.old.table
  
  #ENSEMBL
  ensembl.table <- annotation.table[,c("GID","ENSEMBL")]
  ensembl.table <- ensembl.table[!is.na(ensembl.table$ENSEMBL),]
  ensembl.table <- as.data.frame(ensembl.table %>% separate_rows(all_of(c("ENSEMBL")),sep="\\|"))
  build.package.args$ensembl <- ensembl.table
  
  #ENSEMBL OLD
  ensembl.old.table <- annotation.table[,c("GID","ENSEMBLOLD")]
  ensembl.old.table <- ensembl.old.table[!is.na(ensembl.old.table$ENSEMBLOLD),]
  ensembl.old.table <- as.data.frame(ensembl.old.table %>% separate_rows(all_of(c("ENSEMBLOLD")),sep="\\|"))
  build.package.args$ensemblOld <- ensembl.old.table
  
  #SPECIAL DATA
  if(!is.na(speciesdb.name)){
    speciesdb.table <- annotation.table[,c("GID",speciesdb.name)]
    speciesdb.table <- speciesdb.table[!is.na(speciesdb.table[[speciesdb.name]]),]
    speciesdb.table <- as.data.frame(speciesdb.table %>% separate_rows(all_of(c(speciesdb.name)),sep="\\|"))
    if(nrow(speciesdb.table)>0){
      build.package.args$speciesId <- speciesdb.table
    }
  }
  
  #UNIPROT
  uniprot.table <- annotation.table[,c("GID","UNIPROT")]
  uniprot.table <- uniprot.table[!is.na(uniprot.table$UNIPROT),]
  uniprot.table <- as.data.frame(uniprot.table %>% separate_rows(all_of(c("UNIPROT")),sep="\\|"))
  if(nrow(uniprot.table)>0){
    build.package.args$uniprot <- uniprot.table
  }
  
  #ORTHOLOGS
  for(other.org in species.metadata$species){
    if(org!=other.org){
      ortho.data <- annotation.table[,c("GID",paste0("ORTHO",toupper(other.org)))]
      ortho.data <- ortho.data[!is.na(ortho.data[[paste0("ORTHO",toupper(other.org))]]),]
      ortho.data <- as.data.frame(ortho.data %>% separate_rows(all_of(c(paste0("ORTHO",toupper(other.org)))),sep="\\|"))
      if(nrow(ortho.data)>0){
        build.package.args[[paste0("ortho",other.org)]] <- ortho.data
      }
    }
  }
  
  #GO
  go.table <- annotation.table[,c("GID","GO","GOEVIDENCE","GONAME","GOTYPE")]
  go.table <- go.table[!is.na(go.table$GO),]
  go.table <- as.data.frame(go.table %>% separate_rows(all_of(c("GO","GOEVIDENCE","GONAME","GOTYPE")),sep="\\|"))
  if(nrow(go.table)>0){
    build.package.args$go <- go.table
  }
  
  #REACTOME PATHWAY
  reactome.table <- annotation.table[,c("GID","REACTOMEPATH","REACTOMEPATHNAME")]
  reactome.table <- reactome.table[!is.na(reactome.table$REACTOMEPATH),]
  reactome.table <- as.data.frame(reactome.table %>% separate_rows(all_of(c("REACTOMEPATH","REACTOMEPATHNAME")),sep="\\|"))
  reactome.table <- unique(reactome.table)
  if(nrow(reactome.table)>0){
    build.package.args$reactome <- reactome.table
  }
  
  #WIKIPATHWAYS PATHWAY
  wikipath.table <- annotation.table[,c("GID","WIKIPATH","WIKIPATHNAME")]
  wikipath.table <- wikipath.table[!is.na(wikipath.table$WIKIPATH),]
  wikipath.table <- as.data.frame(wikipath.table %>% separate_rows(all_of(c("WIKIPATH","WIKIPATHNAME")),sep="\\|"))
  if(nrow(wikipath.table)>0){
    build.package.args$wiki <- wikipath.table
  }
  
  #Set additional parameters for makeOrgPackage function
  build.package.args$version <- "1.0"
  build.package.args$maintainer <- "Giovanni Micale <giovanni.micale@unict.it>"
  build.package.args$author <- "Giovanni Micale <giovanni.micale@unict.it>"
  build.package.args$outputDir <- "Output"
  species.taxid <- species.metadata[species.metadata$species==org,"taxid"][[1]]
  build.package.args$tax_id <- species.taxid[length(species.taxid)]
  species.tax <- strsplit(species.metadata[species.metadata$species==org,"official_name"]," ")[[1]]
  build.package.args$genus <- species.tax[1]
  build.package.args$species <- species.tax[2]
  
  #Create annotation database
  do.call(makeOrgPackage,build.package.args)
  
  db.name      <- paste0("org.", substring(species.tax[1],1,1), species.tax[2], ".eg")
  pkg.folder   <- file.path("Output", paste0(db.name, ".db"))
  built.sqlite <- file.path(pkg.folder, "inst", "extdata", paste0(db.name, ".sqlite"))
  final.sqlite <- file.path("Output", paste0(db.name, ".sqlite"))
  
  #Keep a copy of the sqlite database before removing the package folder
  if(file.exists(built.sqlite)){
    file.copy(built.sqlite, final.sqlite, overwrite = TRUE)
  }
  
  #Delete the built package folder (only the .sqlite is wanted)
  unlink(pkg.folder, recursive = TRUE, force = TRUE)
  
}
