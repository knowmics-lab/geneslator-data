#Get main script path
if (!requireNamespace("this.path", quietly = TRUE)) {
  install.packages("this.path", repos = "https://cloud.r-project.org")
}
base.dir <- dirname(this.path::this.path())
script.dir <- file.path(base.dir, "Script")

#Load required scripts
source(file.path(script.dir, "check_packages.R"))
source(file.path(script.dir, "download_functions.R"))
source(file.path(script.dir, "process_functions.R"))
source(file.path(script.dir, "mapping_functions.R"))
source(file.path(script.dir, "merge_functions.R"))

#Load required packages
ensure.packages(
  cran.pkgs = c("data.table","jsonlite","curl","zen4R","rvest","R.utils","ontologyIndex",
                "rentrez","xml2","purrr","dplyr","tidyr","optparse","DBI","RSQLite"),
  bioc.pkgs = c("rtracklayer","rWikiPathways","AnnotationForge")
)

#Read input parameters
option.list <- list(
  make_option(c("--org"), type = "character", default = "all",
              help = "Organisms [default: %default]"),
  make_option(c("--folder"), type = "character", default = "Output",
              help = "Output folder [default: %default]"),
  make_option(c("--saveTxt"), action = "store_true", default = FALSE,
              help = "Also save output files as .txt")
)
opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)
org.set <- opt$org
org.set <- strsplit(org.set,",")[[1]]
output.folder <- file.path(base.dir, opt$folder)
dir.create(output.folder, showWarnings = FALSE, recursive = TRUE)
save.text <- opt$saveTxt

#Load list of URLs
list.urls <- fromJSON(file.path(base.dir, "sources.json"), simplifyVector = TRUE)
global.urls <- list.urls$global_sources
species.url <- list.urls$species_sources
if(org.set!="all"){
  species.url <- species.url[species.url$official_name %in% org.set,]
}

#Prepare annotation data
annotation.data.list <- list()

#Create a taxonomy table
taxonomy.table <- species.url[,c("species","official_name","taxid")]

##-----------DOWNLOAD COMMON ANNOTATION DATA-----------

#NCBI archive data
print("Download NCBI archive data...")
global.ncbi.archive.data <- download.tabular.data(global.urls$ncbi_archive)
global.ncbi.archive.data <- global.ncbi.archive.data[global.ncbi.archive.data$`#tax_id` %in% unlist(species.url$taxid),]

#NCBI data about discontinued ids
print("Get annotations for NCBI discontinued ids...")
if(nrow(global.ncbi.archive.data)>0){
  global.ncbi.discontinued.data <- query.ncbi.discontinued.data(global.ncbi.archive.data)
} else {
  global.ncbi.discontinued.data <- data.frame()
}

#Get list of Ensembl folders
print("Get list of Ensembl folders required to download annotation data...")
ensembl.species <- download.tabular.data(paste0(global.urls$ensembl,"current/species_EnsemblVertebrates.txt"))
ensembl.folders <- paste0("gff3/",ensembl.species$`#name`)
ensembl.species <- gsub("_"," ",gsub("_core_.*$","",ensembl.species$other_alignments))
ensembl.species <- paste0(toupper(substr(ensembl.species,1,1)),substr(ensembl.species,2,nchar(ensembl.species)))
names(ensembl.folders) <- ensembl.species

#Get list of EnsemblGenomes folders
print("Get list of EnsemblGenomes folders required to download annotation data...")
ensembl.genome.species <- download.tabular.data(paste0(global.urls$ensemblGenomes,"current/species.txt"))
ens.divisions <- tolower(gsub("Ensembl","",ensembl.genome.species$species))
ens.collections <- gsub("_core_.*$","",ensembl.genome.species$other_alignments)
ensembl.genome.folders <- ifelse(grepl("_collection$",ens.collections),paste0(ens.divisions,"/gff3/",ens.collections,"/",ensembl.genome.species$`#name`),
                                 paste0(ens.divisions,"/gff3/",ensembl.genome.species$`#name`))
ensembl.genome.species <- gsub("_"," ",gsub("_gca_?[0-9]+.*$","",ensembl.genome.species$`#name`))
ensembl.genome.species <- paste0(toupper(substr(ensembl.genome.species,1,1)),substr(ensembl.genome.species,2,nchar(ensembl.genome.species)))
names(ensembl.genome.folders) <- ensembl.genome.species

#Read UNIPROT proteome species data
print("Get list of UNIPROT proteome identifiers required to download annotation data...")
uniprot.species <- download.tabular.data(paste0(global.urls$uniprot,"README"),skip="Proteome_ID\t",
  alternative.url=paste0(global.urls$uniprot_alternative,"README"))
uniprot.species <- uniprot.species[uniprot.species$Tax_ID %in% unlist(species.url$taxid),]
uniprot.species <- uniprot.species[!duplicated(uniprot.species$Tax_ID),]

#NCBI orthologs
print("Download NCBI orthologs data...")
global.ncbi.orthologs <- download.tabular.data(global.urls$ncbi_orthologs)
global.ncbi.orthologs <- global.ncbi.orthologs[global.ncbi.orthologs$`#tax_id` %in% unlist(species.url$taxid) & 
      global.ncbi.orthologs$`Other_tax_id` %in% unlist(species.url$taxid),]
#Alliance orthologs
print("Download Alliance of Genome Resources data...")
global.alliance.orthologs <- download.tabular.data(global.urls$alliance_genome_orthologs)
global.alliance.orthologs <- global.alliance.orthologs[global.alliance.orthologs$Gene1SpeciesTaxonID %in% paste0("NCBITaxon:",unlist(species.url$taxid)) & 
      global.alliance.orthologs$Gene2SpeciesTaxonID %in% paste0("NCBITaxon:",unlist(species.url$taxid)),]

#GO data
print("Download GO dictionary...")
go.dictionary <- download.go.dictionary(global.urls$go_dictionary)
go.dictionary <- process.go.dictionary(go.dictionary)

#Reactome data
print("Download Reactome data...")
reactome.ncbi <- download.tabular.data(paste0(global.urls$reactome,"NCBI2Reactome_All_Levels.txt"),header=F)
reactome.ensembl <- download.tabular.data(paste0(global.urls$reactome,"Ensembl2Reactome_All_Levels.txt"),header=F)
print("Download Reactome Plant data...")
reactome.plant.ncbi <- download.tabular.data(paste0(global.urls$reactomePlant,"NCBI2PlantReactome_All_Levels.txt"),header=F)
reactome.plant.ensembl <- download.tabular.data(paste0(global.urls$reactomePlant,"Ensembl2PlantReactome_All_Levels.txt"),header=F)

#HCOP data
print("Download HCOP data...")
hcop.data <- download.tabular.data(list.urls$species_sources[list.urls$species_sources$species=="Human","speciesdb_orthologs"])


##-----------BUILD GENERAL INFO ANNOTATION DATA-----------

list.species <- species.url$species
for(species in list.species){
  
  print(paste0("------Build general info annotation data for ",species,"---------"))
  
  #Get species-specific URLs of resources to download
  list.species.urls <- species.url[species.url$species==species,]
  species.scientific.name <- species.url[species.url$species==species,"official_name"]
  species.taxid <- species.url[species.url$species==species,"taxid"][[1]]
  speciesdb.name <- species.url[species.url$species==species,"speciesdb_name"]
  
  #----------NCBI data--------------
  #Process current data
  print("Download NCBI current data...")
  ncbi.data <- download.tabular.data(list.species.urls$ncbi_current)
  print("Process NCBI current data...")
  ncbi.data <- ncbi.data[ncbi.data$`#tax_id` %in% species.taxid,]
  ncbi.data <- process.ncbi.data(ncbi.data,speciesdb.name)
  #Process archive data
  print("Process archive data...")
  ncbi.archive.data <- global.ncbi.archive.data[global.ncbi.archive.data$`#tax_id` %in% species.taxid,]
  ncbi.archive.data <- process.ncbi.archive.data(ncbi.archive.data,global.ncbi.discontinued.data,speciesdb.name)
  ncbi.replaced.data <- ncbi.archive.data[[1]]
  ncbi.discontinued.data <- ncbi.archive.data[[2]]
  #Merge current and archive data
  print("Merge current and archive NCBI data...")
  ncbi.data <- merge.ncbi.data(ncbi.data,ncbi.replaced.data,ncbi.discontinued.data,speciesdb.name)
  
  #----------SPECIES DB data-------------
  if(!is.na(list.species.urls$speciesdb)){
    print("Download species-specific DB data...")
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
      speciesdb.file <- filter.remote.links.html(list.species.urls$speciesdb,"fbgn_annotation_ID")
      speciesdb.data <- download.delim.data(paste0(list.species.urls$speciesdb,speciesdb.file),
                                            header=F,comment.character="#")
    } else if(species=="Worm"){
      speciesdb.data <- download.json.data(list.species.urls$speciesdb,"")
    } else if(species=="Arabidopsis"){
      speciesdb.file <- filter.remote.links.json(list.species.urls$speciesdb,"^Araport11_GFF3_genes_transposons.*\\.gff\\.gz$")
      list.tags <- c("symbol","locus_type","computational_description","ID")
      list.filters <- "gene"
      download.path <- gsub("list\\?dir","download?filePath",list.species.urls$speciesdb)
      speciesdb.data <- download.gff.data(paste0(download.path,speciesdb.file),list.tags,list.filters)
    } else if(species=="AfricanClawedFrog"){
      speciesdb.data <- download.delim.data(list.species.urls$speciesdb,header=F,comment.character="!")
      speciesdb.data <- speciesdb.data[speciesdb.data$V7==paste0("taxon:",species.taxid),]
    } else if(species=="Macaque"){
      speciesdb.data <- download.tabular.data(list.species.urls$speciesdb)
      speciesdb.data <- speciesdb.data[speciesdb.data$taxon_id %in% species.taxid,]
    }
    print("Process species-specific DB data...")
    speciesdb.data <- process.speciesdb.data(speciesdb.data,species)
  } else {
    speciesdb.data <- data.frame()
  }
  
  #-------------ENSEMBL data------------------
  #Process current data
  print("Download Ensembl current data...")
  if(species.scientific.name %in% names(ensembl.folders)){
    ensembl.url <- global.urls$ensembl
    ensembl.folder <- ensembl.folders[species.scientific.name]
  } else {
    ensembl.url <- global.urls$ensemblGenomes
    ensembl.folder <- ensembl.genome.folders[species.scientific.name]
  }
  if(!is.na(ensembl.folder)){
    list.tags.ens <- c("Name","biotype","gene_id","description")
    list.filters.ens <- c("C_gene_segment","gene","J_gene_segment","lincRNA_gene","miRNA_gene",
                        "mt_gene","processed_transcript","pseudogene","RNA","rRNA_gene","snoRNA_gene","snRNA_gene",
                        "V_gene_segment","VD_gene_segment","ncRNA_gene")
    ensembl.path <- paste0(ensembl.url,"current/",ensembl.folder)
    link.content <- curl_fetch_memory(paste0(ensembl.path,"/?C=S;O=D"))
    list.links <- read_html(link.content$content) %>% html_elements("a") %>% html_text(trim=T)
    ensembl.file <- list.links[grep("gff3",list.links)[1]]
    ensembl.data <- download.gff.data(paste0(ensembl.path,"/",ensembl.file),list.tags.ens,list.filters.ens)
    print("Process Ensembl current data...")
    ensembl.data <- process.ensembl.data(ensembl.data,speciesdb.name,ncbi.data,
      if(species == "Zebrafish") hcop.data else NULL,species.taxid,is.archive=F)
    #Process archive data
    ensembl.data[["ENSEMBLOLD ENSEMBL"]] <- NA
    ensembl.archive.data <- list()
    print("Download Ensembl archive data...")
    ensembl.links <- filter.remote.links.html(ensembl.url,"release-")
    for(archive.link in ensembl.links){
      ensembl.version <- strsplit(archive.link,"-|/")[[1]][2]
      release.folder <- paste0(ensembl.url,archive.link,ensembl.folder)
      ensembl.archive <- download.ensembl.archive.data(release.folder, ensembl.version,
        list.tags.ens, list.filters.ens, species.scientific.name)
      if(!is.null(ensembl.archive)){
        ensembl.archive.data[[ensembl.version]] <- ensembl.archive
      }
    }
    print("Merge current and archive Ensembl data...")
    ensembl.data <- merge.with.ensembl.archive.data(ensembl.data,ensembl.archive.data,
      species,speciesdb.name,ncbi.data,if(species == "Zebrafish") hcop.data else NULL,species.taxid)
    #Process GRCH37 data
    if(!is.na(list.species.urls$ensembl_grch37)){
      print("Download GRCh37 Ensembl data...")
      list.tags <- c("gene_id","gene_biotype","gene_name")
      list.filters <- "gene"
      ensembl.grch37.data <- download.gff.data(list.species.urls$ensembl_grch37,list.tags,list.filters)
      ensembl.grch37.data <- process.ensembl.grch37.data(ensembl.grch37.data,speciesdb.name)
      print("Merge GRCh37 data with Ensembl current and archive data...")
      ensembl.data <- merge.with.ensembl.grch37.data(ensembl.data,ensembl.grch37.data,ncbi.data,speciesdb.name)
    }
  } else {
    ensembl.data <- data.frame()
  }
  
  #--------------UNIPROT data------------------
  print("Download Uniprot data...")
  uniprot.class <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"SUPERREGNUM"]
  uniprot.class <- paste0(toupper(substring(uniprot.class,1,1)),substring(uniprot.class,2,nchar(uniprot.class)))
  uniprot.code <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"Proteome_ID"]
  uniprot.tax <- uniprot.species[uniprot.species$Tax_ID %in% species.taxid,"Tax_ID"]
  uniprot.url <- paste0(global.urls$uniprot,uniprot.class,"/",uniprot.code,"/",uniprot.code,"_",uniprot.tax,".idmapping.gz")
  uniprot.url.alternative <- paste0(global.urls$uniprot_alternative,uniprot.class,"/",uniprot.code,"/",uniprot.code,"_",uniprot.tax,".idmapping.gz") 
  uniprot.data <- download.tabular.data(uniprot.url,header=F,alternative.url=uniprot.url.alternative)
  print("Process Uniprot data...")
  uniprot.data <- process.uniprot.data(uniprot.data,species)
  
  #----------Merge NCBI, Ensembl, species-db and Uniprot data-----------
  print("Merge NCBI, Ensembl, species-specific db and Uniprot data...")
  annotation.data <- merge.databases(ncbi.data,ensembl.data,uniprot.data,speciesdb.data,speciesdb.name,species)
  annotation.data.list[[species]] <- annotation.data
  
}

for(species in list.species){
  
  print(paste0("------Integrate functional annotation data for ",species,"---------"))
  
  #Get species-specific URLs of resources to download
  list.species.urls <- species.url[species.url$species==species,]
  species.scientific.name <- species.url[species.url$species==species,"official_name"]
  species.taxid <- species.url[species.url$species==species,"taxid"][[1]]
  speciesdb.name <- species.url[species.url$species==species,"speciesdb_name"]
  
  #-----------NCBI orthologs data----------
  print("Process NCBI orthologs data...")
  if(nrow(global.ncbi.orthologs)>0){
    ncbi.orthologs <- process.ncbi.orthologs.data(global.ncbi.orthologs,species.taxid,
                      taxonomy.table,annotation.data.list)
  } else {
    ncbi.orthologs <- data.frame()
  }
  
  #-----------ENSEMBL orthologs data------------
  if(species.scientific.name %in% names(ensembl.folders)){
    ensembl.url <- global.urls$ensembl
    ensembl.folder <- ensembl.folders[species.scientific.name]
  } else {
    ensembl.url <- global.urls$ensemblGenomes
    ensembl.folder <- ensembl.genome.folders[species.scientific.name]
  }
  if(!is.na(ensembl.folder)){
    print("Download Ensembl orthologs data...")
    ensembl.orthologs.folder <- gsub("gff3/","json/",ensembl.folder)
    ensembl.orthologs.species <- strsplit(ensembl.orthologs.folder,"/")[[1]]
    ensembl.orthologs.species <- ensembl.orthologs.species[length(ensembl.orthologs.species)]
    list.species.taxid <- unlist(species.url$taxid)
    ortho.species.taxid <- list.species.taxid[!list.species.taxid %in% species.taxid]
    if(length(ortho.species.taxid)>0){
      filter.string.ortho.json <- paste0('[.genes[] | .id as $gene_id | .homologues[] | select(.taxonomy_id | IN(',
        paste0(ortho.species.taxid,collapse = ","),
        ')) | {ID: $gene_id, TAXID: .taxonomy_id, ENSEMBL: .stable_id}]')
      ensembl.orthologs <- download.json.data(paste0(ensembl.url,"current/",
        ensembl.orthologs.folder,"/",ensembl.orthologs.species,".json"),filter.string.ortho.json)
      print("Process Ensembl orthologs data...")
      if(class(ensembl.orthologs)=="list"){
        ensembl.orthologs <- data.frame(matrix(NA, nrow = 0, ncol = 3))
      }
      ensembl.orthologs <- process.ensembl.orthologs.data(ensembl.orthologs,
        species.taxid,taxonomy.table,annotation.data.list)
    } else {
      ensembl.orthologs <- data.frame()
    }
  } else {
    ensembl.orthologs <- data.frame()
  }
  
  #----------ALLIANCE orthologs data----------
  print("Process Alliance of Genome orthologs data...")
  if(nrow(global.alliance.orthologs)>0){
    alliance.orthologs <- process.alliance.orthologs.data(global.alliance.orthologs,
      species.taxid,taxonomy.table)
  } else {
    alliance.orthologs <- data.frame()
  }
  
  #----------SPECIESDB orthologs data-------------
  if(!is.na(list.species.urls$speciesdb_orthologs)){
    print("Download species-specific DB orthologs data...")
    speciesdb.orthologs <- download.tabular.data(list.species.urls$speciesdb_orthologs)
    speciesdb.orthologs <- speciesdb.orthologs[speciesdb.orthologs$ortholog_species %in% species.url$taxid,]
    print("Process species-specific DB orthologs data...")
    speciesdb.orthologs <- process.speciesdb.orthologs.data(speciesdb.orthologs,
      species.taxid,taxonomy.table,annotation.data.list)
  } else {
    speciesdb.orthologs <- data.frame()
  }
  
  #-----------Merge orthologs data--------------
  print("Merge orthologs data...")
  orthologs.data <- merge.ortho.databases(ncbi.orthologs,ensembl.orthologs,alliance.orthologs,
    speciesdb.orthologs,species.taxid,taxonomy.table)
  
  #-----------GO data-----------------
  print("Download GO data...")
  go.data <- download.delim.data(list.species.urls$go,header=F,comment.character="!")
  print("Process GO data...")
  go.data <- process.go.data(go.data,go.dictionary,annotation.data.list,species)
  
  #-----------REACTOME data----------------
  print("Download Reactome data...")
  if(species.scientific.name %in% unique(reactome.ncbi$V6) || species.scientific.name %in% unique(reactome.ensembl$V6)){
    reactome.ncbi.data <- reactome.ncbi[reactome.ncbi$V6==species.scientific.name,]
    reactome.ensembl.data <- reactome.ensembl[reactome.ensembl$V6==species.scientific.name,]
    reactome.data <- list(ncbi=reactome.ncbi.data,ensembl=reactome.ensembl.data)
  } else if(species.scientific.name %in% unique(reactome.plant.ncbi$V6) || species.scientific.name %in% unique(reactome.plant.ensembl$V6)) {
    reactome.ncbi.data <- reactome.plant.ncbi[reactome.plant.ncbi$V6==species.scientific.name,]
    reactome.ensembl.data <- reactome.plant.ensembl[reactome.plant.ensembl$V6==species.scientific.name,]
    reactome.data <- list(ncbi=reactome.ncbi.data,ensembl=reactome.ensembl.data)
  } else {
    reactome.data <- list()
  }
  print("Process Reactome data...")
  reactome.data <- process.reactome.data(reactome.data,annotation.data.list,species)
  
  #-----------WIKIPATHWAYS data--------------
  print("Download Wikipathways data...")
  wikipathways.data <- download.wikipathways.data(global.urls$wikipathways,species.scientific.name)
  print("Process Wikipathways data...")
  wikipathways.data <- process.wikipathways.data(wikipathways.data,annotation.data.list,species)
  
  #---------Integrate orthologs, GO and pathway data------------
  print("Merge orthologs, GO and pathway data with current annotation data...")
  annotation.data <- annotation.data.list[[species]]
  if(nrow(orthologs.data)>0){
    annotation.data <- merge(annotation.data,orthologs.data,all.x=T)
  }
  if(nrow(go.data)>0){
    annotation.data <- merge(annotation.data,go.data,all.x=T)
  }
  if(nrow(reactome.data)>0){
    annotation.data <- merge(annotation.data,reactome.data,all.x=T)
  }
  if(nrow(wikipathways.data)>0){
    annotation.data <- merge(annotation.data,wikipathways.data,all.x=T)
  }
  annotation.data <- as.data.frame(annotation.data)
  
  #----------Write annotation DB to output file------------
  print(paste0("Save annotation table for ",species," as SQLite file..."))
  #Add GID to annotation table
  annotation.table <- annotation.data
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
  if(nrow(ncbi.table)>0){
    build.package.args$ncbi <- ncbi.table
  }
  #NCBI OLD
  ncbi.old.table <- annotation.table[,c("GID","ENTREZIDOLD")]
  ncbi.old.table <- ncbi.old.table[!is.na(ncbi.old.table$ENTREZIDOLD),]
  ncbi.old.table <- as.data.frame(ncbi.old.table %>% separate_rows(all_of(c("ENTREZIDOLD")),sep="\\|"))
  if(nrow(ncbi.old.table)>0){
    build.package.args$ncbiOld <- ncbi.old.table
  }
  #ENSEMBL
  ensembl.table <- annotation.table[,c("GID","ENSEMBL")]
  ensembl.table <- ensembl.table[!is.na(ensembl.table$ENSEMBL),]
  ensembl.table <- as.data.frame(ensembl.table %>% separate_rows(all_of(c("ENSEMBL")),sep="\\|"))
  if(nrow(ensembl.table)>0){
    build.package.args$ensembl <- ensembl.table
  }
  #ENSEMBL OLD
  ensembl.old.table <- annotation.table[,c("GID","ENSEMBLOLD")]
  ensembl.old.table <- ensembl.old.table[!is.na(ensembl.old.table$ENSEMBLOLD),]
  ensembl.old.table <- as.data.frame(ensembl.old.table %>% separate_rows(all_of(c("ENSEMBLOLD")),sep="\\|"))
  if(nrow(ensembl.old.table)>0){
    build.package.args$ensemblOld <- ensembl.old.table
  }
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
  for(other.species in species.url$species){
    if(species!=other.species){
      ortho.data <- annotation.table[,c("GID",paste0("ORTHO",toupper(other.species)))]
      ortho.data <- ortho.data[!is.na(ortho.data[[paste0("ORTHO",toupper(other.species))]]),]
      ortho.data <- as.data.frame(ortho.data %>% separate_rows(all_of(c(paste0("ORTHO",toupper(other.species)))),sep="\\|"))
      if(nrow(ortho.data)>0){
        build.package.args[[paste0("ortho",other.species)]] <- ortho.data
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
  if("REACTOMEPATH" %in% colnames(annotation.table)){
    reactome.table <- annotation.table[,c("GID","REACTOMEPATH","REACTOMEPATHNAME")]
    reactome.table <- reactome.table[!is.na(reactome.table$REACTOMEPATH),]
    reactome.table <- as.data.frame(reactome.table %>% separate_rows(all_of(c("REACTOMEPATH","REACTOMEPATHNAME")),sep="\\|"))
    reactome.table <- unique(reactome.table)
    if(nrow(reactome.table)>0){
      build.package.args$reactome <- reactome.table
    }
  }
  #WIKIPATHWAYS PATHWAY
  if("WIKIPATH" %in% colnames(annotation.table)){
    wikipath.table <- annotation.table[,c("GID","WIKIPATH","WIKIPATHNAME")]
    wikipath.table <- wikipath.table[!is.na(wikipath.table$WIKIPATH),]
    wikipath.table <- as.data.frame(wikipath.table %>% separate_rows(all_of(c("WIKIPATH","WIKIPATHNAME")),sep="\\|"))
    if(nrow(wikipath.table)>0){
      build.package.args$wiki <- wikipath.table
    }
  }
  #Set additional parameters for makeOrgPackage function
  build.package.args$version <- "1.0"
  build.package.args$maintainer <- "Giovanni Micale <giovanni.micale@unict.it>"
  build.package.args$author <- "Giovanni Micale <giovanni.micale@unict.it>"
  build.package.args$outputDir <- output.folder
  build.package.args$tax_id <- species.taxid[length(species.taxid)]
  species.tax <- strsplit(species.scientific.name," ")[[1]]
  build.package.args$genus <- species.tax[1]
  build.package.args$species <- species.tax[2]
  #Create annotation database
  do.call(makeOrgPackage,build.package.args)
  #Keep a copy of the sqlite database before removing the package folder
  db.name <- paste0("org.", substring(species.tax[1],1,1), species.tax[2], ".eg")
  pkg.folder <- file.path(output.folder, paste0(db.name, ".db"))
  built.sqlite <- file.path(pkg.folder, "inst", "extdata", paste0(db.name, ".sqlite"))
  final.sqlite <- file.path(output.folder, paste0(gsub(".eg$",".db",db.name), ".sqlite"))
  if(file.exists(built.sqlite)){
    file.copy(built.sqlite, final.sqlite, overwrite = TRUE)
  }
  #Delete the built package folder (only the .sqlite is wanted)
  unlink(pkg.folder, recursive = TRUE, force = TRUE)
  if(file.exists(gsub(".db$",".sqlite",pkg.folder))){
    file.remove(gsub(".db$",".sqlite",pkg.folder))
  }
  #Add KEGG code to db metadata
  Sys.chmod(final.sqlite, mode = "0644", use_umask = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), final.sqlite)
  DBI::dbExecute(con, paste0("INSERT INTO metadata (name, value) VALUES ('KEGGCODE', '",list.species.urls$kegg_code,"')"))
  DBI::dbDisconnect(con)
  #Save annotations to text file (if required)
  if(save.text) {
    print(paste0("Save annotation table for ",species," as text file..."))
    write.table(annotation.data,file.path(output.folder,paste0(gsub(".eg$",".db",db.name),".txt")),
                sep="\t",quote=F,row.names = F)
  }
  
}
print("DONE!")
