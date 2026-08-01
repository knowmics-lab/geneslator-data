process.ncbi.data <- function(ncbi.data,taxid,speciesdb.name)
{
  ncbi.data <- ncbi.data[ncbi.data$`#tax_id` %in% taxid,c("GeneID","Symbol","LocusTag","Symbol_from_nomenclature_authority","type_of_gene","description","Full_name_from_nomenclature_authority","Synonyms","dbXrefs")]
  colnames(ncbi.data) <- c("ENTREZID","UNOFFICIAL SYMBOL","LOCUS","SYMBOL","GENETYPE","UNOFFICIAL GENENAME","GENENAME","ALIAS","ENSEMBL")
  ncbi.data[ncbi.data$SYMBOL=="NA" | ncbi.data$SYMBOL=="-","SYMBOL"] <- NA
  ncbi.data[ncbi.data$LOCUS=="NA" | ncbi.data$LOCUS=="-","LOCUS"] <- NA
  ncbi.data[ncbi.data$GENETYPE=="NA" | ncbi.data$GENETYPE=="-","GENETYPE"] <- NA
  ncbi.data[ncbi.data$GENENAME=="NA" | ncbi.data$GENENAME=="-","GENENAME"] <- NA
  ncbi.data[ncbi.data$ALIAS=="NA" | ncbi.data$ALIAS=="-","ALIAS"] <- NA
  ncbi.data[ncbi.data$ENTREZID=="NA" | ncbi.data$ENTREZID=="-","ENTREZID"] <- NA
  ncbi.data[ncbi.data$ENSEMBL=="NA" | ncbi.data$ENSEMBL=="-","ENSEMBL"] <- NA
  ncbi.data$SYMBOL <- apply(ncbi.data,1,function(row){
    ifelse(is.na(row["SYMBOL"]),row["UNOFFICIAL SYMBOL"],row["SYMBOL"])
  })
  ncbi.data$GENENAME <- apply(ncbi.data,1,function(row){
    ifelse(is.na(row["GENENAME"]),row["UNOFFICIAL GENENAME"],row["GENENAME"])
  })
  ncbi.data$`UNOFFICIAL SYMBOL` <- NULL
  ncbi.data$`UNOFFICIAL GENENAME` <- NULL
  if(!is.na(speciesdb.name)){
    ncbi.data[[speciesdb.name]] <- sapply(ncbi.data$ENSEMBL,function(x){
      if(is.na(x)) {
        return(NA)
      } else {
        list.ids <- unlist(strsplit(x,"\\|"))
        res <- paste0(sub(paste0(speciesdb.name,":"),"",list.ids[grep(paste0("^",speciesdb.name,":"),list.ids,ignore.case=T)],
                    ignore.case=T), collapse="|")
        return(ifelse(res=="",NA,res))
      }
    })
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="SGD"){
    ensembl.col <- "FungiDB"
  } else {
    ensembl.col <- "Ensembl"
  }
  if(!is.na(speciesdb.name) && speciesdb.name %in% c("FLYBASE","WORMBASE","TAIR")){
    ncbi.data$ENSEMBL <- ncbi.data[[speciesdb.name]]
  } else {
    ncbi.data$ENSEMBL <- sapply(ncbi.data$ENSEMBL,function(x){
      if(is.na(x)) {
        return(NA)
      } else {
        list.ids <- unlist(strsplit(x,"\\|"))
        res <- paste0(gsub(paste0(ensembl.col,":"),"",list.ids[grep(paste0("^",ensembl.col,":"),list.ids,ignore.case=T)],
                           ignore.case=T), collapse="|")
        return(ifelse(res=="",NA,res))
      }
    })
  }
  ncbi.data$ALIAS <- apply(ncbi.data,1,function(row){
    ifelse(is.na(row["ALIAS"]),row["SYMBOL"],paste0(row["ALIAS"],"|",row["SYMBOL"]))
  })
  ncbi.data$ENTREZID <- as.character(ncbi.data$ENTREZID)
  if(!is.na(speciesdb.name)){
    ncbi.data <- ncbi.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZID","ENSEMBL",speciesdb.name)]
  } else {
    ncbi.data <- ncbi.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZID","ENSEMBL")]
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="TAIR"){
    ncbi.data <- ncbi.data %>% group_by(SYMBOL,GENENAME) %>% summarise(across(everything(),function(x){
      unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
  } else {
    ncbi.data <- ncbi.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
      unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="TAIR"){
    colnames(ncbi.data)[!colnames(ncbi.data) %in% c("SYMBOL","GENENAME")] <- paste0(colnames(ncbi.data)[!colnames(ncbi.data) %in% c("SYMBOL","GENENAME")]," NCBI")
  } else {
    colnames(ncbi.data)[!colnames(ncbi.data) %in% c("SYMBOL")] <- paste0(colnames(ncbi.data)[!colnames(ncbi.data) %in% c("SYMBOL")]," NCBI")
  }
  return(as.data.frame(ncbi.data))
}

process.ncbi.archive.data <- function(global.ncbi.archive.data,global.ncbi.discontinued.data,taxid,speciesdb.name)
{
  #Filter archive data for the organism
  ncbi.archive.data <- global.ncbi.archive.data[global.ncbi.archive.data$`#tax_id` %in% taxid,
                          c("GeneID","Discontinued_GeneID","Discontinued_Symbol")]
  
  #Process data about replaced IDs
  ncbi.replaced.data <- ncbi.archive.data[ncbi.archive.data$GeneID!="-",]
  ncbi.replaced.data$GeneID <- as.integer(ncbi.replaced.data$GeneID)
  colnames(ncbi.replaced.data) <- c("ENTREZID","ENTREZIDOLD","SYMBOL")
  ncbi.replaced.data <- ncbi.replaced.data %>% group_by(ENTREZID) %>% summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
    res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
    ifelse(res=="",NA,res)
  }))
  ncbi.replaced.data$ENTREZID <- as.character(ncbi.replaced.data$ENTREZID)
  ncbi.replaced.data$ENTREZIDOLD <- as.character(ncbi.replaced.data$ENTREZIDOLD)
  colnames(ncbi.replaced.data)[colnames(ncbi.replaced.data)!="ENTREZID"] <- paste0(colnames(ncbi.replaced.data)[colnames(ncbi.replaced.data)!="ENTREZID"]," ARCHIVE")
  
  #Process data about discontinued IDs
  ncbi.discontinued.data <- global.ncbi.discontinued.data[global.ncbi.discontinued.data$ENTREZIDOLD %in% ncbi.archive.data$Discontinued_GeneID,]
  if(!is.na(speciesdb.name)){
    ncbi.discontinued.data[[speciesdb.name]] <- sapply(ncbi.discontinued.data$ExtIds,function(x){
      if(is.na(x)) {
        return(NA)
      } else {
        list.ids <- unlist(strsplit(x,"\\|"))
        res <- paste0(sub(paste0(speciesdb.name,":"),"",list.ids[grep(paste0("^",speciesdb.name,":"),list.ids,ignore.case=T)],
                      ignore.case=T), collapse="|")
        return(ifelse(res=="",NA,res))
      }
    })
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="SGD"){
    ensembl.col <- "FungiDB"
  } else {
    ensembl.col <- "Ensembl"
  }
  if(!is.na(speciesdb.name) && speciesdb.name %in% c("FLYBASE","WORMBASE","TAIR")){
    ncbi.discontinued.data$ENSEMBL <- ncbi.discontinued.data[[speciesdb.name]]
  } else {
    ncbi.discontinued.data$ENSEMBL <- sapply(ncbi.discontinued.data$ExtIds,function(x){
      if(is.na(x)) {
        return(NA)
      } else {
        list.ids <- unlist(strsplit(x,"\\|"))
        res <- paste0(gsub(paste0(ensembl.col,":"),"",list.ids[grep(paste0("^",ensembl.col,":"),list.ids,ignore.case=T)],
                          ignore.case=T), collapse="|")
        return(ifelse(res=="",NA,res))
      }
    })
  }
  if(nrow(ncbi.discontinued.data)>0){
    ncbi.discontinued.data$ALIAS <- apply(ncbi.discontinued.data,1,function(row){
      if(is.na(row[["ALIAS"]])){
        row[["SYMBOL"]]
      } else {
        list.ids <- unique(c(strsplit(row[["ALIAS"]],"\\|")[[1]],row[["SYMBOL"]]))
        paste0(list.ids,collapse = "|")
      }
    })
  }
  if(!is.na(speciesdb.name)){
    ncbi.discontinued.data <- ncbi.discontinued.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZIDOLD","ENSEMBL",speciesdb.name)]
  } else {
    ncbi.discontinued.data <- ncbi.discontinued.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZIDOLD","ENSEMBL")]
  }
  if(species=="Arabidopsis"){
    ncbi.discontinued.data <- ncbi.discontinued.data %>% group_by(SYMBOL,GENENAME) %>% summarise(across(everything(),function(x){
      unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
  } else {
    ncbi.discontinued.data <- ncbi.discontinued.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
      unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
  }
  ncbi.discontinued.data$ENTREZIDOLD <- as.character(ncbi.discontinued.data$ENTREZIDOLD)
  colnames(ncbi.discontinued.data)[colnames(ncbi.discontinued.data)!="SYMBOL"] <- paste0(colnames(ncbi.discontinued.data)[colnames(ncbi.discontinued.data)!="SYMBOL"]," ARCHIVE")
  
  list.res <- list(as.data.frame(ncbi.replaced.data),as.data.frame(ncbi.discontinued.data))
  return(list.res)
}

process.ensembl.data <- function(ensembl.data,speciesdb.name,ncbi.data,speciesdb.data,hcop.data,taxid,is.archive)
{
  #Pre-process data
  colnames(ensembl.data) <- c("SYMBOL","GENETYPE","ENSEMBL","GENENAME")
  ensembl.data$ALIAS <- ensembl.data$SYMBOL
  ensembl.data <- ensembl.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ENSEMBL")]
  
  #Fix symbols from HGNC in zebrafish
  if(!is.na(speciesdb.name) && speciesdb.name=="ZFIN"){
    ensembl.data[grep("HGNC",ensembl.data$GENENAME),"SYMBOL"] <- tolower(ensembl.data[grep("HGNC",ensembl.data$GENENAME),"SYMBOL"])
  }
  
  #Double check with NCBI entrez ids (only in zebrafish)
  if(!is.na(speciesdb.name) && speciesdb.name=="ZFIN"){
    #Get Entrez Ids from HCOP
    hcop.data <- unique(hcop.data[hcop.data$ortholog_species==taxid,c("ortholog_species_entrez_gene","ortholog_species_ensembl_gene")])
    colnames(hcop.data) <- c("ENTREZID HCOP","ENSEMBL")
    hcop.data <- hcop.data[hcop.data$`ENTREZID HCOP` !="-" & hcop.data$ENSEMBL!="-",]
    ensembl.data <- merge(ensembl.data,hcop.data,all.x=T)
    #Get Entrez Ids from NCBI
    ensembl.data[["ENTREZID ENSEMBL"]] <- unname(sapply(ensembl.data$GENENAME,function(x){
      if(grepl("Source:NCBI",x)){
        gsub("\\]","",strsplit(x,"Acc:")[[1]][2])
      } else {
        NA
      }
    }))
    #Get Entrez Ids from ZFIN
    ensembl.data[["ZFIN"]] <- unname(sapply(ensembl.data$GENENAME,function(x){
      if(grepl(paste0("Source:",speciesdb.name),x)){
        gsub("\\]","",strsplit(x,"Acc:")[[1]][2])
      } else {
        NA
      }
    }))
    ncbi.data.special <- ncbi.data[,c("ENTREZID NCBI",paste0(speciesdb.name," NCBI"))] %>% separate_rows(all_of(c(paste0(speciesdb.name," NCBI"))),sep="\\|")
    colnames(ncbi.data.special) <- c("ENTREZID NCBI",speciesdb.name)
    ncbi.data.special <- ncbi.data.special[!is.na(ncbi.data.special[[speciesdb.name]]),]
    ensembl.data <- merge(ensembl.data,ncbi.data.special,all.x=T)
    ensembl.data[[speciesdb.name]] <- NULL
    #Merge Entrez Ids columns from HCOP, NCBI and ZFIN
    ensembl.data$ENTREZID <- merge.columns(ensembl.data,"ENTREZID",single.val=T)
    ensembl.data$`ENTREZID HCOP` <- NULL
    ensembl.data$`ENTREZID ENSEMBL` <- NULL
    ensembl.data$`ENTREZID NCBI` <- NULL
    #Perform double check
    ensembl.data <- double.check.symbols(ensembl.data,ncbi.data,"ENTREZID","NCBI",is.archive)
    colnames(ensembl.data)[colnames(ensembl.data)=="ENTREZID"] <- "ENTREZIDOLD"
    ensembl.data <- double.check.symbols(ensembl.data,ncbi.data,"ENTREZIDOLD","NCBI",is.archive)
    ensembl.data$ENTREZIDOLD <- NULL
  }
  
  list.desc <- strsplit(as.character(ensembl.data$GENENAME)," \\[Source")
  ensembl.data$GENENAME <- unlist(lapply(list.desc,function(x)x[1]))
  list.ncbi.symbols <- ncbi.data$SYMBOL
  if(!is.na(speciesdb.name) && speciesdb.name=="SGD"){
    ensembl.data[is.na(ensembl.data$SYMBOL) & (ensembl.data$ENSEMBL %in% list.ncbi.symbols | startsWith(ensembl.data$ENSEMBL,"Y")),"SYMBOL"] <- ensembl.data[is.na(ensembl.data$SYMBOL) & (ensembl.data$ENSEMBL %in% list.ncbi.symbols | startsWith(ensembl.data$ENSEMBL,"Y")),"ENSEMBL"]
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="FLYBASE"){
    ensembl.data <- ensembl.data[ensembl.data$GENETYPE!="transposable_element",]
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="HGNC"){
    ensembl.data <- ensembl.data[!startsWith(ensembl.data$ENSEMBL,"LRG_"),]
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="TAIR"){
    ensembl.data[is.na(ensembl.data$SYMBOL) & ensembl.data$ENSEMBL %in% list.ncbi.symbols,"SYMBOL"] <- ensembl.data[is.na(ensembl.data$SYMBOL) & ensembl.data$ENSEMBL %in% list.ncbi.symbols,"ENSEMBL"]
    ensembl.data[!is.na(ensembl.data$SYMBOL) & startsWith(ensembl.data$SYMBOL,"ath-"),"SYMBOL"] <- gsub("ath-","",ensembl.data[!is.na(ensembl.data$SYMBOL) & startsWith(ensembl.data$SYMBOL,"ath-"),"SYMBOL"])
    ensembl.data[startsWith(ensembl.data$ENSEMBL,"ENSRNA"),"ENSEMBL"] <- NA
    ensembl.data <- ensembl.data[!is.na(ensembl.data$ENSEMBL),]
  }
  if(is.na(speciesdb.name)){
    ensembl.data$ENSEMBL <- gsub("gene-","",ensembl.data$ENSEMBL)
  }
  
  #Fill symbol and alias information, if needed
  ensembl.data <- as.data.frame(ensembl.data %>% separate_rows(all_of(c("SYMBOL")),sep="\\|"))
  ensembl.data <- ensembl.data[!is.na(ensembl.data$SYMBOL),]
  if(nrow(ensembl.data)>0){
    if(!is.na(speciesdb.name) && speciesdb.name=="ZFIN"){
      ensembl.data[grepl("^loc[0-9]+",ensembl.data$SYMBOL),"SYMBOL"] <- toupper(ensembl.data[grepl("^loc[0-9]+",ensembl.data$SYMBOL),"SYMBOL"])
    }
    ensembl.data$ALIAS <- apply(ensembl.data,1,function(row){
      if(is.na(row[["ALIAS"]])){
        row[["SYMBOL"]]
      } else {
        list.ids <- unique(c(strsplit(row[["ALIAS"]],"\\|")[[1]],row[["SYMBOL"]]))
        paste0(list.ids,collapse = "|")
      }
    })
  
    if(!is.na(speciesdb.name) && speciesdb.name %in% c("FLYBASE","WORMBASE","TAIR") & !is.archive){
      ensembl.data[[speciesdb.name]] <- ensembl.data$ENSEMBL
    }
  
    if(!is.na(speciesdb.name) && speciesdb.name=="TAIR"){
      ensembl.data <- ensembl.data %>% group_by(SYMBOL,GENENAME) %>% summarise(across(everything(),function(x){
        unique.info <- unique(unlist(strsplit(x,"\\|")))
        res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
        ifelse(res=="",NA,res)
      }))
      if(!is.archive){
        ensembl.data <- ensembl.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ENSEMBL",speciesdb.name)]
      } else {
        ensembl.data <- ensembl.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ENSEMBL")]
      }
      colnames(ensembl.data)[colnames(ensembl.data)!="SYMBOL" & colnames(ensembl.data)!="GENENAME"] <- paste0(colnames(ensembl.data)[colnames(ensembl.data)!="SYMBOL" & colnames(ensembl.data)!="GENENAME"]," ENSEMBL")
    } else {
      ensembl.data <- ensembl.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
        unique.info <- unique(unlist(strsplit(x,"\\|")))
        res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
        ifelse(res=="",NA,res)
      }))
      if(!is.na(speciesdb.name) && speciesdb.name %in% c("FLYBASE","WORMBASE") & !is.archive){
        ensembl.data <- ensembl.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ENSEMBL",speciesdb.name)]
      } else {
        ensembl.data <- ensembl.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ENSEMBL")]
      }
      colnames(ensembl.data)[colnames(ensembl.data)!="SYMBOL"] <- paste0(colnames(ensembl.data)[colnames(ensembl.data)!="SYMBOL"]," ENSEMBL")
    }
  }
  return(as.data.frame(ensembl.data))
}

process.ensembl.grch37.data <- function(ensembl.archive.data,speciesdb.name)
{
  colnames(ensembl.archive.data) <- c("ENSEMBLOLD ARCHIVE","GENETYPE ARCHIVE","SYMBOL")
  ensembl.archive.data <- ensembl.archive.data[!is.na(ensembl.archive.data$SYMBOL),]
  ensembl.archive.data <- ensembl.archive.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(x,"\\|")))
    res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
    ifelse(res=="",NA,res)
  }))
  ensembl.archive.data[["ALIAS ARCHIVE"]] <- ensembl.archive.data$SYMBOL
  ensembl.archive.data[["GENENAME ARCHIVE"]] <- NA
  ensembl.archive.data <- ensembl.archive.data[,c("SYMBOL","ALIAS ARCHIVE","GENETYPE ARCHIVE","GENENAME ARCHIVE","ENSEMBLOLD ARCHIVE")]
  return(as.data.frame(ensembl.archive.data))
}

process.uniprot.data <- function(uniprot.data,speciesdb.name,species)
{
  uniprot.data <- uniprot.data[uniprot.data$V2 %in% c("Gene_Name","GeneID","EnsemblGenome"),c("V1","V2","V3")]
  uniprot.data <- as.data.frame(pivot_wider(uniprot.data, names_from = V2, values_from = V3, values_fn = first))
  colnames(uniprot.data)[colnames(uniprot.data)=="V1"] <- "UNIPROT"
  colnames(uniprot.data)[colnames(uniprot.data)=="EnsemblGenome"] <- "ENSEMBL"
  colnames(uniprot.data)[colnames(uniprot.data)=="GeneID"] <- "ENTREZID"
  colnames(uniprot.data)[colnames(uniprot.data)=="Gene_Name"] <- "SYMBOL"
  uniprot.data$ALIAS <- uniprot.data$SYMBOL
  if(species %in% c("Rapeseed")){
    uniprot.data[!is.na(uniprot.data$SYMBOL) & startsWith(uniprot.data$SYMBOL,"Bna"),"SYMBOL"] <- toupper(
      uniprot.data[!is.na(uniprot.data$SYMBOL) & startsWith(uniprot.data$SYMBOL,"Bna"),"SYMBOL"])
    uniprot.data$ALIAS <- apply(uniprot.data,1,function(row){
      if(is.na(row[["ALIAS"]])){
        row[["SYMBOL"]]
      } else {
        list.ids <- unique(c(strsplit(row[["ALIAS"]],"\\|")[[1]],row[["SYMBOL"]]))
        paste0(list.ids,collapse = "|")
      }
    })
  }
  if(!is.na(speciesdb.name) && speciesdb.name=="ZFIN"){
    uniprot.data[grepl("^loc[0-9]+",uniprot.data$SYMBOL),"SYMBOL"] <- toupper(uniprot.data[grepl("^loc[0-9]+",uniprot.data$SYMBOL),"SYMBOL"])
  }
  if(!species %in% c("Tomato","Cabbage","Rapeseed")){
    uniprot.data <- uniprot.data[!is.na(uniprot.data$SYMBOL),]
    uniprot.data <- uniprot.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
      unique.info <- unique(x)
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
    uniprot.data <- uniprot.data[,c("SYMBOL","UNIPROT")]
    colnames(uniprot.data) <- c("SYMBOL","UNIPROT UNIPROT")
  }
  return(uniprot.data)
}

process.speciesdb.data <- function(speciesdb.data,species)
{
  if(species=="Human"){
    speciesdb.data <- process.hgnc.data(speciesdb.data)
  } else if(species=="Mouse") {
    speciesdb.data <- process.mgi.data(speciesdb.data)
  } else if(species=="Rat") {
    speciesdb.data <- process.rgd.data(speciesdb.data)
  } else if(species=="Yeast") {
    speciesdb.data <- process.sgd.data(speciesdb.data)
  } else if(species=="Worm") {
    speciesdb.data <- process.wormbase.data(speciesdb.data)
  } else if(species=="Fly") {
    speciesdb.data <- process.flybase.data(speciesdb.data)
  } else if(species=="Zebrafish") {
    speciesdb.data <- process.zfin.data(speciesdb.data)
  } else if(species=="Arabidopsis") {
    speciesdb.data <- process.tair.data(speciesdb.data)
  }
  return(speciesdb.data)
}

process.hgnc.data <- function(speciesdb.data)
{
  speciesdb.data <- speciesdb.data[,c("symbol","locus_group","name","alias_symbol","hgnc_id")]
  colnames(speciesdb.data) <- c("SYMBOL","GENETYPE HGNC","GENENAME HGNC","ALIAS HGNC","HGNC HGNC")
  speciesdb.data$`ALIAS HGNC` <- ifelse(is.na(speciesdb.data$`ALIAS HGNC`),speciesdb.data$SYMBOL,paste0(speciesdb.data$`ALIAS HGNC`,"|",speciesdb.data$SYMBOL))
  return(speciesdb.data)
}

process.mgi.data <- function(speciesdb.data)
{
  speciesdb.data <- speciesdb.data[speciesdb.data$`Marker Type` %in% c("Gene","Pseudogene"),]
  speciesdb.data <- speciesdb.data[,c("Marker Symbol","Marker Synonyms (pipe-separated)","Feature Type","Marker Name","MGI Accession ID")]
  colnames(speciesdb.data) <- c("SYMBOL","ALIAS MGI","GENETYPE MGI","GENENAME MGI","MGI MGI")
  speciesdb.data[speciesdb.data==""] <- NA
  speciesdb.data <- speciesdb.data[!is.na(speciesdb.data$SYMBOL),]
  speciesdb.data$`ALIAS MGI` <-ifelse(is.na(speciesdb.data$`ALIAS MGI`),speciesdb.data$SYMBOL,paste0(speciesdb.data$`ALIAS MGI`,"|",speciesdb.data$SYMBOL))
  speciesdb.data <- speciesdb.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    unique.info <- unique(x)
    res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
    ifelse(res=="",NA,res)
  }))
  return(as.data.frame(speciesdb.data))
}

process.rgd.data <- function(speciesdb.data)
{
  speciesdb.data <- speciesdb.data[,c("SYMBOL","OLD_SYMBOL","GENE_TYPE","NAME","GENE_RGD_ID")]
  colnames(speciesdb.data) <- c("SYMBOL","ALIAS SPECIAL","GENETYPE SPECIAL","GENENAME SPECIAL","RGD SPECIAL")
  speciesdb.data <- speciesdb.data[!startsWith(speciesdb.data$SYMBOL,"ENSRNOG"),]
  speciesdb.data <- speciesdb.data %>% separate_rows("ALIAS SPECIAL", sep = ";")
  speciesdb.data <- speciesdb.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    unique.info <- unique(x)
    res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
    ifelse(res=="",NA,res)
  }))
  speciesdb.data$`GENETYPE SPECIAL` <- apply(speciesdb.data,1,function(row){
    if(is.na(row[["GENETYPE SPECIAL"]])){
      return(NA)
    } else {
      unique.info <- strsplit(row[["GENETYPE SPECIAL"]],"\\|")[[1]]
      return(unique.info[length(unique.info)])
    }
  })
  speciesdb.data$`ALIAS SPECIAL` <- apply(speciesdb.data,1,function(row){
    ifelse(is.na(row["ALIAS SPECIAL"]),row["SYMBOL"],paste0(row["ALIAS SPECIAL"],"|",row["SYMBOL"]))
  })
  return(as.data.frame(speciesdb.data))
}

process.sgd.data <- function(speciesdb.data)
{
  speciesdb.data$dbxref <- gsub("SGD:","",speciesdb.data$dbxref)
  speciesdb.data[!is.na(speciesdb.data$so_term_name) & speciesdb.data$so_term_name=="protein_coding_gene","so_term_name"] <- "protein-coding"
  speciesdb.data$GENENAME <- unlist(lapply(speciesdb.data$Alias,function(alias){
    res <- alias[grepl("^[[:lower:]]+$|( )+", alias)]
    ifelse(length(res)==0,NA,res[which.max(nchar(res))])
  }))
  speciesdb.data$ALIAS <- unlist(lapply(speciesdb.data$Alias,function(alias){
    res <- alias[!grepl("^[[:lower:]]+$|( )+", alias)]
    ifelse(length(res)==0,NA,paste0(res,collapse="|"))
  }))
  speciesdb.data$ALIAS <- apply(speciesdb.data,1,function(row){
    ifelse(is.na(row[["ALIAS"]]),row[["Name"]],paste0(row[["ALIAS"]],"|",row[["Name"]]))
  })
  speciesdb.data$gene <- ifelse(is.na(speciesdb.data$gene),speciesdb.data$Name,speciesdb.data$gene)
  speciesdb.data <- speciesdb.data[,c("gene","ALIAS","so_term_name","GENENAME","dbxref")]
  colnames(speciesdb.data) <- c("SYMBOL","ALIAS SPECIAL","GENETYPE SPECIAL","GENENAME SPECIAL","SGD SPECIAL")
  return(as.data.frame(speciesdb.data))
}

process.wormbase.data <- function(speciesdb.data)
{
  speciesdb.data <- speciesdb.data[,c(3,4,6,2)]
  colnames(speciesdb.data) <- c("SYMBOL","ALIAS SPECIAL","GENETYPE SPECIAL","WORMBASE SPECIAL")
  speciesdb.data <- speciesdb.data[!is.na(speciesdb.data$SYMBOL),]
  speciesdb.data$`ALIAS SPECIAL` <- ifelse(is.na(speciesdb.data$`ALIAS SPECIAL`),speciesdb.data$SYMBOL,
                paste0(speciesdb.data$`ALIAS SPECIAL`,"|",speciesdb.data$SYMBOL))
  speciesdb.data[["ENSEMBL SPECIAL"]] <- speciesdb.data[["WORMBASE SPECIAL"]]
  return(as.data.frame(speciesdb.data))
}

process.flybase.data <- function(speciesdb.data)
{
  speciesdb.data <- speciesdb.data[,c("V1","V3","V4","V5","V6")]
  colnames(speciesdb.data) <- c("SYMBOL","FLYBASE SPECIAL","FLYBASE ALT SPECIAL","ALIAS SPECIAL","ALIAS ALT SPECIAL")
  speciesdb.data$`FlyBase SPECIAL` <- apply(speciesdb.data,1,function(row){
    if(is.na(row[["FLYBASE ALT SPECIAL"]])){
      row[["FLYBASE SPECIAL"]]
    } else {
      paste0(paste0(strsplit(row[["FLYBASE ALT SPECIAL"]],",")[[1]],collapse="|"),"|",row[["FLYBASE SPECIAL"]])
    }
  })
  speciesdb.data$`ALIAS SPECIAL` <- apply(speciesdb.data,1,function(row){
    list.alias <- row[["ALIAS SPECIAL"]]
    if(!is.na(row[["ALIAS ALT SPECIAL"]])){
      list.alias <- c(list.alias,strsplit(row[["ALIAS ALT SPECIAL"]],",")[[1]])
    }
    list.alias <- unique(c(list.alias,row[["SYMBOL"]]))
    paste0(list.alias,collapse="|")
  })
  speciesdb.data[["ENSEMBL SPECIAL"]] <- speciesdb.data[["FLYBASE SPECIAL"]]
  speciesdb.data <- speciesdb.data[,c("SYMBOL","ALIAS SPECIAL","ENSEMBL SPECIAL","FLYBASE SPECIAL")]
  return(as.data.frame(speciesdb.data))
}

process.zfin.data <- function(speciesdb.data)
{
  speciesdb.data$ALIAS <- unlist(lapply(speciesdb.data$Alias,function(alias){
    alias.names <- alias[!startsWith(alias,"ENSDARG") & !startsWith(alias,"ENSDARP") & !startsWith(alias,"NM_") & !startsWith(alias,"NR_")
                         & !startsWith(alias,"NP_") & !startsWith(alias,"XM_") & !startsWith(alias,"XR_") & !startsWith(alias,"XP_")]
    alias.names <- alias.names[!grepl("[ ]+",alias.names)]
    ifelse(length(alias.names)==0,NA,paste0(alias.names,collapse="|"))
  }))
  speciesdb.data$ALIAS <- ifelse(is.na(speciesdb.data$ALIAS),speciesdb.data$Name,paste0(speciesdb.data$ALIAS,"|",speciesdb.data$Name))
  speciesdb.data$GENENAME <- unlist(lapply(speciesdb.data$full_name,function(gene.name){
    ifelse(length(gene.name)==0,NA,paste0(gene.name,collapse=","))
  }))
  speciesdb.data$SYMBOL <- speciesdb.data$Name
  speciesdb.data[grepl("^loc[0-9]+",speciesdb.data$SYMBOL),"SYMBOL"] <- toupper(speciesdb.data[grepl("^loc[0-9]+",speciesdb.data$SYMBOL),"SYMBOL"])
  speciesdb.data$GENETYPE <- speciesdb.data$so_term_name
  speciesdb.data$ENSEMBL <- unlist(lapply(speciesdb.data$Alias,function(alias){
    alias.names <- alias[startsWith(alias,"ENSDARG")]
    ifelse(length(alias.names)==0,NA,paste0(alias.names,collapse="|"))
  }))
  speciesdb.data$ZFIN <- unlist(lapply(speciesdb.data$secondaryIds,function(zfin.ids){
    zfin.ids <- zfin.ids[startsWith(zfin.ids,"ZDB-GENE")]
    ifelse(length(zfin.ids)==0,NA,paste0(zfin.ids,collapse="|"))
  }))
  speciesdb.data$ZFIN <- ifelse(is.na(speciesdb.data$ZFIN),speciesdb.data$ID,paste0(speciesdb.data$ID,"|",speciesdb.data$ZFIN))
  speciesdb.data <- speciesdb.data[,c("SYMBOL","ALIAS","GENETYPE","GENENAME","ZFIN","ENSEMBL")]
  colnames(speciesdb.data) <- c("SYMBOL","ALIAS SPECIAL","GENETYPE SPECIAL","GENENAME SPECIAL","ZFIN SPECIAL","ENSEMBL SPECIAL")
  return(as.data.frame(speciesdb.data))
}

process.tair.data <- function(speciesdb.data)
{
  speciesdb.data$GENENAME <- unlist(lapply(speciesdb.data$computational_description,function(gene.name){
    ifelse(length(gene.name)==0,NA,paste0(gene.name,collapse=","))
  }))
  speciesdb.data <- speciesdb.data[!is.na(speciesdb.data$symbol),]
  colnames(speciesdb.data)[colnames(speciesdb.data)=="locus_type"] <- "GENETYPE SPECIAL"
  speciesdb.data <- speciesdb.data[is.na(speciesdb.data$`GENETYPE SPECIAL`) | !speciesdb.data$`GENETYPE SPECIAL` %in% c("pre_trna","novel_transcribed_region"),]
  colnames(speciesdb.data)[colnames(speciesdb.data)=="symbol"] <- "SYMBOL"
  colnames(speciesdb.data)[colnames(speciesdb.data)=="ID"] <- "TAIR SPECIAL"
  speciesdb.data[["ENSEMBL SPECIAL"]] <- speciesdb.data[["TAIR SPECIAL"]]
  speciesdb.data <- speciesdb.data[,c("SYMBOL","GENETYPE SPECIAL","GENENAME","ENSEMBL SPECIAL","TAIR SPECIAL")]
  speciesdb.data <- speciesdb.data %>% group_by(SYMBOL,GENENAME) %>% summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(x,"\\|")))
    res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
    ifelse(res=="",NA,res)
  }))
  return(as.data.frame(speciesdb.data))
}

process.ncbi.orthologs.data <- function(ncbi.orthologs,taxid,taxonomy.table,annotation.data.list)
{
  #Filter data about the requested species
  ncbi.orthologs.ref <- ncbi.orthologs[ncbi.orthologs$`#tax_id` %in% taxid,]
  ncbi.orthologs.other <- ncbi.orthologs[ncbi.orthologs$`Other_tax_id` %in% taxid,]
  ncbi.orthologs <- data.frame("GeneID"=c(ncbi.orthologs.ref$GeneID,ncbi.orthologs.other$Other_GeneID),
                               "Other_tax_id"=c(ncbi.orthologs.ref$Other_tax_id,ncbi.orthologs.other$`#tax_id`),
                               "Other_GeneID"=c(ncbi.orthologs.ref$Other_GeneID,ncbi.orthologs.other$GeneID))
  
  #Map ids to gene symbols
  ncbi.orthologs <- map.ortho.ids.to.symbol(annotation.data.list,ncbi.orthologs,taxid,taxonomy.table,"ENTREZID")
  
  #Reshape orthology data with one column for each species
  ncbi.orthologs <- reshape.orthology.data(ncbi.orthologs,taxid,taxonomy.table)
  colnames(ncbi.orthologs)[colnames(ncbi.orthologs)!="SYMBOL"] <- paste0(colnames(ncbi.orthologs)[colnames(ncbi.orthologs)!="SYMBOL"]," NCBI")
  
  return(ncbi.orthologs)
}

process.ensembl.orthologs.data <- function(ensembl.orthologs,taxid,taxonomy.table,annotation.data.list)
{
  #Rename columns
  colnames(ensembl.orthologs) <- c("GeneID","Other_tax_id","Other_GeneID")
  #Check ensembl ids
  ensembl.orthologs$GeneID <- gsub("gene-","",ensembl.orthologs$GeneID)
  ensembl.orthologs$Other_GeneID <- gsub("gene-","",ensembl.orthologs$Other_GeneID)
  #Map ids to gene symbols
  ensembl.orthologs <- map.ortho.ids.to.symbol(annotation.data.list,ensembl.orthologs,taxid,taxonomy.table,"ENSEMBL")
  #Reshape orthology data with one column for each species
  ensembl.orthologs <- reshape.orthology.data(ensembl.orthologs,taxid,taxonomy.table)
  colnames(ensembl.orthologs)[colnames(ensembl.orthologs)!="SYMBOL"] <- paste0(colnames(ensembl.orthologs)[colnames(ensembl.orthologs)!="SYMBOL"]," ENSEMBL")
  return(ensembl.orthologs)
}

process.alliance.orthologs.data <- function(alliance.orthologs,taxid,taxonomy.table)
{
  #Filter data about the requested species and the species present in the geneslator DB
  alliance.orthologs$Gene1SpeciesTaxonID <- gsub("NCBITaxon:","",alliance.orthologs$Gene1SpeciesTaxonID)
  alliance.orthologs$Gene2SpeciesTaxonID <- gsub("NCBITaxon:","",alliance.orthologs$Gene2SpeciesTaxonID)
  alliance.orthologs <- alliance.orthologs[alliance.orthologs$Gene1SpeciesTaxonID %in% taxid | alliance.orthologs$Gene2SpeciesTaxonID %in% taxid,]
  alliance.orthologs.ref <- alliance.orthologs[alliance.orthologs$Gene1SpeciesTaxonID %in% taxid,]
  alliance.orthologs.other <- alliance.orthologs[alliance.orthologs$Gene2SpeciesTaxonID %in% taxid,]
  alliance.orthologs <- data.frame("Ortho_TAXID"=c(alliance.orthologs.ref$Gene2SpeciesTaxonID,alliance.orthologs.other$Gene1SpeciesTaxonID),
                                   "SYMBOL"=c(alliance.orthologs.ref$Gene1Symbol,alliance.orthologs.other$Gene2Symbol),
                                   "Ortho_SYMBOL"=c(alliance.orthologs.ref$Gene2Symbol,alliance.orthologs.other$Gene1Symbol))
  alliance.orthologs <- unique(alliance.orthologs)
  
  #Reshape orthology data with one column for each species
  alliance.orthologs <- reshape.orthology.data(alliance.orthologs,taxid,taxonomy.table)
  colnames(alliance.orthologs)[colnames(alliance.orthologs)!="SYMBOL"] <- paste0(colnames(alliance.orthologs)[colnames(alliance.orthologs)!="SYMBOL"]," ALLIANCE")
 
  return(alliance.orthologs) 
}

process.speciesdb.orthologs.data <- function(speciesdb.orthologs,taxid,taxonomy.table,annotation.data.list)
{
  if(species=="Human"){
    speciesdb.orthologs <- process.hgnc.orthologs.data(speciesdb.orthologs,taxid,taxonomy.table,annotation.data.list)
  }
  return(speciesdb.orthologs)
}

process.hgnc.orthologs.data <- function(speciesdb.orthologs,taxid,taxonomy.table,annotation.data.list)
{
  #Filter out rows where ids of genes are missing
  speciesdb.orthologs <- speciesdb.orthologs[speciesdb.orthologs$human_entrez_gene!="-" | speciesdb.orthologs$human_ensembl_gene!="-",]
  speciesdb.orthologs <- speciesdb.orthologs[speciesdb.orthologs$ortholog_species_entrez_gene!="-" | speciesdb.orthologs$ortholog_species_ensembl_gene!="-",]
  #Select ref gene ids for human and ortho species
  speciesdb.orthologs$human_entrez_gene <- ifelse(speciesdb.orthologs$human_entrez_gene=="-",speciesdb.orthologs$human_ensembl_gene,speciesdb.orthologs$human_entrez_gene)
  speciesdb.orthologs$ortholog_species_entrez_gene <- ifelse(speciesdb.orthologs$ortholog_species_entrez_gene=="-",speciesdb.orthologs$ortholog_species_ensembl_gene,speciesdb.orthologs$ortholog_species_entrez_gene)
  #Take only relevant columns for the analysis
  speciesdb.orthologs <- speciesdb.orthologs[,c("human_entrez_gene","ortholog_species","ortholog_species_entrez_gene")]
  colnames(speciesdb.orthologs) <- c("GeneID","Other_tax_id","Other_GeneID")
  
  #Get species info
  valid.pos <- sapply(taxonomy.table$taxid, function(x) any(x %in% taxid))
  species.info <- strsplit(taxonomy.table[valid.pos,"official_name"]," ")[[1]]
  species.name <- taxonomy.table[valid.pos,"species"]
  
  #Map ref species ids to gene symbols
  map.ref.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species.name]],speciesdb.orthologs$GeneID,"SYMBOL","ENTREZID")
  speciesdb.orthologs <- merge(speciesdb.orthologs,map.ref.ids.to.symbol,by.x="GeneID",by.y="ENTREZID",all.x=T)
  colnames(speciesdb.orthologs)[colnames(speciesdb.orthologs)=="SYMBOL"] <- "SYM_ENTREZ"
  map.ref.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species.name]],speciesdb.orthologs$GeneID,"SYMBOL","ENSEMBL")
  speciesdb.orthologs <- merge(speciesdb.orthologs,map.ref.ids.to.symbol,by.x="GeneID",by.y="ENSEMBL",all.x=T)
  speciesdb.orthologs[["Ref_Symbol"]] <- ifelse(is.na(speciesdb.orthologs$SYM_ENTREZ),speciesdb.orthologs$SYMBOL,speciesdb.orthologs$SYM_ENTREZ)
  speciesdb.orthologs <- as.data.frame(speciesdb.orthologs[!is.na(speciesdb.orthologs$Ref_Symbol),c("Other_tax_id","Other_GeneID","Ref_Symbol")])
  
  #Map other species ids to gene symbols
  map.other.ids.to.symbol <- data.frame(ID=character(),SYMBOL=character())
  for(i in 1:nrow(taxonomy.table)) {
    if(all(unlist(taxonomy.table[i,"taxid"])!=taxid)){
      list.gene.ids <- speciesdb.orthologs[speciesdb.orthologs$`Other_tax_id` %in% unlist(taxonomy.table[i,"taxid"]),"Other_GeneID"]
      other.species.info <- strsplit(taxonomy.table[i,"official_name"],"_")[[1]]
      other.species.name <- taxonomy.table[i,"species"]
      res.table.entrez <- map.keys.to.values(annotation.data.list[[other.species.name]],list.gene.ids,"SYMBOL","ENTREZID")
      colnames(res.table.entrez)[colnames(res.table.entrez)=="SYMBOL"] <- "SYM_ENTREZ"
      res.table.ensembl <- map.keys.to.values(annotation.data.list[[other.species.name]],list.gene.ids,"SYMBOL","ENSEMBL")
      res.table <- merge(res.table.entrez,res.table.ensembl,by.x="ENTREZID",by.y="ENSEMBL")
      res.table$SYMBOL <- ifelse(is.na(res.table$SYM_ENTREZ),res.table$SYMBOL,res.table$SYM_ENTREZ)
      res.table <- res.table[,c("ENTREZID","SYMBOL")]
      colnames(res.table) <- c("ID","SYMBOL")
      map.other.ids.to.symbol <- rbind(map.other.ids.to.symbol,res.table)
    }
  }
  speciesdb.orthologs <- merge(speciesdb.orthologs,map.other.ids.to.symbol,by.x="Other_GeneID",by.y="ID",all.x=T)
  speciesdb.orthologs <- speciesdb.orthologs[,c("Other_tax_id","Ref_Symbol","SYMBOL")]
  colnames(speciesdb.orthologs) <- c("Ortho_TAXID","SYMBOL","Ortho_SYMBOL")
  
  #Reshape orthology data with one column for each species
  speciesdb.orthologs <- reshape.orthology.data(speciesdb.orthologs,taxid,taxonomy.table)
  colnames(speciesdb.orthologs)[colnames(speciesdb.orthologs)!="SYMBOL"] <- paste0(colnames(speciesdb.orthologs)[colnames(speciesdb.orthologs)!="SYMBOL"]," SPECIAL")
  
  return(speciesdb.orthologs)
}

reshape.orthology.data <- function(orthologs.data,taxId,taxonomy.table)
{
  for(i in 1:nrow(taxonomy.table)) {
    if(all(unlist(taxonomy.table[i,"taxid"])!=taxId)){
      orthologs.data[[paste0("ORTHO",toupper(taxonomy.table[i,"species"]))]] <- ifelse(orthologs.data$Ortho_TAXID %in% unlist(taxonomy.table[i,"taxid"]),orthologs.data$Ortho_SYMBOL,NA)
    }
  }
  orthologs.data <- orthologs.data[,!startsWith(colnames(orthologs.data),"Ortho_")]
  orthologs.data <- orthologs.data %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
    res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
    ifelse(res=="",NA,res)
  }))
  return(orthologs.data)
}

process.go.dictionary <- function(go.dictionary)
{
  go.dictionary <- go.dictionary[,c("id","alt_id","name","namespace")]
  go.dictionary$id <- apply(go.dictionary,1,function(row){
    if(row[["alt_id"]]==""){
      return(row[["id"]])
    } else {
      return(paste0(c(row[["id"]],strsplit(row[["alt_id"]],"; ")[[1]]),collapse = "|"))
    }
  })
  go.dictionary <- go.dictionary[,c("id","name","namespace")]
  colnames(go.dictionary) <- c("GO","GONAME","GOTYPE")
  go.dictionary <- go.dictionary[go.dictionary$GOTYPE!="external",]
  go.dictionary[go.dictionary$GOTYPE=="biological_process","GOTYPE"] <- "BP"
  go.dictionary[go.dictionary$GOTYPE=="cellular_component","GOTYPE"] <- "CC"
  go.dictionary[go.dictionary$GOTYPE=="molecular_function","GOTYPE"] <- "MF"
  go.dictionary <- as.data.frame(go.dictionary %>% separate_rows(GO,sep="\\|"))
  return(go.dictionary)
}

process.go.data <- function(go.data,go.dictionary,annotation.data.list,species)
{
  go.data <- go.data[go.data$V1=="UniProtKB",c(2,5,7)]
  colnames(go.data) <- c("UNIPROT","GO","GOEVIDENCE")
  #if(speciesdb.name=="TAIR"){
  #  go.data[grep("^At[0-9]g",go.data$SYMBOL),"SYMBOL"] <- toupper(go.data[grep("^At[0-9]g",go.data$SYMBOL),"SYMBOL"])
  #}
  go.data <- unique(go.data)
  go.data <- merge(go.data,go.dictionary,all.x=T)
  go.data <- go.data %>% group_by(UNIPROT) %>% summarise(across(everything(),function(x){
    unique.info <- unlist(strsplit(as.character(x),"\\|"))
    res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
    ifelse(res=="",NA,res)
  }))
  map.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species]],
    go.data$UNIPROT,"SYMBOL","UNIPROT")
  go.data <- merge(go.data,map.ids.to.symbol,all.x=T)
  go.data <- go.data[!is.na(go.data$SYMBOL),c("SYMBOL","GO","GOEVIDENCE","GONAME","GOTYPE")]
  return(go.data)
}

process.reactome.data <- function(reactome.data,annotation.data.list,species)
{
  #Process reactome data for NCBI entrez IDs
  reactome.data.ncbi <- reactome.data$ncbi
  reactome.data.ncbi <- reactome.data.ncbi[grep("^[0-9]+$",reactome.data.ncbi$V1),c(1,2,4)]
  reactome.data.ncbi <- unique(reactome.data.ncbi)
  colnames(reactome.data.ncbi) <- c("ENTREZID","REACTOMEPATH NCBI","REACTOMEPATHNAME NCBI")
  reactome.data.ncbi <- reactome.data.ncbi %>% group_by(ENTREZID) %>% summarise(across(everything(),function(x){
    res <- paste0(as.character(x),collapse="|")
    ifelse(res=="",NA,res)
  }))
  map.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species]],reactome.data.ncbi$ENTREZID,"SYMBOL","ENTREZID")
  reactome.data.ncbi <- merge(reactome.data.ncbi,map.ids.to.symbol,all.x=T)
  reactome.data.ncbi <- reactome.data.ncbi[!is.na(reactome.data.ncbi$SYMBOL),c("SYMBOL","REACTOMEPATH NCBI","REACTOMEPATHNAME NCBI")]
  
  #Process reactome data for Ensembl IDs
  reactome.data.ens <- reactome.data$ensembl
  reactome.data.ens <- reactome.data.ens[,c(1,2,4)]
  reactome.data.ens <- unique(reactome.data.ens)
  colnames(reactome.data.ens) <- c("ENSEMBL","REACTOMEPATH ENSEMBL","REACTOMEPATHNAME ENSEMBL")
  reactome.data.ens <- reactome.data.ens %>% group_by(ENSEMBL) %>% summarise(across(everything(),function(x){
    res <- paste0(as.character(x),collapse="|")
    ifelse(res=="",NA,res)
  }))
  map.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species]],reactome.data.ens$ENSEMBL,"SYMBOL","ENSEMBL")
  reactome.data.ens <- merge(reactome.data.ens,map.ids.to.symbol,all.x=T)
  reactome.data.ens <- reactome.data.ens[!is.na(reactome.data.ens$SYMBOL),c("SYMBOL","REACTOMEPATH ENSEMBL","REACTOMEPATHNAME ENSEMBL")]
  
  #Merge data
  reactome.data <- merge(reactome.data.ncbi,reactome.data.ens,all=T)
  reactome.data$REACTOMEPATH <- apply(reactome.data,1,function(row){
    if(!is.na(row[["REACTOMEPATH NCBI"]]) && !is.na(row[["REACTOMEPATH ENSEMBL"]])){
      paste0(row[["REACTOMEPATH NCBI"]],"|",row[["REACTOMEPATH ENSEMBL"]])
    } else if(is.na(row[["REACTOMEPATH NCBI"]]) && !is.na(row[["REACTOMEPATH ENSEMBL"]])){
      row[["REACTOMEPATH ENSEMBL"]]
    } else if(is.na(row[["REACTOMEPATH ENSEMBL"]]) && !is.na(row[["REACTOMEPATH NCBI"]])){
      row[["REACTOMEPATH NCBI"]]
    } else {
      NA
    }
  })
  reactome.data$REACTOMEPATHNAME <- apply(reactome.data,1,function(row){
    if(!is.na(row[["REACTOMEPATHNAME NCBI"]]) && !is.na(row[["REACTOMEPATHNAME ENSEMBL"]])){
      paste0(row[["REACTOMEPATHNAME NCBI"]],"|",row[["REACTOMEPATHNAME ENSEMBL"]])
    } else if(is.na(row[["REACTOMEPATHNAME NCBI"]]) && !is.na(row[["REACTOMEPATHNAME ENSEMBL"]])){
      row[["REACTOMEPATHNAME ENSEMBL"]]
    } else if(is.na(row[["REACTOMEPATHNAME ENSEMBL"]]) && !is.na(row[["REACTOMEPATHNAME NCBI"]])){
      row[["REACTOMEPATHNAME NCBI"]]
    } else {
      NA
    }
  })
  reactome.data <- reactome.data[,c("SYMBOL","REACTOMEPATH","REACTOMEPATHNAME")]
  reactome.data <- reactome.data %>% separate_rows(all_of(c("REACTOMEPATH","REACTOMEPATHNAME")),sep="\\|")
  reactome.data <- unique(reactome.data)
  reactome.data <- reactome.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    res <- paste0(as.character(x),collapse="|")
    ifelse(res=="",NA,res)
  }))
  reactome.data <- as.data.frame(reactome.data)
  
  return(reactome.data)
}

process.wikipathways.data <- function(wikipathways.data,annotation.data.list,species)
{
  wikipathways.data <- unique(wikipathways.data[,c("gene","wpid","name")])
  colnames(wikipathways.data) <- c("ENTREZID","WIKIPATH","WIKIPATHNAME")
  wikipathways.data <- wikipathways.data %>% group_by(ENTREZID) %>% summarise(across(everything(),function(x){
    res <- paste0(as.character(x),collapse="|")
    ifelse(res=="",NA,res)
  }))
  map.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species]],wikipathways.data$ENTREZID,"SYMBOL","ENTREZID")
  wikipathways.data <- merge(wikipathways.data,map.ids.to.symbol,all.x=T)
  wikipathways.data <- wikipathways.data[!is.na(wikipathways.data$SYMBOL),c("SYMBOL","WIKIPATH","WIKIPATHNAME")]
  return(wikipathways.data)
}
