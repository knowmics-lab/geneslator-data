merge.ncbi.data <- function(ncbi.data,ncbi.replaced.data,ncbi.discontinued.data,speciesdb.name)
{
  #Merge with replaced data
  ncbi.data <- merge(ncbi.data,ncbi.replaced.data,all.x=T,by.x="ENTREZID NCBI",by.y="ENTREZID")
  ncbi.data$`ENTREZIDOLD NCBI` <- merge.columns(ncbi.data,"ENTREZIDOLD",single.val=F)
  ncbi.data$`ALIAS NCBI` <- apply(ncbi.data,1,function(row){
    if(is.na(row[["SYMBOL ARCHIVE"]])){
      row[["ALIAS NCBI"]]
    } else {
      list.ids <- unique(c(strsplit(row[["SYMBOL ARCHIVE"]],"\\|")[[1]],strsplit(row[["ALIAS NCBI"]],"\\|")[[1]]))
      paste0(list.ids,collapse = "|")
    }
  })
  ncbi.data$`ENTREZIDOLD ARCHIVE` <- NULL
  ncbi.data$`SYMBOL ARCHIVE` <- NULL
  
  #Merge with discontinued data
  colnames(ncbi.discontinued.data)[colnames(ncbi.discontinued.data)=="ENSEMBL ARCHIVE"] <- "ENSEMBLOLD ARCHIVE"
  ncbi.data <- merge(ncbi.data,ncbi.discontinued.data,all=T)
  ncbi.data$`ALIAS NCBI` <- merge.columns(ncbi.data,"ALIAS",single.val=F)
  if("LOCUS NCBI" %in% colnames(ncbi.data)){
    ncbi.data$`LOCUS NCBI` <- merge.columns(ncbi.data,"LOCUS",single.val=F)
  }
  ncbi.data$`GENETYPE NCBI` <- merge.columns(ncbi.data,"GENETYPE",single.val=T)
  ncbi.data$`GENENAME NCBI` <- merge.columns(ncbi.data,"GENENAME",single.val=T)
  ncbi.data$`ENSEMBLOLD NCBI` <- merge.columns(ncbi.data,"ENSEMBLOLD",single.val=F)
  if(!is.na(speciesdb.name)){
    ncbi.data[[paste0(speciesdb.name," NCBI")]] <- merge.columns(ncbi.data,speciesdb.name,single.val=F)
  }
  ncbi.data$`ENTREZIDOLD NCBI` <- merge.columns(ncbi.data,"ENTREZIDOLD",single.val=F)
  ncbi.data <- ncbi.data[,!endsWith(colnames(ncbi.data),"ARCHIVE")]
  
  return(ncbi.data)
}

merge.columns <- function(annotation.data,column,single.val)
{
  ref.cols <- colnames(annotation.data)[startsWith(colnames(annotation.data),paste0(column," "))]
  return(apply(annotation.data,1,function(row){
    unique.info <- unname(unlist(sapply(ref.cols,function(col){
      strsplit(as.character(row[col]),"\\|")
    })))
    unique.info <- unique(unique.info[!is.na(unique.info)])
    if(single.val){
      res <- unique.info[1]
    } else {
      res <- paste0(unique.info,collapse="|")
    }
    ifelse(res=="",NA,res)
  }))
}

rename.genetypes <- function(annotation.data)
{
  annotation.data[is.na(annotation.data$GENETYPE),"GENETYPE"] <- "unknown"
  annotation.data[annotation.data$GENETYPE %in% c("unclassified gene","gene","gene segment","heritable phenotypic marker","processed_transcript"),"GENETYPE"] <- "other"
  annotation.data[annotation.data$GENETYPE %in% c("processed_pseudogene","unprocessed_pseudogene","transcribed_processed_pseudogene","polymorphic pseudogene",
                                                  "transcribed_unprocessed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","transcribed_unitary_pseudogene",
                                                  "unitary_pseudogene","IG_V_pseudogene","IG_C_pseudogene","translated_processed_pseudogene","rRNA_pseudogene","IG_J_pseudogene",
                                                  "IG_pseudogene","pseudogene","translated_unprocessed_pseudogene","IG_D_pseudogene","pseudogenic gene segment","polymorphic_pseudogene"),"GENETYPE"] <- "pseudo"
  annotation.data[annotation.data$GENETYPE %in% c("protein_coding","protein coding gene","protein_coding_gene"),"GENETYPE"] <- "protein-coding"
  annotation.data[annotation.data$GENETYPE %in% c("ribozyme","ribozyme gene","unclassified non-coding RNA gene","non-coding RNA gene","non-coding RNA","ncrna","sense_intronic","sense_overlapping","antisense","y_rna","Y_RNA",
                                                  "3prime_overlapping_ncrna"),"GENETYPE"] <- "ncRNA"
  annotation.data[annotation.data$GENETYPE %in% c("antisense lncRNA gene","lincRNA gene","sense overlapping lncRNA gene","sense intronic lncRNA gene","lncRNA gene","lincrna","lncrna","lincRNA"),"GENETYPE"] <- "lncRNA"
  annotation.data[annotation.data$GENETYPE %in% c("miRNA gene","mirna"),"GENETYPE"] <- "miRNA"
  annotation.data[annotation.data$GENETYPE %in% c("snoRNA gene","snorna"),"GENETYPE"] <- "snoRNA"
  annotation.data[annotation.data$GENETYPE %in% c("Mt_rRNA","rRNA gene","rrna","mt_rrna"),"GENETYPE"] <- "rRNA"
  annotation.data[annotation.data$GENETYPE %in% c("snRNA gene","snrna"),"GENETYPE"] <- "snRNA"
  annotation.data[annotation.data$GENETYPE %in% c("IG_V_gene","IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","ig_v_gene"),"GENETYPE"] <- "IG-gene"
  annotation.data[annotation.data$GENETYPE %in% c("TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene","tr_v_gene","tr_c_gene","tr_j_gene"),"GENETYPE"] <- "TR-gene"
  annotation.data[annotation.data$GENETYPE %in% c("misc_RNA","misc_rna"),"GENETYPE"] <- "miscRNA"
  annotation.data[annotation.data$GENETYPE %in% c("Mt_tRNA","mt_trna","trna"),"GENETYPE"] <- "tRNA"
  annotation.data[annotation.data$GENETYPE %in% c("vault_RNA"),"GENETYPE"] <- "vaultRNA"
  annotation.data[annotation.data$GENETYPE %in% c("tec"),"GENETYPE"] <- "TEC"
  annotation.data[annotation.data$GENETYPE %in% c("scarna"),"GENETYPE"] <- "scaRNA"
  return(annotation.data)
}

merge.databases <- function(ncbi.data,ensembl.data,uniprot.data,speciesdb.data,speciesdb.name,species)
{
  if(species %in% c("Tomato","Cabbage","Rapeseed")){
    
    #Join NCBI database with UNIPROT data
    uniprot.data.ncbi <- uniprot.data[!is.na(uniprot.data$ENTREZID),c("UNIPROT","ENTREZID","ENSEMBL","ALIAS")]
    uniprot.data.ncbi <- uniprot.data.ncbi %>% group_by(ENTREZID) %>% summarise(across(everything(),function(x){
      unique.info <- unique(x)
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
    colnames(uniprot.data.ncbi)[colnames(uniprot.data.ncbi)=="ENSEMBL"] <- "ENSEMBL UNIPROT NCBI"
    colnames(uniprot.data.ncbi)[colnames(uniprot.data.ncbi)=="ALIAS"] <- "ALIAS UNIPROT NCBI"
    colnames(ncbi.data)[colnames(ncbi.data)=="SYMBOL"] <- "SYMBOL NCBI"
    colnames(ncbi.data)[colnames(ncbi.data)=="ENTREZID NCBI"] <- "ENTREZID"
    annotation.data.ncbi <- merge(ncbi.data,uniprot.data.ncbi,all.x=T)
    colnames(annotation.data.ncbi)[colnames(annotation.data.ncbi)=="SYMBOL NCBI"] <- "SYMBOL"
    colnames(annotation.data.ncbi)[colnames(annotation.data.ncbi)=="UNIPROT"] <- "UNIPROT NCBI"
    colnames(annotation.data.ncbi)[colnames(annotation.data.ncbi)=="ENTREZID"] <- "ENTREZID NCBI"
    colnames(annotation.data.ncbi)[colnames(annotation.data.ncbi)=="ENSEMBL"] <- "ENSEMBL NCBI"
    
    #Join Ensembl database with UNIPROT data
    uniprot.data.ens <- uniprot.data[!is.na(uniprot.data$ENSEMBL),]
    uniprot.data.ens <- uniprot.data.ens %>% group_by(ENSEMBL) %>% summarise(across(everything(),function(x){
      unique.info <- unique(x)
      res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
      ifelse(res=="",NA,res)
    }))
    colnames(uniprot.data.ens)[colnames(uniprot.data.ens)=="ENTREZID"] <- "ENTREZID UNIPROT ENS"
    colnames(uniprot.data.ens)[colnames(uniprot.data.ens)=="SYMBOL"] <- "SYMBOL UNIPROT ENS"
    colnames(uniprot.data.ens)[colnames(uniprot.data.ens)=="ALIAS"] <- "ALIAS UNIPROT ENS"
    colnames(ensembl.data)[colnames(ensembl.data)=="SYMBOL"] <- "SYMBOL ENSEMBL"
    colnames(ensembl.data)[colnames(ensembl.data)=="ENSEMBL ENSEMBL"] <- "ENSEMBL"
    annotation.data.ens <- merge(ensembl.data,uniprot.data.ens,all.x=T)
    annotation.data.ens$SYMBOL <- merge.columns(annotation.data.ens,"SYMBOL",single.val=T)
    annotation.data.ens$`SYMBOL ENSEMBL` <- NULL
    annotation.data.ens$`SYMBOL UNIPROT ENS` <- NULL
    annotation.data.ens$ALIAS <- merge.columns(annotation.data.ens,"ALIAS",single.val=T)
    annotation.data.ens$`ALIAS ENSEMBL` <- NULL
    annotation.data.ens$`ALIAS UNIPROT ENS` <- NULL
    colnames(annotation.data.ens)[colnames(annotation.data.ens)=="ENSEMBL"] <- "ENSEMBL ENSEMBL"
    colnames(annotation.data.ens)[colnames(annotation.data.ens)=="UNIPROT"] <- "UNIPROT ENSEMBL"
    
    #Merge data
    annotation.data <- merge(annotation.data.ncbi,annotation.data.ens,all=T)
    
  } else {
    #Join databases on SYMBOL column
    annotation.data <- merge(ncbi.data,ensembl.data,all=T)
    annotation.data <- merge(annotation.data,uniprot.data,all.x=T)
    annotation.data <- merge(annotation.data,speciesdb.data,all.x=T)
  }
  
  #Merge columns coming from different datasets
  annotation.data$ALIAS <- merge.columns(annotation.data,"ALIAS",single.val=F)
  if("LOCUS NCBI" %in% colnames(annotation.data)){
    annotation.data$LOCUS <- merge.columns(annotation.data,"LOCUS",single.val=F)
  }
  annotation.data$GENETYPE <- merge.columns(annotation.data,"GENETYPE",single.val=T)
  if(species!="Arabidopsis"){
    annotation.data$GENENAME <- merge.columns(annotation.data,"GENENAME",single.val=T)
  }
  annotation.data$ENTREZID <- merge.columns(annotation.data,"ENTREZID",single.val=F)
  annotation.data$ENSEMBL <- merge.columns(annotation.data,"ENSEMBL",single.val=F)
  if(species %in% c("Tomato","Cabbage","Rapeseed")){
    annotation.data <- annotation.data[is.na(annotation.data$ENSEMBL) | !duplicated(annotation.data$ENSEMBL),]
  }
  if(!is.na(speciesdb.name)){
    annotation.data[[speciesdb.name]] <- merge.columns(annotation.data,speciesdb.name,single.val=F)
  }
  annotation.data$UNIPROT <- merge.columns(annotation.data,"UNIPROT",single.val=F)
  annotation.data$ENTREZIDOLD <- merge.columns(annotation.data,"ENTREZIDOLD",single.val=F)
  annotation.data$ENSEMBLOLD <- merge.columns(annotation.data,"ENSEMBLOLD",single.val=F)
  
  #Select only joined columns
  if(!is.na(speciesdb.name)){
    annotation.data <- annotation.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZID","ENSEMBL",speciesdb.name,
                                          "UNIPROT","ENTREZIDOLD","ENSEMBLOLD")]
  } else {
    annotation.data <- annotation.data[,c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZID","ENSEMBL",
                                          "UNIPROT","ENTREZIDOLD","ENSEMBLOLD")]
  }
  
  #Fix GENETYPE categories
  annotation.data <- rename.genetypes(annotation.data)
  
  #Fill ALIAS information if missing
  annotation.data$ALIAS <- ifelse(is.na(annotation.data$ALIAS),annotation.data$SYMBOL,annotation.data$ALIAS)
  
  #Fill ENSEMBL information if gene symbol correspond to ensembl id
  if(!is.na(speciesdb.name) && speciesdb.name=="TAIR"){
    annotation.data[is.na(annotation.data$ENSEMBL) & is.na(annotation.data$ENSEMBLOLD) & 
      grepl("^AT[0-9]G",toupper(annotation.data$SYMBOL)),"ENSEMBLOLD"] <- toupper(annotation.data[is.na(annotation.data$ENSEMBL) & 
      is.na(annotation.data$ENSEMBLOLD) & grepl("^AT[0-9]G",toupper(annotation.data$SYMBOL)),"SYMBOL"])
  }
  
  #Fill locus data if ENSEMBL id is present
  if("LOCUS" %in% colnames(annotation.data) && !is.na(speciesdb.name) && speciesdb.name %in% c("SGD","TAIR")){
    annotation.data$LOCUS <- ifelse(is.na(annotation.data$LOCUS),annotation.data$ENSEMBL,annotation.data$LOCUS)
    annotation.data$LOCUS <- ifelse(is.na(annotation.data$LOCUS),annotation.data$ENSEMBLOLD,annotation.data$LOCUS)
  }
  
  #Check and remove columns that have all NA
  #annotation.data <- annotation.data[,colSums(!is.na(annotation.data))>0]
  
  #Order annotation data by symbol
  annotation.data <- annotation.data[order(annotation.data$SYMBOL),]
  annotation.data <- as.data.frame(annotation.data)
  row.names(annotation.data) <- 1:nrow(annotation.data)
  
  return(annotation.data)
}

merge.ortho.databases <- function(ncbi.orthologs,ensembl.orthologs,alliance.orthologs,speciesdb.orthologs,taxid,taxonomy.table)
{
  #Join databases on SYMBOL column
  if(nrow(ncbi.orthologs)>0){
    orthologs.data <- merge(ncbi.orthologs,ensembl.orthologs,all=T)
  } else {
    orthologs.data <- ensembl.orthologs
  }
  if(nrow(alliance.orthologs)>0){
    orthologs.data <- merge(orthologs.data,alliance.orthologs,all=T)
  }
  if(nrow(speciesdb.orthologs)>0){
    orthologs.data <- merge(orthologs.data,speciesdb.orthologs,all=T)
  }
  
  #Merge columns coming from different datasets
  for(i in 1:nrow(taxonomy.table)) {
    if(all(unlist(taxonomy.table[i,"taxid"])!=taxid)){
      orthologs.data[[paste0("ORTHO",toupper(taxonomy.table[i,"species"]))]] <- merge.columns(orthologs.data,
                                          paste0("ORTHO",toupper(taxonomy.table[i,"species"])),single.val=F)
    }
  }
  
  #Select only joined columns
  valid.pos <- sapply(taxonomy.table$taxid, function(x) all(!x %in% taxid))
  orthologs.data <- orthologs.data[,c("SYMBOL",paste0("ORTHO",toupper(taxonomy.table[valid.pos,"species"])))]
  
  return(as.data.frame(orthologs.data))
}

merge.with.ensembl.archive.data <- function(ensembl.data,ensembl.archive.data,species,speciesdb.name,
                                            ncbi.data,speciesdb.data,hcop.data,taxid)
{
  list.versions <- sort(as.numeric(names(ensembl.archive.data)),decreasing=T)
  for(i in list.versions){
    ens.data <- ensembl.archive.data[[as.character(i)]]
    #print(paste0("Merging Ensembl version ",i,"..."))
    if(nrow(ens.data)>0){
      
      #Get current list of ensembl ids (OLD and NEW)
      list.ncbi.ens.ids <- (ncbi.data %>% separate_rows(all_of(c("ENSEMBL NCBI")),sep="\\|"))$`ENSEMBL NCBI`
      list.ncbi.ens.ids <- unique(list.ncbi.ens.ids[!is.na(list.ncbi.ens.ids)])
      list.ensembl.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBL ENSEMBL")),sep="\\|"))$`ENSEMBL ENSEMBL`
      list.ensembl.ids <- unique(list.ensembl.ids[!is.na(list.ensembl.ids)])
      list.ensembl.old.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBLOLD ENSEMBL")),sep="\\|"))$`ENSEMBLOLD ENSEMBL`
      list.ensembl.old.ids <- unique(list.ensembl.old.ids[!is.na(list.ensembl.old.ids)])
      list.curr.ens.ids <- unique(c(list.ncbi.ens.ids,list.ensembl.ids,list.ensembl.old.ids))
    
      #Process archive data
      ensembl.archive <- process.ensembl.data(ens.data,speciesdb.name,ncbi.data,speciesdb.data,hcop.data,taxid,is.archive=T)
      if(species=="Arabidopsis"){
        colnames(ensembl.archive) <- c("SYMBOL","GENENAME","ALIAS ARCHIVE","GENETYPE ARCHIVE","ENSEMBL ARCHIVE")
      } else {
        colnames(ensembl.archive) <- c("SYMBOL","ALIAS ARCHIVE","GENETYPE ARCHIVE","GENENAME ARCHIVE","ENSEMBL ARCHIVE")
      }
    
      #Filter out ensembl ids that are already present in current annotations
      colnames(ensembl.archive)[colnames(ensembl.archive)=="ENSEMBL ARCHIVE"] <- "ENSEMBLOLD ARCHIVE"
      ensembl.archive.extended <- ensembl.archive %>% separate_rows(all_of("ENSEMBLOLD ARCHIVE"),sep="\\|")
      ensembl.archive.extended <- ensembl.archive.extended[!ensembl.archive.extended$`ENSEMBLOLD ARCHIVE` %in% list.curr.ens.ids,]
      if(species=="Arabidopsis"){
        ensembl.archive <- ensembl.archive.extended %>% group_by(SYMBOL,GENENAME) %>% summarise(across(everything(),function(x){
          unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
          res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
          ifelse(res=="",NA,res)
        }))
      } else {
        ensembl.archive <- ensembl.archive.extended %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
          unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
          res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
          ifelse(res=="",NA,res)
        }))
      }
    
      #Integrate with current annotation db
      if(nrow(ensembl.archive)>0){
        ensembl.data <- merge(ensembl.data,ensembl.archive,all=T)
        ensembl.data$`ALIAS ENSEMBL` <- merge.columns(ensembl.data,"ALIAS",single.val=F)
        ensembl.data$`GENETYPE ENSEMBL` <- merge.columns(ensembl.data,"GENETYPE",single.val=T)
        if(!is.na(speciesdb.name) && speciesdb.name!="TAIR"){
          ensembl.data$`GENENAME ENSEMBL` <- merge.columns(ensembl.data,"GENENAME",single.val=T)
        }
        ensembl.data$`ENSEMBLOLD ENSEMBL` <- merge.columns(ensembl.data,"ENSEMBLOLD",single.val=F)
        ensembl.data <- ensembl.data[,!endsWith(colnames(ensembl.data),"ARCHIVE")]
        print(paste0("Integrated ",nrow(ensembl.archive)," rows from Ensembl version ",i,"..."))
      }
    }
  }
  if(species %in% c("Rapeseed")){
    ensembl.data[!is.na(ensembl.data$SYMBOL) & startsWith(ensembl.data$SYMBOL,"Bna"),"SYMBOL"] <- toupper(
      ensembl.data[!is.na(ensembl.data$SYMBOL) & startsWith(ensembl.data$SYMBOL,"Bna"),"SYMBOL"])
    ensembl.data$`ALIAS ENSEMBL` <- apply(ensembl.data,1,function(row){
      if(is.na(row[["ALIAS ENSEMBL"]])){
        row[["SYMBOL"]]
      } else {
        list.ids <- unique(c(strsplit(row[["ALIAS ENSEMBL"]],"\\|")[[1]],row[["SYMBOL"]]))
        paste0(list.ids,collapse = "|")
      }
    })
  }
  return(ensembl.data)
}

merge.with.ensembl.grch37.data <- function(ensembl.data,ensembl.grch37.data,ncbi.data,speciesdb.name)
{
  #Get current list of ensembl ids (OLD and NEW)
  list.ncbi.ens.ids <- (ncbi.data %>% separate_rows(all_of(c("ENSEMBL NCBI")),sep="\\|"))$`ENSEMBL NCBI`
  list.ncbi.ens.ids <- unique(list.ncbi.ens.ids[!is.na(list.ncbi.ens.ids)])
  list.ensembl.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBL ENSEMBL")),sep="\\|"))$`ENSEMBL ENSEMBL`
  list.ensembl.ids <- unique(list.ensembl.ids[!is.na(list.ensembl.ids)])
  list.ensembl.old.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBLOLD ENSEMBL")),sep="\\|"))$`ENSEMBLOLD ENSEMBL`
  list.ensembl.old.ids <- unique(list.ensembl.old.ids[!is.na(list.ensembl.old.ids)])
  list.curr.ens.ids <- unique(c(list.ncbi.ens.ids,list.ensembl.ids,list.ensembl.old.ids))
  
  #Integrate data
  ensembl.archive.extended <- ensembl.grch37.data %>% separate_rows(all_of("ENSEMBLOLD ARCHIVE"),sep="\\|")
  ensembl.archive.extended <- ensembl.archive.extended[!ensembl.archive.extended$`ENSEMBLOLD ARCHIVE` %in% list.curr.ens.ids,]
  ensembl.archive <- ensembl.archive.extended %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
    res <- paste0(sort(unique.info[!is.na(unique.info)]),collapse="|")
    ifelse(res=="",NA,res)
  }))
  if(nrow(ensembl.archive)>0){
    ensembl.data <- merge(ensembl.data,ensembl.archive,all=T)
    ensembl.data$`ALIAS ENSEMBL` <- merge.columns(ensembl.data,"ALIAS",single.val=F)
    ensembl.data$`GENETYPE ENSEMBL` <- merge.columns(ensembl.data,"GENETYPE",single.val=T)
    if(speciesdb.name!="TAIR"){
      ensembl.data$`GENENAME ENSEMBL` <- merge.columns(ensembl.data,"GENENAME",single.val=T)
    }
    ensembl.data$`ENSEMBLOLD ENSEMBL` <- merge.columns(ensembl.data,"ENSEMBLOLD",single.val=F)
    ensembl.data <- ensembl.data[,!endsWith(colnames(ensembl.data),"ARCHIVE")]
  }
  
  return(ensembl.data)
}
