merge.ncbi.data <- function(ncbi.data,ncbi.replaced.data,ncbi.discontinued.data,speciesdb.name)
{
  #Merge with replaced data
  ncbi.data <- merge(ncbi.data,ncbi.replaced.data,all.x=T,by.x="ENTREZID NCBI",by.y="ENTREZID")
  ncbi.data$`ENTREZIDOLD NCBI` <- merge.single.column(ncbi.data,"ENTREZIDOLD",single.val=F)
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
  ncbi.data$`ALIAS NCBI` <- merge.single.column(ncbi.data,"ALIAS",single.val=F)
  if("LOCUS NCBI" %in% colnames(ncbi.data)){
    ncbi.data$`LOCUS NCBI` <- merge.single.column(ncbi.data,"LOCUS",single.val=F)
  }
  ncbi.data$`GENETYPE NCBI` <- merge.single.column(ncbi.data,"GENETYPE",single.val=T)
  ncbi.data$`GENENAME NCBI` <- merge.single.column(ncbi.data,"GENENAME",single.val=T)
  ncbi.data$`ENSEMBLOLD NCBI` <- merge.single.column(ncbi.data,"ENSEMBLOLD",single.val=F)
  if(!is.na(speciesdb.name)){
    ncbi.data[[paste0(speciesdb.name," NCBI")]] <- merge.single.column(ncbi.data,speciesdb.name,single.val=F)
  }
  ncbi.data$`ENTREZIDOLD NCBI` <- merge.single.column(ncbi.data,"ENTREZIDOLD",single.val=F)
  ncbi.data <- ncbi.data[,!endsWith(colnames(ncbi.data),"ARCHIVE")]
  
  return(ncbi.data)
}

merge.single.column <- function(annotation.data,column,single.val)
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

merge.columns <- function(annotation.data){
  #Merge columns coming from different datasets
  annotation.data$ALIAS <- merge.single.column(annotation.data,"ALIAS",single.val=F)
  if(any(startsWith(colnames(annotation.data),"LOCUS"))){
    annotation.data$LOCUS <- merge.single.column(annotation.data,"LOCUS",single.val=F)
  }
  annotation.data$GENETYPE <- merge.single.column(annotation.data,"GENETYPE",single.val=T)
  #if(species!="Arabidopsis"){
    annotation.data$GENENAME <- merge.single.column(annotation.data,"GENENAME",single.val=T)
  #}
  annotation.data$ENTREZID <- merge.single.column(annotation.data,"ENTREZID",single.val=F)
  annotation.data$ENSEMBL <- merge.single.column(annotation.data,"ENSEMBL",single.val=F)
  #if(species %in% c("Tomato","Cabbage","Rapeseed")){
    #annotation.data <- annotation.data[is.na(annotation.data$ENSEMBL) | !duplicated(annotation.data$ENSEMBL),]
  #}
  if(!is.na(speciesdb.name)){
    annotation.data[[speciesdb.name]] <- merge.single.column(annotation.data,speciesdb.name,single.val=F)
  }
  if(any(startsWith(colnames(annotation.data),"UNIPROT"))){
    annotation.data$UNIPROT <- merge.single.column(annotation.data,"UNIPROT",single.val=F)
  }
  annotation.data$ENTREZIDOLD <- merge.single.column(annotation.data,"ENTREZIDOLD",single.val=F)
  annotation.data$ENSEMBLOLD <- merge.single.column(annotation.data,"ENSEMBLOLD",single.val=F)
  #Select only joined columns
  columns.to.select <- c("SYMBOL","ALIAS","LOCUS","GENETYPE","GENENAME","ENTREZID","ENSEMBL",
                         "ENTREZIDOLD","ENSEMBLOLD")
  if("UNIPROT" %in% colnames(annotation.data)){
    columns.to.select <- c(columns.to.select,"UNIPROT")
  }
  if(!is.na(speciesdb.name)){
    columns.to.select <- c(columns.to.select,speciesdb.name)
  }
  annotation.data <- annotation.data[,columns.to.select]
  return(annotation.data)
}

merge.databases <- function(ncbi.data,ensembl.data,uniprot.data,speciesdb.data,speciesdb.name,species)
{
  #Join NCBI and Ensembl data on SYMBOL column
  if(nrow(ensembl.data)>0){
    annotation.data <- merge(ncbi.data,ensembl.data,all=T)
  } else {
    annotation.data <- ncbi.data
  }
  annotation.data <- merge.columns(annotation.data)
  colnames(annotation.data)[colnames(annotation.data)!="SYMBOL"] <- paste0(colnames(annotation.data)[colnames(annotation.data)!="SYMBOL"]," NCBIENS")
  
  #Join current annotation data with Uniprot data
  map.ncbi.ids.to.symbol <- map.keys.to.values(annotation.data,uniprot.data$ENTREZID,"SYMBOL","ENTREZID NCBIENS")
  colnames(map.ncbi.ids.to.symbol) <- c("ENTREZID","SYMBOL UNIPROT NCBI")
  uniprot.data <- merge(uniprot.data,map.ncbi.ids.to.symbol,all.x=T)
  map.ens.ids.to.symbol <- map.keys.to.values(annotation.data,uniprot.data$ENSEMBL,"SYMBOL","ENSEMBL NCBIENS")
  colnames(map.ens.ids.to.symbol) <- c("ENSEMBL","SYMBOL UNIPROT ENS")
  uniprot.data <- merge(uniprot.data,map.ens.ids.to.symbol,all.x=T)
  uniprot.data$SYMBOL <- ifelse(is.na(uniprot.data$`SYMBOL UNIPROT NCBI`),
    ifelse(is.na(uniprot.data$`SYMBOL UNIPROT ENS`),uniprot.data$`SYMBOL UNIPROT`,
           uniprot.data$`SYMBOL UNIPROT ENS`),uniprot.data$`SYMBOL UNIPROT NCBI`)
  ref.cols.alias <- c("SYMBOL UNIPROT","SYMBOL UNIPROT NCBI","SYMBOL UNIPROT ENS")
  uniprot.data$ALIAS <- apply(uniprot.data,1,function(row){
    unique.info <- unname(unlist(sapply(ref.cols.alias,function(col){
      strsplit(as.character(row[col]),"\\|")
    })))
    unique.info <- unique(unique.info[!is.na(unique.info)])
    res <- paste0(unique.info,collapse="|")
    ifelse(res=="",NA,res)
  })
  uniprot.data$`SYMBOL UNIPROT` <- NULL
  uniprot.data$`SYMBOL UNIPROT NCBI` <- NULL
  uniprot.data$`SYMBOL UNIPROT ENS` <- NULL
  uniprot.data <- uniprot.data[!is.na(uniprot.data$SYMBOL),]
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
  if(species=="ZFIN"){
    uniprot.data[grepl("^loc[0-9]+",uniprot.data$SYMBOL),"SYMBOL"] <- toupper(uniprot.data[grepl("^loc[0-9]+",uniprot.data$SYMBOL),"SYMBOL"])
  }
  uniprot.data <- uniprot.data %>% group_by(SYMBOL) %>% summarise(across(everything(),function(x){
    unique.info <- unique(unlist(strsplit(as.character(x),"\\|")))
    res <- paste0(unique.info[!is.na(unique.info)],collapse="|")
    ifelse(res=="",NA,res)
  }))
  colnames(uniprot.data)[colnames(uniprot.data)!="SYMBOL"] <- paste0(colnames(uniprot.data)[colnames(uniprot.data)!="SYMBOL"]," UNIPROT")
  annotation.data <- merge(annotation.data,uniprot.data,all.x=T)
  annotation.data$`ENTREZID NCBIENSUNI` <- ifelse(is.na(annotation.data$`ENTREZID NCBIENS`),
    annotation.data$`ENTREZID UNIPROT`,annotation.data$`ENTREZID NCBIENS`)
  annotation.data$`ENTREZID NCBIENS` <- NULL
  annotation.data$`ENTREZID UNIPROT` <- NULL
  annotation.data$`ENSEMBL NCBIENSUNI` <- ifelse(is.na(annotation.data$`ENSEMBL NCBIENS`),
    annotation.data$`ENSEMBL UNIPROT`,annotation.data$`ENSEMBL NCBIENS`)
  annotation.data$`ENSEMBL NCBIENS` <- NULL
  annotation.data$`ENSEMBL UNIPROT` <- NULL
  annotation.data$`ALIAS NCBIENSUNI` <- merge.single.column(annotation.data,"ALIAS",single.val=F)
  annotation.data$`ALIAS NCBIENS` <- NULL
  annotation.data$`ALIAS UNIPROT` <- NULL
  
  #Join current data with species-specific DB data
  if(nrow(speciesdb.data)>0){
    annotation.data <- merge(annotation.data,speciesdb.data,all.x=T)
  }
  annotation.data <- merge.columns(annotation.data)
  
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

merge.ortho.databases <- function(ncbi.orthologs,ensembl.orthologs,alliance.orthologs,speciesdb.orthologs,species.taxid,taxonomy.table)
{
  #Join databases on SYMBOL column
  if(nrow(ncbi.orthologs)>0){
    orthologs.data <- merge(ncbi.orthologs,ensembl.orthologs,all=T)
  } else {
    orthologs.data <- ensembl.orthologs
  }
  if(nrow(orthologs.data)>0){
    if(nrow(alliance.orthologs)>0){
      orthologs.data <- merge(orthologs.data,alliance.orthologs,all=T)
    }
  } else {
    orthologs.data <- alliance.orthologs
  }
  if(nrow(orthologs.data)>0){
    if(nrow(speciesdb.orthologs)>0){
      orthologs.data <- merge(orthologs.data,speciesdb.orthologs,all=T)
    }
  } else {
    orthologs.data <- speciesdb.orthologs
  }
  if(nrow(orthologs.data)>0){
    #Merge columns coming from different datasets
    for(i in 1:nrow(taxonomy.table)) {
      if(all(unlist(taxonomy.table[i,"taxid"])!=species.taxid)){
        orthologs.data[[paste0("ORTHO",toupper(taxonomy.table[i,"species"]))]] <- merge.single.column(orthologs.data,
                                          paste0("ORTHO",toupper(taxonomy.table[i,"species"])),single.val=F)
      }
    }
    #Select only joined columns
    valid.pos <- sapply(taxonomy.table$taxid, function(x) all(!x %in% species.taxid))
    orthologs.data <- orthologs.data[,c("SYMBOL",paste0("ORTHO",toupper(taxonomy.table[valid.pos,"species"])))]
  }
  return(as.data.frame(orthologs.data))
}

merge.with.ensembl.archive.data <- function(ensembl.data,ensembl.archive.data,species,speciesdb.name,
                                            ncbi.data,hcop.data,species.taxid)
{
  list.versions <- sort(as.numeric(names(ensembl.archive.data)),decreasing=T)
  for(i in list.versions){
    ens.data <- ensembl.archive.data[[as.character(i)]]
    ens.data <- ens.data[!is.na(ens.data$gene_id),]
    print(paste0("Merging Ensembl version ",i,"..."))
    if(nrow(ens.data)>0){
      
      #Get current list of ensembl ids (OLD and NEW)
      list.ncbi.ens.ids <- (ncbi.data %>% separate_rows(all_of(c("ENSEMBL NCBI")),sep="\\|"))$`ENSEMBL NCBI`
      list.ncbi.ens.ids <- unique(list.ncbi.ens.ids[!is.na(list.ncbi.ens.ids)])
      list.ensembl.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBL ENSEMBL")),sep="\\|"))$`ENSEMBL ENSEMBL`
      list.ensembl.ids <- unique(list.ensembl.ids[!is.na(list.ensembl.ids)])
      list.ensembl.old.ids <- (ensembl.data %>% separate_rows(all_of(c("ENSEMBLOLD ENSEMBL")),sep="\\|"))$`ENSEMBLOLD ENSEMBL`
      list.ensembl.old.ids <- unique(list.ensembl.old.ids[!is.na(list.ensembl.old.ids)])
      list.curr.ens.ids <- unique(c(list.ncbi.ens.ids,list.ensembl.ids,list.ensembl.old.ids))
      #Filter out ensembl ids that are already present in current annotations
      ens.data <- ens.data[!ens.data$gene_id %in% list.curr.ens.ids,]
      #Filter out RNAs
      ens.data[startsWith(ens.data$gene_id,"ENSRNA"),"gene_id"] <- NA
      ens.data <- ens.data[!is.na(ens.data$gene_id),]
      
      if(nrow(ens.data)>0){
        #Process archive data
        ensembl.archive <- process.ensembl.data(ens.data,speciesdb.name,ncbi.data,hcop.data,species.taxid,is.archive=T)
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
          ensembl.data$`ALIAS ENSEMBL` <- merge.single.column(ensembl.data,"ALIAS",single.val=F)
          ensembl.data$`GENETYPE ENSEMBL` <- merge.single.column(ensembl.data,"GENETYPE",single.val=T)
          if(!is.na(speciesdb.name) && speciesdb.name!="TAIR"){
            ensembl.data$`GENENAME ENSEMBL` <- merge.single.column(ensembl.data,"GENENAME",single.val=T)
          }
          ensembl.data$`ENSEMBLOLD ENSEMBL` <- merge.single.column(ensembl.data,"ENSEMBLOLD",single.val=F)
          ensembl.data <- ensembl.data[,!endsWith(colnames(ensembl.data),"ARCHIVE")]
          print(paste0("Integrated ",nrow(ensembl.archive)," rows from Ensembl version ",i,"..."))
        }
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
    ensembl.data$`ALIAS ENSEMBL` <- merge.single.column(ensembl.data,"ALIAS",single.val=F)
    ensembl.data$`GENETYPE ENSEMBL` <- merge.single.column(ensembl.data,"GENETYPE",single.val=T)
    if(speciesdb.name!="TAIR"){
      ensembl.data$`GENENAME ENSEMBL` <- merge.single.column(ensembl.data,"GENENAME",single.val=T)
    }
    ensembl.data$`ENSEMBLOLD ENSEMBL` <- merge.single.column(ensembl.data,"ENSEMBLOLD",single.val=F)
    ensembl.data <- ensembl.data[,!endsWith(colnames(ensembl.data),"ARCHIVE")]
  }
  
  return(ensembl.data)
}
