map.ortho.ids.to.symbol <- function(annotation.data.list,orthologs.data,taxid,taxonomy.table,id.type)
{
  valid.pos <- sapply(taxonomy.table$taxid, function(x) any(x %in% taxid))
  species.info <- strsplit(taxonomy.table[valid.pos,"official_name"]," ")[[1]]
  species.name <- taxonomy.table[valid.pos,"species"]
  
  #Map ref species ids to gene symbols
  map.ref.ids.to.symbol <- map.keys.to.values(annotation.data.list[[species.name]],orthologs.data$GeneID,"SYMBOL",id.type)
  orthologs.data <- merge(orthologs.data,map.ref.ids.to.symbol,by.x="GeneID",by.y=id.type,all.x=T)
  orthologs.data <- orthologs.data[!is.na(orthologs.data$SYMBOL),c("Other_tax_id","Other_GeneID","SYMBOL")]
  colnames(orthologs.data)[colnames(orthologs.data)=="SYMBOL"] <- "Ref_Symbol"
  
  #Map other species ids to gene symbols
  map.other.ids.to.symbol <- data.frame(ID=character(),SYMBOL=character())
  colnames(map.other.ids.to.symbol) <- c(id.type,"SYMBOL")
  for(i in 1:nrow(taxonomy.table)) {
    if(all(unlist(taxonomy.table[i,"taxid"])!=taxid)){
      list.gene.ids <- orthologs.data[orthologs.data$`Other_tax_id` %in% unlist(taxonomy.table[i,"taxid"]),"Other_GeneID"]
      other.species.info <- strsplit(taxonomy.table[i,"official_name"]," ")[[1]]
      other.species.name <- taxonomy.table[i,"species"]
      res.table <- map.keys.to.values(annotation.data.list[[other.species.name]],list.gene.ids,"SYMBOL",id.type)
      map.other.ids.to.symbol <- rbind(map.other.ids.to.symbol,res.table)
    }
  }
  orthologs.data <- merge(orthologs.data,map.other.ids.to.symbol,by.x="Other_GeneID",by.y=id.type,all.x=T)
  orthologs.data <- orthologs.data[,c("Other_tax_id","Ref_Symbol","SYMBOL")]
  colnames(orthologs.data) <- c("Ortho_TAXID","SYMBOL","Ortho_SYMBOL")
  return(orthologs.data)
}

double.check.symbols <- function(annot.data, check.data, pivot.col, check.data.type, is.archive)
{
  rev.mapping <- check.data[!is.na(check.data[[paste0(pivot.col," ",check.data.type)]]),c("SYMBOL",paste0(pivot.col," ",check.data.type))]
  colnames(rev.mapping) <- c(paste0("SYMBOL ",check.data.type),pivot.col)
  rev.mapping <- rev.mapping %>% separate_rows(all_of(c(pivot.col)),sep="\\|")
  rev.mapping <- rev.mapping %>% group_by_at(pivot.col) %>% summarise(across(everything(),function(x){
    res <- paste0(as.character(x),collapse="|")
    ifelse(res=="",NA,res)
  }))
  annot.data <- merge(annot.data,rev.mapping,all.x=T)
  annot.data$SYMBOL <- ifelse(is.na(annot.data[[paste0("SYMBOL ",check.data.type)]]),annot.data$SYMBOL,
                              ifelse(is.na(annot.data$SYMBOL),annot.data[[paste0("SYMBOL ",check.data.type)]],
                                     ifelse(annot.data$SYMBOL!=annot.data[[paste0("SYMBOL ",check.data.type)]] & !is.archive,annot.data[[paste0("SYMBOL ",check.data.type)]],annot.data$SYMBOL)))
  annot.data[[paste0("SYMBOL ",check.data.type)]] <- NULL
  return(annot.data)
}

map.keys.to.values <- function(annotation.data,list.keys,value.type,key.type){
  
  #Build list of keys and values to check
  key.types <- key.type
  if(key.type %in% c("ENTREZID","ENSEMBL")){
    key.types <- c(key.types,paste0(key.type,"OLD"))
  }
  if(key.type=="SYMBOL"){
    key.types <- c(key.types,"ALIAS")
  }
  value.types <- value.type
  if(value.type %in% c("ENTREZID","ENSEMBL")){
    value.types <- c(value.types,paste0(value.type,"OLD"))
  }
  
  #Run queries
  res.mapping <- data.frame(ID=list.keys)
  i <- 1
  for(kt in key.types){
    annot.data.ext <- annotation.data %>% separate_rows(all_of(kt),sep="\\|")
    annot.data.ext <- annot.data.ext[!is.na(annot.data.ext[[kt]]),]
    for(vt in value.types){
      annot.data.ext.sub <- annot.data.ext[,c(kt,vt)]
      res.mapping <- merge(res.mapping,annot.data.ext.sub,by.x="ID",by.y=kt,all.x=T)
      colnames(res.mapping)[colnames(res.mapping)==vt] <- paste0(value.type," ",LETTERS[i])
      i <- i+1
    }
  }
  colnames(res.mapping)[colnames(res.mapping)=="ID"] <- key.type
  
  #Merge value columns
  res.mapping[[value.type]] <- merge.columns(res.mapping,value.type,single.val=T)
  res.mapping <- res.mapping[,c(key.type,value.type)]
  res.mapping <- res.mapping %>% separate_rows(all_of(key.type),sep="\\|")
  res.mapping <- unique(res.mapping)
  
  return(res.mapping)
}

