options(timeout = max(3600, getOption("timeout")))

download.annot.file <- function(url,max.retries = 5, retry.delay = 10,alternative.url=NULL){
  temp.file <- tempfile(fileext = paste0(".",tools::file_ext(url)))
  attempt <- 0
  success <- FALSE
  while(!success && attempt < max.retries) {
    attempt <- attempt + 1
    tryCatch({
      download.file(url, temp.file, method = "libcurl")
      success <- TRUE
    }, error = function(e) {
      message(paste0("Error during download: ", e$message))
      Sys.sleep(10)
    }, warning = function(w) {
      # Catch "downloaded length != reported length" warning
      if (grepl("lunghezza scaricata|downloaded length|Failure when receiving", w$message)) {
        message(paste0("Incomplete download detected: ", w$message,
                       ". Retrying in ", retry.delay, "s..."))
        Sys.sleep(10)
      } else {
        warning(w)  # Propagate unrelated warnings
      }
    })
    if(!success && !is.null(alternative.url)){
      #Try alternative URL for download
      tryCatch({
        download.file(alternative.url, temp.file, method = "libcurl")
        success <- TRUE
      }, error = function(e) {
        message(paste0("Error during download: ", e$message))
        Sys.sleep(10)
      }, warning = function(w) {
        # Catch "downloaded length != reported length" warning
        if (grepl("lunghezza scaricata|downloaded length|Failure when receiving", w$message)) {
          message(paste0("Incomplete download detected: ", w$message,
                         ". Retrying in ", retry.delay, "s..."))
          Sys.sleep(10)
        } else {
          warning(w)  # Propagate unrelated warnings
        }
      })
    }
  }
  return(temp.file)
}

download.tabular.data <- function(url, header="auto",skip=0,alternative.url=NULL)
{
  file <- download.annot.file(url,alternative.url=alternative.url)
  tabular.data <- fread(file,header=header,data.table=F,showProgress=F,skip=skip)
  tabular.data[tabular.data==""] <- NA
  file.remove(file)
  return(tabular.data)
}

download.json.data <- function(url, filter.string)
{
  file <- download.annot.file(url)
  if(filter.string==""){
    json.data <- fromJSON(file,simplifyDataFrame = FALSE)
  } else {
    filtered.json <- system(paste0('jq ',shQuote(filter.string),' ',file),intern=T)
    json.data <- fromJSON(filtered.json)
    json.data[json.data==""] <- NA
  }
  file.remove(file)
  return(json.data)
}

download.gff.data <- function(url, list.tags, list.filters)
{
  file <- download.annot.file(url)
  gff.data <- as.data.frame(readGFF(file,columns=character(0),tags=list.tags,
              filter=list(type=list.filters)))
  gff.data[gff.data==""] <- NA
  file.remove(file)
  return(gff.data)
}

download.delim.data <- function(url, header, comment.character)
{
  file <- download.annot.file(url)
  delim.data <- read.delim(file, header=header, comment.char=comment.character, quote="")
  delim.data[delim.data==""] <- NA
  file.remove(file)
  return(delim.data)
}

download.go.dictionary <- function(url)
{
  file.go.dictionary <- download.annot.file(url)
  go.dictionary <- as.data.frame(get_ontology(file.go.dictionary,extract_tags="everything"))
  file.remove(file.go.dictionary)
  return(go.dictionary)
}

download.wikipathways.data <- function(url, species)
{
  wiki.content <- curl_fetch_memory(url)
  wiki.files <- read_html(wiki.content$content) %>% html_elements("a") %>% html_attr("href")
  wiki.file <- wiki.files[grep(gsub(" ","_",species),wiki.files)]
  file <- download.annot.file(paste0(url,wiki.file))
  wikipath.data <- readPathwayGMT(file)
  file.remove(file)
  return(wikipath.data)
}

filter.remote.links.html <- function(url, file.filter.string)
{
  folder.content <- curl_fetch_memory(url)
  list.links <- read_html(folder.content$content) %>% html_elements("a") %>% html_text(trim=T)
  list.links <- list.links[grep(file.filter.string,list.links)]
  return(list.links)
}

filter.remote.links.json <- function(url, file.filter.string) {
  folder.content <- curl_fetch_memory(url)
  listing <- fromJSON(rawToChar(folder.content$content))
  #Keep only entries of type "file" at root level
  is_file <- sapply(listing, function(x) is.list(x) && !is.null(x$type) && x$type == "file")
  files_root <- listing[is_file]
  #Filter by requested pattern
  list.links <- names(files_root)[grepl(file.filter.string, names(files_root))]
  return(list.links)
}

download.ensembl.archive.data <- function(archive.folder, ensembl.version, list.tags.ens, 
                                          list.filters.ens, species.scientific.name)
{
  link.content <- curl_fetch_memory(paste0(archive.folder,"/?C=S;O=D"))
  if(link.content$status_code!=200){
    archive.folder <- gsub("_gca_?[0-9]+.*$","",archive.folder)
    link.content <- curl_fetch_memory(paste0(archive.folder,"/?C=S;O=D"))
  }
  if(link.content$status_code==200){
    list.links <- read_html(link.content$content) %>% html_elements("a") %>% html_text(trim=T)
    ensembl.file <- list.links[grep("gff3",list.links)[1]]
    archive.data <- download.gff.data(paste0(archive.folder,"/",ensembl.file),list.tags.ens,
                                      list.filters.ens)
    archive.data <- archive.data[!is.na(archive.data$gene_id),]
    #if(nrow(archive.data)>0){
      return(archive.data)
    #}
  }
  #return(NULL)
}

query.ncbi.discontinued.data <- function(ncbi.archive.data)
{
  #Filter discontinued ids and symbols
  list.discontinued.ids <- ncbi.archive.data[ncbi.archive.data$GeneID=="-","Discontinued_GeneID"]
  list.discontinued.symbols <- ncbi.archive.data[ncbi.archive.data$GeneID=="-","Discontinued_Symbol"]
  #Fill information about discontinued ids
  set_entrez_key("b32dfed41d4aab5892c80e46e26d10341108")
  Sys.getenv("ENTREZ_KEY")
  num.requests <- 300
  for(i in seq(1,length(list.discontinued.ids),num.requests)){
    end.index <- min(length(list.discontinued.ids),i+num.requests-1)
    current.ids <- list.discontinued.ids[i:end.index]
    current.symbols <- list.discontinued.symbols[i:end.index]
    print(paste0("Processing NCBI discontinued ids: ",i,"-",end.index," of ",length(list.discontinued.ids)," ..."))
    while(TRUE) {
      res.xml <- tryCatch({
        entrez_fetch(db="gene",id=current.ids,rettype="xml")
      }, error = function(e) {
        if (grepl("Timeout was reached", e$message) || grepl("Connection timed out", e$message)) {
          message("Timeout, retry...")
        }
        return(NULL)
      })
      if(!is.null(res.xml)){
        batch.data <- tryCatch({
          parsed <- parse.entrez.xml(res.xml)
          # Check row count matches
          if(nrow(parsed) != length(current.ids)){
            message(paste0("Row mismatch: expected ", length(current.ids), 
                           " rows but got ", nrow(parsed), ". Retrying..."))
            NULL
          } else {
            parsed$ENTREZIDOLD <- current.ids
            parsed$SYMBOL <- current.symbols
            parsed
          }
        }, error = function(e) {
          message(paste0("Error assigning batch data: ", e$message, ". Retrying..."))
          NULL
        })
        if(!is.null(batch.data)){
          break
        }
      }
      #Sys.sleep(2)
    }
    if(i==1){
      ncbi.discontinued.data <- batch.data
    } else {
      ncbi.discontinued.data <- rbind(ncbi.discontinued.data,batch.data)
    }
  }
  ncbi.discontinued.data[ncbi.discontinued.data==""] <- NA
  #Return final data
  return(ncbi.discontinued.data)
}

parse.entrez.xml <- function(res.xml)
{
  xml.data <- read_xml(res.xml) 
  genes <- xml_find_all(xml.data, "//Entrezgene")
  gene.info <- map_df(genes, function(gene) {
    
    # 1. ID Esterni (Percorso mirato)
    db_nodes <- xml_find_all(gene, "./Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag")
    ext_ids <- map_chr(db_nodes, function(node) {
      db <- xml_text(xml_find_first(node, "Dbtag_db"))
      id <- xml_text(xml_find_first(node, ".//*[starts-with(name(), 'Object-id_')]"))
      paste0(db, ":", id)
    }) %>% unique() %>% paste(collapse = "|")
    
    # 2. Alias + Extra Terms
    synonyms <- xml_text(xml_find_all(gene, ".//Gene-ref_syn_E"))
    extra_terms <- xml_text(xml_find_all(gene, ".//Entrezgene_xtra-index-terms_E"))
    all_aliases <- unique(c(extra_terms, synonyms)) %>% .[!is.na(.)] %>% paste(collapse = "|")
    
    # 3. Estrazione Locus Tag
    locus_tag <- xml_text(xml_find_first(gene, ".//Gene-ref_locus-tag"))
    
    # 4. Risultato finale
    list(
      GENENAME = xml_text(xml_find_first(gene, ".//Gene-nomenclature_name")),
      GENETYPE     = xml_attr(xml_find_first(gene, "Entrezgene_type"), "value"),
      LOCUS     = locus_tag,
      ExtIds  = ext_ids,
      ALIAS       = all_aliases
    )
  })
  return(gene.info)
}

