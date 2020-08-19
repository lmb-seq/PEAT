#library(prodlim) # for row.match - turns out this scales poorly
library(httr)       #InterMineGetAllIdentifiers
library(jsonlite)   #InterMineGetAllIdentifiers
library(InterMineR) #getAnnot
#library(gProfileR) #getAnnot, gProfile enrichment, orthologs

#Given:
#Returns:
#Requires:


#Takes a vector of the names of things you want and an object (vector/matrix/dataframe) of the things you have. Returns only those parts of have that are in want and separately the parts that aren't
#Given: vector of things to match against (want),
#  object (vector/matrix/dataframe) of the things you have, where the Ncolumn column (if has columns) is to match using (have);
#  column to match against in have (Ncolumn), defaults to 1;
#Returns: list of the object modified to only contain those matching entries (matches) - NA if empty,
#  the discarded parts of the object (nonmatches_have) - NA if empty,
#  those elements of have not found in want (nonmatches_want) - NA if empty,
##  vector of what any given row of have refers to in matches (trans_matches),
##  vector of what any given row of have refers to in nonmatches (trans_nonmatches)
validateFilter <- function(want,have,Ncolumn=1){
  if(is.vector(have)){
    positions = match(have,want)
    matches = have[!is.na(positions)]
    nonmatches = have[is.na(positions)]
    nonmatches_want = want[is.na(match(want,have))]
    
    #trans_matches =  match(have,matches)
    #trans_nonmatches = match(have,nonmatches)
    
    if(length(matches) < 1){
      matches = NA
    }
    if(length(nonmatches) < 1){
      nonmatches = NA
    }
  }else{
    positions = match(have[,Ncolumn],want)
    matches = have[!is.na(positions),]
    nonmatches = have[is.na(positions),]
    nonmatches_want = want[is.na(match(want,have[,Ncolumn]))]

    if(dim(matches)[2] < 1){
      matches = NA
      #trans_matches = rep(NA,dim(have)[1])
    }else{
      #trans_matches = row.match(have,matches)
    }
    if(dim(nonmatches)[2] < 1){
      nonmatches = NA
      #trans_nonmatches = rep(NA,dim(have)[1])
    }else{
      #trans_nonmatches = row.match(have,nonmatches)
    }
  }
  if(length(nonmatches_want) < 1){
    nonmatches_want = NA
  }
  return(list('matches' = matches, 'nonmatches_have' = nonmatches, 'nonmatches_want' = nonmatches_want))#, 
              #'trans_matches' = trans_matches, 'trans_nonmatches' = trans_nonmatches))
}



#Function to return information about species - names, functions to use, annotations available etc.
#Given: vector of species (species)
#Returns: list of values for that species used in app
#Given species return relavant bits of text, e.g. Database name
speciesNames <- function(species){
  if (is.null(species)){
    return(NULL)
  }else{
    #Same for all species - may be added to in a per species manner, only altered in "Other" species
    built_in_list = c("brief_description_checkbox_intermine", "GO_checkbox_intermine") #can add to this list per species, must write retrieval function in annotateBuiltIn - brief description MUST be first as called specifically as first in a few places
    built_in_list_text = c("Add brief description","Add GO terms") #the text for checkboxes on the UI
    
    #Initialise
    species_specific = FALSE
    #if TRUE, must have species_specific_function that defines a lookup for built_in_list[1] 
    # at a minimum
    species_specific_per_gene_data = NULL
    
    if(species == "LIST"){
      names = list("H. sapiens" = "Human", 
                  "C. elegans" = "Worm", 
                  "M. musculus" = "Mouse", 
                  "D. melanogaster" = "Fly", 
                  "S. cerevisiae" = "BYeast",
                  "Other" = "Other")
      return(names)
    }else if (species == 'Worm'){ #some of this may be shared with other species/functionalised
      #Names
      database = "WormBase"
      mine = "WormMine"
      gene.lookup = "WormBase ID"
      gene.lookup.long = "WormBase gene identifier"
      other.lookup.long = "Sequence name"
      latin.species = "Caenorhabditis elegans"
      short.species = "celegans"
      
      #InterMine tags
      #primaryIdentifier - WormBase gene ID
      #secondaryIdentifier - Sequence name
      #geneSymbol - Gene symbol
      lookup_use = "primary"
      #WormMine
      mine_link = "http://intermine.wormbase.org/tools/wormmine/service/query/results"
      
      #Wormbase links
      #Uses primaryIdentifier for webpage addresses
      link_root_gene = "http://www.wormbase.org/species/c_elegans/gene/"
      link_root_locus = "http://www.wormbase.org/tools/genome/gbrowse/c_elegans_PRJNA13758?name="
      link_root_locus_jb = "http://www.wormbase.org/tools/genome/jbrowse-simple/?tracks=Curated_Genes%2CClassical_alleles%2CYACs_Fosmids_Cosmids&loc="
      link_tail_locus_jb = "&data=data%2Fc_elegans_PRJNA13758"
      
      #Species_specific method - species_specific_function must exist as a function like worm_lookup that 
      # takes a gene and returns data depending on input - it must be defined for built_in_list[1]
      species_specific = TRUE
      species_specific_function = "worm_lookup"
      species_specific_per_gene_data = "per_gene_wormbase"
      species_link_test = "http://api.wormbase.org/"
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 database = database, mine = mine, gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long, other.lookup.long = other.lookup.long, latin.species = latin.species, 
                 short.species = short.species, lookup_use = lookup_use, mine_link = mine_link,
                 link_root_gene = link_root_gene, link_root_locus = link_root_locus, link_root_locus_jb = link_root_locus_jb, link_tail_locus_jb = link_tail_locus_jb,
                 species_specific = species_specific, species_specific_function = species_specific_function, species_specific_per_gene_data = species_specific_per_gene_data, species_link_test = species_link_test)
    }else if(species == 'Human'){
      #Names
      database = "Ensembl"
      mine = "HumanMine"
      gene.lookup = "Ensembl gene ID"
      gene.lookup.long = "Ensembl gene identifier"
      other.lookup.long = "primaryIdentifier"
      latin.species = "Homo sapiens"
      short.species = "hsapiens"
      
      #InterMine tags
      #primaryIdentifier - Entrez ID
      #secondaryIdentifier - Ensembl gene ID
      #geneSymbol - Gene symbol
      lookup_use = "secondary"
      #HumanMine
      mine_link = "http://www.humanmine.org/humanmine/service/query/results"
      
      #Ensembl links
      #Uses secondaryIdentifier for webpage addresses
      link_root_gene = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
      link_root_locus = NULL
      link_root_locus_jb = NULL
      link_tail_locus_jb = NULL
      
      #Species_specific method - species_specific_function must exist as a function like worm_lookup that 
      # takes a gene and returns data depending on input - it must be defined for built_in_list[1]
      species_specific = TRUE
      species_specific_function = "ensembl_lookup"
      species_specific_per_gene_data = NULL
      species_link_test = "http://rest.ensembl.org/"
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 database = database, mine = mine, gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long, other.lookup.long = other.lookup.long, latin.species = latin.species, 
                 short.species = short.species, lookup_use = lookup_use, mine_link = mine_link,
                 link_root_gene = link_root_gene, link_root_locus = link_root_locus, link_root_locus_jb = link_root_locus_jb, link_tail_locus_jb = link_tail_locus_jb,
                 species_specific = species_specific, species_specific_function = species_specific_function, species_specific_per_gene_data = species_specific_per_gene_data, species_link_test = species_link_test)
    }else if(species == 'Fly'){
      #Names
      database = "FlyMine"
      mine = "FlyMine"
      gene.lookup = "Flybase gene ID"
      gene.lookup.long = "Flybase gene identifier"
      other.lookup.long = "Flybase gene identifier"
      latin.species = "Drosophila melanogaster"
      short.species = "dmelanogaster"
      
      #InterMine tags
      #primaryIdentifier - Flybase gene ID
      #secondaryIdentifier - CG/CR ID
      #geneSymbol - Gene symbol
      lookup_use = "primary"
      #HumanMine
      mine_link = "http://www.flymine.org/flymine/service/query/results"
      
      #FlyBase links
      #Uses primaryIdentifier for webpage addresses
      link_root_gene = "http://flybase.org/reports/"
      link_root_locus = NULL
      link_root_locus_jb = NULL
      link_tail_locus_jb = NULL
      
      #Species_specific method - species_specific_function must exist as a function like worm_lookup that 
      # takes a gene and returns data depending on input - it must be defined for built_in_list[1]
      species_specific = FALSE
      species_specific_function = NULL
      species_specific_per_gene_data = NULL
      species_link_test = NULL
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 database = database, mine = mine, gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long, other.lookup.long = other.lookup.long, latin.species = latin.species, 
                 short.species = short.species, lookup_use = lookup_use, mine_link = mine_link,
                 link_root_gene = link_root_gene, link_root_locus = link_root_locus, link_root_locus_jb = link_root_locus_jb, link_tail_locus_jb = link_tail_locus_jb,
                 species_specific = species_specific, species_specific_function = species_specific_function, species_specific_per_gene_data = species_specific_per_gene_data, species_link_test = species_link_test)
    }else if(species == 'Mouse'){
      #Names
      database = "MGI"
      mine = "MouseMine"
      gene.lookup = "MGI gene ID"
      gene.lookup.long = "MGI gene identifier"
      other.lookup.long = "none"
      latin.species = "Mus musculus"
      short.species = "mmusculus"
      
      #InterMine tags
      #primaryIdentifier - Mix of MGI ID and Ensembl and others! - mostly MGI
      #secondaryIdentifier - Non existant
      #geneSymbol - Gene symbol
      lookup_use = "primary"
      #HumanMine
      mine_link = "http://www.mousemine.org/mousemine/service/query/results"
      
      #MGI links
      #Uses primaryIdentifier for webpage addresses where it is MGI
      link_root_gene = "http://www.informatics.jax.org/marker/"
      link_root_locus = NULL
      link_root_locus_jb = NULL
      link_tail_locus_jb = NULL
      
      #Species_specific method - species_specific_function must exist as a function like worm_lookup that 
      # takes a gene and returns data depending on input - it must be defined for built_in_list[1]
      species_specific = FALSE
      species_specific_function = NULL
      species_specific_per_gene_data = NULL
      species_link_test = NULL
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 database = database, mine = mine, gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long, other.lookup.long = other.lookup.long, latin.species = latin.species, 
                 short.species = short.species, lookup_use = lookup_use, mine_link = mine_link,
                 link_root_gene = link_root_gene, link_root_locus = link_root_locus, link_root_locus_jb = link_root_locus_jb, link_tail_locus_jb = link_tail_locus_jb,
                 species_specific = species_specific, species_specific_function = species_specific_function, species_specific_per_gene_data = species_specific_per_gene_data, species_link_test = species_link_test)
    }else if (species == 'BYeast'){ #some of this may be shared with other species/functionalised
      #Names
      database = "Saccharomyces Genome Database"
      mine = "YeastMine"
      gene.lookup = "SGD DBID"
      gene.lookup.long = "SGD DBID"
      other.lookup.long = "Gene Systematic Name"
      latin.species = "Saccharomyces cerevisiae"
      short.species = "scerevisiae"
      
      #InterMine tags
      #primaryIdentifier - SGD DBID
      #secondaryIdentifier - Gene Systematic Name
      #geneSymbol - Gene Standard Name
      lookup_use = "primary"
      #YeastMine
      mine_link = "http://yeastmine.yeastgenome.org:8080/yeastmine/service/query/results"
      
      #Wormbase links
      #Uses primaryIdentifier for webpage addresses
      link_root_gene = "http://www.yeastgenome.org/locus/"
      link_root_locus = NULL
      link_root_locus_jb = NULL
      link_tail_locus_jb = NULL
      
      #Species_specific method - species_specific_function must exist as a function like worm_lookup that 
      # takes a gene and returns data depending on input - it must be defined for built_in_list[1]
      species_specific = FALSE
      species_specific_function = NULL
      species_specific_per_gene_data = NULL
      species_link_test = NULL
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 database = database, mine = mine, gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long, other.lookup.long = other.lookup.long, latin.species = latin.species, 
                 short.species = short.species, lookup_use = lookup_use, mine_link = mine_link,
                 link_root_gene = link_root_gene, link_root_locus = link_root_locus, link_root_locus_jb = link_root_locus_jb, link_tail_locus_jb = link_tail_locus_jb,
                 species_specific = species_specific, species_specific_function = species_specific_function, species_specific_per_gene_data = species_specific_per_gene_data, species_link_test = species_link_test)
    }else if(species == 'Other'){
      built_in_list = NULL
      built_in_list_text = NULL
      
      gene.lookup = "Gene ID"
      gene.lookup.long = "Gene ID"
      
      lookup_use = NULL
      
      
      out = list(built_in_list = built_in_list, built_in_list_text = built_in_list_text,
                 gene.lookup = gene.lookup, gene.lookup.long = gene.lookup.long,  
                 lookup_use = lookup_use,
                 species_specific = species_specific)
    }else{
      out = NULL
    }
    return(out)
  }
}




##Web related functions
#Given: vector of a web address (link);
#  vector of the text associated with the link (text)
#Returns: The html for the link
createLink <- function(link,text) {
  sprintf('<a href="%s" target="_blank" class="btn btn-primary">%s</a>',link,text)
}

#Given: vector of the text with the link (text)
#Returns: The original text used to create the link in createLink
removeLink <- function(text) {
  a = gsub(pattern='</a>',replacement="",text)
  a = gsub(pattern='<.+>',replacement="",a)
  return(a)
}

#Given: web address (addr)
#Returns: list of
#  TRUE if can connect to it, else FALSE (connection)
#  date of the time connection was tested (time)
testConnection = function(addr){
  options = "--head -m 2"
  test = system(paste0("curl ",options," ",addr),intern=T,ignore.stderr = TRUE)
  last_test = Sys.time()
  if (length(test) > 0){
    return(list(connection = TRUE, time = last_test))
  }else{
    return(list(connection = FALSE, time = last_test))
  }
}
#Version of this for running on a windows machine with curl installed - requires updating with curl location, then rename to replace the generic version
winTestConnection = function(addr){
  options = " --head -m 2"
  curl_loc = "c:/curl-7.53.1/src/curl.exe"
  test = system(paste0(curl_loc,options," ",addr),intern=T,ignore.stderr = TRUE)
  last_test = Sys.time()
  if (length(test) > 0){
    return(list(connection = TRUE, time = last_test))
  }else{
    return(list(connection = FALSE, time = last_test))
  }
}


##



##Get annotations
#Get annotations from file or from Web
#Given: vector of what annotations want (type), possibles below;
# vector of species name (species), as defined in speciesNames();
# vector of genes to get annotation for, defaults to NA in which case can only use 'bulk' annotation sources;
# boolean TRUE where only want to retrieve local copy of annotations (offline);
# (local copy of annotations if connection_error = TRUE or offline = TRUE);
#  vector of where to get annotations from (source), defaults to InterMine, possibles below;
#  boolean of whether to save local copy (save), defaults to TRUE
#Returns: NULL if offline and no local version of annotation, else list of
#  matrix of annotations (table);
#  boolean of whether was a problem connecting to source (connection_error);
#  date of the time connection was tested (time);
#  (local copy of annotation if current one older than n days);
#  filename of the annotation used (file_ref)
#Requires: speciesNames(), testConnection()

#Possible types
#trans_table
#brief_description_checkbox_intermine
#GO_checkbox_intermine
#GO_list_intermine - stores pre GO term info, not per gene like GO_checkbox_intermine

##Possible sources
#intermine    'bulk'    currently only 'bulk' annotation source is intermine edit code when more than one
#species_specific specified in speciesNames() - only works for brief_description_checkbox_intermine, only takes a single gene

getAnnot = function(type,species,genes=NA,offline=FALSE,source="intermine",save=TRUE){

   if(as.logical(sum(is.na(genes)))){
    #use stored data if data newer than n_days if looking up 'bulk' data
    today = Sys.Date()
    filename = paste0(species,"/",species,"_",type,"_",today,".tsv")   
    n_days = 14
    v_days = seq.Date(today - (n_days - 1),today,1)
    filenames = paste0("data/",species,"/",species,"_",type,"_",v_days,".tsv")
    # if (sum(file.exists(filenames)) > 0){
      offline = TRUE
      file_ref = filename
      offline2 = TRUE
    # }else{
    #   offline2 = FALSE
    # }
   }else{
     offline2 = FALSE
   }
  
  if (offline){
    connection_error = TRUE
    test_time = Sys.Date()
  }else{
    species_data = speciesNames(species)
    
    if(as.logical(sum(is.na(genes)))){ #currently only 'bulk' annotation source is intermine
      source="intermine"
    }
    
    if (source=="intermine"){
      latin.species = species_data$latin.species
      mine = species_data$mine
      lookup_use = species_data$lookup_use
      
      if(is.na(genes)){
        constraint = setConstraints(
          paths = "Gene.organism.name",
          operators = "=",
          values = list(latin.species)
        )
      }else{
        constraint = setConstraints(
          paths = c("Gene.organism.name",
                    paste0("Gene.",lookup_use,"Identifier")),
          operators = rep("=", 2),
          values = list(latin.species,genes)
        )
      }
      if (type == "trans_table"){
        query = setQuery(
          select = c("Gene.primaryIdentifier",
                     "Gene.secondaryIdentifier",
                     "Gene.symbol"),
          #orderBy = list(c(Gene.primaryIdentifier = "ASC")),
          where = constraint
        )
      }else if (type == "brief_description_checkbox_intermine"){
        query = setQuery(
          select = c(paste0("Gene.",lookup_use,"Identifier"),
                     "Gene.briefDescription"),
          #orderBy = list(c(Gene.briefDescription = "ASC")),
          where = constraint
        )
        if(species == "Mouse"){
          query = setQuery(
            select = c(paste0("Gene.",lookup_use,"Identifier"),
                       "Gene.description"),
            #orderBy = list(c(Gene.briefDescription = "ASC")),
            where = constraint
          )
        }
      }else if (type == "GO_checkbox_intermine"){
        query = setQuery(
          select = c(paste0("Gene.",lookup_use,"Identifier"),
                     "Gene.goAnnotation.ontologyTerm.identifier",
                     "Gene.goAnnotation.ontologyTerm.name"),
          #orderBy = list(c(Gene.briefDescription = "ASC")),
          where = constraint
        )
      }else if (type == "GO_list_intermine"){
        query = setQuery(
          select = c(paste0("Gene.",lookup_use,"Identifier"),
                     "Gene.goAnnotation.ontologyTerm.identifier",
                     "Gene.goAnnotation.ontologyTerm.name"),
          #orderBy = list(c(Gene.briefDescription = "ASC")),
          where = constraint
        )
      }
      
      #connection error
      im <- initInterMine(mine=listMines()[mine])
      test = testConnection(im$mine)
      test_time = test$time
      connection_error = !(test$connection)
      
      if (!connection_error){
        #get data
        data <- runQuery(im, query)
        
        #data processing
        if(is.null(data)){
          return(NULL)
        }
        
        data = as.matrix(data)
        if (type == "trans_table"){
          colnames(data) = c("primary","secondary","symbol")
        }else if (type == "brief_description_checkbox_intermine"){
          colnames(data) = c('gene_id','brief.description')
        }else if (type == "GO_checkbox_intermine"){
          #Join GO per gene
          combined_GO_ID = unlist(lapply((split(data[,2],factor(data[,1]))),paste,collapse=", "))
          combined_GO_names = unlist(lapply((split(data[,3],factor(data[,1]))),paste,collapse=", "))
          combined_GO = cbind(levels(factor(data[,1])),combined_GO_ID,combined_GO_names)    
          
          ##return for given genes
          #data = combined_GO[(match(genes,combined_GO[,1])),]
          data = combined_GO
          if(!is.null(dim(data))){
            if(dim(data)[1] > 0){
              colnames(data) = c('gene_id','GO.ID','GO.names')
            }
          }else{
            data = c(NA,NA,NA)
          }
          
        }else if (type == "GO_list_intermine"){
          # #Join GO per gene
          # combined_GO_ID = unlist(lapply((split(data[,2],factor(data[,1]))),paste,collapse=", "))
          # combined_GO_names = unlist(lapply((split(data[,3],factor(data[,1]))),paste,collapse=", "))
          # combined_GO = cbind(levels(factor(data[,1])),combined_GO_ID,combined_GO_names)    
          # 
          # ##return for given genes
          # #data = combined_GO[(match(genes,combined_GO[,1])),]
          # data = combined_GO
          if(!is.null(dim(data))){
            if(dim(data)[1] > 0){
              colnames(data) = c('gene_id','GO.ID','GO.name')
            }
          }else{
            data = c(NA,NA,NA)
          }
          
        }
      }
    }else{ #add sources here
      if((species_data$species_specific)&(length(genes) == 1)){
        ss_function = match.fun(species_data$species_specific_function)
        species_link_test = species_data$species_link_test
        possible_names = species_data$built_in_list
        
        test_data = testConnection(species_link_test)
        test = test_data$connection
        time = test_data$time
        file_ref = NULL
        
        if (test){
          if (type == possible_names[1]){
            ss_data = ss_function(genes,possible_names[1],species)
            table = ss_data$table
            connection_error = ss_data$connection_error
            last_test = ss_data$time
            
            #fall back to intermine if connection_error - to not overwrite ss_connection_error 
            #we ignore im connection_error here
            if (connection_error){
              ss_data = getAnnot(type,species,genes)
              table = ss_data$table
              #connection_error = ss_data$connection_error
              #last_test = ss_data$time
            }
            
            #return(list(table = table, connection_error = connection_error, time = last_test))
            return(list(table = table, connection_error = connection_error, time = last_test, file_ref = file_ref))
          }else{
            return(NULL)
          }
        }else{
          table = NULL     
          connection_error = TRUE  
          
          #fall back to intermine if connection_error - to not overwrite ss_connection_error we ignore im connection_error here - or move this out of function
          ss_data = getAnnot(type,species,genes)
          if (is.null(ss_data)){ 
            return(NULL)
          }else{
            table = ss_data$table
          }
          
          #return(list(table = table, connection_error = connection_error, time = time))
          return(list(table = table, connection_error = connection_error, time = time, file_ref = file_ref))
        }
      }else{
        return(NULL)
      }
    }
  }
  
  #for local backup/load of bulk data
  today = Sys.Date()
  filename = paste0(species,"/",species,"_",type,"_",today,".tsv")   
  if(connection_error){
    #load trans_table from file
    possibles = list.files(paste0("data/",species,"/"),pattern=paste0(species,"_",type,"_"))
    if(length(possibles) > 0){
      #get most recent - this will be last due to numeric sorting
      filename = paste0(species,"/",possibles[length(possibles)])
      
      new_table = read.delim(paste0("data/",filename),stringsAsFactors=FALSE) #ideally use as.is per type #changed from read.table, check doesn't break
      data = as.matrix(new_table)
      data[,1] = new_table[,1]
      
      # pat <- ".+_([0-9]+-[0-9]+-[0-9]+).tsv"
      # file_ref = sub(pat, "\\1", filename[grepl(pat, filename)])
      file_ref = filename
    }else{
      data = NULL
      new_column = NULL
      lookup = NULL
      rows_not = NULL
      
      return(NULL)
    }
  }else{
    #backup data if backup older than n_days
    n_days = 28
    v_days = seq.Date(today - (n_days - 1),today,1)
    filenames = paste0("data/",species,"/",species,"_",type,"_",v_days,".tsv")
    if (sum(file.exists(filenames)) < 1){
      if (save){
        write.table(data,paste0("data/",filename),sep="\t")       
      }
    }
    
    file_ref = filename #today
  }
  
  if(offline2){
    connection_error = FALSE
  }
  
  return(list(table = data, connection_error = connection_error, time = test_time, file_ref = file_ref))
}



#ENSEMBL
#Given: vector of gene (gene);
#  vector of what you want to lookup (type)    
#Returns: list of data from ensembl
request_ensembl <- function(gene,type) {
  #bit to check what goes with type
  search_by = "id"
  
  fromJSON(getURL(paste0("http://rest.ensembl.org/",type,"/",search_by,"/",gene), httpheader=list(`Content-Type` = "application/json")))
} 

#Given: vector of gene (gene);
#  vector of what you want to lookup (method);
#  vector of which species (species)
#Returns: NULL if method not recognized, else list of
# matrix of lookup gene identifier and annotation (table);
# boolean TRUE if connection fails - then table = NULL, else FALSE (connection_error);
# date of the time connection was tested (time)
#Requires: request_ensembl()
ensembl_lookup <- function(gene,method,species){
  species_data = speciesNames(species)
  possible_names = species_data$built_in_list
  
  last_test = Sys.time()
  
  #try to get data from Ensembl
  if (method == possible_names[1]){
    type = "lookup"
    response=request_ensembl(gene,type)
    if(is.null(response)) { #Ensembl offline - might not actually look like this, check
      table = NULL
      connection_error = TRUE  
    }else{ 
      if (is.null(response$description)){ #occurs when no description exists
        table = "No description avaliable" #- gene may no longer exist - this isn't true for the Human genes - it appears description is not a required field
      }else{
        table = response$description
      }
      connection_error = FALSE  
    }       
    return(list(table = table, connection_error = connection_error, time = last_test))
  }else{
    return(NULL)
  }
} 




#WORM
#Given: vector of webaddress (url)
#Returns: 
#WormBase API
request_wormbase <- function(url) {
  lookup = (getURL(url, httpheader=list(`Content-Type` = "application/json")))
  if(length(grep("<!DOCTYPE",lookup)) < 1){ #handle page not found
    fromJSON(lookup)
  }else{
    return(NULL)
  }
}    

#Given: vector of gene(s) (gene);
#  vector of what you want to lookup (method);
#  vector of which species (species) - ignored only included for cross-species compatibility
#Returns: NULL if method not recognized, else list of
# matrix of lookup gene identifier and annotation (table);
# boolean TRUE if connection fails - then table = NULL, else FALSE (connection_error);
# date of the time connection was tested (time)
#Requires: request_wormbase(), speciesNames()
worm_lookup <- function(gene,method,species){
  species_data = speciesNames("Worm")
  possible_names = species_data$built_in_list
  
  last_test = Sys.time()
  
  #try to get data from WormBase
  if (method == possible_names[1]){
    response=request_wormbase(paste0('http://api.wormbase.org/rest/field/gene/',gene,'/concise_description'))
    if(is.null(response)) { #WormBase offline - might not actually look like this, check
      table = NULL
      connection_error = TRUE  
    }else{ 
      if (is.null(response$concise_description$data$text[1])){ #occurs when no description exisst
        table = "No description avaliable - gene may no longer exist"
      }else{
        table = response$concise_description$data$text[1]  
      }
      connection_error = FALSE  
    }       
    return(list(table = table, connection_error = connection_error, time = last_test))
  }else if (method == "wormbase_blastp_table"){
    response=request_wormbase(paste0("http://api.wormbase.org/rest/field/gene/",gene,"/best_blastp_matches"))
    if(is.null(response)) { #WormBase offline - might not actually look like this, check
      data3 = NULL
      connection_error = TRUE  
    }else{ 
      data = response$best_blastp_matches$data$hits
      connection_error = FALSE
      
      if (is.null(data)){ #this table doesn't exist for this gene
        data3 = NULL
      }else{
        #data2 = as.data.frame(data)
        data2 = matrix(unlist(data),nrow = dim(data)[1])
        #merge species/genus
        data2[,3] = paste0(data2[,3],'.',data2[,4])
        data2 = data2[,-4]
        #select columns
        data3 = cbind(data2[,1],data2[,3],data2[,5],data2[,7],data2[,2])
        colnames(data3) = c("BLAST e-value","Species","Hit","Description","% Length")
      }
    }
    return(list(table = data3, connection_error = connection_error, time = last_test))
  }else if (method == "wormbase_ortholog_table"){
    response=request_wormbase(paste0("http://api.wormbase.org/rest/field/gene/",gene,"/nematode_orthologs"))
    if(is.null(response)) { #WormBase offline - might not actually look like this, check
      new_table2 = NULL
      connection_error = TRUE  
    }else{ 
      connection_error = FALSE
      new_table = response$nematode_orthologs$data
      if (is.null(new_table)){ #this table doesn't exist for this gene
        new_table2 = "No data available"
      }else{
        new_table2 = cbind(new_table$species$genus,new_table$species$species,new_table$ortholog$label)
        new_table2[,1] = paste0(new_table2[,1],'.',new_table2[,2])
        new_table2 = new_table2[,-2]
        new_table2 = cbind(sort(new_table2[,1]),new_table2[order(new_table2[,1]),2])
        
        colnames(new_table2) = c("Species","Ortholog")
      }
    }
    return(list(table = new_table2, connection_error = connection_error, time = last_test))
  }else{
    return(NULL)
  }
}



##Analysis
#Given: vector of gene_names that are expected lookup in speciesNames(),
#  vector of species from speciesNames,
#  boolean of whether to use Benjamini & Hochberg multiple testing correction (correction)
#Returns: matrix of per GO term enrichment
#Requires: getAnnot()
go_enrichment <- function(gene_names,species,correction=TRUE){
  #my version
  data3 = getAnnot("GO_list_intermine",species)
  go_table = data3$table
  d4 = subset(go_table, go_table[,1] %in% gene_names)
  
  go_term_counts_gene_names = (summary(factor(d4[,2]),maxsum=1000000))
  go_term_counts = (summary(factor(go_table[,2]),maxsum=1000000))
  
  all_gene_names_n = sum(!is.na(match(gene_names,go_table[,1])))
  all_genes_n = length(summary(factor(go_table[,1]),maxsum=1000000))
  
  
  table = matrix(NA,nrow = length(go_term_counts_gene_names), ncol = 5)
  colnames(table) = c("GO ID","description","p-value","count in selected genes","count in all genes")
  for (i in 1:length(go_term_counts_gene_names)){
    prob = phyper((go_term_counts_gene_names[i] - 1),
                  go_term_counts[(match(names(go_term_counts_gene_names[i]),names(go_term_counts)))],
                  (all_genes_n - go_term_counts[(match(names(go_term_counts_gene_names[i]),names(go_term_counts)))]),
                  all_gene_names_n,lower.tail=F)
    table[i,] = c(names(go_term_counts_gene_names[i]),go_table[match(names(go_term_counts_gene_names[i]),go_table[,2]),3],
                  prob,go_term_counts_gene_names[i],go_term_counts[(match(names(go_term_counts_gene_names[i]),names(go_term_counts)))])
  }
  table = table[(order(table[,3])),]
  
  if(correction){
    #Benjamini Hochberg correction
    new_p = p.adjust(table[,3],method="BH")
    table[,3] = new_p
    table = table[(order(table[,3])),]
  }
  
  return(table)
}










##Depreciated

#Get translation table matrix of InterMine primaryID, secondaryID, gene symbol for all genes from file or from Web
#Given: vector of species name (species), as defined in speciesNames();
# boolean TRUE where only want to retrieve local copy of transtable (offline);
# (local copy of transtable if InterMineGetAllIdentifiers()$connection_error = TRUE or offline = TRUE)
#Returns: NULL if offline and no local version of trans_table, else list of
#  matrix from InterMine of primaryID, secondaryID, gene symbol (table);
#  boolean of whether was a problem connecting in InterMineGetAllIdentifiers() (connection_error);
#  date of the time connection was tested (time);
#  (local copy of transtable if current one older than n days);
#  filename of the translation_table used (file_ref)
#Requires: InterMineGetAllIdentifiers()
retrieveTranslation = function(species,offline=FALSE){
  if (!offline){
    #Retrieve translation matrix
    data = InterMineGetAllIdentifiers(species)
    trans_table = data$table
    test_time = data$time
    connection_error = data$connection_error
  }else{
    connection_error = TRUE
    test_time = Sys.Date()
  }
  
  #for local backup/load
  today = Sys.Date()
  filename = paste0(species,"_trans_table_",today,".tsv")   
  
  if(connection_error){
    #load trans_table from file
    possibles = list.files("data/",pattern=paste0(species,"_trans_table_"))
    if(length(possibles) > 0){
      #get most recent - this will be last due to numeric sorting
      filename = possibles[length(possibles)]
      
      new_table = read.table(paste0("data/",filename),as.is=c(1,2,3))
      trans_table = as.matrix(new_table)
      trans_table[,1] = new_table[,1]
      
      # pat <- ".+_([0-9]+-[0-9]+-[0-9]+).tsv"
      # file_ref = sub(pat, "\\1", filename[grepl(pat, filename)])
      file_ref = filename
    }else{
      table = NULL
      new_column = NULL
      lookup = NULL
      rows_not = NULL
      
      return(NULL)
    }
  }else{
    #backup trans_table if backup older than n_days
    n_days = 7
    v_days = seq.Date(today - (n_days - 1),today,1)
    filenames = paste0("data/",species,"_trans_table_",v_days,".tsv")
    if (sum(file.exists(filenames)) < 1){
      write.table(trans_table,paste0("data/",filename),sep="\t")       
    }
    
    file_ref = filename #today
  }
  return(list(table = trans_table, connection_error = connection_error, time = test_time, file_ref = file_ref))
}


#InterMine
#Name conversions
#Given: vector of species name (species)
#Returns: NULL is species not in speciesNames(), else list of
# matrix of InterMine primaryIdentifier, secondaryIdentifier, gene symbol (table);
# boolean TRUE if connection fails - then table = NULL, else FALSE (connection_error);
# date of the time connection was tested (time)
#Requires: testConnection(), speciesNames() 
InterMineGetAllIdentifiers = function(species){
  species_data = speciesNames(species)
  mine_link = species_data$mine_link
  latin.species = species_data$latin.species
  if (is.null(species_data)){
    return(NULL)
  }else{
    test_data = testConnection(mine_link)
    test = test_data$connection
    time = test_data$time
    if (test){
      query=paste0('<query model="genomic" view="Gene.primaryIdentifier Gene.secondaryIdentifier Gene.symbol" ><constraint path="Gene.organism" op="LOOKUP" value="',latin.species,'" /></query>')
      ret=POST(mine_link,
               body=list(query=query, format='json'),
               encode='form')
      response=fromJSON(content(ret,as='text'))
      if(response$wasSuccessful) {
        if (length(response$results) > 0){
          table = response$results
          connection_error = FALSE
          colnames(table) = c("primary","secondary","symbol")
        }else{
          table = NULL
          connection_error = TRUE   
        }
      }
    }else{
      table = NULL
      connection_error = TRUE          
    }
    return(list(table = table, connection_error = connection_error, time = time))
  }
}  