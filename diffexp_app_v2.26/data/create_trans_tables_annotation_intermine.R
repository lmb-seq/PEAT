library(InterMineR) 

names = list("H. sapiens" = "Human", 
             "C. elegans" = "Worm", 
             #"M. musculus" = "Mouse",  #Mouse doesn't work, gets stuck on the download, unclear why
             "D. melanogaster" = "Fly", 
             "S. cerevisiae" = "BYeast")

types = c("trans_table","brief_description_checkbox_intermine","GO_checkbox_intermine","GO_list_intermine")

#Generate all files
for (i in 1:length(names)){
  for (j in 1:length(types)){
    A = getAnnotInterMine(names[[i]],types[j])
  }
}

##FUNCTION
getAnnotInterMine = function(species,type){
  source="intermine"
  genes = NA
  
  today = Sys.Date()
  filename = paste0(species,"_",type,"_",today,".tsv")   
  
  species_data = speciesNames(species) #load from crispAnnotationCollection.R - un-updated version below
  
  
  # if (source=="intermine"){
  latin.species = species_data$latin.species
  mine = species_data$mine
  lookup_use = species_data$lookup_use
  
  # if(is.na(genes)){
    constraint = setConstraints(
      paths = "Gene.organism.name",
      operators = "=",
      values = list(latin.species)
    )
  # }else{
  #   constraint = setConstraints(
  #     paths = c("Gene.organism.name",
  #               paste0("Gene.",lookup_use,"Identifier")),
  #     operators = rep("=", 2),
  #     values = list(latin.species,genes)
  #   )
  # }
    
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
  test = testConnection(im$mine) #load from crispAnnotationCollection.R - un-updated version below
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
    write.table(data,paste0(filename),sep="\t")   
  }
  
}


##OLD FUNCTIONS
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


