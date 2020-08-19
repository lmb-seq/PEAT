#Requires crispAnnotationCollection.R,
# gene_lists_UI+functions.R

library(gProfileR) #getAnnot, gProfile enrichment, orthologs
library(ComplexHeatmap) #annotation heatmap

##REQUIRES PER PROGRAM FUNCTIONS  
# return_primary_gene_ID() #a reactive function returning a list of 
#raw - vector of gene ids for enrichment testing
#primary - vector of primary gene ids (for genes in raw) as defined by intermine for intermine widgets

# return_primary_gene_ID_orthologs() as above but for the orthologs tab

##MUST BE CREATED IN PROGRAM
#input$species #vector of current species
#data_cache$connectionError #NULL if untested, TRUE if no connection possible, FALSE if have made a connection
#data_cache$data_loaded #TRUE if data loaded, else NULL

##Downloads section
# output$downloadplot6 <- downloadHandler(
#   filename = function() {
#     name = paste0("gProfile_enrichment_data_figure-", Sys.Date(), ".svg")
#     return(name)
#   },
#   content = function(file) {
#     if (is.null(make_gprofile_enrichment_table())){
#       return(NULL)
#     }else{
#       #ensg_gene_names = make_gprofile_enrichment_table()$genes_used
#       update = plot6_plot()
#       max_genes = data_cache$max_genes# length(unique(ensg_gene_names))
#       vert_size = 0.15 * max_genes
#       vert_size = ceiling(vert_size) + 8
#       # if (max_genes > 50){
#       #   svg(file,height=28,width=10)
#       # }else{
#       #   svg(file,height=14,width=10)
#       # }
#       svg(file,height=vert_size,width=10)
#       draw(plot6_plot())
#       dev.off()
#     }
#   }
# )

##Extra UI bits to add on pages
# output$analysis_extra <- renderUI({
#   return(NULL)
# })
# output$orthologs_extra <- renderUI({
#   return(NULL)
# })



##Tabs
# tabPanel("Enrichment analysis", uiOutput("analysis"))
# tabPanel("Orthologs", uiOutput("orthologs"))


#Include in createGeneListInput() to be able to create relevant gene lists
# if (user_input$gene_list_genes_source == "from_enrich"){
#   gene_names = make_gprofile_enrichment_table()$genes_used
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_annotations")
# }else if (user_input$gene_list_genes_source == "from_enrich_selected"){
#   gene_names = make_gprofile_enrichment_table()$genes_used
#   gene_names = gene_names[input$gprofile_enrichment_table_rows_selected]
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_annotations")
# }else if (user_input$gene_list_genes_source == "from_enrich_ipd"){
#   gene_names = make_publication_enrichment_table()$genes_used
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_publications")
# }else if (user_input$gene_list_genes_source == "from_enrich_selected_ipd"){
#   gene_names = make_publication_enrichment_table()$genes_used
#   gene_names = gene_names[input$protein_domain_enrichment_table_rows_selected]
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_publications")
# }else if (user_input$gene_list_genes_source == "from_ortholog"){
#   gene_names = make_ortholog_table()$genes_used
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_orthologs")
# }else if (user_input$gene_list_genes_source == "from_ortholog_selected"){
#   gene_names = make_ortholog_table()$genes_used
#   gene_names = gene_names[input$ortholog_table_rows_selected]
#   gene_names = unlist(strsplit(as.character(gene_names), ","))
#   gene_names = unique(as.vector(gene_names))
#   #restore gene names - per species
#   gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
#   
#   list_name = paste0(list_name,"_orthologs")
# }


##Other variables
# data_cache$max_genes #set in plot6_plot()








#Text
output$need_file_text = renderText({"Select valid file for analysis..."})
output$enrichment_text <- renderText({"Enrichment calculation is slow. We recommend hiding all but one enrichment before changing filters."})
output$enrichment_text2 <- renderText({"All enrichments use Benjamini Hochberg multiple testing correction for their probabilities."})
output$enrichment_text3 <- renderText({"Please note that query size for enrichment may be (much) smaller than the selected number of genes, as query size is only the number of genes with that type of annotation."}) 
output$no_enrichment_text <- renderText({"No enrichment found."}) 
output$no_enrich_text <- renderText({"No enrichment analysis possible for this species or when offline."}) 


##Enrichment analysis TAB - UI + GENERATE & RENDER OUTPUT
#Final tab output
#Given:
#Returns: UI element for the 'Figure' tab
#Requires UI elements: output$need_file_text, output$volcano_plot
#Requires reactive functions: modifyData()
#Requires reactive variables: input$species
#Requires functions: speciesNames()
output$analysis <- renderUI({
  out = list(br())
  
  if (is.null(data_cache$data_loaded)){ #Prior to file loading
    out[[2]] = textOutput("need_file_text")
  }else{
    if (!is.null(uiOutput("analysis_extra"))){
      out2 = list(uiOutput("analysis_extra"))
      out = c(out,out2)
    }
    
    out2 = list(textOutput("enrichment_text"))
    out2[[2]] = textOutput("enrichment_text2")
    out2[[3]] = textOutput("enrichment_text3")
    out2[[4]] = HTML("<br>")
    out = c(out,out2) 
    
    #Handle connection_error
    if (is.null(data_cache$connectionError)){
      species_data = speciesNames(input$species)
      mine = species_data$mine
      im <- initInterMine(mine=listMines()[mine])
      if(is.na(im$mine)){
        widgets = NULL
      }else{
        widgets = as.data.frame(getWidgets(im))
      }
    }else if (!data_cache$connectionError){
      species_data = speciesNames(input$species)
      mine = species_data$mine
      im <- initInterMine(mine=listMines()[mine])
      if(is.na(im$mine)){
        widgets = NULL
      }else{
        widgets = as.data.frame(getWidgets(im))
      }
    }else{
      widgets = NULL
    }
    
    #if(!is.na(match("GO terms enriched for items in this list.",widgets$description))){
    # out2 = list(HTML(paste0("<h4>GO enrichment</h4>")))
    # out2[[2]] = tags$div(title=paste0("Click to show/hide the GO enrichment analysis of currently chosen data"),
    #                     actionButton("show_go_enrichment","Show GO enrichment"))
    # out2[[3]] = HTML("<br>")
    # out2[[4]] = uiOutput("GO_enrichment")
    # 
    # out2[[5]] = HTML("<br>")
    # 
    # out = c(out,out2)
    #}
    # if(!is.na(match("Publications enriched for genes in this list.",widgets$description))){
    #   out2 = list(HTML(paste0("<h4>Publication enrichment</h4>")))
    #   out2[[2]] = tags$div(title=paste0("Click to show/hide the Publication enrichment analysis of currently chosen data"),
    #                        actionButton("show_publication_enrichment","Show Publication enrichment"))
    #   out2[[3]] = HTML("<br>")
    #   out2[[4]] = uiOutput("Publication_enrichment")
    #   
    #   out2[[5]] = HTML("<br>")
    #   
    #   out = c(out,out2)
    # }
    
    if(input$species == "Human"){
      g_profile_choices = list("Gene Ontology" = "GO",                                            
                               "Kyoto Encyclopedia of Genes and Genomes" = "KEGG",               
                               "Reactome" = "REAC",                                               
                               "TRANSFAC TFBS" = "TF",                                           
                               "miRBase microRNAs" = "MI",                                        
                               "COmprehensive Resource of Mammalian protein complexes" = "CORUM",
                               "Human Phenotype Ontology" = "HP",                                 
                               "Human Protein Atlas" = "HPA",                                    
                               "Online Mendelian Inheritance in Man" = "OMIM")
    }else{
      g_profile_choices = list("Gene Ontology" = "GO",                                            
                               "Kyoto Encyclopedia of Genes and Genomes" = "KEGG",               
                               "Reactome" = "REAC",                                               
                               "TRANSFAC TFBS" = "TF",                                           
                               "miRBase microRNAs" = "MI",                                        
                               "COmprehensive Resource of Mammalian protein complexes" = "CORUM")
    }
    
    out2 = list(HTML(paste0("<h4>gProfile enrichment</h4>")))
    out2[[2]] = (tags$div(title="Select which data source to look for enrichment in",
                          selectInput("gProfile_data_source", label = ("Select data source"), 
                                      choices = g_profile_choices)
    ))
    out2[[3]] = tags$div(title=paste0("Click to show/hide the gProfile enrichment analysis of currently chosen data"),
                         actionButton("show_gProfile_enrichment","Show gProfile enrichment"))
    out2[[4]] = HTML("<br>")
    out2[[5]] = uiOutput("gProfile_enrichment")
    
    out2[[6]] = HTML("<br>")
    
    #Handle connection_error
    if (is.null(data_cache$connectionError)){
      out = c(out,out2)
    }else if (!data_cache$connectionError){
      out = c(out,out2)
    }
    
    if(!is.na(match("Protein Domains enriched for items in this list.",widgets$description))){
      out2 = list(HTML(paste0("<h4>InterMine Protein Domain enrichment</h4>")))
      out2[[2]] = tags$div(title=paste0("Click to show/hide the Protein Domain enrichment analysis of currently chosen data"),
                           actionButton("show_protein_domain_enrichment","Show Protein Domain enrichment"))
      out2[[3]] = HTML("<br>")
      out2[[4]] = uiOutput("Protein_domain_enrichment")
      
      out2[[5]] = HTML("<br>")
      
      out = c(out,out2)
    }
    
    
    if (length(out) < 2){
      out = list(br())
      out[[2]] = textOutput("no_enrich_text")
    }
    
  }
  return(out)
})

#UI parts per enrichment type
#Publication enrichment
#Given: Number of time button clicked
#Returns: Publication enrichment table
#Requires UI elements: 
#Requires reactive functions: showPublicationenrichment()
output$Publication_enrichment <- renderUI({
  n_clicks = showPublicationenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      return(NULL)
    }else{
      if (is.null(dataTableOutput("publication_enrichment_table"))){
        out = list(textOutput("no_enrichment_text"))
      }else{
        out = list(dataTableOutput("publication_enrichment_table"))
      }
      return(out)
    }
  }else{
    return(NULL)
  }
})
#Functions to control show/hide button  
#Given: Clicks from input$show_publication_enrichment
#Returns: The number of clicks on input$show_publication_enrichment
#Requires reactive variables: input$show_publication_enrichment
showPublicationenrichment <- eventReactive(input$show_publication_enrichment, {
  return(input$show_publication_enrichment)
})  
#Given: The number of clicks on button
#Returns: Updates the name on the button
#Requires reactive functions: showPublicationenrichment()
#Requires reactive variables: input$show_publication_enrichment
observe({ #Updates label depending on whether homology shown or not
  req(input$show_publication_enrichment)
  n_clicks = showPublicationenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      updateActionButton(session, "show_publication_enrichment",
                         label = "Show Publication enrichment")
    }else{
      updateActionButton(session, "show_publication_enrichment",
                         label = "Hide Publication enrichment")        
    }
  }
})   

#Protein domain enrichment
#Given: Number of time button clicked
#Returns: Protein domain enrichment table
#Requires UI elements: 
#Requires reactive functions: showProteindomainenrichment()
output$Protein_domain_enrichment <- renderUI({
  n_clicks = showProteindomainenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      return(NULL)
    }else{
      if (is.null(dataTableOutput("protein_domain_enrichment_table"))){
        out = list(textOutput("no_enrichment_text"))
      }else{
        out = list(dataTableOutput("protein_domain_enrichment_table"))
        #out[[2]] = uiOutput("protein_domain_enrichment_list_buttons") #make_proteindomain_enrichment_table() doesn't currently return the genes used so can't make list from them
      }
      return(out)
    }
  }else{
    return(NULL)
  }
})
output$protein_domain_enrichment_list_buttons <- renderUI({
  out= list(HTML("<br>"))
  out[[2]] = (tags$div(title="Create gene list from filtered genes that have annotation. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_wa_ipd_button", "Create gene list")
  ))       
  out[[3]] = HTML("<br>")
  out[[4]] = (tags$div(title="Create gene list from above selected genes. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_wa_ipd_selected_button", "Create gene list from selected")
  ))   
  return(out)
})
#Functions to control show/hide button  
#Given: Clicks from input$show_protein_domain_enrichment
#Returns: The number of clicks on input$show_protein_domain_enrichment
#Requires reactive variables: input$show_protein_domain_enrichment
showProteindomainenrichment <- eventReactive(input$show_protein_domain_enrichment, {
  return(input$show_protein_domain_enrichment)
})  
#Given: The number of clicks on button
#Returns: Updates the name on the button
#Requires reactive functions: showProteindomainenrichment()
#Requires reactive variables: input$show_protein_domain_enrichment
observe({ #Updates label depending on whether homology shown or not
  req(input$show_protein_domain_enrichment)
  n_clicks = showProteindomainenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      updateActionButton(session, "show_protein_domain_enrichment",
                         label = "Show Protein Domain enrichment")
    }else{
      updateActionButton(session, "show_protein_domain_enrichment",
                         label = "Hide Protein Domain enrichment")        
    }
  }
})   

#gProfile enrichment
#Given: Number of time button clicked
#Returns: gProfile enrichment table
#Requires UI elements: 
#Requires reactive functions: showgProfilerenrichment()
output$gProfile_enrichment <- renderUI({
  n_clicks = showgProfilerenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      return(NULL)
    }else{
      if (is.null(dataTableOutput("gprofile_enrichment_table"))){
        out = list(textOutput("no_enrichment_text"))
      }else{
        out = list(dataTableOutput("gprofile_enrichment_table"))
        out[[2]] = uiOutput("gProfile_enrichment_list_buttons")
        out[[3]] = uiOutput("gProfile_enrichment_figure")
      }
      return(out)
    }
  }else{
    return(NULL)
  }
})
output$gProfile_enrichment_list_buttons <- renderUI({
  out= list(HTML("<br>"))
  out[[2]] = (tags$div(title="Create gene list from filtered genes that have annotation. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_wa_button", "Create gene list")
  ))       
  out[[3]] = HTML("<br>")
  out[[4]] = (tags$div(title="Create gene list from above selected genes. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_wa_selected_button", "Create gene list from selected")
  ))   
  return(out)
})
output$gProfile_enrichment_figure <- renderUI({
  out= list(HTML("<br>"))
  #ensg_gene_names = make_gprofile_enrichment_table()$genes_used
  update = plot6_plot()
  if(is.null(update)){
    return(out)
  }
  max_genes = data_cache$max_genes #length(unique(ensg_gene_names))
  vert_size = (10 * max_genes) + 600
  #vert_size = ceiling(vert_size)
  vert_size = paste0(vert_size,"px")
  # if(max_genes > 50){
  #   out[[2]] = plotOutput("plot6", width = "100%", height = "1600px")
  # }else{
  #   out[[2]] = plotOutput("plot6", width = "100%", height = "800px")
  # }
  out[[2]] = plotOutput("plot6", width = "100%", height = vert_size)
  return(out)
})
#Functions to control show/hide button  
#Given: Clicks from input$show_gProfile_enrichment
#Returns: The number of clicks on input$show_gProfile_enrichment
#Requires reactive variables: input$show_gProfile_enrichment
showgProfilerenrichment <- eventReactive(input$show_gProfile_enrichment, {
  return(input$show_gProfile_enrichment)
})  
#Given: The number of clicks on button
#Returns: Updates the name on the button
#Requires reactive functions: showgProfilerenrichment()
#Requires reactive variables: input$show_gProfile_enrichment
observe({ #Updates label depending on whether homology shown or not
  req(input$show_gProfile_enrichment)
  n_clicks = showgProfilerenrichment() #records number of times button clicked
  if (n_clicks > 0){
    if (n_clicks %% 2 == 0){
      updateActionButton(session, "show_gProfile_enrichment",
                         label = "Show gProfile enrichment")
    }else{
      updateActionButton(session, "show_gProfile_enrichment",
                         label = "Hide gProfile enrichment")        
    }
  }
})   


#Generate data for display
#Given: Gene list
#Returns: Table of Publications enriched for those genes
#Requires reactive functions: modifyData()
#Requires reactive variables: input$species
#Requires functions: speciesNames()
make_publication_enrichment_table <- reactive({
  if (is.null(data_cache$data_loaded)){
    return(NULL)
  }else{
    gene_names = return_primary_gene_ID()$primary
    if(is.null(gene_names)){
      return(NULL)
    }
    species_data = speciesNames(input$species)
    
    mine = species_data$mine
    im <- initInterMine(mine=listMines()[mine])
    tryCatch(publication_enrichResult = doEnrichment(
      im = im,
      ids = gene_names,
      widget = "publication_enrichment",
      correction = "Benjamini Hochberg"
      # organism = "Homo sapiens" # optional if results from more than one organism are retrieved
    ),error=function(err){return(NULL)})
    if(!exists("publication_enrichResult")){
      return(as.matrix("No enrichment found."))
    }
    if(is.null(publication_enrichResult)){
      return(as.matrix("No enrichment found."))
    }
    if(is.null(publication_enrichResult$data)){
      return(as.matrix("No enrichment found."))
    }
    if (dim(publication_enrichResult$data)[1] > 0){
      colnames(publication_enrichResult$data) = c("pubmed ID","description","adjusted p-value","count in selected genes","count in all genes")
    }else{
      return(as.matrix("No enrichment found."))
    }
    
    return(publication_enrichResult$data)
  }
})
output$publication_enrichment_table <- DT::renderDataTable({
  data = make_publication_enrichment_table()
  data = as.data.frame(data)
  if(dim(data)[2] > 1){
    data[,1] = createLink(paste0("https://www.ncbi.nlm.nih.gov/pubmed/",data[,1]),data[,1])
  }
  
},rownames = FALSE, escape = FALSE,options = list(scrollX = TRUE))  

#Given: Gene list
#Returns: list of data.frame of protein domains enriched for those genes (table);
#  vector of gene names with annotations (genes_used)
#Requires reactive functions: modifyData()
#Requires reactive variables: input$species
#Requires functions: speciesNames()
make_proteindomain_enrichment_table <- reactive({
  if (is.null(data_cache$data_loaded)){
    return(NULL)
  }else{
    gene_names = return_primary_gene_ID()$primary
    if(is.null(gene_names)){
      return(NULL)
    }
    species_data = speciesNames(input$species)
    
    mine = species_data$mine
    im <- initInterMine(mine=listMines()[mine])
    prot_dom_enrichResult = doEnrichment(
      im = im,
      ids = gene_names,
      widget = "prot_dom_enrichment_for_gene",
      correction = "Benjamini Hochberg"
      # organism = "Homo sapiens" # optional if results from more than one organism are retrieved
    )
    if(!exists("prot_dom_enrichResult")){
      return(as.matrix("No enrichment found."))
    }
    if(is.null(prot_dom_enrichResult)){
      return(as.matrix("No enrichment found."))
    }
    if(is.null(prot_dom_enrichResult$data)){
      return(as.matrix("No enrichment found."))
    }
    if (dim(prot_dom_enrichResult$data)[1] > 0){
      #colnames(prot_dom_enrichResult$data) = c("Domain ID","description","adjusted p-value","count in selected genes","count in all genes")
    }else{
      return(as.matrix("No enrichment found."))
    }
    
    genes_wa = NULL
    return(list("table" = prot_dom_enrichResult$data,"genes_used" = genes_wa))
  }
})
output$protein_domain_enrichment_table <- DT::renderDataTable({
  data = make_proteindomain_enrichment_table()
},rownames = FALSE,options = list(scrollX = TRUE))  

#Given: Gene list
#Returns: list of data.frame of various enriched for those genes (table);
#  vector of gene names with annotations (genes_used)
#Requires reactive functions: return_primary_gene_ID(), 
#Requires reactive variables: input$species,
#Requires functions: speciesNames
make_gprofile_enrichment_table <- reactive({
  if (is.null(data_cache$data_loaded)){
    return(NULL)
  }else{
    # data.object = filteredData()
    # 
    # species_data = speciesNames(input$species)
    # short.species = species_data$short.species
    # 
    # format = data.object$data_format
    # format_details = formatLookup(format)
    # gene_id_position = format_details$gene_id_position
    # correct_colnames = format_details$correct_colnames
    # gene_id_name = correct_colnames[gene_id_position]
    # 
    # data = data.object$filtered_data
    # 
    # gene_names = (data[,gene_id_name])
    # 
    # #convert to ENSG
    # #ensg_gene_names = gconvert(gene_names,organism=short.species,target="ENSG")
    # # human needs ENSG
    # # fly fine with FBgn
    # # worm fine with WBGene
    # # mouse fine with MGI
    # ensg_gene_names = as.character(gene_names)
    ensg_gene_names = return_primary_gene_ID()$raw
    if(is.null(ensg_gene_names)){
      return(NULL)
    }
    
    species_data = speciesNames(input$species)
    short.species = species_data$short.species
    
    #lookup enrichment
    table = gprofiler(unique(ensg_gene_names),organism=short.species,significant=F,correction_method = "fdr",src_filter=input$gProfile_data_source)
    
    #order by p-value
    table = table[order(table[,3]),]
    
    #just columns I want
    #id, description, p-value, count in chosen, count in all
    table = table[,c(9,12,3,6,4,14)]
    colnames(table) = c("Term ID","description","adjusted p-value","count in selected genes","count in all genes","genes with annotation")
    
    genes_wa = table[,6]
    
    return(list("table" = table, "genes_used" = genes_wa))
  }
})
output$gprofile_enrichment_table <- DT::renderDataTable({
  data = make_gprofile_enrichment_table()$table
  data = as.data.frame(data)
  if(dim(data)[2] > 1){
    data[,3] = as.numeric(as.character(data[,3]))
    data[,4] = as.numeric(as.character(data[,4]))
    data[,5] = as.numeric(as.character(data[,5]))
    
    types = c("GO","KEGG","REAC","TF","MI","CORUM","HP","HPA","OMIM")
    
    if(input$gProfile_data_source == types[1]){
      data[,1] = createLink(paste0("http://amigo.geneontology.org/amigo/term/",data[,1]),data[,1])
    }else if(input$gProfile_data_source == types[2]){
      short_data = sub("KEGG:","",data[,1])
      data[,1] = createLink(paste0("http://www.genome.jp/dbget-bin/www_bget?path:map",short_data),data[,1])
    }else if(input$gProfile_data_source == types[3]){
      short_data = sub("REAC:","",data[,1])
      data[,1] = createLink(paste0("https://reactome.org/content/detail/",short_data),data[,1])
    }else if(input$gProfile_data_source == types[4]){
      #requires login - everything seems to be in professional version only
      short_data = sub("TF:","",data[,1])
      short_data = sub("_\\d","",short_data)
      data[,1] = createLink(paste0("http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=",short_data),data[,1])
    }else if(input$gProfile_data_source == types[5]){
      short_data = sub("MI:","",data[,1])
      short_data = sub("\\*","",short_data)
      data[,1] = createLink(paste0("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=",short_data),data[,1])
    }else if(input$gProfile_data_source == types[6]){
      # https://mips.helmholtz-muenchen.de/corum/ #doesn't have per search links
      # short_data = sub("OMIM:","",data[,1])
      data[,1] = createLink(paste0("https://mips.helmholtz-muenchen.de/corum/"),data[,1])
    }else if(input$gProfile_data_source == types[7]){
      short_data = sub("HP:","",data[,1])
      short_data = as.numeric(short_data) #remove leading zeros
      data[,1] = createLink(paste0("https://mseqdr.org/hpo_browser.php?",short_data),data[,1])
    }else if(input$gProfile_data_source == types[8]){
      # short_data = sub("OMIM:","",data[,1])
      # data[,1] = createLink(paste0("https://www.omim.org/entry/",short_data),data[,1])
    }else if(input$gProfile_data_source == types[9]){
      short_data = sub("OMIM:","",data[,1])
      data[,1] = createLink(paste0("https://www.omim.org/entry/",short_data),data[,1])
    }
  }
  return(data)
},rownames = FALSE, escape = FALSE,options = list(scrollX = TRUE))  


#Heatmap for annotations of selected genes
#Given: matrix from genes_for_heatmap
#Returns: Heatmap for annotations of selected genes
#Requires reactive functions: make_gprofile_enrichment_table()
#Requires reactive variables: 
plot6_plot <- function(){
  if (is.null(make_gprofile_enrichment_table())){
    return(NULL)
  }
  ensg_gene_names = make_gprofile_enrichment_table()$genes_used
  if(length(ensg_gene_names) < 1){
    return(NULL)
  }
  table = make_gprofile_enrichment_table()$table
  
  max_col = 10
  n_col = min(max_col, dim(table)[1])
  
  if(length(unique(ensg_gene_names)) == 1){ #handle cases with one gene
    dot_plot_data = matrix(0,ncol = n_col, nrow = (1+length(unique(ensg_gene_names))))
    rownames(dot_plot_data) = c(unique(ensg_gene_names),"")
  }else{
    dot_plot_data = matrix(0,ncol = n_col, nrow = length(unique(ensg_gene_names)))
    rownames(dot_plot_data) = unique(ensg_gene_names)
  }
  colnames(dot_plot_data) = paste0(table[1:n_col,1]," - ",table[1:n_col,2])
  #colnames(dot_plot_data) = table[1:n_col,1]
  
  #unique genes with annotations
  genes_wa = table[1:n_col,6]
  
  #restore gene names - per species
  genes_wa = gsub("WBGENE","WBGene",genes_wa) #Worm
  
  genes_wa = (strsplit(as.character(genes_wa), ","))
  
  for (i in 1:length(genes_wa)){
    dot_plot_data[match(genes_wa[[i]],rownames(dot_plot_data)),i] = 1
  }
  
  #cut it down to just the genes with annotations
  dot_plot_data2 = dot_plot_data[(rowSums(dot_plot_data) > 0),]
  
  #handle cases where this reduces us to one row
  if(is.null(dim(dot_plot_data2))){
    dot_plot_data = rbind(dot_plot_data,rep(1,dim(dot_plot_data)[2]))
    dot_plot_data = dot_plot_data[(rowSums(dot_plot_data) > 0),]
    dot_plot_data[dim(dot_plot_data)[1],] = rep(0,dim(dot_plot_data)[2])
    data_cache$max_genes = 1
  }else{
    dot_plot_data = dot_plot_data[(rowSums(dot_plot_data) > 0),]
    data_cache$max_genes = dim(dot_plot_data)[1]
  }
  
  if(sum(dot_plot_data) != (dim(dot_plot_data)[1]*dim(dot_plot_data)[2])){ #can't handle all genes having all annotations
    ComplexHeatmap::Heatmap(dot_plot_data,col=c("white","black"),cluster_columns = F, show_column_dend = F, row_names_side="left",column_names_side="top",show_heatmap_legend=F,
                            show_row_dend = F,cluster_rows = T,
                            column_names_gp = gpar(fontsize = 10),rect_gp = gpar(col= "white"),column_names_max_height= unit(1000,"mm"),
                            row_names_gp = gpar(fontsize = 10))
  }else{
    return(NULL)
  }
  # heatmap.2(dot_plot_data,col=c("white","black"),dendrogram = "none",trace="none",key = F,
  #           sepwidth=c(0.1,0.1),
  #           sepcolor="white",
  #           colsep=1:ncol(dot_plot_data),
  #           rowsep=1:nrow(dot_plot_data),cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8))
}
output$plot6 <- renderPlot({
  draw(plot6_plot())
})



#Orthologs
#Given: table from output$ortholog_table
#Returns: this as a datatable
#Requires reactive functions: filteredData()
#Requires reactive variables: input$species,
#Requires functions: speciesNames
output$orthologs <- renderUI({
  if (is.null(data_cache$data_loaded)){ #Prior to file loading
    return(NULL)
  }else{
    ortholog_species = list("Homo sapiens" = "hsapiens","Mus musculus" = "mmusculus","Saccharomyces cerevisiae" = "scerevisiae","Rattus norvegicus" = "rnorvegicus","Drosophila melanogaster" = "dmelanogaster","Caenorhabditis elegans" = "celegans","Danio rerio" = "drerio","Ailuropoda melanoleuca" = "amelanoleuca","Algerian mouse" = "mspreteij","Anas platyrhynchos" = "aplatyrhynchos","Anolis carolinensis" = "acarolinensis","Astyanax mexicanus" = "amexicanus","Bos taurus" = "btaurus","Brazilian guinea pig" = "caperea","Callithrix jacchus" = "cjacchus","Canis familiaris" = "cfamiliaris","Cavia porcellus" = "cporcellus","Chinese hamster CHOK1GS" = "cchok1gshd","Chinese hamster CriGri" = "ccrigri","Chlorocebus sabaeus" = "csabaeus","Choloepus hoffmanni" = "choffmanni","Ciona intestinalis" = "cintestinalis","Ciona savignyi" = "csavignyi","Damara mole rat" = "fdamarensis","Dasypus novemcinctus" = "dnovemcinctus","Degu" = "odegus","Dipodomys ordii" = "dordii","Echinops telfairi" = "etelfairi","Equus caballus" = "ecaballus","Erinaceus europaeus" = "eeuropaeus","Felis catus" = "fcatus","Ficedula albicollis" = "falbicollis","Gadus morhua" = "gmorhua","Gallus gallus" = "ggallus","Gasterosteus aculeatus" = "gaculeatus","Golden Hamster" = "mauratus","Gorilla gorilla" = "ggorilla","Ictidomys tridecemlineatus" = "itridecemlineatus","Latimeria chalumnae" = "lchalumnae","Lepisosteus oculatus" = "loculatus","Lesser Egyptian jerboa" = "jjaculus","Long-tailed chinchilla" = "clanigera","Loxodonta africana" = "lafricana","Macaca mulatta" = "mmulatta","Meleagris gallopavo" = "mgallopavo","Microcebus murinus" = "mmurinus","Monodelphis domestica" = "mdomestica","Mustela putorius furo" = "mfuro","Myotis lucifugus" = "mlucifugus","Naked mole-rat female" = "hfemale","Naked mole-rat male" = "hmale","Nomascus leucogenys" = "nleucogenys","Northern American deer mouse" = "pbairdii","Ochotona princeps" = "oprinceps","Oreochromis niloticus" = "oniloticus","Ornithorhynchus anatinus" = "oanatinus","Oryctolagus cuniculus" = "ocuniculus","Oryzias latipes" = "olatipes","Otolemur garnettii" = "ogarnettii","Ovis aries" = "oaries","Pan troglodytes" = "ptroglodytes","Papio anubis" = "panubis","Pelodiscus sinensis" = "psinensis","Petromyzon marinus" = "pmarinus","Poecilia formosa" = "pformosa","Pongo abelii" = "pabelii","Prairie vole" = "mochrogaster","Procavia capensis" = "pcapensis","Pteropus vampyrus" = "pvampyrus","Ryukyu mouse" = "mcaroli","Sarcophilus harrisii" = "sharrisii","Shrew mouse" = "mpahari","Sorex araneus" = "saraneus","Sus scrofa" = "sscrofa","Taeniopygia guttata" = "tguttata","Takifugu rubripes" = "trubripes","Tarsier" = "csyrichta","Tetraodon nigroviridis" = "tnigroviridis","Tupaia belangeri" = "tbelangeri","Tursiops truncatus" = "ttruncatus","Upper Galilee mountains blind mole rat" = "ngalili","Vicugna pacos" = "vpacos","Wallaby" = "neugenii","Xenopus tropicalis" = "xtropicalis","Xiphophorus maculatus" = "xmaculatus")
    
    names = as.list(c("None",unlist(ortholog_species))) #add to beginning of list
    
    #remove current species
    species_data = speciesNames(input$species)
    short.species = species_data$short.species
    names[[match(short.species,names)]] = NULL
    
    out = list(br())
    if (!is.null(uiOutput("orthologs_extra"))){
      out2 = list(uiOutput("orthologs_extra"))
      out = c(out,out2)
    }
    
    out2 = list(tags$div(title=paste0("Select a species to find orthologs in."),
                         selectInput("select_orth_org", label = ("Orthologs in:"), 
                                     choices = names, 
                                     selected = 1)))
    out2[[2]] = (dataTableOutput("ortholog_table"))
    out2[[3]] = uiOutput("ortholog_enrichment_list_buttons")
    out = c(out,out2) 
    
    return(out)
  }
})
output$ortholog_enrichment_list_buttons <- renderUI({
  out= list(HTML("<br>"))
  out[[2]] = (tags$div(title="Create gene list from filtered genes that have orthologs WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_ortholog_button", "Create gene list")
  ))       
  out[[3]] = HTML("<br>")
  out[[4]] = (tags$div(title="Create gene list from above selected genes. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_ortholog_selected_button", "Create gene list from selected")
  ))   
  return(out)
})

#Given: Gene list
#Returns: list of data.frame of various enriched for those genes (table);
#  vector of gene names with orthologs (genes_used)
#Requires reactive functions: return_primary_gene_ID_orthologs(), 
#Requires reactive variables: input$species,
#Requires functions: speciesNames
make_ortholog_table <- reactive({
  if (is.null(data_cache$data_loaded)){
    return(NULL)
  }else{
    if(input$select_orth_org == "None"){
      return(NULL)
    }
    # data.object = filteredData()
    # 
    # species_data = speciesNames(input$species)
    # short.species = species_data$short.species
    # 
    # format = data.object$data_format
    # format_details = formatLookup(format)
    # gene_id_position = format_details$gene_id_position
    # correct_colnames = format_details$correct_colnames
    # gene_id_name = correct_colnames[gene_id_position]
    # 
    # data = data.object$filtered_data
    # 
    # gene_names = (data[,gene_id_name])
    # 
    # #convert to ENSG
    # #ensg_gene_names = gconvert(gene_names,organism=short.species,target="ENSG")
    # # human needs ENSG
    # # fly fine with FBgn
    # # worm fine with WBGene
    # # mouse fine with MGI
    # ensg_gene_names = as.character(gene_names)
    ensg_gene_names = return_primary_gene_ID_orthologs()$raw
    if(is.null(ensg_gene_names)){
      return(NULL)
    }
    
    species_data = speciesNames(input$species)
    short.species = species_data$short.species
    
    #lookup enrichment
    table = gorth(unique(ensg_gene_names),source_organism=short.species,target_organism = input$select_orth_org, filter_na = T)
    
    
    if(dim(table)[1] > 0){
      #just columns I want
      table = table[,c(2,5,6,7)]
      colnames(table) = c("Gene ID","Target ID","Target name","Target description")
      
      genes_wa = table[,1]
    }else{
      return(NULL)
    }
    
    return(list("table" = table, "genes_used" = genes_wa))
  }
})
output$ortholog_table <- DT::renderDataTable({
  data = make_ortholog_table()
  if(is.null(data)){
    data = "No orthologs found"
    data = as.matrix(data)
  }else{
    data = data$table
    if(dim(data)[1] > 0){
      data[,2] = createLink(paste0("http://www.ensembl.org/Gene/Summary?g=",data[,2]),data[,2])
    }else{
      data = "No orthologs found"
      data = as.matrix(data)
    }
  }
  return(data)
},rownames = FALSE, escape = FALSE,options = list(scrollX = TRUE))


#Gene list observes
observeEvent(input$create_gene_list_wa_button, {
  user_input$gene_list_genes_source = "from_enrich"
  do = createGeneList()
})  
observeEvent(input$create_gene_list_wa_selected_button, {
  user_input$gene_list_genes_source = "from_enrich_selected"
  do = createGeneList()
})  
observeEvent(input$create_gene_list_ortholog_button, {
  user_input$gene_list_genes_source = "from_ortholog"
  do = createGeneList()
})  
observeEvent(input$create_gene_list_ortholog_selected_button, {
  user_input$gene_list_genes_source = "from_ortholog_selected"
  do = createGeneList()
})  
observeEvent(input$create_gene_list_wa_ipd_button, {
  user_input$gene_list_genes_source = "from_enrich_ipd"
  do = createGeneList()
})  
observeEvent(input$create_gene_list_wa_ipd_selected_button, {
  user_input$gene_list_genes_source = "from_enrich_selected_ipd"
  do = createGeneList()
})  
