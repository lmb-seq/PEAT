
##ADD TO UI TO USE
#To sidebar
# menuItem("Gene lists", tabName = "cufflinks_gene_level", icon = icon("briefcase"),
#          uiOutput("gene_lists")
# ),

#To tabs
# tabPanel("Gene lists", uiOutput("gene_list_table_ui"))
# tabPanel("Gene lists genes", uiOutput("gene_list_genes_table"))


##REQUIRES PER PROGRAM FUNCTIONS  
# filterFileData() #a reactive function returning a list of 
#table - the data to display about the genes in a particular gene list, could just be the genes
#lookup - NULL when only first two columns of table are non-numeric, else first three columns are non-numeric

# hyper_test_input_background() #a reactive function returning a vector of the total population for the hypergeometric comparison

# createGeneListInput() #a reactive function returning a list of 
#gene_names - vector of gene IDs
#metadata - list of species, length of gene_names, input file name, ui$s1_fpkm_num,ui$s2_fpkm_num,
#ui$log2FC_num,ui$sig_checkbox,select2InputFoundIn,select2InputNotFoundIn,user_input$select_gene_list,description,Other loaded files

##PER PROGRAM UI
#gene_list_table_ui_extra - containing anything extra you want to put on "Gene lists" tab


##MUST BE CREATED IN PROGRAM
#input$species #vector of current species
#data_cache$gene_list_table_columns #vector of which columns of metadata should be displayed in this program



#General storage variables that will be/need to be created - I think they don't need anymore than is found in the below functions
# data_cache$gene_lists_metadata
# data_cache$gene_lists
# 
# data_cache$gene_lists_DT #made by make_gene_list_table()


# #Created by observeEvent(input$newValuepHyper
# user_input$background_gene_n 
# user_input$background_gene_with_n 
# user_input$query_gene_n 
# user_input$query_gene_with_n 


##MENU
output$gene_lists <- renderUI({
  species = input$species
  
  species_data = speciesNames(species)
  if (is.null(species_data)){
    gene.lookup.long = "Gene ID"
  }else{
    gene.lookup.long = species_data$gene.lookup.long
  }
  
  choice_gene_name = "Gene_list"
  existing_names = make_gene_list_table()
  if(!is.null(existing_names)){
    existing_names = existing_names[,1]
    i = 2
    temp_name = choice_gene_name
    while(!is.na(match(temp_name,existing_names))){
      temp_name = paste0(choice_gene_name,"_",i)
      i = i+1
    }
    choice_gene_name = temp_name
  }
  
  out = list(HTML("<br>"))
  out[[2]] = (HTML(paste0("<h4>Create lists</h4>")))
  out[[3]] = (tags$div(title="Name for gene list",textInput("create_gene_list_name","Name gene list", value = choice_gene_name)))
  out[[4]] = (tags$div(title="Create gene list from current selection, with the above name. WARNING this deselects current gene list filtering.",
                       actionButton("create_gene_list_button", "Create gene list")
  ))
  out[[5]] = tags$div(title=paste0("Create gene list from file; can be a list created in this app (uses saved name) or a list of ",gene.lookup.long,"s which will be given the above name. WARNING this deselects current gene list filtering."),
                      fileInput("gene_list_file", label = "Load gene list")
  )
  out[[6]] = (HTML(paste0("<h4>Compare lists</h4>")))
  out[[7]] = (tags$div(title="Create gene list containing all genes from selected gene lists. Must select at least two lists.",
                       actionButton("union_gene_list_button", "Union of selected gene lists")
  ))
  out[[8]] = (tags$div(title="Create gene list containing only genes in all selected gene lists. Must select at least two lists.",
                       actionButton("intersect_gene_list_button", "Intersect of selected gene lists")
  ))
  out[[9]] = (tags$div(title="Create gene list containing genes in first selected gene list, but not in subsequently selected gene lists. Must select at least two lists. WARNING order of selection of lists matters.",
                       actionButton("setdiff_gene_list_button",  HTML("Set difference of <br/> selected gene lists"))
  ))
  out[[10]] = HTML("<br>")
  
  return(out)
})

  #create a new union gene list
  observeEvent(input$union_gene_list_button, {
    metadata = (data_cache$gene_lists_metadata)
    gene_lists = (data_cache$gene_lists)
    if (is.null(metadata)){
      return(NULL)
    }else{
      rows = input$gene_list_table_rows_selected
      if(!is.null(rows)){
        if(length(rows) > 1){
          new_name = paste0("Union of ",metadata[[rows[1]]]$gene_list_name)
          for (i in 2:(length(rows))){
            new_name = paste0(new_name,",",metadata[[rows[i]]]$gene_list_name)
          }
          
          new_genes = union(gene_lists[[rows[1]]],gene_lists[[rows[2]]])
          new_metadata = metadata[[rows[1]]]
          for (i in 1:length(new_metadata)){
            if(new_metadata[[i]] != metadata[[rows[2]]][[i]]){
              new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[2]]][[i]])
            }
          }
          if(length(rows) > 2){
            for (j in 3:(length(rows))){
              new_genes = union(new_genes,gene_lists[[rows[j]]])
              for (i in 1:length(new_metadata)){
                if(new_metadata[[i]] != metadata[[rows[j]]][[i]]){
                  new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[j]]][[i]])
                }
              }
            }
          }
          new_metadata$gene_list_name = new_name
          new_metadata[[2]] = length(new_genes)
          
          gene_lists[[(length(gene_lists)+1)]] = new_genes
          metadata[[(length(metadata)+1)]] = new_metadata
          
          data_cache$gene_lists_metadata = metadata
          data_cache$gene_lists = gene_lists
        }
      }
    }
  })  
  #create a new intersect gene list
  observeEvent(input$intersect_gene_list_button, {
    metadata = (data_cache$gene_lists_metadata)
    gene_lists = (data_cache$gene_lists)
    if (is.null(metadata)){
      return(NULL)
    }else{
      rows = input$gene_list_table_rows_selected
      if(!is.null(rows)){
        if(length(rows) > 1){
          new_name = paste0("Intersect of ",metadata[[rows[1]]]$gene_list_name)
          for (i in 2:(length(rows))){
            new_name = paste0(new_name,",",metadata[[rows[i]]]$gene_list_name)
          }
          
          new_genes = intersect(gene_lists[[rows[1]]],gene_lists[[rows[2]]])
          new_metadata = metadata[[rows[1]]]
          for (i in 1:length(new_metadata)){
            if(new_metadata[[i]] != metadata[[rows[2]]][[i]]){
              new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[2]]][[i]])
            }
          }
          if(length(rows) > 2){
            for (j in 3:(length(rows))){
              new_genes = intersect(new_genes,gene_lists[[rows[j]]])
              for (i in 1:length(new_metadata)){
                if(new_metadata[[i]] != metadata[[rows[j]]][[i]]){
                  new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[j]]][[i]])
                }
              }
            }
          }
          new_metadata$gene_list_name = new_name
          new_metadata[[2]] = length(new_genes)
          
          gene_lists[[(length(gene_lists)+1)]] = new_genes
          metadata[[(length(metadata)+1)]] = new_metadata
          
          data_cache$gene_lists_metadata = metadata
          data_cache$gene_lists = gene_lists
        }
      }
    }
  })  
  #create a new setdiff gene list
  observeEvent(input$setdiff_gene_list_button, {
    metadata = (data_cache$gene_lists_metadata)
    gene_lists = (data_cache$gene_lists)
    if (is.null(metadata)){
      return(NULL)
    }else{
      rows = input$gene_list_table_rows_selected
      if(!is.null(rows)){
        if(length(rows) > 1){
          new_name = paste0("Set difference of ",metadata[[rows[1]]]$gene_list_name, " from ",metadata[[rows[2]]]$gene_list_name)
          
          new_genes = setdiff(gene_lists[[rows[1]]],gene_lists[[rows[2]]])
          new_metadata = metadata[[rows[1]]]
          for (i in 1:length(new_metadata)){
            if(new_metadata[[i]] != metadata[[rows[2]]][[i]]){
              new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[2]]][[i]])
            }
          }
          if(length(rows) > 2){
            for (j in 3:(length(rows))){
              new_genes = setdiff(new_genes,gene_lists[[rows[j]]])
              new_name = paste0(new_name,",",metadata[[rows[j]]]$gene_list_name)
              for (i in 1:length(new_metadata)){
                if(new_metadata[[i]] != metadata[[rows[j]]][[i]]){
                  new_metadata[[i]] = paste0(new_metadata[[i]],", ",metadata[[rows[j]]][[i]])
                }
              }
            }
          }
          new_metadata$gene_list_name = new_name
          new_metadata[[2]] = length(new_genes)
          
          gene_lists[[(length(gene_lists)+1)]] = new_genes
          metadata[[(length(metadata)+1)]] = new_metadata
          
          data_cache$gene_lists_metadata = metadata
          data_cache$gene_lists = gene_lists
        }
      }
    }
  })  
  
  #create gene list observes - additional ones with diff source should be created in this style and createGeneListInput() modified to use them
  #create a new gene list from selected data
  observeEvent(input$create_gene_list_button, {
    user_input$gene_list_genes_source = "standard"
    do = createGeneList()
  })  

##TABS
output$gene_list_table_ui <- renderUI({
  #dataTableOutput("gene_list_table")
  gene_lists = data_cache$gene_lists
  if (is.null(gene_lists)){
    return(NULL)
  }else{
    out = list(box(width=14,
        column(12,dataTableOutput("gene_list_table")),
        tags$script("$(document).on('click', '#gene_list_table button', function () {
                        Shiny.onInputChange('lastClickId',this.id);
                        Shiny.onInputChange('lastClick', Math.random())
          });"),
        column(6,offset = 6,
               HTML('<div class="btn-group" role="group" aria-label="Basic example">'),
               #actionButton(inputId = "Add_row_head",label = "Add a new row"),
               actionButton(inputId = "Del_row_head",label = "Delete selected gene list(s)"),
               actionButton(inputId = "Modify_row_head",label = "Modify selected row"),
               (tags$div(title="Save selected gene list - saves the last selected gene list.",
                         downloadButton("save_gene_list_button", "Download selected gene list"))),
               (tags$div(title="Perform hypergeometric test to calculate enrichment - Select exactly two gene lists, the first selected should be the list of all genes with a particular feature, the second the list of genes you wish to see if are enriched for this feature. WARNING - do not use this feature if you do not understand it, the authors take no responsibility for a misuse of statistics.",
                         actionButton(inputId = "hyper_test_gene",label = "Perform enrichment test"))),
               # (tags$div(title="Click to produce Venn diagram from selected lists",
               #           actionButton(inputId = "venn_gene",label = "Draw Venn diagram"))),
               HTML('</div>')
        )
    ))
    out[[2]] = uiOutput("gene_list_table_ui_extra")
    return(out)
  }
})
  
  output$gene_list_table<-DT::renderDataTable({
    update = make_gene_list_table()
    DT=data_cache$gene_lists_DT
    
    
    if(!is.null(data_cache$gene_list_table_columns)){
      cnames = colnames(DT)
      cnames = cnames[data_cache$gene_list_table_columns]
      DT = subset(DT,select=cnames)
    }
    
    #data_cache$gene_lists_DT2 = subset(DT,select=cnames)
    
    if(is.null(update)){
      return(NULL)
    }
    datatable(DT,escape=F,rownames = FALSE,extensions = 'Buttons', options = list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE))
  })
  
  ##Managing in row modification
  #Display modification popup on click
  observeEvent(input$Modify_row_head,{
    showModal(modal_modify)
  })
  
  #Delete row
  observeEvent(input$Del_row_head,{
    gene_lists = (data_cache$gene_lists)
    metadata = (data_cache$gene_lists_metadata)
    rows = input$gene_list_table_rows_selected
    if(!is.null(rows)){
      rows = sort(rows)
      for(i in length(rows):1){
        gene_lists[[rows[i]]] = NULL
        metadata[[rows[i]]] = NULL
      }
      if(length(gene_lists) < 1){
        gene_lists = NULL
        metadata = NULL
      }
      data_cache$gene_lists = gene_lists
      data_cache$gene_lists_metadata = metadata
    }
  })
  
  modal_modify <- modalDialog(
    fluidPage(
      h3(strong("Row modification"),align="center"),
      hr(),
      dataTableOutput('row_modif'),
      actionButton("save_changes","Save changes"),
      
      tags$script(HTML("$(document).on('click', '#save_changes', function () {
                       var list_value=[]
                       for (i = 0; i < $( '.new_input' ).length; i++)
                       {
                       list_value.push($( '.new_input' )[i].value)
                       }
                       Shiny.onInputChange('newValue', list_value)
                       });"))
      )
  ,size="l",easyClose=TRUE)
  
  #Row modification table
  output$row_modif<-renderDataTable({
    selected_row = input$gene_list_table_rows_selected 
    selected_row = selected_row[length(selected_row)]
    
    if(is.null(selected_row)){
      return(NULL)
    }
    
    old_row=data_cache$gene_lists_DT[selected_row]
    row_change=list()
    
    #subset to only change desired columns
    old_row = old_row[,c(1,13)]
    
    for (i in colnames(old_row)){
      if (is.numeric(data_cache$gene_lists_DT[[i]])){
        row_change[[i]]<-paste0('<input class="new_input" type="number" id=new_',i,', value=',data_cache$gene_lists_DT[selected_row][[i]],'><br>')
      }else{
        row_change[[i]]<-paste0("<input class='new_input' type='text' id=new_",i,", value=",gsub('([[:space:]])|\\s+','_',data_cache$gene_lists_DT[selected_row][[i]]),"><br>")
      }
    }
    row_change=as.data.table(row_change)
    if(length(row_change) > 0){
      setnames(row_change,colnames(old_row))
    }
    DT=rbind(old_row,row_change)
    if(dim(DT)[1] == 2){
      rownames(DT)<-c("Current values","New values")
    }
    return(DT)
    
  },escape=F,options=list(dom='t',ordering=F),selection="none")
  
  #Update with the new values
  observeEvent(input$newValue,{
    removeModal()
    newValue=lapply(input$newValue, function(col) {
      if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
        as.numeric(as.character(col))
      } else {
        col
      }
    })
    DF=data.frame(lapply(newValue, function(x) t(data.frame(x))))
    
    selected_row = input$gene_list_table_rows_selected 
    selected_row = selected_row[length(selected_row)]
    
    metadata = data_cache$gene_lists_metadata[[selected_row]]
    
    #gene list name must be unique
    existing_names = make_gene_list_table()
    if(!is.null(existing_names)){
      existing_names = existing_names[,1]
      i = 2
      temp_name = as.character(DF[,1])
      while(!is.na(match(temp_name,existing_names))){
        temp_name = paste0(as.character(DF[,1]),"_",i)
        i = i+1
      }
      metadata$gene_list_name = temp_name
    }
    
    metadata[[12]] = as.character(DF[,2])
    data_cache$gene_lists_metadata[[selected_row]] = metadata
    
    #add back in the other columns
    old_row=data_cache$gene_lists_DT[selected_row]
    old_row = old_row[,c(2:12,14)]
    DF=cbind(DF[,1],old_row[,1:11],DF[,2],old_row[,12])
    
    colnames(DF)=colnames(data_cache$gene_lists_DT)
    data_cache$gene_lists_DT[selected_row]<-DF
  })
  
  
  
  ##Hypergeometric enrichment test
  modal_hyper_test_gene <- modalDialog(
    fluidPage(
      h3(strong("Hypergeometric enrichment test"),align="center"),
      dataTableOutput('hyper_test_gene_table'),
      hr(),
      actionButton("calculate_phyper","Calculate enrichment"),
      tags$script(HTML("$(document).on('click', '#calculate_phyper', function () {
                       var list_value=[]
                       for (i = 0; i < $( '.new_input' ).length; i++)
                       {
                       list_value.push($( '.new_input' )[i].value)
                       }
                       Shiny.onInputChange('newValuepHyper', list_value)
                       });")),
        hr(),
      uiOutput('hyper_test_gene_table_results')
      )
  ,size="l",easyClose=TRUE)
  
  modal_hyper_test_gene_wrong_n <- modalDialog(
    fluidPage(HTML("Select exactly two gene lists, the first selected should be the list of all genes with a particular feature, the second the list of genes you wish to see if are enriched for this feature.")),size="s",easyClose=TRUE
  )
  
  #Inputs for hyper test
  output$hyper_test_gene_table<-renderDataTable({
    selected_row = input$gene_list_table_rows_selected
    # DT[[1]] = numericInput("background_gene_n", label = "Population size", value = data_cache$gene_lists_DT[(selected_row[1])][[3]])
    # DT[[2]] = numericInput("background_gene_with_n", label = "Number of successes in population", value = 0)
    # DT[[3]] = numericInput("query_gene_n", label = "Sample size", value = 0)
    # DT[[4]] = numericInput("query_gene_with_n", label = "Number of successes in sample (x)", value = 0)
    
    # data.object = modifyData()
    # 
    # format = data.object$data_format #user_input$format
    # format_details = formatLookup(format)
    # gene_id_position = format_details$gene_id_position
    # correct_colnames = format_details$correct_colnames
    # gene_id_name = correct_colnames[gene_id_position]
    # 
    # data = data.object$modified_data
    # 
    # background = length(unique(data[gene_id_name][[1]],incomparables=c(NA,NULL,"-")))
    
    background = hyper_test_input_background()
    
    #calculate number of successes
    sucess_genes = data_cache$gene_lists[[(selected_row[1])]]
    query_genes = data_cache$gene_lists[[(selected_row[2])]]
    query_sucess_genes = intersect(query_genes,sucess_genes)
    query_successes = length(unique(query_sucess_genes))
    
    DT=list()
    DT[[1]] = paste0('<input class="new_input" type="number" id=new_',1,', value=',background,'><br>')
    DT[[2]] = paste0('<input class="new_input" type="number" id=new_',2,', value=',data_cache$gene_lists_DT[(selected_row[1])][[3]],'><br>')
    DT[[3]] = paste0('<input class="new_input" type="number" id=new_',3,', value=',data_cache$gene_lists_DT[(selected_row[2])][[3]],'><br>')
    DT[[4]] = paste0('<input class="new_input" type="number" id=new_',4,', value=',query_successes,'><br>')
    
    DT = rbind(DT) #matrix(NA,ncol=4)
    colnames(DT) = c("Population size","Number of successes in population","Sample size","Number of successes in sample (x)")
    
    return(DT)
  },escape=F,selection="none",rownames = FALSE,options=list(dom='t',ordering=F))
  
  #Results of hyper test
  output$hyper_test_gene_table_results<-renderText({
    background_gene_n = user_input$background_gene_n
    background_gene_with_n = user_input$background_gene_with_n
    query_gene_n = user_input$query_gene_n
    query_gene_with_n = user_input$query_gene_with_n
    
    if (is.numeric(c(background_gene_n,background_gene_with_n,query_gene_n,query_gene_with_n))){
      if(query_gene_with_n > query_gene_n){
        return("Number of successes in sample cannot be greater than Sample size")
      }
      if(background_gene_with_n > background_gene_n){
        return("Number of successes in population cannot be greater than Population size")
      }
      if(query_gene_n > background_gene_n){
        return("Sample size cannot be greater than Population size")
      }
      
      results = 1.0-phyper(query_gene_with_n-1, background_gene_with_n, background_gene_n-background_gene_with_n, query_gene_n)
      results = paste0("p-value: ",results)
    }else{
      results = NULL
    }
    return(results)
  }) #,escape=F,selection="none",rownames = FALSE,options=list(dom='t',ordering=F)
  
  #Show hyper enrichment popup
  observeEvent(input$hyper_test_gene,{
    selected_row = input$gene_list_table_rows_selected 
    if (!is.null(selected_row)){
      if(length(selected_row) == 2){
        showModal(modal_hyper_test_gene)
      }else{
        showModal(modal_hyper_test_gene_wrong_n)
      }
    }
  })
  
  #Update with the new values
  observeEvent(input$newValuepHyper,{
    newValue=lapply(input$newValuepHyper, function(col) {
      if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
        as.numeric(as.character(col))
      } else {
        col
      }
    })
    DF=data.frame(lapply(newValue, function(x) t(data.frame(x))))
    
    user_input$background_gene_n = DF[,1]
    user_input$background_gene_with_n = DF[,2]
    user_input$query_gene_n = DF[,3]
    user_input$query_gene_with_n = DF[,4]
  })
  
  #Given: 
  #Returns: 
  #Requires reactive variables: input$gene_list_table_rows_selected, data_cache$gene_lists, data_cache$gene_lists_metadata
  output$save_gene_list_button <- downloadHandler(
    filename = function() {
      gene_lists = (data_cache$gene_lists)
      metadata = (data_cache$gene_lists_metadata)
      if (is.null(gene_lists)|is.null(input$gene_list_table_rows_selected)){
        return(NULL)
      }else{
        rows = input$gene_list_table_rows_selected
        rows = rows[length(rows)]
        rootname = metadata[[rows]]$gene_list_name
        
        name = paste0(rootname,"-", Sys.Date(), ".tsv")
        return(name)
      }
    },
    content = function(file) {
      gene_lists = (data_cache$gene_lists)
      if (is.null(gene_lists)){
        return(NULL)
      }else if (is.null(input$gene_list_table_rows_selected)){
        return(NULL)
      }else{
        rows = input$gene_list_table_rows_selected
        rows = rows[length(rows)]
        
        metadata = (make_gene_list_table())[rows,]
        gene_list = gene_lists[[rows]]
        data = as.matrix(c("METADATA",metadata,"GENELIST",gene_list))
        
        write.table(data, file,sep="\t",row.names=F,quote=F,col.names = F)
      }
    },contentType = "text/csv"
  )

  #Venn diagram stuff
  modal_venn_gene <- modalDialog(
    fluidPage(
      h3(strong("Venn diagram"),align="center"),
      
      # dataTableOutput('hyper_test_gene_table'),
      # hr(),
      # actionButton("calculate_phyper","Calculate enrichment"),
      # tags$script(HTML("$(document).on('click', '#calculate_phyper', function () {
      #                  var list_value=[]
      #                  for (i = 0; i < $( '.new_input' ).length; i++)
      #                  {
      #                  list_value.push($( '.new_input' )[i].value)
      #                  }
      #                  Shiny.onInputChange('newValuepHyper', list_value)
      #                  });")),
      # hr(),
      uiOutput('venn_plot_ui')
    )
    ,size="l",easyClose=TRUE)
  
  output$venn_plot_ui <- renderUI({
    out = list(plotOutput("plotVenn"))
    out[[2]] = (fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Venn diagram in svg format",
                                                       downloadButton("downloadVenn", "Download Venn diagram")))))
    return(out)
  })
  output$plotVenn <- renderPlot({
    plotVenn_plot()
  })
  plotVenn_plot <- function(){
    # library(VennDiagram)
    # library(ascii)
    # library(eulerr)
    # 
    # 
    # # list_names = c("Name1","N2","n3")#letters[1:n_lists]
    # # a = c(1:3,8,10)
    # # b = c(2:5,9,10)
    # # c = c(1,4:7,10)
    # # all_together = list(a,b,c)
    
    #load data
    selected_rows = input$gene_list_table_rows_selected
    all_together = list(as.character(data_cache$gene_lists[[(selected_rows[1])]]))
    for (i in 2:length(selected_rows)){
      all_together[[i]] = as.character(data_cache$gene_lists[[(selected_rows[i])]])
    }
    all_names = geneListNames()
    list_names = all_names[selected_rows]
    names(all_together) = list_names
    
    # data_cache$test = selected_rows
    # data_cache$test2 = all_together
    # 
    # overlap = calculate.overlap(all_together)
    # sizes = sapply(overlap, length)
    # 
    # posnam = paste(paste0(combn(list_names,length(list_names))),collapse="&")
    # for (i in (length(list_names)-1):1){
    #   posnam = c(posnam,paste.matrix(combn(list_names,i),collapse="&",byrow=TRUE))
    # }
    # names(sizes) = posnam
    # 
    # vd = euler(sizes)
    # plot(vd,col="black")
    
    venn(all_together,names = list_names)
  }
  
  output$downloadVenn <- downloadHandler(
    filename = function() {
      name = paste0("Venn_diagram-", Sys.Date(), ".svg")
      return(name)
    },
    content = function(file) {
      # if (is.null(make_gprofile_enrichment_table())){
      #   return(NULL)
      # }else{
      #ensg_gene_names = make_gprofile_enrichment_table()$genes_used
      update = plotVenn_plot()
      # max_genes = data_cache$max_genes# length(unique(ensg_gene_names))
      # vert_size = 0.15 * max_genes
      # vert_size = ceiling(vert_size) + 8
      # # if (max_genes > 50){
      # #   svg(file,height=28,width=10)
      # # }else{
      # #   svg(file,height=14,width=10)
      # # }
      svg(file)
      (plotVenn_plot())
      dev.off()
      #   }
    }
  )
  
  
  #Show Venn diagram popup
  observeEvent(input$venn_gene,{
    selected_row = input$gene_list_table_rows_selected
    if (!is.null(selected_row)){
      if(length(selected_row) > 1){
        showModal(modal_venn_gene)
        #   }else{
        #     showModal(modal_hyper_test_gene_wrong_n)
      }
    }
  })
  
  
  

  
  
#Gene list genes tab
#Given: table from output$genelistgenestable
#Returns: this as a datatable
#Requires reactive functions: geneListNames()
#Requires reactive variables: 
output$gene_list_genes_table <- renderUI({
  names = geneListNames()
  if (is.null(names)){
    return(NULL)
  }else{
    names = c("None",names)
  }
  
  out = list(tags$div(title=paste0("Select a gene list to display data for those genes; gene lists are created from the Gene lists menu."),
                      selectInput("select_gene_list2", label = ("Show data for genes in this list"), 
                                  choices = names, 
                                  selected = 1)))
  out[[2]] = (dataTableOutput("genelistgenestable"))
  
  return(out)
})  

# Display gene list genes  
#Given: table and lookup from filterFileData()
#Returns: this as a datatable
#Requires reactive functions: filterFileData()
#Requires reactive variables: user_input$select_gene_list2
output$genelistgenestable <- DT::renderDataTable({
  if(input$select_gene_list2 == "None"){ #if no genelist or nonselected
    return(NULL)
  }else{
    data = filterFileData()
    lookup = data$lookup
    data = data.frame(x=data$table)
    
    #set expression values to numbers
    if (is.null(lookup)){
      for (i in 3:dim(data)[2]){
        data[,i] = as.numeric(as.character(data[,i])) 
      }
    }else{
      for (i in 4:dim(data)[2]){
        data[,i] = as.numeric(as.character(data[,i])) 
      }
    }
    
    return(data)
  }
}, escape = FALSE, extensions = 'FixedHeader',
options = list(pageLength = 25, order = list(list(2, 'desc'), list(1, 'asc'))
               , orderClasses = TRUE,searchHighlight = TRUE,
               fixedHeader = TRUE, scrollX = TRUE),selection = 'single',rownames = FALSE)






##FUNCTIONS
#Load gene_list file
#Given: file (input$gene_list_file), input$create_gene_list_name
#Returns reactive variables: data_cache$gene_lists, data_cache$gene_lists_metadata
#Returns NULL if no file, else list of 
#  list of vectors of metadata (metadata),
#  vector of genes (new_genes)
geneListFile <- reactive({
  inFile <- input$gene_list_file
  if (is.null(inFile)){ #No file selected
    return(NULL)
  }
  data = read.csv(inFile$datapath,as.is=1,header=F,sep="\t")
  
  #metadata?
  if((data)[1,] == "METADATA"){
    last = match("GENELIST",data[,1])
    metadata = as.list(data[3:(last-1),])
    
    #backwards compatibility
    if(length(metadata) == 12){
      metadata[[13]] = NA
    }
    
    metadata$gene_list_name = data[2,] #must be last element added for output$gene_list_table to work
    
    new_genes = data[((last+1):dim(data)[1]),1]
    new_genes = unique(new_genes)
  }else{
    metadata = as.list(rep("NULL",13))
    metadata$gene_list_name = isolate(input$create_gene_list_name) #must be last element added for output$gene_list_table to work
    
    new_genes = data[,1]
    new_genes = unique(new_genes)
    
    #update metadata
    metadata[[2]] = length(new_genes)
    metadata[[3]] = input$gene_list_file$name
  }
  
  return(list(new_genes = new_genes, metadata = metadata))
})

#Update/Retrieve gene lists  
#Given: list of
#  list of metadata, vector of genes from geneListFile()
#Returns: NULL if geneListFile is NULL or list of
#  list of vector of genes as above with new one added (gene_lists);
#  list of list of vectors of metadata (gene_lists_metadata)
#Returns reactive variables: data_cache$gene_lists, data_cache$gene_lists_metadata
#Requires reactive functions: geneListFile()
#Requires reactive variables: data_cache$gene_lists, data_cache$gene_lists_metadata
queryGeneLists <- reactive({
  if (is.null(geneListFile())){ #No file selected
    return(NULL) 
  }else{
    input = geneListFile()
    gene_lists = isolate(data_cache$gene_lists) #isolate necessary to prevent recursion
    gene_lists_metadata = isolate(data_cache$gene_lists_metadata)
    
    #former combineAnnotationSets
    if(is.null(gene_lists)){ #nothing to add
      gene_lists = list(input$new_genes)
      gene_lists_metadata = list(input$metadata)
    }else{
      gene_lists = append(gene_lists,list(input$new_genes))
      gene_lists_metadata = append(gene_lists_metadata,list(input$metadata))        
    }
    
    data_cache$gene_lists = gene_lists
    data_cache$gene_lists_metadata = gene_lists_metadata
    
    output = list('gene_lists' = gene_lists, 'gene_lists_metadata' = gene_lists_metadata)
    return(output)
  }
})

#Make table of gene lists
#Requires reactive variables: data_cache$gene_lists, data_cache$gene_lists_metadata
#Returns reactive variables: data_cache$gene_lists_DT
make_gene_list_table <- reactive ({
  update = queryGeneLists()
  gene_lists = data_cache$gene_lists
  metadata = data_cache$gene_lists_metadata
  if (is.null(gene_lists)){
    return(NULL)
  }else{
    data = matrix(NA,ncol=length(metadata[[1]]))
    
    for (i in 1:length(data_cache$gene_lists)){
      row = metadata[[i]]$gene_list_name
      for (m in 1:(length(metadata[[i]])-1)){
        row = c(row,(metadata[[i]])[[m]])
      }
      data = rbind(data,row)
    }
    data = data[-1,]
    if (is.null(dim(data))){
      data = t(as.matrix(data))
    }
    colnames(data) = c("List name","Species","n genes","source file name","Comparison","Sample 1 Exp >=","Sample 2 Exp >=",
                       "Fold change threshold","Sig. DE only","In other comparisons","Not in other comparisons","Gene list filter","Description","Other loaded files")
    

    data_cache$gene_lists_DT<-data.table(data)
    return(data)
  }
})

#Gene lists names
#Given reactive variables: data_cache$gene_lists_metadata()
#Returns: if data_cache$gene_lists_metadata = NULL returns NULL, else vector of gene list names
geneListNames <- reactive({
  metadata = (data_cache$gene_lists_metadata)
  if (is.null(metadata)){
    return(NULL)
  }else{
    names = 0
    for (i in 1:length(metadata)){
      names = c(names,metadata[[i]]$gene_list_name)
    }
    names = names[-1]
    return(names)
  }
})

#Create a gene list from current settings
#Returns: NULL
#Returns reactive variables: data_cache$gene_lists, data_cache$gene_lists_metadata
#Requires reactive functions: createGeneListInput()
createGeneList <- reactive({
  if (is.null(createGeneListInput())){
    return(NULL)
  }else{
    data = createGeneListInput()
    gene_names = data$gene_names
    metadata = data$metadata
    
    if (is.null(data_cache$gene_lists)){
      data_cache$gene_lists = list(gene_names)
      data_cache$gene_lists_metadata = list(metadata)
    }else{
      data_cache$gene_lists[[(length(data_cache$gene_lists) + 1)]] = (gene_names)
      data_cache$gene_lists_metadata[[(length(data_cache$gene_lists_metadata) + 1)]] = (metadata)
    }
    return(NULL)
  }
})


