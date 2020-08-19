library(shiny)
library(shinyjs) #error messages
#Appearance
library(shinydashboard)
library(DT)
library(data.table) #editable datatable
library(plotly)
library(gplots) #expression heatmap
library(svglite) #save isoform svg
library(ComplexHeatmap) #annotation heatmap
#Web annotation
library(httr)
library(jsonlite)
library(RCurl)
library(InterMineR) #getAnnot
library(gProfileR) #getAnnot, gProfile enrichment, orthologs
#Other
#library(prodlim) # for row.match - turns out this scales poorly
library(eulerr) #for venn diagrams

#My functions
source("crispAnnotationCollection.R")

#increase max filesize
options(shiny.maxRequestSize=400*2024^2)

#sort graph plotting font
options(bitmapType="cairo")

function(input, output, session) {
  #My functions specific to session
  source("gene_lists_UI+functions.R",local=T)
  source("various_enrichment.R",local=T)
  
  #Test outputs
  output$value <- renderPrint({ data_cache$test })
  output$value2 <- renderPrint({ data_cache$test2 })
  output$value3 <- renderPrint({ input$select_comp })
  output$value4 <- renderPrint({ user_input$select_comp })
  
  #Collapsible sidebar
  runjs({'
        var el2 = document.querySelector(".skin-blue");
        el2.className = "skin-blue sidebar-mini";
        var clicker = document.querySelector(".sidebar-toggle");
        clicker.id = "switchState";
    '})
  onclick('switchState', runjs({'
        var title = document.querySelector(".logo")
        if (title.style.visibility == "hidden") {
          title.style.visibility = "visible";
        } else {
          title.style.visibility = "hidden";
        }
  '}))
  
  
  ##REACTIVE VARIABLES
  #Cache data - includes webdata (to save reloading)... 
  data_cache = reactiveValues()
  #Collect of user inputs sanitised/altered by other functions
  user_input = reactiveValues()  
  #Store annotations
  global_annotation = reactiveValues()
  
  ##UI
  #USER INPUT INITIALISATION - except annotations
  #Collect User input and sanitise
  #Given: user input
  #Returns: user input as reactive variables dependent upon various calculations
  #Returns reactive variables: most parts of user_input
  #Requires reactive variables: data_cache$lookup, data_cache$names_list
  #Requires functions: speciesNames()
  collectUser <- reactive({
    #Initialise file format 
    user_input$format = input$format
    
    # #Reset user_input$table_rows_all if made any changes to filters when not on table tab - otherwise summary/download not accurate till switch back to table tab
    # reset = resetTableRowsAll(1)
    
    #Initialise and sanitise numerical inputs 
    #- prevents errors when you delete the values when changing them
    user_input$s1_fpkm_num = input$s1_fpkm_num
    user_input$s2_fpkm_num = input$s2_fpkm_num
    
    if (is.null(input$s1_fpkm_num)){
      user_input$s1_fpkm_num = 0;
    }else{
      if (is.na(input$s1_fpkm_num)){
        user_input$s1_fpkm_num = 0;
      }
    }
    if (is.null(input$s2_fpkm_num)){
      user_input$s2_fpkm_num = 0;
    }else{
      if (is.na(input$s2_fpkm_num)){
        user_input$s2_fpkm_num = 0;
      }
    }
    if (is.null(input$log2FC_num)){
      user_input$log2FC_num = 0;
    }else{
      if (is.na(input$log2FC_num)){
        user_input$log2FC_num = 0;
      }else{
        user_input$log2FC_num = abs(input$log2FC_num)
      }
    }
    
    #initialise and sanitise checkboxes
    if (is.null(input$sig_checkbox)){
      user_input$sig_checkbox = FALSE
    }else{
      user_input$sig_checkbox = input$sig_checkbox
    }
    if (is.null(input$altsplice_checkbox)){
      user_input$altsplice_checkbox = FALSE
    }else{
      user_input$altsplice_checkbox = input$altsplice_checkbox
    }
    if (is.null(input$log2FC_up_checkbox)){
      user_input$log2FC_up_checkbox = TRUE
    }else{
      user_input$log2FC_up_checkbox = input$log2FC_up_checkbox
    }
    if (is.null(input$log2FC_down_checkbox)){
      user_input$log2FC_down_checkbox = TRUE
    }else{
      user_input$log2FC_down_checkbox = input$log2FC_down_checkbox
    }
    if (is.null(input$invertFCcomp_checkbox)){
      user_input$invertFCcomp_checkbox = FALSE
    }else{
      user_input$invertFCcomp_checkbox = input$invertFCcomp_checkbox
    }
    
    if (is.null(input$select_gene_list)){
      user_input$select_gene_list = "None"
    }else{
      if (input$select_gene_list == "Create a gene list to filter by"){
        user_input$select_gene_list = "None"
      }else{
        user_input$select_gene_list = input$select_gene_list
      }
    }
    if (is.null(input$select_gene_list2)){
      user_input$select_gene_list2 = "None"
    }else{
      user_input$select_gene_list2 = input$select_gene_list2
    }
    
    if (is.null(input$invertfilter_file_checkbox)){
      user_input$invertfilter_file_checkbox = FALSE
    }else{
      user_input$invertfilter_file_checkbox = input$invertfilter_file_checkbox
    }
    if (is.null(input$clearannot_file_checkbox)){
      user_input$clearannot_file_checkbox = FALSE
    }else{
      user_input$clearannot_file_checkbox = input$clearannot_file_checkbox
    }
    
    user_input$select2InputFoundIn = input$select2InputFoundIn
    user_input$select2InputNotFoundIn = input$select2InputNotFoundIn
    
    #initialise select comparison
    if (is.null(input$select_comp)){
      choices = comparisonNames() #not run before now in this case
      if (length(choices) > 0){
        user_input$select_comp = choices[1]
      }
    }else if(nchar(input$select_comp) < 1){ #occurs when Filter menu closed
      choices = comparisonNames() #not run before now in this case
      if (length(choices) > 0){
        user_input$select_comp = choices[1]
      }     
    }else{
      user_input$select_comp = input$select_comp
    }
    
    user_input$species = input$species
    
    
    #set up 'built-in' annotation checkboxes
    species_data = speciesNames(user_input$species)
    if (is.null(species_data)){
      data_cache$built_in_checkboxes_values = NULL
    }else{
      names = species_data$built_in_list
      if (is.null(names)){
        data_cache$built_in_checkboxes_values = NULL
      }else{
        #vectorise the following
        checkboxes_values = eval(parse(text = paste0("input$",names[1])))
        if (length(names) > 1){
          for (i in 2:length(names)){
            checkboxes_values = c(checkboxes_values,eval(parse(text = paste0("input$",names[i]))))
          }
        }
        #data_cache$built_in_checkboxes_values = checkboxes_values
        
        #initialise where annotation menu tab not open
        if(is.null(checkboxes_values)){
          checkboxes_values = rep(FALSE,length(names))
          #names(checkboxes_values) = names
          #checkboxes_values = as.list(checkboxes_values)
        }
        data_cache$built_in_checkboxes_values = checkboxes_values
      }
    }
    
    
    
    #For and selected on per gene plot
    user_input$select_DE = input$select_DE
    user_input$remove_max = input$remove_max
    if(!is.null(input$rotate_sample)){
      user_input$rotate_sample = input$rotate_sample
    }
    #user_input$rotate_sample = input$rotate_sample
    user_input$gene_exp_log_scale = input$gene_exp_log_scale
    user_input$gene_exp_log_scale_lt_zero = input$gene_exp_log_scale_lt_zero
    user_input$pergene_cexx = input$pergene_cexx
    if(!is.null(input$heatmap_log_scale)){
      user_input$heatmap_log_scale = input$heatmap_log_scale
    }
    if(!is.null(input$heatmap_no_rearrange)){
      user_input$heatmap_no_rearrange = input$heatmap_no_rearrange
    }
    if(!is.null(input$heatmap_all_genes)){
      user_input$heatmap_all_genes = input$heatmap_all_genes
    }
    if(!is.null(input$heatmap_cexCol)){
      user_input$heatmap_cexCol = input$heatmap_cexCol
    }
    if(!is.null(input$heatmap_cexRow)){
      user_input$heatmap_cexRow = input$heatmap_cexRow
    }
    if(!is.null(input$volcano_plot_grey_data)){
      user_input$volcano_plot_grey_data = input$volcano_plot_grey_data
    }
    if(!is.null(input$volcano_plot_no_rearrange)){
      user_input$volcano_plot_no_rearrange = input$volcano_plot_no_rearrange
    }
    if(!is.null(input$scatter_plot_log)){
      user_input$scatter_plot_log = input$scatter_plot_log
    }
    if(!is.null(input$scatter_plot_line)){
      user_input$scatter_plot_line = input$scatter_plot_line
    }
    if(!is.null(input$pairwise_plot_original)){
      user_input$pairwise_plot_original = input$pairwise_plot_original
    }
    if(!is.null(input$pairwise_plot_line)){
      user_input$pairwise_plot_line = input$pairwise_plot_line
    }
    
    #Initialise
    if (is.null(input$pergene_text)){
    }else{
      if (nchar((input$pergene_text)[1]) > 0){
        names_list = data_cache$names_list
        trans_names_list = data_cache$trans_names_list
        user_input$pergene_text = trans_names_list[match(input$pergene_text,names_list)]
        
      }else{
        if (is.null(user_input$pergene_text)){
          user_input$pergene_text = ""
        }
      }
    }
    
    if (!is.null(input$select2InputSampleOrder)){
      user_input$select2InputSampleOrder = input$select2InputSampleOrder
    }
    if (!is.null(input$select2InputSampleOrder_heatmap)){
      user_input$select2InputSampleOrder_heatmap = input$select2InputSampleOrder_heatmap
    }
    
    #Update the column checkboxes if they exist/are shown
    if (is.null(input$col1_checkbox)){}else{
      user_input$col1_checkbox = input$col1_checkbox
    }
    if (is.null(input$col2_checkbox)){}else{
      user_input$col2_checkbox = input$col2_checkbox
    }
    if (is.null(input$col3_checkbox)){}else{
      user_input$col3_checkbox = input$col3_checkbox
    }
    if(is.null(data_cache$lookup)){
    }else{
      #if(!(data_cache$lookup)){
      if (is.null(input$col4_checkbox)){}else{
        user_input$col4_checkbox = input$col4_checkbox
      }
      #}
    }
    if (is.null(input$col5_checkbox)){}else{
      user_input$col5_checkbox = input$col5_checkbox
    }
    if (is.null(input$col6_checkbox)){}else{
      user_input$col6_checkbox = input$col6_checkbox
    }
    if (is.null(input$col7_checkbox)){}else{
      user_input$col7_checkbox = input$col7_checkbox
    }
    if (is.null(input$col8_checkbox)){}else{
      user_input$col8_checkbox = input$col8_checkbox
    }
    if (is.null(input$col9_checkbox)){}else{
      user_input$col9_checkbox = input$col9_checkbox
    }
    if (is.null(input$col10_checkbox)){}else{
      user_input$col10_checkbox = input$col10_checkbox
    }
    if (is.null(input$col11_checkbox)){}else{
      user_input$col11_checkbox = input$col11_checkbox
    }
    if (is.null(input$col12_checkbox)){}else{
      user_input$col12_checkbox = input$col12_checkbox
    }
    if (is.null(input$col13_checkbox)){}else{
      user_input$col13_checkbox = input$col13_checkbox
    }
    if (is.null(input$col14_checkbox)){}else{
      user_input$col14_checkbox = input$col14_checkbox
    }
    if (is.null(input$col15_checkbox)){}else{
      user_input$col15_checkbox = input$col15_checkbox
    }
    
    return(user_input)
  })  
  
  #Initialise species
  user_input$species = "Human"
  
  #Initialise column checkboxes - need this because of ability to hide the options for selecting them - it starts hidden
  user_input$col1_checkbox = TRUE
  user_input$col2_checkbox = TRUE
  user_input$col3_checkbox = TRUE
  user_input$col4_checkbox = TRUE
  user_input$col5_checkbox = TRUE
  user_input$col6_checkbox = TRUE
  user_input$col7_checkbox = TRUE
  user_input$col8_checkbox = TRUE
  user_input$col9_checkbox = TRUE
  user_input$col10_checkbox = TRUE
  user_input$col11_checkbox = TRUE
  user_input$col12_checkbox = TRUE
  user_input$col13_checkbox = TRUE
  user_input$col14_checkbox = TRUE
  user_input$col15_checkbox = TRUE
  
  #Initialise checkboxes+values in per gene plot+figures, need this because not shown till plot generated,
  #  but wanted for plot generation
  user_input$remove_max = FALSE
  user_input$rotate_sample = FALSE
  user_input$gene_exp_log_scale = FALSE
  user_input$gene_exp_log_scale_lt_zero = FALSE
  user_input$heatmap_log_scale = FALSE
  user_input$heatmap_no_rearrange = TRUE
  user_input$heatmap_all_genes = FALSE
  user_input$heatmap_cexCol = 1
  user_input$heatmap_cexRow = 0.8
  user_input$pergene_cexx = 1 #10 if transcript
  user_input$volcano_plot_grey_data = FALSE
  user_input$volcano_plot_no_rearrange = FALSE
  user_input$scatter_plot_log = TRUE
  user_input$scatter_plot_line = TRUE
  user_input$pairwise_plot_original = FALSE
  user_input$pairwise_plot_line = TRUE
  
  #For gene lists
  user_input$from_enrich = FALSE
  user_input$from_enrich_selected = FALSE
  user_input$from_enrich_ipd = FALSE
  user_input$from_enrich_selected_ipd = FALSE
  user_input$from_ortholog = FALSE
  user_input$from_ortholog_selected = FALSE
  
  user_input$background_gene_n = NA
  user_input$background_gene_with_n = NA
  user_input$query_gene_n = NA
  user_input$query_gene_with_n = NA
  
  
  #Reused text
  output$need_file_text = renderText({"Select valid file for analysis..."})
  
  #Various page text
  output$gene_not_found <- renderText({"Gene not found in dataset..."}) 
  output$plot_text <- renderText({paste0("Points with adjusted p-value < ",pAdjThreshold()," are coloured red, points with log2FoldChange > 1 are coloured orange and points with both are coloured green")})
  output$plot2_text <- renderText({
    if(user_input$select_DE == ""){
      text = "Select at least one sample to display graph."
    }else{
      text = paste0("Conditions significantly DE w.r.t. ",user_input$select_DE," are indicated with an asterix.")
    }
    return(text)
  })
  output$plot2_text2 <- renderText({
    if(user_input$select_DE == ""){
      return(NULL)
    }
    text = paste0("Conditions not compared to ",user_input$select_DE," are coloured grey.")
    return(text)
  })
  output$plot2_text2_isoform <- renderText({paste0("Warning not all conditions may have been compared to ",user_input$select_DE)})
  output$plot3_text <- renderText({paste0("Genes/transcripts from the currently selected data that are differently expressed and whose log2FoldChange is not infinite are clustered by their log2FoldChange.")})
  output$plot4_text <- renderText({paste0("Genes/transcripts from the currently selected data are clustered by their expression in all samples. Where more than 50 genes are selected, the 50 genes with the highest log2FoldChange in the current comparison are shown.")})
  output$noplot3_text <- renderText({"Not enough data to plot."})
  output$noplot4_text <- renderText({"Not enough data to plot."})
  output$plot7_text <- renderText({"R value calculated using Spearman correlation coefficient."})
  
  output$intermine_not_found <- renderText({"Could not connect to InterMine server..."}) #alter this to be species specific
  output$choose_gene <- renderText({"Select a gene from the above list."}) 
  output$species_specific_not_found <- renderText({
    species = user_input$species
    data = speciesNames(species)
    site = data$database
    
    out = paste0("Could not connect to ",site," server...")
    return(out)
  }) 
  
  #Put things in here that you want to happen no matter what menu option is selected,
  #or if none are selected
  #Initialise and update user_input 
  output$user_update <- renderUI({
    ui = collectUser() #needed to initialise and update user_input
    
    out = list(NULL)
    return(out)
  })
  
  #Menu options
  output$input_selection <- renderUI({
    out = list(tags$div(title="Select which species data is from",
                        selectInput("species", label = h4("Select species"), 
                                    choices = speciesNames("LIST"),
                                    selected = isolate(user_input$species))
    ))
    out[[2]] = tags$div(title="Select file format",
                        selectInput("format", label = h4("Select file format"), 
                                    choices = list("Cufflinks gene_exp.diff" = "cufflinks_gene",
                                                   "Cufflinks isoform_exp.diff" = "cufflinks_transcript",
                                                   "PRAGUI DESeq results gene" = "pragui_gene",
                                                   "PRAGUI DESeq results isoform" = "pragui_transcript"), 
                                    selected = isolate(user_input$species))
    )
    out[[3]] = tags$div(title="Click here to select a file for analysis",
                        fileInput("file", label = h4("Input file"))
    )
    out[[4]] = br()
    
    return(out)
  })
  
  output$filter_selection <- renderUI({
    format = user_input$format
    format_details = formatLookup(format)
    exp_unit = format_details$exp_unit
    
    out = list(tags$div(title="Select a comparison to display in Summary and Table tabs",
                        uiOutput("select_comparisons")
    ))
    out[[2]] = uiOutput("significance_selection")
    out[[3]] = uiOutput("isoform_selection")
    out[[4]] = numericInput("s1_fpkm_num", label = paste0("Sample 1 ",exp_unit," at least"), value = 0)
    out[[5]] = numericInput("s2_fpkm_num", label = paste0("Sample 2 ",exp_unit," at least"), value = 0)
    out[[6]] = tags$div(title="Takes the absolute value of this (i.e. -1 = 1); Following checkboxes control which direction(s) this filter is applied in",
                        numericInput("log2FC_num", label = "Log2 fold-change minimum", value = 0)
    )
    out[[7]] = fluidRow(
      column(3, 
             tags$div(title="Up regulated in Sample 2 (Log2 fold-change >=)",
                      checkboxInput("log2FC_up_checkbox", label = "Up", value = TRUE)
             )
      ),
      column(3, 
             tags$div(title="Down regulated in Sample 2 (Log2 fold-change <=)",
                      checkboxInput("log2FC_down_checkbox", label = "Down", value = TRUE)
             )
      )
    )
    out[[8]] = br()
    out[[9]] = fluidRow(column(12, tags$div(title="Filter so displayed genes are also found in other comparisons - where multiple comparisons are chosen the displayed genes must be found in ALL comparisons; delete to remove",
                                            uiOutput("comparison_found_in")
    )))
    out[[10]] = tags$div(title="Change the filter so displayed genes must have a fold-change in the opposite direction",
                        checkboxInput("invertFCcomp_checkbox", label = "Fold-change in opposite direction in comparison", value = FALSE)
    )
    out[[11]] = fluidRow(column(12, tags$div(title="Filter so displayed genes must not be found in any of the comparisons listed - where multiple comparisons are chosen the displayed genes must NOT be found in ANY of the comparisons; delete to remove",
                                             uiOutput("comparison_not_found_in")
    )))
    out[[12]] = br()
    out[[13]] = uiOutput("filter_file_selection")
    out[[14]] = tags$div(title="When selected displayed genes must not be found in the chosen list",
                         checkboxInput("invertfilter_file_checkbox", label = "Not in this list", value = FALSE)
    )
    out[[15]] = br()
    
    return(out)
  })
    #List of possible comparisons in dataset
    #Given: data object from modifyData()
    #Returns: if modifyData() = NULL returns NULL, else vector of comparisons possible in 
    #         modified_data from modifyData
    #Returns reactive variables: user_input$allComparisonNames
    #Requires reactive functions: modifyData()
    comparisonNames <- reactive({
      if (is.null(modifyData())){ #Prior to loading a dataset
        return(NULL)
      }else{
        data.object = modifyData()
        temp_raw = data.object$modified_data
        
        format = data.object$data_format
        format_details = formatLookup(format)
        condition_1_position = format_details$condition_1_position
        condition_2_position = format_details$condition_2_position
        
        if(is.na(condition_1_position)){ #does not give condition names, assume that all data from one comparison
          comps_names = "sample_1 vs sample_2"
        }else{
          correct_colnames = format_details$correct_colnames
          condition_1_name = correct_colnames[condition_1_position]
          condition_2_name = correct_colnames[condition_2_position]
          
          comparisons = unique(temp_raw[c(condition_1_name,condition_2_name)])
          comps_names = paste(comparisons[condition_1_name][[1]],"vs",comparisons[condition_2_name][[1]],sep=" ")
        }
        
        user_input$allComparisonNames = comps_names
        
        return(comps_names)
      }
    })
    #Make list of possible comparisons in dataset into UI element
    #Given: vector of comparisons from comparisonNames()
    #Returns: UI element selectbox with first comparison selected
    #Requires reactive functions: comparisonNames()
    output$select_comparisons <- renderUI({
      out = list(selectInput("select_comp", label = h4("Select comparison"), 
                             choices = comparisonNames(), 
                             selected = 1))
      return(out)
    })
    #What adjusted p-value threshold was used for significance  
    #Given: data.object from modifyData()
    #Returns: 0.05 if modifyData() = NULL, else integer of the estimated q-value threshold used
    #Requires reactive functions: modifyData()
    pAdjThreshold <- reactive ({
      if (is.null(modifyData())){ #Prior to loading a dataset
        return(0.05)
      }else{
        format = user_input$format
        format_details = formatLookup(format)
        correct_colnames = format_details$correct_colnames
        
        adjusted_p_position = format_details$adjusted_p_position
        if (is.na(adjusted_p_position)){
          return(0.05)
        }else{
          sig_position = format_details$sig_position
          sig_name = correct_colnames[sig_position]
          adjusted_p_name = correct_colnames[adjusted_p_position]
          
          data.object = modifyData()
          data = data.object$modified_data
          if (dim((data[(data[sig_name] == 'yes'),])[adjusted_p_name])[1] > 0){
            maximum_value_sig = max((data[(data[sig_name] == 'yes'),])[adjusted_p_name],na.rm=T)
            digits = 1 #how many decimal places to display
            display_value = round(maximum_value_sig,digits)
            while (display_value == 0){ #if using very low threshold then need to increase number of decimal places to avoid rounding to zero
              digits = digits + 1
              display_value = round(maximum_value_sig,digits)
            }
          }else{
            return(0.05)
          }
          return(display_value)
        }
      }
    })
    #Given: integer of adjusted p-value threshold from pAdjThreshold()
    #Returns: UI element of checkbox to filter on threshold
    #Requires reactive functions: pAdjThreshold()
    output$significance_selection <- renderUI({
      pAdj = pAdjThreshold()
      
      out = list(tags$div(title=paste0("Display only those genes which are significantly differentially expressed (adjusted p-value < ",pAdj,")"),
                          checkboxInput("sig_checkbox", label = "Significantly DE", value = FALSE)
      ))
      
      return(out)
    })
    #Given: integer of adjusted p-value threshold from pAdjThreshold(), vector of data format (user_input$format)
    #Returns: if data format includes isoform IDs UI element of checkbox to only show genes with multiple isoforms
    #Requires reactive functions: pAdjThreshold()
    output$isoform_selection <- renderUI({
      pAdj = pAdjThreshold()
      format = user_input$format
      format_details = formatLookup(format)
      
      isoform_id_position = format_details$isoform_id_position
      if (!is.null(isoform_id_position)){
        out = list(tags$div(title=paste0("Display only those genes with more than one isoform in the selection, useful for finding those in which splicing patterns have been changed"),
                            checkboxInput("altsplice_checkbox", label = "Multiple isoforms only", value = FALSE)
        ))
      }else{
        out = NULL
      }
      
      return(out)
    })
    #Given: vector of comparisons possible in data.frame table from modifyData (user_input$allComparisonNames), 
    #  vector of current comparison (user_input$select_comp)
    #Returns: UI element allowing the selection of all comparisons apart from current one, with preserved previous selection
    output$comparison_found_in <- renderUI({
      all = user_input$allComparisonNames
      current = user_input$select_comp
      chosen = isolate(user_input$select2InputFoundIn)
      not = match(current,all)
      all_but_one = all[-not]
      
      out = list(selectInput("select2InputFoundIn", label = "In other comparisons", 
                             choices = all_but_one, 
                             selected = chosen,multiple=TRUE))
      return(out)
    })
    #Given: vector of comparisons possible in data.frame table from modifyData (user_input$allComparisonNames), 
    #  vector of current comparison (user_input$select_comp),
    #  vector of the comparisons selected in select2InputFoundIn
    #Returns: list of all possible comparisons that are not the current one or the one(s) selected in select2InputFoundIn (choices),
    #  vector of the comparisons selected in select2InputNotFoundIn
    comparison_not_found_inChoices = reactive({
      all = user_input$allComparisonNames
      current = user_input$select_comp
      not = match(current,all)
      all_but_one = all[-not]
      
      in_use = (user_input$select2InputFoundIn)
      if (length(in_use) > 0){
        not = match(in_use,all_but_one)
        all_but_one = all_but_one[-not]
      }
      chosen = isolate(user_input$select2InputNotFoundIn)
      
      return(list(choices = all_but_one, chosen = chosen))
    })
    #Given: list of all possible comparisons that are not the current one or the one(s) selected in select2InputFoundIn (choices),
    #  vector of the comparisons selected in select2InputNotFoundIn
    #Returns: UI element allowing the selection of all comparisons apart from current one or ones selected in select2InputFoundIn, with preserved previous selection
    #Requires reactive functions: comparison_not_found_inChoices()
    output$comparison_not_found_in <- renderUI({
      data = comparison_not_found_inChoices()
      choices = data$choices
      chosen = data$chosen
      
      out = list(selectInput("select2InputNotFoundIn", label = "Not in other comparisons", 
                             choices = choices, 
                             selected = chosen,multiple=TRUE))
      return(out)
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
    #Given: vector of species name from collectUser()
    #Returns: UI element of file select for filtering genes by
    #Requires reactive functions: collectUser()
    #Requires functions: speciesNames()
    output$filter_file_selection <- renderUI({
      names = geneListNames()
      if (is.null(names)){
        names = "Create a gene list to filter by"
      }else{
        names = c("None",names)
      }
      
      out = list(tags$div(title=paste0("Select a gene list to filter by; gene lists are created from the Gene lists menu."),
                          selectInput("select_gene_list", label = ("Genes in this list"), 
                             choices = names, 
                             selected = 1)))
      return(out)
    })  
  
  output$annotation_selection <- renderUI({
    out = list(uiOutput("built_in"))
    out[[2]] = uiOutput("annotation_file_selection")
    out[[3]] = tags$div(title="Temporarily hide all annotation columns",
                        checkboxInput("clearannot_file_checkbox", label = "Hide annotation(s) from file(s)", value = FALSE)
    )
    out[[4]] = tags$div(title="Each button removes the columns added by the relevant file",
                        uiOutput("annotation_buttons")
    )
    out[[5]] = br()
    
    return(out)
  })
    #Which built-in annotations  
    #Given: vector of species name (user_input$species)
    #Returns: NULL if species not found in species_data or 'Other', 
    #  else UI element of one checkbox for each
    #  built-in annotation specified for species and a return at the end
    #Requires reactive variables: user_input$species
    #Requires functions: speciesNames()
    output$built_in <- renderUI({
      species = user_input$species
      species_data = speciesNames(species)
      if (is.null(species_data)){
        out = NULL
      }else{
        names = species_data$built_in_list
        labels = species_data$built_in_list_text
        
        if (is.null(names)){
          out = NULL
        }else{
          out = list(checkboxInput(names[1], label = labels[1], value = FALSE))
          if (length(names) > 1){
            for (i in 2:length(names)){
              out[[i]] = checkboxInput(names[i], label = labels[i], value = FALSE)
            }
          }
          out[[(length(names) + 1)]] = br()
        }
      }
      
      #Handle connection_error
      if (is.null(data_cache$connectionError)){
        
      }else if (data_cache$connectionError){
        out = NULL
      }
      
      return(out)
    })
    #Given: vector of species name from user_input$species
    #Returns: UI element of file select for annotating genes with
    #Requires reactive variables: user_input$species
    #Requires functions: speciesNames()
    output$annotation_file_selection <- renderUI({
      species = user_input$species
      
      species_data = speciesNames(species)
      if (is.null(species_data)){
        gene.lookup.long = "Gene ID"
      }else{
        gene.lookup.long = species_data$gene.lookup.long
      }
      
      tags$div(title=paste0("A file with annotations you wish to add to the displayed genes; file must have a minimum of two columns, the first with ",gene.lookup.long,"s and subsequent ones with annotation to add, the first row must contain names for the columns; multiple files may be added one after the other"),
               fileInput("annot_file", label = "Annotation file")
      )
    })
    #Generate removal buttons for annotations and render as UI  
    #Given: global_annotation$annotation_sets_names
    #Returns: a UI element (list) with a button corresponding to each element in the list 
    #  global_annotation$annotation_sets_names (can be an empty list)
    #Returns reactive variables: annotation_removal_buttons$buttons
    #Requires reactive variables: global_annotation$annotation_sets_names
    annotation_removal_buttons <- reactiveValues()
    output$annotation_buttons <- renderUI({
      annotation_sets_names = (global_annotation$annotation_sets_names)
      if ((is.null(annotation_sets_names))|(length(annotation_sets_names) < 1)){ #No file selected / no annotations left in set
        temp_list = list(buttons = list(NULL))
        annotation_removal_buttons$buttons = temp_list$buttons
      }else if (length(annotation_sets_names) == 1){ #One annotation
        temp_list = list(buttons = list(actionButton(inputId = paste0("button",1),label = paste("Remove",annotation_sets_names[[1]]))))
        annotation_removal_buttons$buttons = temp_list$buttons
      }else{ #More than one annotation
        temp_list = list(buttons = list(actionButton(inputId = paste0("button",1),label = paste("Remove",annotation_sets_names[[1]]))))
        for (i in 2:length(annotation_sets_names)){
          temp_list$buttons[[i]] = actionButton(inputId = paste0("button",i),label = paste("Remove",annotation_sets_names[[i]]))
        }
        annotation_removal_buttons$buttons = temp_list$buttons
      }
      return(annotation_removal_buttons$buttons)
    })


  #UI checkboxes to select columns in table  
  #Given: integer of number of times the show_column_selection button has been pressed from 
  #       showColumnChecks(), data on species for names of columns, and data on lookup for number 
  #        of columns
  #Returns: NULL if showColumnChecks() = 0 or an exact multiple of 2, else returns UI element
  #         of checkboxes for showing/hiding each column
  #Requires reactive functions: showColumnChecks() 
  #Requires reactive variables: user_input$species, data_cache$lookup
  #Requires functions: speciesNames()
  output$column_selection <- renderUI({
    if (is.null(rawData())){ #Prior to loading a dataset
      out = list(br())
      out[[2]] = fluidRow(column(12, uiOutput("need_file_text")))
      return(out)
    }else{
      #n_clicks = showColumnChecks() #records number of times button clicked
      n_clicks = 1
      species = user_input$species
      species_data = speciesNames(species)
      gene.lookup = species_data$gene.lookup
      other.lookup.long = species_data$other.lookup.long
      
      if (is.null(gene.lookup)){
        gene.lookup = "Gene ID"
      }
      
      lookup = data_cache$lookup
      
      format = user_input$format
      format_details = formatLookup(format)
      long_names = format_details$long_colnames
      n_columns = format_details$n_columns
      gene_id_position = format_details$gene_id_position
      extra_position = format_details$extra_position
      extra_colname = format_details$extra_colname
      extra_long_colname = format_details$extra_long_colname
      gene_symbol_position = format_details$gene_symbol_position
      correct_colnames = format_details$correct_colnames
      
      #Check to see if had to add gene symbol column
      if (is.na(gene_symbol_position)){
        gene_symbol_name = "Gene Symbol"
        
        n_columns = n_columns + 1
        long_names = c(long_names,gene_symbol_name)
      }
      
      # #Check to see if had to add binary significance column
      if (is.na(match("significant",correct_colnames))){
        sig_name = "DE significant"
        
        n_columns = n_columns + 1
        long_names = c(long_names,sig_name)
      }
      
      #Generate for n_columns, with an extra column inserted at extra_position if lookup!=NULL
      temp_list = list(checks = list(checkboxInput("col1_checkbox", label = long_names[1], value = user_input$col1_checkbox)))
      
      if (is.null(lookup)){
        past_ep = 0
        for (i in 2:n_columns){
          if (i == extra_position){
            j = i + 1
            temp_list$checks[[i]] = checkboxInput(paste0("col",j,"_checkbox"), label = long_names[i], value = eval(parse(text = paste0("user_input$col",j,"_checkbox"))))
            past_ep = 1
          }else if (past_ep > 0){
            j = i + 1
            temp_list$checks[[i]] = checkboxInput(paste0("col",j,"_checkbox"), label = long_names[i], value = eval(parse(text = paste0("user_input$col",j,"_checkbox"))))
          }else{  
            temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = long_names[i], value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
          }
        }
      }else{
        n_columns = n_columns + 1
        past_ep = 0
        for (i in 2:n_columns){
          if (i == extra_position){
            temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = extra_long_colname, value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
            past_ep = 1
          }else if (past_ep > 0){
            j = i - 1
            temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = long_names[j], value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
          }else{  
            temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = long_names[i], value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
          }
        }
      }
      
      #Set position-specific species-specific special names
      if (is.null(lookup)){
        i = gene_id_position
        temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = gene.lookup, value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
      }else{
        if (lookup){
          i = gene_id_position
          temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = gene.lookup, value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
          j = extra_position
          temp_list$checks[[j]] = checkboxInput(paste0("col",j,"_checkbox"), label = other.lookup.long, value = eval(parse(text = paste0("user_input$col",j,"_checkbox"))))
        }else{
          i = gene_id_position
          temp_list$checks[[i]] = checkboxInput(paste0("col",i,"_checkbox"), label = other.lookup.long, value = eval(parse(text = paste0("user_input$col",i,"_checkbox"))))
          j = extra_position
          temp_list$checks[[j]] = checkboxInput(paste0("col",j,"_checkbox"), label = gene.lookup, value = eval(parse(text = paste0("user_input$col",j,"_checkbox"))))
        }
      }
      
      return(temp_list$checks)
    }
  })
    #UI show column selection  
    #Returns: integer of number of times the show_column_selection button has been pressed
    showColumnChecks <- eventReactive(input$show_column_selection, {
      return(input$show_column_selection)
    })
    observe({ #Updates label depending on whether checkboxes shown or not
      req(input$show_column_selection)
      n_clicks = showColumnChecks() #records number of times button clicked
      if (n_clicks > 0){
        if (n_clicks %% 2 == 0){
          updateActionButton(session, "show_column_selection",
                             label = "Show select columns")
        }else{
          updateActionButton(session, "show_column_selection",
                             label = "Hide select columns")        
        }
      }
    })  
  
  output$download_menu <- renderUI({
    out = list(HTML(paste0("<h4>Tables</h4>")))
    if (is.null(filteredData())){
      #return(NULL)
    }else{
      out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the filtered version of the data in tsv format",
                                                         downloadButton("downloadData", "Download table")))))
      out2[[2]] = br()
      out = c(out,out2)
    }
    if(user_input$select_gene_list2 == "None"){ #if no genelist or nonselected
      #return(NULL)
    }else{
      out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the table from the Gene lists genes tab in tsv format",
                                                         downloadButton("downloadfilterData", HTML("Download Gene lists <br/> genes table"))))))
      out2[[2]] = br()
      out = c(out,out2)
    }
    if (!is.null(genes_for_heatmap())){
      out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the table of expression of genes used to draw the heatmap in tsv format",
                                                         downloadButton("downloadHeatmapdata", "Download Heatmap table")))))
      out2[[2]] = br()
      out = c(out,out2)
    }
    if(!is.null(input$show_gProfile_enrichment)){
      n_clicks = input$show_gProfile_enrichment[1]
      if (n_clicks > 0){
        out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the gProfile enrichment data in tsv format",
                                                           downloadButton("downloadgProfileenrichmentdata", HTML("Download gProfile <br/> enrichment data"))))))
        out2[[2]] = br()
        out = c(out,out2)
      }
    }
    if(!is.null(input$show_protein_domain_enrichment)){
      n_clicks = input$show_protein_domain_enrichment[1]
      if (n_clicks > 0){
        out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Protein domain enrichment data in tsv format",
                                                           downloadButton("downloadProteindomainenrichmentdata", HTML("Download Protein domain <br/> enrichment data"))))))
        out2[[2]] = br()
        out = c(out,out2)
      }
    }
    if (!is.null(input$select_orth_org)){
      if (!is.null(make_ortholog_table())){
        out2 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the ortholog data in tsv format",
                                                           downloadButton("downloadorthologdata", "Download ortholog data")))))
        #out2[[2]] = br()
        out = c(out,out2)
      }
    }
    ##CODED BUT NOT CURRENTLY USED WIDGETS
    # out[[6]] = br()
    # out[[7]] = (fluidRow(column(3,offset = 1, tags$div(title="Click here to download the GO enrichment data in tsv format",
    #                                                   downloadButton("downloadGOenrichmentdata", "Download GO enrichment data")
    # ))))
    # out[[8]] = br()
    # out[[9]] = (fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Publication enrichment data in tsv format",
    #                                                    downloadButton("downloadPublicationenrichmentdata", "Download Publication enrichment data")
    # ))))
    
    out2 = list(HTML(paste0("<h4>Figures</h4>")))
    out = c(out,out2)
    
    
    if (!is.null(user_input$pergene_text)){
      if ((user_input$pergene_text) != ""){
        out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Download Per gene tab plot in svg format",
                                                           downloadButton("downloadplot2", "Download Per gene plot")))))
        out3[[2]] = br()
        out = c(out,out3)
      }
    }
    if (!is.null(filteredData())){
      out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Volcano plot in svg format",
                                                         downloadButton("downloadplot", "Download Volcano plot")))))
      out3[[2]] = br()
      out = c(out,out3)
    }
    if (!is.null(filteredData())){
      out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Scatter plot in svg format",
                                                         downloadButton("downloadplot5", "Download Scatter plot")))))
      out3[[2]] = br()
      out = c(out,out3)
    }
    if(!is.null(input$show_H_Clust)){
      n_clicks = input$show_H_Clust[1]
      if (n_clicks > 0){
          out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Hierarchical clustering in svg format",
                                                             downloadButton("downloadplot3", HTML("Download Hierarchical <br/> clustering"))))))
          out3[[2]] = br()
          out = c(out,out3)
      }
    }
    if(!is.null(input$show_filter_heatmap)){
      n_clicks = input$show_filter_heatmap[1]
      if (n_clicks > 0){
        out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the Heatmap in svg format",
                                                           downloadButton("downloadplot4", "Download Heatmap")))))
        out3[[2]] = br()
        out = c(out,out3)
      }
    }
    if(!is.null(input$show_gProfile_enrichment)){
      n_clicks = input$show_gProfile_enrichment[1]
      if (n_clicks > 0){
        out3 = list(fluidRow(column(3,offset = 1, tags$div(title="Click here to download the gProfile enrichment data figure in svg format",
                                                           downloadButton("downloadplot6", HTML("Download gProfile <br/> enrichment figure"))))))
        #out3[[2]] = br()
        out = c(out,out3)
      }
    }
    
    #remove Figures header if no figures
    if(!exists("out3")){
      out[[length(out)]] = NULL
    }
    #Remove all text if no buttons
    if(length(out) < 2){
      return(NULL)
    }
    return(out)
  })
    ##Generate downloads
    #Tables
    #Given: data.object from filteredData()
    #Returns: csv of the data.object, with default filename determined from comparison in collectUser()
    #Returns reactive variables:
    #Requires reactive functions: collectUser(), filteredData()
    #Requires reactive variables:
    #Requires functions: 
    output$downloadData <- downloadHandler(
      filename = function() {
        ui = collectUser()
        if (is.null(ui$select_comp)){
          name = paste0("DEdata-", Sys.Date(), ".tsv")
        }else{
          name = paste0(ui$select_comp,"-", Sys.Date(), ".tsv")
        }
        return(name)
      },
      content = function(file) {
        if (is.null(filteredData())){
          return(NULL)
        }else{
          data.object = filteredData()
          linkless_data = data.object$filtered_data
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(linkless_data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: data.object from filterFileData()
    #Returns: csv of the data.object
    #Returns reactive variables:
    #Requires reactive functions: filterFileData()
    #Requires reactive variables: user_input$select_gene_list2
    #Requires functions: 
    output$downloadfilterData <- downloadHandler(
      filename = function() {
        name = paste0("Genelistdata-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if(user_input$select_gene_list2 == "None"){ #if no genelist or nonselected
          return(NULL)
        }else{
          data = filterFileData()
          data = data.frame(x=data$table)
          lookup = data_cache$lookup
          
          #set FPKM values to numbers
          if (is.null(lookup)){
            for (i in 3:dim(data)[2]){
              data[,i] = as.numeric(as.character(data[,i])) #levels(data[,i])[as.numeric(data[,i])]
            }
          }else{
            for (i in 4:dim(data)[2]){
              data[,i] = as.numeric(as.character(data[,i])) #levels(data[,i])[as.numeric(data[,i])]
            }
          }
          write.table(data, file,sep="\t",row.names=F,quote=F)
          #return(data)
        }
        
      },contentType = "text/csv"
    )
    
    #Given: matrix from genes_for_heatmap()
    #Returns: csv of the matrix
    #Requires reactive functions: genes_for_heatmap()
    output$downloadHeatmapdata <- downloadHandler(
      filename = function() {
        name = paste0("Heatmaptable-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(genes_for_heatmap())){
          return(NULL)
        }else{
          data = genes_for_heatmap()
          data = data$unmod
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: matrix from make_go_enrichment_table()
    #Returns: csv of the matrix
    #Requires reactive functions: make_go_enrichment_table()
    output$downloadGOenrichmentdata <- downloadHandler(
      filename = function() {
        name = paste0("GOenrichment-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(make_go_enrichment_table())){
          return(NULL)
        }else{
          data = make_go_enrichment_table()
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: matrix from make_publication_enrichment_table()
    #Returns: csv of the matrix
    #Requires reactive functions: make_publication_enrichment_table()
    output$downloadPublicationenrichmentdata <- downloadHandler(
      filename = function() {
        name = paste0("Publicationenrichment-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(make_publication_enrichment_table())){
          return(NULL)
        }else{
          data = make_publication_enrichment_table()
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: matrix from make_proteindomain_enrichment_table()
    #Returns: csv of the matrix
    #Requires reactive functions: make_proteindomain_enrichment_table()
    output$downloadProteindomainenrichmentdata <- downloadHandler(
      filename = function() {
        name = paste0("Proteindomainenrichment-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(make_proteindomain_enrichment_table())){
          return(NULL)
        }else{
          data = make_proteindomain_enrichment_table()
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: matrix from make_gprofile_enrichment_table()
    #Returns: csv of the matrix
    #Requires reactive functions: make_gprofile_enrichment_table()
    output$downloadgProfileenrichmentdata <- downloadHandler(
      filename = function() {
        name = paste0("gProfileenrichment-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(make_gprofile_enrichment_table())){
          return(NULL)
        }else{
          data = make_gprofile_enrichment_table()$table
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    #Given: matrix from make_ortholog_table()
    #Returns: csv of the matrix
    #Requires reactive functions: make_ortholog_table()
    output$downloadorthologdata <- downloadHandler(
      filename = function() {
        name = paste0("orthologs-", Sys.Date(), ".tsv")
        return(name)
      },
      content = function(file) {
        if (is.null(make_ortholog_table())){
          return(NULL)
        }else{
          data = make_ortholog_table()$table
          
          # #filter the data based on search boxes on Table tab
          # rows = table_rows() #user_input$table_rows_all
          # if (!is.null(rows)){
          #   linkless_data = linkless_data[rows,]
          # }
          
          write.table(data, file,sep="\t",row.names=F,quote=F)
        }
      },contentType = "text/csv"
    )
    
    ##Figures
    #Given: plot from plot2_plot()
    #Returns: svg
    #Requires reactive variables: user_input$pergene_text
    #Requires functions: plot2_plot()
    output$downloadplot2 <- downloadHandler(
      filename = function() {
        name = paste0("expression_",user_input$pergene_text,"-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(user_input$pergene_text)){
          return(NULL)
        }else if ((user_input$pergene_text) == ""){
          return(NULL)
        }else{
          data.object = modifyData()
          format = data.object$data_format
          format_details = formatLookup(format)
          isoform_id_position = format_details$isoform_id_position
          
          if(is.null(isoform_id_position)){
            svg(file,width=10)
            plot2_plot()
            dev.off()
          }else{
            ggsave(file, plot = plot2_plot(), device = "svg",width=10)
          }
        }
      }
    )
    
    output$downloadplot <- downloadHandler(
      filename = function() {
        name = paste0("Volcano_plot-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(filteredData())){
          return(NULL)
        }else{
          svg(file)
          plot_plot()
          dev.off()
        }
      }
    )
    
    output$downloadplot3 <- downloadHandler(
      filename = function() {
        name = paste0("Hierarchical_clusters-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(genes_for_h_clust())){
          return(NULL)
        }else{
          svg(file,height=7,width=15)
          plot3_plot()
          dev.off()
        }
      }
    )
    
    output$downloadplot4 <- downloadHandler(
      filename = function() {
        name = paste0("Heatmap-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(genes_for_heatmap())){
          return(NULL)
        }else{
          svg(file,height=14,width=10)
          plot4_plot()
          dev.off()
        }
      }
    )
    
    output$downloadplot5 <- downloadHandler(
      filename = function() {
        name = paste0("Scatter_plot-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(filteredData())){
          return(NULL)
        }else{
          svg(file)
          plot5_plot()
          dev.off()
        }
      }
    )
    
    output$downloadplot6 <- downloadHandler(
      filename = function() {
        name = paste0("gProfile_enrichment_data_figure-", Sys.Date(), ".svg")
        return(name)
      },
      content = function(file) {
        if (is.null(make_gprofile_enrichment_table())){
          return(NULL)
        }else{
          #ensg_gene_names = make_gprofile_enrichment_table()$genes_used
          update = plot6_plot()
          max_genes = data_cache$max_genes# length(unique(ensg_gene_names))
          vert_size = 0.15 * max_genes
          vert_size = ceiling(vert_size) + 8
          # if (max_genes > 50){
          #   svg(file,height=28,width=10)
          # }else{
          #   svg(file,height=14,width=10)
          # }
          svg(file,height=vert_size,width=10)
          draw(plot6_plot())
          dev.off()
        }
      }
    )
  
  output$bottom_left_text <- renderUI({
    data.object = modifyData()
    file_ref = data.object$file_ref
    file_ref = as.character(file_ref)
    if (is.null(data.object)){
      out = list(br())
    }else{
      if (file_ref == "OFFLINE MODE"){
        out = list(fluidRow(column(12,offset = 1, file_ref)))
        out[[2]] = br()
      }else{
        #convert file_ref to date
        pat <- ".+_([0-9]+-[0-9]+-[0-9]+).tsv"
        file_ref = sub(pat, "\\1", file_ref[grepl(pat, file_ref)])
        
        header = "Gene ID conversion table from:"
        out = list(fluidRow(column(12,offset = 1, header)))
        out[[2]] = (fluidRow(column(12,offset = 1, file_ref)))
        out[[3]] = br()
      }
    }
    
    return(out)
  })

  
  
  ##Tabs 
  output$tab_panel <- renderUI({
    out = list(tabBox(id = "tabset1", width = 11, selected = "Summary",#user_input$which_tab, #temp disabled due to large genelist problems
                      tabPanel("Summary", tableOutput("summary")), 
                      tabPanel("Table", div(dataTableOutput("table"),style = "font-size:90%")),
                      tabPanel("Per gene", uiOutput("per_gene")),
                      tabPanel("Figures", uiOutput("graphs")),
                      tabPanel("Enrichment analysis", uiOutput("analysis")),
                      tabPanel("Orthologs", uiOutput("orthologs")),
                      tabPanel("Gene lists", uiOutput("gene_list_table_ui")),
                      tabPanel("Gene lists genes", uiOutput("gene_list_genes_table"))
    ))
    return(out)
  })

  #Summary  
  #Given: data.object from filteredData()
  #Returns: Summary of data.object as defined in mySummaryTable()
  #Requires reactive functions: filteredData()
  #Requires reactive variables: user_input$select_comp
  #Requires functions:mySummaryTable()
  output$summary <- renderTable({
    if (is.null(filteredData())){
      return("Select valid file for analysis...")
    }else{
      data.object = filteredData()
      if (!is.null(data.object)){
        format = data.object$data_format
        Summary = data.object$filtered_data
        
        # #filter the data based on search boxes on Table tab
        #  rows = table_rows() #user_input$table_rows_all
        #  if (!is.null(rows)){
        #    Summary = Summary[rows,]
        #  }
        
        if (is.null(dim(Summary))){ #this isn't working to catch these cases - fix
          return("Invalid comparison")
        }else{
          comparison = user_input$select_comp
          s1_name = strsplit(as.character(comparison)," vs ")[[1]][1] #Sample 1 name
          s2_name = strsplit(as.character(comparison)," vs ")[[1]][2] #Sample 2 name
          sum = mySummaryTable(Summary,s1_name,s2_name,format = format) #summarise data
          return(sum)
        }
      }else{
        return("Invalid comparison")
      }
    }
  },colnames=FALSE,digits=2)
  
  #Table
  # Generate an datatable view of the data  
  #Given: data.object from filteredData()
  #Returns: datatable of data.object$filtered_data from filteredData(), containing only the columns from whichColumns()
  #Returns reactive variables: user_input$table_rows_all
  #Requires reactive functions: filteredData(), whichColumns(), collectUser()
  #Requires functions: formatLookup()
  output$table <- DT::renderDataTable({
    if (is.null(filteredData())){
      return(NULL)
    }else{
      user_input$table_rows_all = NULL
      
      data.object = filteredData()
      format = data.object$data_format
      lookup = data.object$lookup
      
      format_details = formatLookup(format)
      numeric_columns = format_details$numeric_columns
      extra_colname = format_details$extra_colname
      
      ui = collectUser()
      species = ui$species
      
      #add links
      #make version of modified data with links
      if (is.null(lookup)){
        table = data.object$filtered_data
      }else{
        new_column = data.object$filtered_data[extra_colname][[1]] #data_cache$new_column
        data = createTableLinks(species,data.object$filtered_data,data.object$no_link_rows,lookup,new_column,format)
        table = data$table
      }
      data = data.frame(x=table)
      
      #make correct columns numeric
      for (i in 1:length(numeric_columns)){
        data[,numeric_columns[i]] = as.numeric(as.character(data[,numeric_columns[i]]))
      }
      
      #choose wanted columns
      wanted_columns = whichColumns()
      data = data[,wanted_columns]
      
      #user_input$nrows4 = dim(data)[1]
      
      return(data)
    }
  }, escape = FALSE, extensions = 'FixedHeader', 
  options = list(pageLength = 25, fixedHeader = TRUE, orderClasses = TRUE,
                 searchHighlight = TRUE, searchDelay = 1000,scrollX = TRUE),
  #columnDefs = columnTitles()),
  rownames = FALSE,selection = 'single', filter = 'bottom')
    #Show per gene tab for selected gene when selected on main table
    observeEvent(input$table_rows_selected,{
      if (is.null(input$table_rows_selected)){
      }else{
        if (length(input$table_rows_selected) > 0){
          user_input$table_rows_selected = input$table_rows_selected
        }
        data.object = filteredData()
        format = data.object$data_format
        #format = user_input$format
        format_details = formatLookup(format)
        gene_id_position = format_details$gene_id_position
        correct_colnames = format_details$correct_colnames
        gene_id_name = correct_colnames[gene_id_position]
        
        
        table = data.object$filtered_data 
        row = user_input$table_rows_selected
        
        isoform_id_position = format_details$isoform_id_position
        if (is.null(isoform_id_position)){
          gene = table[row,gene_id_name]
          gene = removeLink(gene)
        }else{
          isoform_id_name = correct_colnames[isoform_id_position]
          gene = table[row,isoform_id_name]
        }
        user_input$pergene_text = gene
        #user_input$which_tab = "Per gene"
        updateTabsetPanel(session,"tabset1",selected="Per gene")
      }
    })
    #Which columns to output
    #Given: number of columns in data from modifyData()
    #  booleans of UI checkboxes on whether columns should be shown from collectUser();
    #  boolean of whether the "gene_id" column in original data was the lookup_identifier;
    #  number of columns in data from filteredData() - i.e. the number of annotation columns added
    #Returns: vector of column numbers of columns that should be shown
    #Requires reactive functions: modifyData(), collectUser(), filteredData()
    whichColumns <- reactive({
      data.object = modifyData()
      data = data.object$modified_data
      n_data_columns = dim(data)[2]
      ui = collectUser()
      lookup = data.object$lookup #data_cache$lookup
      
      format = data.object$data_format
      format_details = formatLookup(format)
      n_columns = format_details$n_columns
      extra_position = format_details$extra_position
      gene_symbol_position = format_details$gene_symbol_position
      correct_colnames = format_details$correct_colnames
      
      #Check to see if had to add gene symbol column
      if (is.na(gene_symbol_position)){
        n_columns = n_columns + 1
      }
      
      #Check to see if had to add binary significance column
      if (is.na(match("significant",correct_colnames))){
        n_columns = n_columns + 1
      }
      
      checkbox_columns = ui$col1_checkbox
      if (is.null(lookup)){
        past_ep = 0
        for (i in 2:n_columns){
          if (i == extra_position){
            past_ep = 1
            j = i + 1
            checkbox_columns = c(checkbox_columns,eval(parse(text = paste0("ui$col",j,"_checkbox"))))
          }else if (past_ep > 0){
            j = i + 1
            checkbox_columns = c(checkbox_columns,eval(parse(text = paste0("ui$col",j,"_checkbox"))))
          }else{
            checkbox_columns = c(checkbox_columns,eval(parse(text = paste0("ui$col",i,"_checkbox"))))
          }
        }
      }else{
        n_columns = n_columns + 1
        for (i in 2:n_columns){
          checkbox_columns = c(checkbox_columns,eval(parse(text = paste0("ui$col",i,"_checkbox"))))
        }
      }
      
      wanted_columns = (1:n_data_columns)[checkbox_columns]
      
      #add columns for annotations
      data.object = filteredData()
      data = data.object$filtered_data
      n_data_columns_f = dim(data)[2]
      
      extra_columns = setdiff((1:n_data_columns_f),(1:n_data_columns))
      wanted_columns = c(wanted_columns,extra_columns)
      
      return(wanted_columns)
    })
  
  
  #Per gene
  ##PER GENE TAB - UI + GENERATE & RENDER OUTPUT
  #Final tab output  
  #Given:
  #Returns: UI element for the 'Per gene' tab containing description of gene, plot of gene and possibly other info
  #Requires UI elements: output$per_gene_plot, output$per_gene_description
  #Requires reactive functions: modifyData(), geneNamesUpdate()
  #Requires functions: speciesNames(), formatLookup()
  output$per_gene <- renderUI({
    data.object = modifyData()
    species = data.object$species #species = user_input$species
    species_data = speciesNames(species)
    gene.lookup.long = species_data$gene.lookup.long
    
    lookup = data.object$lookup
    
    out = list(br())
    
    if (is.null(modifyData())){ #Prior to file loading
      out[[2]] = textOutput("need_file_text")
    }else{
      out[[2]] = tags$div(title=paste0("Enter ",gene.lookup.long," to display data about that gene"),
                          selectizeInput("pergene_text", label = "Choose gene",
                                         choices = NULL,
                                         selected = "",multiple=FALSE)
      )
      
      if (species == 'Other'){
        out[[3]] = uiOutput("per_gene_plot")
      }else if (is.null(lookup)){ #can't work out what are identifers so have done the equivalent of species = other
        out[[3]] = uiOutput("per_gene_plot")
      }else{
        out[[3]] = uiOutput("per_gene_description")
        out[[4]] = uiOutput("per_gene_plot")
      }
      
      out[[5]] = uiOutput("DE_versus")
      out[[6]] = HTML("<br>")
      
      #Include species-specific data
      if (species_data$species_specific){
        if (is.null(species_data$species_specific_per_gene_data)){}else{
          if(!is.null(data_cache$geneIDsUnknown)){
            if (!data_cache$geneIDsUnknown){ #have not fallen back to 'Other' without updating the species name
              out[[7]] = uiOutput(species_data$species_specific_per_gene_data)
            }
          }else{
            out[[7]] = uiOutput(species_data$species_specific_per_gene_data)
          }
        }
      }
      
    }
    update = geneNamesUpdate() #updates the list of gene names - much quicker to generate the selectizeInput and then update than to generate with choices in place
    return(out)
  })
    #Fill in input$pergene_text with possible gene names from data_cache$names_list
    geneNamesUpdate <- reactive({
      names_list = data_cache$names_list
      
      updateSelectizeInput(session, 'pergene_text', choices = names_list, server = TRUE, selected=input$pergene_text)
      
      return(TRUE)
    })
    
    #Separate plot part to allow for gene not in dataset error - constantly resets if input$pergene_text is used in same function that generates UI for it  
    #Given: 
    #Returns:
    #Returns reactive variables:
    #Requires UI elements: output$plot2
    #Requires reactive functions: modifyData()
    #Requires reactive variables: user_input$pergene_text
    #Requires functions: formatLookup(), findGene(), findTranscript()
    output$per_gene_plot <- renderUI({
      data.object = modifyData()
      raw_gene_names = data.object$modified_data
      species = data.object$species
      format = data.object$data_format
      lookup = data.object$lookup
      format_details = formatLookup(format)
      gene = user_input$pergene_text
      isoform_id_position = format_details$isoform_id_position
      
      #Check gene in table - allow different versions of name
      isoform_id_position = format_details$isoform_id_position
      if (is.null(isoform_id_position)){
        geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
      }else{
        geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
        transcript = gene
        gene = geneLookupID
        if(is.null(gene)){
          geneLookupID = NULL
        }else{
          geneLookupID = lookup
        }
      }
      
      if(is.null(geneLookupID)){ #not in table
        out = list(textOutput("gene_not_found"))
        out[[2]] = HTML("<br>")
      }else{
        if (is.null(plotOutput("plot2"))){ #this doesn't work, should be changed
          out = list(HTML("<br>"))
        }else{
          out = list(plotOutput("plot2"))
        }
        out[[2]] = textOutput("plot2_text")
        if (is.null(isoform_id_position)){
          out[[3]] = textOutput("plot2_text2")
        }else{
          out[[3]] = textOutput("plot2_text2_isoform")
        }
        out[[4]] = HTML("<br>")  
        out[[5]] = uiOutput("sampleorder")
        out[[6]] = uiOutput("logscale")
        out[[7]] = uiOutput("removebar")
        out[[8]] = uiOutput("rotatename")
        out[[9]] = tags$div(title="Set the sample name text size, default 1 for genes, 10 for transcripts",numericInput("pergene_cexx", label = "Sample text size", value = user_input$pergene_cexx))
        out[[10]] = HTML("<br>")  
      }
      
      return(out)
    })
      #Separate part to allow log scale on graph
      output$logscale <- renderUI({
        out = list(tags$div(title="Use log10 scale on y-axis",
                            checkboxInput("gene_exp_log_scale", label = "Use log scale", value = user_input$gene_exp_log_scale)))
        out[[2]] = (tags$div(title="Remove values less than zero, may be helpful when using log scale",
                             checkboxInput("gene_exp_log_scale_lt_zero", label = "Remove values < 0", value = user_input$gene_exp_log_scale_lt_zero)))
        
        return(out)
      })
      #Separate part to allow removal of largest bar  
      output$removebar <- renderUI({
        list(tags$div(title="Remove the largest bar on the graph, useful where one bar is much largest than the rest",
                      checkboxInput("remove_max", label = "Remove largest bar", value = user_input$remove_max)))
      })
      #Separate part to allow removal of largest bar  
      output$rotatename <- renderUI({
        list(tags$div(title="Rotate sample names.",
                      checkboxInput("rotate_sample", label = "Rotate sample names", value = user_input$rotate_sample)))
      })
      #Separate part to allow setting order of samples on graph
      output$sampleorder <- renderUI({
        all = sampleNames() #user_input$allSampleNames
        chosen = isolate(user_input$select2InputSampleOrder) #all #isolate(user_input$allSampleNames)
        
        out = list(tags$div(title="Select which samples to display and the order to display them in",
                            selectInput("select2InputSampleOrder", label = "Sample order", 
                                        choices = all, 
                                        selected = chosen,multiple=TRUE)
        ))
        
        return(out)
      })
    
    #Separate parts to allow lookup using pergene_text - this can't be done in same renderUI as the 
    # input that generates pergene_text
    #Given: 
    #Returns: UI element consisting of a description of the gene from the web and a table of gene name variants
    #Returns reactive variables:
    #Requires reactive functions: modifyData()
    #Requires reactive variables: user_input$pergene_text
    #  data_cache$connectionErrorSS etc., data_cache$connectionError etc.
    #Requires functions: formatLookup(), speciesNames(), findGene(), findTranscript(), getAnnot()
    output$per_gene_description <- renderUI({
      data.object = modifyData()
      species = data.object$species
      format = data.object$data_format
      lookup = data.object$lookup
      raw_gene_names = data.object$modified_data
      format_details = formatLookup(format)
      gene = user_input$pergene_text
      species_data = speciesNames(species)
      
      name = (species_data$built_in_list)[1]
      
      #Check gene in table
      isoform_id_position = format_details$isoform_id_position
      if (is.null(isoform_id_position)){
        geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
      }else{
        geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
        transcript = gene
        gene = geneLookupID
        if(is.null(gene)){
          geneLookupID = NULL
        }else{
          geneLookupID = lookup
        }
      }
      
      if(is.null(geneLookupID)){ #not in table
        out = list(HTML(""))
      }else{
        if (species_data$species_specific){ #has species specific method
          if (data_cache$connectionErrorSS){
            time_since_test = as.numeric(difftime(Sys.time(),data_cache$connectionErrorTimeSS,units="mins"))
            
            if (time_since_test >= data_cache$connectionRetryTime){
              attempt = getAnnot(name,species,gene,source="species_specific",save=FALSE) #annotateBuiltIn(gene,species,name,mine=FALSE,single=TRUE)
              connection_error = attempt$connection_error
              if (connection_error){ 
                if (is.null(attempt$table)){ #mine=FALSE allows fallback to IM
                  response = NULL
                  data_cache$connectionErrorTimeSS = attempt$time
                }else{
                  response = attempt$table
                  #data_cache$connectionErrorTimeSS = attempt$time #if we set this would prevent fallback
                }
                data_cache$connectionErrorSS = TRUE
              }else{
                response = attempt$table
              }
            }else{
              response = NULL
            }
          }else{
            response = (getAnnot(name,species,gene,source="species_specific",save=FALSE))$table #(annotateBuiltIn(gene,species,name,mine=FALSE,single=TRUE))$table
          }
        }else{
          if (data_cache$connectionError){
            time_since_test = as.numeric(difftime(Sys.time(),data_cache$connectionErrorTime,units="mins"))
            
            if (time_since_test >= data_cache$connectionRetryTime){
              attempt = getAnnot(name,species,gene,save=FALSE) #annotateBuiltIn(gene,species,name,single=TRUE)
              connection_error = attempt$connection_error
              if (connection_error){
                response = NULL
                data_cache$connectionErrorTime = attempt$time
              }else{
                response = attempt$table
              }
            }else{
              response = NULL
            }
          }else{
            response = (getAnnot(name,species,gene,save=FALSE))$table #(annotateBuiltIn(gene,species,name,single=TRUE))$table
          }
        }
        
        if (is.null(response)){
          #pre looking up a gene
          if (length(data_cache$empty_graph) < 1){
            out = list(textOutput("choose_gene"))
            out[[2]] = HTML("<br>")          
          }else{
            out = list(textOutput("intermine_not_found"))
            out[[2]] = HTML("<br>")          
          }
        }else if (is.null(lookup)){ #happens when switching species to Other with a description already loaded
          out = list(HTML("<br>"))                
        }else{
          if (is.na(response)[1]){
            response = "No description found"
          }
          out = list(tableOutput("gene_names_table"))
          out[[2]] = HTML(response)
          out[[3]] = HTML("<br>")
        }        
      }
      
      return(out)
    })
      #Gene name versions  
      #Given: gene name (user_input$pergene_text), data.object containing all data on genes (modifyData())
      #Returns: table of gene name variants
      #Requires reactive functions: modifyData()
      #Requires reactive variables: user_input$pergene_text
      #Requires functions: formatLookup(), speciesNames()
      output$gene_names_table <- renderTable({
        data.object = modifyData()
        raw_gene_names = data.object$modified_data
        gene = user_input$pergene_text
        species = data.object$species
        format = data.object$data_format
        
        format_details = formatLookup(format)
        gene_id_position = format_details$gene_id_position
        correct_colnames = format_details$correct_colnames
        gene_id_name = correct_colnames[gene_id_position]
        extra_colname = format_details$extra_colname
        gene_symbol_position = format_details$gene_symbol_position
        if (!is.na(gene_symbol_position)){
          gene_symbol_name = correct_colnames[gene_symbol_position]
        }else{ 
          gene_symbol_name = "gene_symbol"
        }
        lookup = data.object$lookup #data_cache$lookup
        
        #convert transcript to gene name
        isoform_id_position = format_details$isoform_id_position
        if (!is.null(isoform_id_position)){
          transcript = gene
          gene = findTranscript(transcript,raw_gene_names,format = format,lookup)
        }
        
        if (is.null(lookup)){
          if (!is.na(gene_symbol_position)){
            symbol = raw_gene_names[match(gene,raw_gene_names[gene_id_name][[1]]),gene_symbol_name] #needs [[1]] because gene_id isn't factorised
            symbol = as.character(symbol)
          }else{
            symbol = NA
          }
        }else{
          symbol = raw_gene_names[match(gene,raw_gene_names[gene_id_name][[1]]),gene_symbol_name] #needs [[1]] because gene_id isn't factorised
          symbol = as.character(symbol)
        }
        
        
        
        if(is.na(match(extra_colname,names(raw_gene_names)))){ #species = Other or can't do translation - shouldn't reach here
          seq_name = NA
        }else{
          seq_name = raw_gene_names[match(gene,raw_gene_names[gene_id_name][[1]]),extra_colname] #needs [[1]] because gene_id isn't factorised
          seq_name = as.character(seq_name)
        }
        
        species_data = speciesNames(species)
        gene.lookup = species_data$gene.lookup
        other.lookup.long = species_data$other.lookup.long
        
        
        
        if (is.na(seq_name)){
          data = c(gene.lookup,"Gene symbol","",gene,symbol,"")
          if (length(data) < 6){
            data = c(gene.lookup,"Gene symbol","",NA,NA,"")
          }
          
          out = matrix(data,3,2)
        }else{
          data = c(gene.lookup,"Gene symbol",other.lookup.long,"",gene,symbol,seq_name,"")
          if (length(data) < 8){
            data = c(gene.lookup,"Gene symbol",other.lookup.long,"",NA,NA,NA,"")
          }
          
          out = matrix(data,4,2)
        }
        
        return(out)
      },colnames=FALSE)
    
    #Separate part to allow for output$sampleorder calculation
    output$DE_versus <- renderUI({
      #all_conditions = sampleNames()
      conditions = input$select2InputSampleOrder
      
      out = list(tags$div(title="Choose which condition to compare others to, this bar will be coloured black; those significantly DE versus this conditions will be indicated on the graph with *, this will only be done for those comparisons included in the dataset - samples where comparisons were not made will be coloured grey",
                          selectInput("select_DE", label = ("Select DE versus"),
                                      choices = conditions,
                                      selected = 'WT')))
      return(out)
    })
    
    #SPECIES SPECIFIC PER GENE TAB
    #All wormbase bits  
    #Given: gene name and table of gene names
    #Returns: UI element of button to show/hide WormBase annotation if gene in table
    #Requires reactive functions: modifyData()
    #Requires reactive variables: user_input$pergene_text
    #Requires functions: formatLookup(), speciesNames(), findGene(), findTranscript()
    output$per_gene_wormbase <- renderUI({
      data.object = modifyData()
      raw_gene_names = data.object$modified_data
      species = data.object$species
      format = data.object$data_format
      lookup = data.object$lookup
      format_details = formatLookup(format)
      gene = user_input$pergene_text
      
      species_data = speciesNames(species)
      database = species_data$database
      
      #Check gene in table - allow different versions of name
      isoform_id_position = format_details$isoform_id_position
      if (is.null(isoform_id_position)){
        geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
      }else{
        geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
        transcript = gene
        gene = geneLookupID
        if(is.null(gene)){
          geneLookupID = NULL
        }else{
          geneLookupID = lookup
        }
      }
      
      if(is.null(geneLookupID)){
        out = list(HTML(""))
      }else{
        out = list(HTML(paste0("<h4>Data from ",database,"</h4>")))
        out[[2]] = tags$div(title=paste0("Click to show/hide the ",database," tables on homology"),
                            actionButton("show_wormbase_homology","Show Homology"))
        out[[3]] = uiOutput("wormbase_homology")       
      }
      
      return(out)
    })
      #Homology
      #Create wormbase homology table - these include links that need to change if the API changes  
      #Display or not  
      #Given: Number of time button clicked and info about whether online
      #Returns: Tables with info about gene homology
      #Requires reactive functions: showWormbaseHomology()
      #Requires reactive variables: user_input$species
      #Requires functions: speciesNames(), testConnection()
      output$wormbase_homology <- renderUI({
        n_clicks = showWormbaseHomology() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            #connection test
            species_data = speciesNames(user_input$species)
            species_link_test = species_data$species_link_test
            
            test_data = testConnection(species_link_test)
            test = test_data$connection
            time = test_data$time
            
            if (test){
              out = list(tableOutput("wormbase_blastp_table"))
              out[[2]] = tableOutput("wormbase_ortholog_table")
            }else{
              out = list(textOutput("species_specific_not_found"))
            }
            
            if (is.null(out)){
              out = list(textOutput("species_specific_not_found"))
            }
            return(out)
          }
        }else{
          return(NULL)
        }
      })
        #Given: Gene name
        #Returns: Table of blastp hits for the gene
        #Requires reactive functions: modifyData()
        #Requires reactive variables: user_input$pergene_text
        #Requires functions: formatLookup(), worm_lookup(), findGene(), findTranscript()
        output$wormbase_blastp_table <- renderTable({
          if (is.null(modifyData())){
            return(NULL)
          }else{
            data.object = modifyData()
            format = data.object$data_format
            lookup = data.object$lookup
            format_details = formatLookup(format)
            raw_gene_names = data.object$modified_data
            gene = user_input$pergene_text
            
            #Check gene in table - allow different versions of name
            isoform_id_position = format_details$isoform_id_position
            if (is.null(isoform_id_position)){
              geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
            }else{
              geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
              transcript = gene
              gene = geneLookupID
              if(is.null(gene)){
                geneLookupID = NULL
              }else{
                geneLookupID = lookup
              }
            }
            
            if(is.null(geneLookupID)){ #not in table - gets changed to lookup identifier elsewhere
              return(NULL)
            }else{
              response = worm_lookup(gene,"wormbase_blastp_table")
              if(is.null(response)){
                return(NULL)
              }else{
                data = response$table
                connection_error = response$connection_error
                time = response$time
                
                return(data)
              }
            }
          }
        })  
        #Given: Gene name
        #Returns: Table of orthologs for the gene
        #Requires reactive functions: modifyData()
        #Requires reactive variables: user_input$pergene_text
        #Requires functions: formatLookup(), worm_lookup(), findGene(), findTranscript()
        output$wormbase_ortholog_table <- renderTable({
          if (is.null(modifyData())){
            return(NULL)
          }else{
            data.object = modifyData()
            format = data.object$data_format
            lookup = data.object$lookup
            format_details = formatLookup(format)
            raw_gene_names = data.object$modified_data
            gene = user_input$pergene_text
            
            #Check gene in table - allow different versions of name
            isoform_id_position = format_details$isoform_id_position
            if (is.null(isoform_id_position)){
              geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
            }else{
              geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
              transcript = gene
              gene = geneLookupID
              if(is.null(gene)){
                geneLookupID = NULL
              }else{
                geneLookupID = lookup
              }
            }
            
            if(is.null(geneLookupID)){ #not in table - gets changes to lookup identifier elsewhere
              return(NULL)
            }else{
              response = worm_lookup(gene,"wormbase_ortholog_table")
              if(is.null(response)){
                return(NULL)
              }else{
                data = response$table
                connection_error = response$connection_error
                time = response$time
                
                return(data)
              }
            }
          }
        })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_wormbase_homology
      #Returns: The number of clicks on input$show_wormbase_homology
      #Requires reactive variables: input$show_wormbase_homology
      showWormbaseHomology <- eventReactive(input$show_wormbase_homology, {
        return(input$show_wormbase_homology)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showWormbaseHomology()
      #Requires reactive variables: input$show_wormbase_homology
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_wormbase_homology)
        n_clicks = showWormbaseHomology() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_wormbase_homology",
                               label = "Show Homology")
          }else{
            updateActionButton(session, "show_wormbase_homology",
                               label = "Hide Homology")        
          }
        }
      })    
    
  
  #Figures
  ##FIGURES TAB - UI + GENERATE & RENDER OUTPUT
  #Final tab output  
  #Given:
  #Returns: UI element for the 'Figure' tab
  #Requires UI elements: output$need_file_text, output$volcano_plot
  #Requires reactive functions: modifyData()
  #Requires functions: 
  output$graphs <- renderUI({
    out = list(br())
    
    if (is.null(modifyData())){ #Prior to file loading
      out[[2]] = textOutput("need_file_text")
    }else{
      out[[2]] = (HTML(paste0("<h4>Volcano plot</h4>")))
      out[[3]] = tags$div(title=paste0("Click to show/hide the volcano plot of currently chosen data"),
                          actionButton("show_volcano_plot","Show volcano plot"))
      out[[4]] = uiOutput("volcano_plot")     
      
      out[[5]] = HTML("<br>")
      
      out[[6]] = (HTML(paste0("<h4>Scatter plot</h4>")))
      out[[7]] = tags$div(title=paste0("Click to show/hide the scatter plot of the expression of the currently chosen data. Currently selected data is shown in red. All data is shown in grey."),
                          actionButton("show_scatter_plot","Show scatter plot"))
      out[[8]] = uiOutput("scatter_plot")     
      
      out[[9]] = HTML("<br>")
      
      out[[10]] = (HTML(paste0("<h4>Pairwise comparison</h4>")))
      out[[11]] = tags$div(title=paste0("Click to show/hide the scatter plot of the log expression in comparison 1 (chosen in Filter Input menu) and comparison 2 (chosen below)."),
                          actionButton("show_pairwise_plot","Show pairwise comparison plot"))
      out[[12]] = uiOutput("pairwise_plot")     
      
      out[[13]] = HTML("<br>")
      
      out[[14]] = (HTML(paste0("<h4>Clustering</h4>")))
      out[[15]] = tags$div(title=paste0("Click to show/hide the hierarchical clustering by log2foldchange of significantly DE gene in currently chosen data"),
                           actionButton("show_H_Clust","Show hierarchical clustering"))
      out[[16]] = uiOutput("h_clust")     
      
      out[[17]] = HTML("<br>")
      
      out[[18]] = (HTML(paste0("<h4>Heatmap</h4>")))
      out[[19]] = tags$div(title=paste0("Click to show/hide the heatmap of gene expression in currently chosen data"),
                           actionButton("show_filter_heatmap","Show heatmap"))
      out[[20]] = HTML("<br>")
      out[[21]] = uiOutput("filter_heatmap")     
      
      out[[22]] = HTML("<br>")
    }
    return(out)
  })
    #Volcano Plot
    #Given: Number of time button clicked
    #Returns: Volcano Plot
    #Requires UI elements: output$plot, output$plot_text
    #Requires reactive functions: showVolcanoPlot()
    output$volcano_plot <- renderUI({
      n_clicks = showVolcanoPlot() #records number of times button clicked
      if (n_clicks > 0){
        if (n_clicks %% 2 == 0){
          return(NULL)
        }else{
          out = list(plotOutput("plot"))
          out[[2]] = textOutput("plot_text")
          out[[3]] = HTML("<br>")
          out[[4]] = tags$div(title=paste0("Click to show/hide the graphical options"),
                              actionButton("show_volcano_plot_options","Show options"))
          out[[5]] = uiOutput("volcano_plot_options")
          
          return(out)
        }
      }else{
        return(NULL)
      }
    })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_volcano_plot
      #Returns: The number of clicks on input$show_volcano_plot
      #Requires reactive variables: input$show_volcano_plot
      showVolcanoPlot <- eventReactive(input$show_volcano_plot, {
        return(input$show_volcano_plot)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showVolcanoPlot()
      #Requires reactive variables: input$show_volcano_plot
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_volcano_plot)
        n_clicks = showVolcanoPlot() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_volcano_plot",
                               label = "Show volcano plot")
          }else{
            updateActionButton(session, "show_volcano_plot",
                               label = "Hide volcano plot")        
          }
        }
      })    
      #Given: Number of time button clicked
      #Returns: Volcano plot options
      #Requires UI elements: 
      #Requires reactive functions: showVolcanoPlotOptions()
      #Requires reactive variables: user_input$volcano_plot_grey_data, user_input$volcano_plot_no_rearrange
      output$volcano_plot_options <- renderUI({
        n_clicks = showVolcanoPlotOptions() #records number of times button clicked
        #n_clicks = 1
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            out = list(tags$div(title="Show all data in grey",
                                checkboxInput("volcano_plot_grey_data", label = "Show all data in grey", value = user_input$volcano_plot_grey_data)))
            out[[2]] = (tags$div(title="Keep the scale used for all data. If 'Show all data in grey' checked does nothing.",
                                 checkboxInput("volcano_plot_no_rearrange", label = "Keep original scale", value = user_input$volcano_plot_no_rearrange)))
            
            return(out)
          }
        }else{
          return(NULL)
        }
      })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_volcano_plot_options
      #Returns: The number of clicks on input$show_volcano_plot_options
      #Requires reactive variables: input$show_volcano_plot_options
      showVolcanoPlotOptions <- eventReactive(input$show_volcano_plot_options, {
        return(input$show_volcano_plot_options)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showVolcanoPlotOptions()
      #Requires reactive variables: input$show_volcano_plot_options
      observe({ #Updates label depending on whether options shown or not
        req(input$show_volcano_plot_options)
        n_clicks = showVolcanoPlotOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_volcano_plot_options",
                               label = "Show options")
          }else{
            updateActionButton(session, "show_volcano_plot_options",
                               label = "Hide options")        
          }
        }
      })    
      
    #Hierarchical clustering
    #Given: Number of time button clicked
    #Returns: Hierarchical clustering
    #Requires UI elements: output$plot3
    #Requires reactive functions: showHClust(), genes_for_h_clust()
    output$h_clust <- renderUI({
      n_clicks = showHClust() #records number of times button clicked
      if (n_clicks > 0){
        if (n_clicks %% 2 == 0){
          return(NULL)
        }else{
          data = genes_for_h_clust()
          if (is.null(data)){
            out = list(HTML("<br>"))
            out[[2]] =(textOutput("noplot3_text"))
          }else{
            out = list(plotOutput("plot3"))
            out[[2]] = HTML("<br>")
            out[[3]] = textOutput("plot3_text")
          }
          
          return(out)
        }
      }else{
        return(NULL)
      }
    })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_H_Clust
      #Returns: The number of clicks on input$show_H_Clust
      #Requires reactive variables: input$show_H_Clust
      showHClust <- eventReactive(input$show_H_Clust, {
        return(input$show_H_Clust)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showHClust()
      #Requires reactive variables: input$show_H_Clust
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_H_Clust)
        n_clicks = showHClust() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_H_Clust",
                               label = "Show hierarchical clustering")
          }else{
            updateActionButton(session, "show_H_Clust",
                               label = "Hide hierarchical clustering")        
          }
        }
      })    
      #Given: data.object from modifyData(), format and ui
      #Returns: Hierarchical clustering of the significantly DE genes
      #Requires reactive functions: modifyData(), collectUser(), geneListNames()
      #Requires reactive variables: user_input$select_gene_list, data_cache$gene_lists
      #Requires functions: formatLookup(), combineComparisons(), applyFilters()
      genes_for_h_clust <- reactive({
        data.object = modifyData()
        ui = collectUser()
        format = data.object$data_format
        
        format_details = formatLookup(format)
        
        gene_symbol_position = format_details$gene_symbol_position
        FC_position = format_details$FC_position
        value_1_position = format_details$value_1_position
        value_2_position = format_details$value_2_position
        adjusted_p_position = format_details$adjusted_p_position
        gene_id_position = format_details$gene_id_position
        
        correct_colnames = format_details$correct_colnames
        value_1_name = correct_colnames[value_1_position]
        value_2_name = correct_colnames[value_2_position]
        gene_symbol_name = correct_colnames[gene_symbol_position]
        FC_name = correct_colnames[FC_position]
        gene_id_name = correct_colnames[gene_id_position]
        
        if (is.na(adjusted_p_position)){
          sig_name = "significant"
        }else{
          sig_position = format_details$sig_position
          sig_name = correct_colnames[sig_position]
        }
        
        isoform_id_position = format_details$isoform_id_position
        
        #Use transcript names if that type of data
        if (is.null(isoform_id_position)){
        }else{
          isoform_id_name = correct_colnames[isoform_id_position]
          gene_symbol_name = isoform_id_name
        }
        
        #Filter based on UI
        #Which comparison and multiple comparisons filter
        if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
          return(NULL)
        }else{
          data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
        }
        data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)

        if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
        }else{
          which = match(user_input$select_gene_list,geneListNames())
          gene_list = data_cache$gene_lists[[which]]
          new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
          if (ui$invertfilter_file_checkbox){
            data.object$filtered_data = new_data$nonmatches_have
            data.object$unused_genes = gene_list
          }else{
            data.object$filtered_data = new_data$matches
            data.object$unused_genes = new_data$nonmatches_want
          }
        }
        
        data = data.object$filtered_data
        
        if (dim(data)[1] > 0){
          
          data2 = as.matrix(data[c(FC_name,value_1_name,value_2_name)])
          rownames(data2) = data[gene_symbol_name][[1]]
          #sig only
          data2 = data2[as.logical(data[sig_name] == "yes"),]
          #remove Inf
          data2 = data2[!is.infinite(data2[,1]),]
          #remove NA
          data2 = data2[!is.na(data2[,1]),]
          
          return(data2)
        }else{
          return(NULL)
        }
      })
      
    #Heatmap
    #Given: Number of time button clicked
    #Returns: Heatmap
    #Requires UI elements: output$plot4
    #Requires reactive functions: showFHMap(), genes_for_heatmap()
    output$filter_heatmap <- renderUI({
      n_clicks = showFHMap() #records number of times button clicked
      if (n_clicks > 0){
        if (n_clicks %% 2 == 0){
          return(NULL)
        }else{
          data = genes_for_heatmap()
          if (is.null(data)){
            out=list(textOutput("noplot4_text"))
          }else{
            out = list(textOutput("plot4_text"))
            out[[2]] = HTML("<br>")
            out[[3]] = (plotOutput("plot4", width = "100%", height = "800px"))
            out[[4]] = tags$div(title=paste0("Click to show/hide the graphical options"),
                                actionButton("show_heatmap_options","Show options"))
            out[[5]] = uiOutput("heatmap_options")
            out[[6]] = HTML("<br>")
            out[[7]] = dataTableOutput("plot4table")
          }
          return(out)
        }
      }else{
        return(NULL)
      }
    })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_filter_heatmap
      #Returns: The number of clicks on input$show_filter_heatmap
      #Requires reactive variables: input$show_filter_heatmap
      showFHMap <- eventReactive(input$show_filter_heatmap, {
        return(input$show_filter_heatmap)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showFHMap()
      #Requires reactive variables: input$show_filter_heatmap
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_filter_heatmap)
        n_clicks = showFHMap() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_filter_heatmap",
                               label = "Show heatmap")
          }else{
            updateActionButton(session, "show_filter_heatmap",
                               label = "Hide heatmap")        
          }
        }
      })   
      #Given: Number of time button clicked
      #Returns: Heatmap options
      #Requires UI elements: 
      #Requires reactive functions: showHeatmapOptions()
      #Requires reactive variables: user_input$heatmap_log_scale, user_input$heatmap_no_rearrange
      output$heatmap_options <- renderUI({
        n_clicks = showHeatmapOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            out = list(tags$div(title="Use log10 expression",
                                checkboxInput("heatmap_log_scale", label = "Use log expression", value = user_input$heatmap_log_scale)))
            out[[2]] = (tags$div(title="Cluster genes, unticked this will prevent the genes from being reordered",
                                 checkboxInput("heatmap_no_rearrange", label = "Cluster genes", value = user_input$heatmap_no_rearrange)))
            out[[3]] = (tags$div(title="Display heatmap for all filtered genes, not just 50 most DE. WARNING may be very slow to calculate.",
                                 checkboxInput("heatmap_all_genes", label = "All genes", value = user_input$heatmap_all_genes)))
            #out[[4]] = uiOutput("heatmap_sampleorder")
            out[[4]] = tags$div(title="Set the sample name text size, default 1",numericInput("heatmap_cexCol", label = "Sample text size", value = user_input$heatmap_cexCol))
            out[[5]] = tags$div(title="Set the sample name text size, default 0.8",numericInput("heatmap_cexRow", label = "Gene text size", value = user_input$heatmap_cexRow))
            return(out)
          }
        }else{
          return(NULL)
        }
      })
        #Separate part to allow setting order of samples on graph
        output$heatmap_sampleorder <- renderUI({
          all = sampleNames() #user_input$allSampleNames
          chosen = all #isolate(user_input$select2InputSampleOrder_heatmap) #all #isolate(user_input$allSampleNames)
          
          out = list(tags$div(title="Select which samples to display and the order to display them in",
                              selectInput("select2InputSampleOrder_heatmap", label = "Sample order", 
                                          choices = all, 
                                          selected = chosen,multiple=TRUE)
          ))
          
          return(out)
        })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_heatmap_options
      #Returns: The number of clicks on input$show_heatmap_options
      #Requires reactive variables: input$show_heatmap_options
      showHeatmapOptions <- eventReactive(input$show_heatmap_options, {
        return(input$show_heatmap_options)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showHeatmapOptions()
      #Requires reactive variables: input$show_heatmap_options
      observe({ #Updates label depending on whether options shown or not
        req(input$show_heatmap_options)
        n_clicks = showHeatmapOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_heatmap_options",
                               label = "Show options")
          }else{
            updateActionButton(session, "show_heatmap_options",
                               label = "Hide options")        
          }
        }
      })    
      #Given: data.object from modifyData(), format and ui
      #Returns: NULL if no data in input else list of matrix of per sample expression of genes to plot the heatmap with (table)
      #  matrix of per sample expression of genes to plot the heatmap with that has no NAs (unmod)
      #Requires reactive functions: modifyData(), collectUser(), geneListNames()
      #Requires reactive variables: user_input$heatmap_log_scale, user_input$select_gene_list, data_cache$gene_lists
      #Requires functions: formatLookup(), combineComparisons(), applyFilters()
      genes_for_heatmap <- reactive({
        data.object = modifyData()
        ui = collectUser()
        format = data.object$data_format
        
        format_details = formatLookup(format)
        
        gene_id_position = format_details$gene_id_position
        gene_symbol_position = format_details$gene_symbol_position
        condition_1_position = format_details$condition_1_position
        condition_2_position = format_details$condition_2_position
        value_1_position = format_details$value_1_position
        value_2_position = format_details$value_2_position
        FC_position = format_details$FC_position
        
        correct_colnames = format_details$correct_colnames
        gene_id_name = correct_colnames[gene_id_position]
        gene_symbol_name = correct_colnames[gene_symbol_position]
        condition_1_name = correct_colnames[condition_1_position]
        condition_2_name = correct_colnames[condition_2_position]
        value_1_name = correct_colnames[value_1_position]
        value_2_name = correct_colnames[value_2_position]
        FC_name = correct_colnames[FC_position]
        
        isoform_id_position = format_details$isoform_id_position
        
        #Use transcript names if that type of data
        if (is.null(isoform_id_position)){
        }else{
          isoform_id_name = correct_colnames[isoform_id_position]
          gene_symbol_name = isoform_id_name
        }
        
        #Filter based on UI
        #Which comparison and multiple comparisons filter
        if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
          return(NULL)
        }else{
          data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
        }
        data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)
        
        if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
        }else{
          which = match(user_input$select_gene_list,geneListNames())
          gene_list = data_cache$gene_lists[[which]]
          new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
          if (ui$invertfilter_file_checkbox){
            data.object$filtered_data = new_data$nonmatches_have
            data.object$unused_genes = gene_list
          }else{
            data.object$filtered_data = new_data$matches
            data.object$unused_genes = new_data$nonmatches_want
          }
        }
        
        data = data.object$filtered_data
        gene_names = (data[gene_symbol_name][[1]])
        
        if (dim(data)[1] > 50){ #limit on size due to computation time issues
          if(!user_input$heatmap_all_genes){
            gene_names = (data[order(abs(data[FC_name][[1]]),decreasing = TRUE),gene_symbol_name])[1:50]
          }
        }
        
        data = data.object$modified_data
        data = data[(data[gene_symbol_name][[1]] %in% gene_names), ]
        
        
        #Get expression values for all comparisons for genes
        allcomparisons = unique(c(levels(factor(data[condition_1_name][[1]])),levels(factor(data[condition_2_name][[1]]))))
        count_conditions = summary(factor(c(as.character(data[condition_1_name][[1]]),as.character(data[condition_2_name][[1]]))))
        if(length(count_conditions) < 1){
          return(NULL)
        }
        comp_to_use = match(max(count_conditions),count_conditions) #the one with most comparisons
        comp_to_use_name = sort(allcomparisons)[comp_to_use]
        allcomparisons = allcomparisons[-match(comp_to_use_name,allcomparisons)]
        
        first_col = data[((data[condition_1_name][[1]] == comp_to_use_name)&(data[condition_2_name][[1]] == allcomparisons[1])),value_1_name]
        names(first_col) = data[((data[condition_1_name][[1]] == comp_to_use_name)&(data[condition_2_name][[1]] == allcomparisons[1])),gene_symbol_name]
        if (length(first_col) < 1){
          first_col = data[((data[condition_2_name][[1]] == comp_to_use_name)&(data[condition_1_name][[1]] == allcomparisons[1])),value_2_name]
          names(first_col) = data[((data[condition_2_name][[1]] == comp_to_use_name)&(data[condition_1_name][[1]] == allcomparisons[1])),gene_symbol_name]
        }
        data3 = list(first_col)
        
        for (i in 1:(length(allcomparisons))){
          first_col = data[((data[condition_1_name][[1]] == comp_to_use_name)&(data[condition_2_name][[1]] == allcomparisons[i])),value_2_name]
          names(first_col) = data[((data[condition_1_name][[1]] == comp_to_use_name)&(data[condition_2_name][[1]] == allcomparisons[i])),gene_symbol_name]
          if (length(first_col) < 1){ #if comparison only exists other way round
            first_col = data[((data[condition_2_name][[1]] == comp_to_use_name)&(data[condition_1_name][[1]] == allcomparisons[i])),value_1_name]
            names(first_col) = data[((data[condition_2_name][[1]] == comp_to_use_name)&(data[condition_1_name][[1]] == allcomparisons[i])),gene_symbol_name]
          }
          data3[[(i+1)]] = (first_col)
        }
        
        
        ## Identify the unique "rownames" for all list items
        Rows <- unique(unlist(lapply(data3, names)))
        
        ## Create a matrix of NA values 
        ##   with appropriate dimensions and dimnames
        data4 <- matrix(NA, nrow = length(Rows), ncol = length(data3), 
                        dimnames = list(Rows, c(comp_to_use_name,allcomparisons)))
        
        ## Use your `for` loop to fill it in
        ##   with the appropriate values from your list
        for (i in seq_along(data3)) {
          data4[names(data3[[i]]), i] <- data3[[i]]
        }
        
        
        if (dim(data4)[1] < 1){
          return(NULL)
        }
        
        if(user_input$heatmap_log_scale){
          data4 = log10(data4)
        }
        
        
        #if (!user_input$heatmap_no_rearrange){
        #attempt to make greater distinction between zeros and non-zeros
        data5 = data4
        data4[(data4 == 0)] = NA
        pos_columns_zero = (1:(dim(data4)[2]))[(colSums(data4,na.rm=T) == 0)]
        def_columns = intersect((1:(dim(data4)[2]))[is.na(colSums(data4))],pos_columns_zero)
        if (length(def_columns) > 0){ #reset if any column becomes completely NA as then hclust fails
          data4 = data5
        }
        #}
        
        #handle NAs + infs
        data4[is.infinite(data4)&(0 > data4)] = 0
        #data4[is.infinite(data4)&(0 < data4)] = max(data4)
        #data4[is.na(data4)] = 0
        
        #remove columns with no comparisons
        counts = apply(apply(data5,2,is.na),2,sum)
        to_remove = counts == dim(data5)[1]
        data4 = data4[,!to_remove]
        data5 = data5[,!to_remove]
        
        return(list(table = data4, unmod = data5))
      })
      output$plot4table<- DT::renderDataTable({
        data = genes_for_heatmap()
        data = data$unmod
        if (is.null(data)){
          return(NULL)
        }else{
          return(data)
        }
      },options = list(scrollX = TRUE))
      
    #Scatter plot of expression
    #Given: Number of time button clicked
    #Returns: Scatter plot
    #Requires UI elements: output$plot5
    #Requires reactive functions: showScatterPlot()
    output$scatter_plot <- renderUI({
      n_clicks = showScatterPlot() #records number of times button clicked
      if (n_clicks > 0){
        if (n_clicks %% 2 == 0){
          return(NULL)
        }else{
          out = list(plotOutput("plot5"))
          out[[2]] = tags$div(title=paste0("Click to show/hide the graphical options"),
                              actionButton("show_scatter_plot_options","Show options"))
          out[[3]] = uiOutput("scatter_plot_options")
          
          return(out)
        }
      }else{
        return(NULL)
      }
    })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_scatter_plot
      #Returns: The number of clicks on input$show_scatter_plot
      #Requires reactive variables: input$show_scatter_plot
      showScatterPlot <- eventReactive(input$show_scatter_plot, {
        return(input$show_scatter_plot)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showScatterPlot()
      #Requires reactive variables: input$show_scatter_plot
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_scatter_plot)
        n_clicks = showScatterPlot() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_scatter_plot",
                               label = "Show scatter plot")
          }else{
            updateActionButton(session, "show_scatter_plot",
                               label = "Hide scatter plot")        
          }
        }
      })    
      #Given: Number of time button clicked
      #Returns: Heatmap options
      #Requires UI elements: 
      #Requires reactive functions: showScatterPlotOptions()
      #Requires reactive variables: user_input$scatter_plot_log, user_input$scatter_plot_line
      output$scatter_plot_options <- renderUI({
        n_clicks = showScatterPlotOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            out = list(tags$div(title="Use log10 expression",
                                checkboxInput("scatter_plot_log", label = "Use log expression", value = user_input$scatter_plot_log)))
            out[[2]] = tags$div(title="Add a dashed line for x=y",
                                checkboxInput("scatter_plot_line", label = "Plot x=y", value = user_input$scatter_plot_line))
            return(out)
          }
        }else{
          return(NULL)
        }
      })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_scatter_plot_options
      #Returns: The number of clicks on input$show_scatter_plot_options
      #Requires reactive variables: input$show_scatter_plot_options
      showScatterPlotOptions <- eventReactive(input$show_scatter_plot_options, {
        return(input$show_scatter_plot_options)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showScatterPlotOptions()
      #Requires reactive variables: input$show_scatter_plot_options
      observe({ #Updates label depending on whether options shown or not
        req(input$show_scatter_plot_options)
        n_clicks = showScatterPlotOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_scatter_plot_options",
                               label = "Show options")
          }else{
            updateActionButton(session, "show_scatter_plot_options",
                               label = "Hide options")        
          }
        }
      })   
      
      
      
      #Pairwise comparison plot of expression
      #Given: Number of time button clicked
      #Returns: Pairwise plot
      #Requires UI elements: output$plot7
      #Requires reactive functions: showPairwisePlot()
      output$pairwise_plot <- renderUI({
        n_clicks = showPairwisePlot() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            out = list(uiOutput("select_comparisons_paiwise"))
            out[[2]] = (plotOutput("plot7"))
            out[[3]] = textOutput("plot7_text")
            out[[4]] = tags$div(title=paste0("Click to show/hide the graphical options"),
                                actionButton("show_pairwise_plot_options","Show options"))
            out[[5]] = uiOutput("pairwise_plot_options")
            
            return(out)
          }
        }else{
          return(NULL)
        }
      })
        output$select_comparisons_paiwise <- renderUI({
          out = list(selectInput("select_comp_pairwise", label = h4("Select comparison 2"), 
                                 choices = comparisonNames(), #user_input$allComparisonNames, 
                                 selected = 1))
          return(out)
        })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_pairwise_plot
      #Returns: The number of clicks on input$show_pairwise_plot
      #Requires reactive variables: input$show_pairwise_plot
      showPairwisePlot <- eventReactive(input$show_pairwise_plot, {
        return(input$show_pairwise_plot)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showPairwisePlot()
      #Requires reactive variables: input$show_pairwise_plot
      observe({ #Updates label depending on whether homology shown or not
        req(input$show_pairwise_plot)
        n_clicks = showPairwisePlot() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_pairwise_plot",
                               label = "Show pairwise comparison plot")
          }else{
            updateActionButton(session, "show_pairwise_plot",
                               label = "Hide pairwise comparison plot")        
          }
        }
      })    
      #Given: Number of time button clicked
      #Returns: Heatmap options
      #Requires UI elements: 
      #Requires reactive functions: showPairwisePlotOptions()
      #Requires reactive variables: 
      output$pairwise_plot_options <- renderUI({
        n_clicks = showPairwisePlotOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            return(NULL)
          }else{
            out = list(tags$div(title="Show all data in grey",
                                checkboxInput("pairwise_plot_original", label = "Show all data", value = user_input$pairwise_plot_original)))
            out[[2]] = tags$div(title="Add a dashed line for x=y",
                                checkboxInput("pairwise_plot_line", label = "Plot x=y", value = user_input$pairwise_plot_line))
            return(out)
          }
        }else{
          return(NULL)
        }
      })
      #Functions to control show/hide button  
      #Given: Clicks from input$show_pairwise_plot_options
      #Returns: The number of clicks on input$show_pairwise_plot_options
      #Requires reactive variables: input$show_pairwise_plot_options
      showPairwisePlotOptions <- eventReactive(input$show_pairwise_plot_options, {
        return(input$show_pairwise_plot_options)
      })  
      #Given: The number of clicks on button
      #Returns: Updates the name on the button
      #Requires reactive functions: showPairwisePlotOptions()
      #Requires reactive variables: input$show_pairwise_plot_options
      observe({ #Updates label depending on whether options shown or not
        req(input$show_pairwise_plot_options)
        n_clicks = showPairwisePlotOptions() #records number of times button clicked
        if (n_clicks > 0){
          if (n_clicks %% 2 == 0){
            updateActionButton(session, "show_pairwise_plot_options",
                               label = "Show options")
          }else{
            updateActionButton(session, "show_pairwise_plot_options",
                               label = "Hide options")        
          }
        }
      })   
  
      
      
  ##Enrichment analysis TAB - UI + GENERATE & RENDER OUTPUT

  #The widgets only take primary ID so it is necessary to convert to it
  #Given: Currently selected data
  #Returns: list of 
  # raw - vector of currently selected gene ids
  # primary - vector of currently selected primary gene ids 
  #Requires reactive functions: filteredData(), 
  #Requires reactive variables: user_input$species
  #Requires functions:speciesNames(), formatLookup(), getAnnot()
  return_primary_gene_ID <- reactive({
    data.object = filteredData()
    species_data = speciesNames(user_input$species)
    
    format = data.object$data_format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    
    data = data.object$filtered_data
    gene_names = as.character(data[gene_id_name][[1]])
    gene_names2 = gene_names
    
    #convert to primary id
    if (species_data$lookup_use == "secondary"){
      if (user_input$species == "Human"){
        data = getAnnot("trans_table","Human",offline = TRUE)
      }else{
        data = getAnnot("trans_table",user_input$species) 
      }
      if (is.null(data)){
        return (NULL)
      }
      trans_table = data$table
      
      gene_names = trans_table[match(gene_names,trans_table[,2]),1]
      gene_names = gene_names[!is.na(gene_names)]
    }
    return(list("raw" = gene_names2, "primary" = gene_names))
  }) 
  return_primary_gene_ID_orthologs <- reactive({
    data.object = filteredData()
    species_data = speciesNames(user_input$species)
    
    format = data.object$data_format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    
    data = data.object$filtered_data
    gene_names = as.character(data[gene_id_name][[1]])
    gene_names2 = gene_names
    
    #convert to primary id
    if (species_data$lookup_use == "secondary"){
      if (user_input$species == "Human"){
        data = getAnnot("trans_table","Human",offline = TRUE)
      }else{
        data = getAnnot("trans_table",user_input$species) 
      }
      if (is.null(data)){
        return (NULL)
      }
      trans_table = data$table
      
      gene_names = trans_table[match(gene_names,trans_table[,2]),1]
      gene_names = gene_names[!is.na(gene_names)]
    }
    return(list("raw" = gene_names2, "primary" = gene_names))
  }) 
  
  output$analysis_extra <- renderUI({
    return(NULL)
  })
  output$orthologs_extra <- renderUI({
    return(NULL)
  })
  
  #Gene lists
  output$gene_list_table_ui_extra <- renderUI({
    return(NULL)
    # out = list(tags$div(title="Click to produce Venn diagram from selected lists",
    #           actionButton(inputId = "venn_gene",label = "Draw Venn diagram")))
    # return(out)
  })
  data_cache$gene_list_table_columns = c(1:14)
  

  #Create lists of genes and metadata to create a gene list from
  #Returns: a list of 
  # gene_names - vector of gene IDs
  # metadata - list of species, length of gene_names, input file name, ui$s1_fpkm_num,ui$s2_fpkm_num,
  # ui$log2FC_num,ui$sig_checkbox,select2InputFoundIn,select2InputNotFoundIn,user_input$select_gene_list,description
  #Requires reactive functions: filteredData(), collectUser()
  #Requires reactive variables: user_input$species, input$gprofile_enrichment_table_rows_selected
  #Requires functions:speciesNames(), formatLookup()
  createGeneListInput <- reactive ({
    if (is.null(modifyData())){
      return(NULL)
    }else{
      data.object = filteredData()
      species_data = speciesNames(user_input$species)
      ui = collectUser()
      
      format = data.object$data_format
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      
      list_name = input$create_gene_list_name
      
      if (user_input$gene_list_genes_source == "from_enrich"){
        gene_names = make_gprofile_enrichment_table()$genes_used
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_annotations")
      }else if (user_input$gene_list_genes_source == "from_enrich_selected"){
        gene_names = make_gprofile_enrichment_table()$genes_used
        gene_names = gene_names[input$gprofile_enrichment_table_rows_selected]
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_annotations")
      }else if (user_input$gene_list_genes_source == "from_enrich_ipd"){
        gene_names = make_publication_enrichment_table()$genes_used
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_publications")
      }else if (user_input$gene_list_genes_source == "from_enrich_selected_ipd"){
        gene_names = make_publication_enrichment_table()$genes_used
        gene_names = gene_names[input$protein_domain_enrichment_table_rows_selected]
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_publications")
      }else if (user_input$gene_list_genes_source == "from_ortholog"){
        gene_names = make_ortholog_table()$genes_used
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_orthologs")
      }else if (user_input$gene_list_genes_source == "from_ortholog_selected"){
        gene_names = make_ortholog_table()$genes_used
        gene_names = gene_names[input$ortholog_table_rows_selected]
        gene_names = unlist(strsplit(as.character(gene_names), ","))
        gene_names = unique(as.vector(gene_names))
        #restore gene names - per species
        gene_names = gsub("WBGENE","WBGene",gene_names) #Worm
        
        list_name = paste0(list_name,"_orthologs")
      }else{ #"standard"
        data = data.object$filtered_data
        gene_names = (data[gene_id_name][[1]])
      }
      
      #sanitise select2Inputs
      if(is.null(ui$select2InputFoundIn)){
        select2InputFoundIn = "NULL"
      }else{
        select2InputFoundIn = paste0(as.character(ui$select2InputFoundIn),collapse=", ")
      }
      if(is.null(ui$select2InputNotFoundIn)){
        select2InputNotFoundIn = "NULL"
      }else{
        select2InputNotFoundIn = paste0(as.character(ui$select2InputNotFoundIn),collapse=", ")
      }
      
      description = "NULL"
      other_files = "NULL"
      
      metadata = list(ui$species, length(gene_names),input$file$name,ui$select_comp,ui$s1_fpkm_num,ui$s2_fpkm_num,
                      ui$log2FC_num,ui$sig_checkbox,select2InputFoundIn,select2InputNotFoundIn,user_input$select_gene_list,description,other_files) #ui$select2InputFoundIn + ui$select2InputNotFoundIn not returning anything
      #FC mod
      if (ui$log2FC_up_checkbox & ui$log2FC_down_checkbox){
        metadata[7] = paste0(metadata[7])
      }else if(ui$log2FC_up_checkbox){
        metadata[7] = paste0(metadata[7]," (Up)")
      }else if(ui$log2FC_down_checkbox){
        metadata[7] = paste0(metadata[7]," (Down)")
      }else{
        metadata[7] = 0
      }
      #Comp in mod
      if (ui$invertFCcomp_checkbox){
        metadata[9] = paste0(metadata[9]," (Opposite FC)")
      }
      
      metadata$gene_list_name = list_name #must be last element added for output$gene_list_table to work
      if (user_input$gene_list_genes_source == "from_enrich"){
        if (is.null(metadata$gene_list_name)){
          metadata$gene_list_name = paste0("Gene_list_with_annotations")
        }else{
          metadata$gene_list_name = paste0(metadata$gene_list_name,"_with_annotations")
        }
      }
      
      #gene list name must be unique
      existing_names = make_gene_list_table()
      if(!is.null(existing_names)){
        existing_names = existing_names[,1]
        i = 2
        temp_name = metadata$gene_list_name
        while(!is.na(match(temp_name,existing_names))){
          temp_name = paste0(metadata$gene_list_name,"_",i)
          i = i+1
        }
        metadata$gene_list_name = temp_name
      }
      
      # if (is.null(data_cache$gene_lists)){
      #   data_cache$gene_lists = list(gene_names)
      #   data_cache$gene_lists_metadata = list(metadata)
      # }else{
      #   data_cache$gene_lists[[(length(data_cache$gene_lists) + 1)]] = (gene_names)
      #   data_cache$gene_lists_metadata[[(length(data_cache$gene_lists_metadata) + 1)]] = (metadata)
      # }
      return(list(gene_names = gene_names, metadata = metadata))
    }
  })
  

  #Generate data about filterfile genes
  #Given: gene list from data_cache$gene_lists, data.object from filteredData()
  #Returns: information about the genes taken from the data.object
  #Requires reactive functions: geneListNames(), filteredData()
  #Requires reactive variables: user_input$select_gene_list2, data_cache$gene_lists
  #Requires functions: speciesNames(), formatLookup()
  filterFileData <- reactive({
    data.object = filteredData()
    raw_gene_names = data.object$modified_data
    unused_genes = data.object$unused_genes
    
    species = data.object$species 
    link_root_gene = (speciesNames(species))$link_root_gene
    lookup = data.object$lookup 
    
    which = match(user_input$select_gene_list2,geneListNames())
    filter_genes = data_cache$gene_lists[[which]]
    filter_genes_link = createLink(paste(link_root_gene,filter_genes,sep=""),filter_genes) #this may not be the correct way to do this in non-Worm species, check
    
    format = data.object$data_format 
    format_details = formatLookup(format)
    exp_unit = format_details$exp_unit
    gene_id_position = format_details$gene_id_position
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    FC_position = format_details$FC_position
    locus_position = format_details$locus_position
    condition_1_position = format_details$condition_1_position
    condition_2_position = format_details$condition_2_position
    gene_symbol_position = format_details$gene_symbol_position
    
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    extra_colname = format_details$extra_colname
    
    if (!is.na(gene_symbol_position)){
      gene_symbol_name = correct_colnames[gene_symbol_position]
      other_genes = paste(raw_gene_names[match(filter_genes,raw_gene_names[,gene_id_name]),gene_symbol_name],raw_gene_names[match(filter_genes,raw_gene_names[,gene_id_name]),extra_colname],sep=", ")
      #remove double NAs
      tr = as.logical(match(other_genes,"NA, NA"))
      tr[is.na(tr)] = FALSE
      other_genes[tr] = NA
    }else{
      other_genes = raw_gene_names[match(filter_genes,raw_gene_names[,gene_id_name]),extra_colname]
    }
    
    #remove NAs from names
    other_genes = gsub(", NA$","",other_genes)
    other_genes = gsub("^NA, ","",other_genes)
    
    #initiate boolean for whether genes found in filtered data or not and first columns of table    
    bool_used = rep(1,length(filter_genes))
    if (is.null(lookup)){
      bool_used[match(unused_genes,filter_genes)] = 0   
      out_table = cbind(filter_genes,as.logical(bool_used))   
    }else{
      bool_used[match(unused_genes,filter_genes)] = 0
      out_table = cbind(filter_genes_link,other_genes,as.logical(bool_used))
    }
    
    adjusted_p_position = format_details$adjusted_p_position
    if (is.na(adjusted_p_position)){
      sig_name = "significant"
    }else{
      sig_position = format_details$sig_position
      sig_name = correct_colnames[sig_position]
    }
    if (!is.na(condition_1_position)){
      condition_1_name = correct_colnames[condition_1_position]
      condition_2_name = correct_colnames[condition_2_position]
      
      #set up output table
      ttable = raw_gene_names[1,c(condition_1_name,condition_2_name,value_1_name,value_2_name,sig_name)] #can do this from a single row because sample_1+2 are factorised
      conditions = unique(c(levels(ttable[condition_1_name][[1]]),levels(ttable[condition_2_name][[1]])))
      if (!is.na((match("WT",conditions)))){
        if (length(which(conditions == "WT")) < 2){
          conditions = c(conditions[(match("WT",conditions))],conditions[-(match("WT",conditions))]) #move WT to first, only if only one WT
        }
      }
    }else{
      conditions = c("sample_1","sample_2")
    }
    
    for (i in 1:length(conditions)){
      out_table = cbind(out_table,rep(NA,length(filter_genes)))
    }
    if (is.null(lookup)){
      colnames(out_table) = c(gene_id_name,"found in filtered data",paste0(conditions," ",exp_unit))
    }else{
      colnames(out_table) = c(gene_id_name,"Other gene IDs","found in filtered data",paste0(conditions," ",exp_unit))
    }

    #get values
    if (length(grep("WT",conditions)) > 0){
      versus = conditions[grep("WT",conditions)[1]]
    }else{
      versus = conditions[1] #should really be whatever current condition is 
    }
    
    if(length(filter_genes) > 100){ #needs testing to see what this value should actually be
      v_pos = match(versus,conditions) + 3
      if (!is.na(condition_1_position)){
        #subset by versus each way
        rgn = subset(raw_gene_names, raw_gene_names[condition_1_name] == versus)
        k = match(rgn[,gene_id_name],filter_genes)
        l = 3 + match(rgn[,condition_2_name],conditions)
        
        for (i in 1:dim(rgn)[1]){
          out_table[k[i],l[i]] = rgn[i,value_2_name]
          out_table[k[i],v_pos] = rgn[i,value_1_name]
        }
        
        rgn = subset(raw_gene_names, raw_gene_names[condition_2_name] == versus)
        k = match(rgn[,gene_id_name],filter_genes)
        l = 3 + match(rgn[,condition_1_name],conditions)
        
        for (i in 1:dim(rgn)[1]){
          out_table[k[i],l[i]] = rgn[i,value_1_name]
          out_table[k[i],v_pos] = rgn[i,value_2_name]
        }
      }else{
        k = match(rgn[,gene_id_name],filter_genes)
        l = 3 + match(rgn[,condition_2_name],conditions)
        
        for (i in 1:dim(rgn)[1]){
          out_table[k[i],l[i]] = rgn[i,value_2_name]
          out_table[k[i],v_pos] = rgn[i,value_1_name]
        }
      }
    }else{
      ##This code is faster for small gene lists but does not scale
      if (!is.na(condition_1_position)){
        for (k in 1:length(filter_genes)){
          gene = filter_genes[k]
          rows = which(gene == raw_gene_names[gene_id_name])
          table = raw_gene_names[rows,c(condition_1_name,condition_2_name,value_1_name,value_2_name,sig_name)]

          y_vals = rep(NA,length(conditions))
          DE = rep(FALSE,length(conditions))
          
          for (i in 1:length(conditions)){ #this assumes all samples normalised together
            if (is.na(match(conditions[i],table[condition_1_name][[1]]))){
              if (is.na(match(conditions[i],table[condition_2_name][[1]]))){
                y_vals[i] = NA
              }else{
                y_vals[i] = table[match(conditions[i],table[condition_2_name][[1]]),value_2_name]
              }
            }else{
              y_vals[i] = table[match(conditions[i],table[condition_1_name][[1]]),value_1_name]
            }
            #is DE vs 'versus'?
            if (conditions[i] != versus){
              row = which(table[,condition_1_name] == versus &(table[,condition_2_name] == conditions[i]))
              if (length(row) > 0){
                if (table[row,sig_name] == 'yes'){
                  DE[i] = TRUE
                }
              }else{
                row = which(table[,condition_2_name] == versus &(table[,condition_1_name] == conditions[i]))
                if (length(row) > 0){
                  if (table[row,sig_name] == 'yes'){
                    DE[i] = TRUE
                  }
                }
              }
            }
          }
          y_vals = round(y_vals,digits=2)
          if (is.null(lookup)){
            out_table[k,3:dim(out_table)[2]] = as.numeric(y_vals)
          }else{
            out_table[k,4:dim(out_table)[2]] = as.numeric(y_vals)
          }
        }
      }else{
        for (k in 1:length(filter_genes)){
          gene = filter_genes[k]
          rows = which(gene == raw_gene_names[gene_id_name]) #should only be one, as only one comparison made
          
          if (length(rows) > 1){ #should really throw an error here
            rows = rows[1]
          }
          
          y_vals = raw_gene_names[rows,c(value_1_name,value_2_name)]
          DE = raw_gene_names[rows,sig_name]
          
          y_vals = round(y_vals,digits=2)
          if (is.null(lookup)){
            out_table[k,3:dim(out_table)[2]] = as.numeric(y_vals) #throws error if this is not as.numeric, don't understand why
          }else{
            out_table[k,4:dim(out_table)[2]] = as.numeric(y_vals) #throws error if this is not as.numeric, don't understand why
          }
        }
      }
    }
    out = list(table = out_table, unused = unused_genes, "lookup" = lookup)
    return(out)
  }) 


  #Show per gene tab for selected gene when selected on filterfile table
  observeEvent(input$genelistgenestable_rows_selected,{
    if (is.null(input$genelistgenestable_rows_selected)){
    }else{
      if (length(input$genelistgenestable_rows_selected) > 0){
        user_input$genelistgenestable_rows_selected = input$genelistgenestable_rows_selected
      }
      format = user_input$format
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      
      data = filterFileData()
      table = data$table 
      row = user_input$genelistgenestable_rows_selected
      
      isoform_id_position = format_details$isoform_id_position
      if (is.null(isoform_id_position)){
        gene = table[row,gene_id_name]
        gene = removeLink(gene)
      }else{
        isoform_id_name = correct_colnames[isoform_id_position]
        gene = table[row,isoform_id_name]
      }
      
      user_input$pergene_text = gene
      updateTabsetPanel(session,"tabset1",selected="Per gene")
    }
  })
    
  
  #Get the background population for hypergeometric test
  
  hyper_test_input_background <- reactive({
    data.object = modifyData()
    
    format = data.object$data_format #user_input$format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    
    data = data.object$modified_data
    
    background = length(unique(data[gene_id_name][[1]],incomparables=c(NA,NULL,"-")))
    return(background)
  })
  
  
  
  ##INITIAL DATA PROCESSING
  #Load dataset
  #Given: file (input$file), user_input$format, user_input$species
  #Returns: NULL if no input file, or input file wrong dimensions, or wrong column names, else
  #  list (refered to as data.object from now on) of $raw_data - data.frame of the input file
  #  $species - vector of species name (user_input$species)
  #  $data_type - vector of 'bulkRNASeq'
  #  $data_subtype - vector of 'gene' or 'transcript'
  #  $data_format - vector of e.g. cufflinks (user_input$format)
  # user_input$pergene_cexx set to 10 if transcript data
  rawData <- reactive({
    inFile <- input$file
    if (is.null(inFile)){ #Prior to loading a dataset
      return(NULL)
    }
    format = user_input$format
    format_details = formatLookup(format)
    n_columns = format_details$n_columns
    correct_colnames = format_details$correct_colnames
    
    data = read.csv(inFile$datapath,sep="\t") #Load for testing if right dimensions
    
    #Basic check of whether right data format was chosen
    #Reset some variables
    #Create data object
    if (dim(data)[2] != n_columns){ #Wrong dimensions
      return(NULL)
    }else{
      # #load properly
      # data = read.csv(inFile$datapath,sep="\t",as.is=3)
      
      #check correct column names for format
      if (sum(colnames(data) == correct_colnames) == length(correct_colnames)){
      }else{
        return(NULL)
      }
      
      #set subtype
      if (is.null(formatLookup(format)$isoform_id_position)){
        data_subtype = 'gene'
      }else{
        data_subtype = 'transcript'
        user_input$pergene_cexx = 10
      }
      
      #Reset incase they are not applicable to new file
      #'Other comparison' select boxes
      user_input$select2InputFoundIn = NULL
      user_input$select2InputNotFoundIn = NULL
      #UI bits
      user_input$pergene_text = NULL
      
      
      data = list('raw_data' = data, 'species' = user_input$species, 'data_type' = 'bulkRNASeq', 'data_subtype' =  data_subtype, 'data_format' = format)
      return(data)
    }
  })
  
  #Get all sample names
  sampleNames <- reactive({
    if (is.null(rawData())){ #Prior to loading a dataset
      return(NULL)
    }else{
      data.object = rawData()
      raw_gene_names = data.object$raw_data
      format = data.object$data_format
      
      format_details = formatLookup(format)
      condition_1_position = format_details$condition_1_position
      condition_2_position = format_details$condition_2_position
      correct_colnames = format_details$correct_colnames
      
      if (!is.na(condition_1_position)){
        condition_1_name = correct_colnames[condition_1_position]
        condition_2_name = correct_colnames[condition_2_position]
        
        conditions = unique(c(levels(factor(raw_gene_names[condition_1_name][[1]])),levels(factor(raw_gene_names[condition_2_name][[1]]))))
      }else{
        condition_1_name = "sample_1"
        condition_2_name = "sample_2"
        conditions = c(condition_1_name,condition_2_name)
      }
      
      if (!is.na((match("WT",conditions)))){
        if (length(which(conditions == "WT")) < 2){
          conditions = c(conditions[(match("WT",conditions))],conditions[-(match("WT",conditions))]) #move WT to first, only if only one WT
        }
      }
      #fix pairing of ilk17OE - replace with generalised pairing fix
      if (conditions[length(conditions)] == "ilk17OE"){
        pair = match("ilk17",conditions)
        if (length(conditions) - 1 != pair){ #they aren't already paired
          conditions = c(conditions[1:pair],conditions[length(conditions)],conditions[(pair+1):(length(conditions)-1)])
        }
      }
      
      return(conditions)
    }
  })
  
  #Check columns and calculate gene names where necessary and cache
  #Given: data object rawData()
  #Returns: NULL if rawData() = NULL, else data object rawData(), with $modified_data containing translation of gene names if possible
  #Returns reactive variables: boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #     was the other of primary/secondary identifier as defined by InterMine - FALSE
  #     was neither so couldn't translate to the lookup identifier - NULL (data_cache$lookup);
  ##  where lookup = FALSE a vector of the lookup_identifier in table row order (data_cache$new_column);
  #  boolean of whether couldn't connect to InterMine version (data_cache$connectionError);
  #  date of the time connection was tested (data_cache$connectionErrorTime);
  #  time before allow a retry if a connection to intermine failed (data_cache$connectionRetryTime);
  #  boolean of whether couldn't connect to species specfic method - starts FALSE (data_cache$connectionErrorSS);
  #  boolean of whether could identify gene IDs - only ever NULL or TRUE (data_cache$geneIDsUnknown);
  #  boolean of whether could offline and no backup Translation table - only ever NULL or TRUE (data_cache$noTransTable)
  #  vector of possible gene names in linkless_table (data_cache$names_list);
  #  vector of primary gene ID for each gene names in data_cache$names_list (data_cache$trans_names_list)
  #Requires reactive functions: rawData()
  #Requires functions: checkNameColumns()
  modifyData <- reactive({
    if (is.null(rawData())){ #Prior to loading a dataset
      return(NULL)
    }else{
      data_cache$data_loaded = TRUE
      
      conditions = sampleNames() 
      if (max(nchar(conditions)) > 8){
        user_input$rotate_sample = TRUE
      }
      user_input$select2InputSampleOrder = conditions
      user_input$select2InputSampleOrder_heatmap = conditions
      
      data.object = rawData()
      raw_data = data.object$raw_data
      species = data.object$species
      data_subtype = data.object$data_subtype
      format = data.object$data_format
      
      format_details = formatLookup(format)
      n_columns = format_details$n_columns
      starting_checks = format_details$starting_checks
      correct_colnames = format_details$correct_colnames
      gene_id_position = format_details$gene_id_position
      gene_id_name = correct_colnames[gene_id_position]

      #Create binary "Is DE significant?" column if doesn't exist
      if (is.na(match("significant",correct_colnames))){
        threshold = 0.05 #default adjusted p_value threshold
        
        #Identify column to threshold
        adjusted_p_position = format_details$adjusted_p_position
        
        #Create new column
        #if is NA assume all DE is significant
        if (is.na(adjusted_p_position)){
          new_column = rep(TRUE,dim(data)[1])
        }else{
          new_column = (data[,adjusted_p_position] <= threshold)
        }
        new_column = factor(new_column,levels=c(FALSE,TRUE),labels=c("no","yes"))
        
        #Add column at end + update n_columns
        new_column_name = "significant"
        new_column_number = n_columns + 1
        starting_checks = c(starting_checks,TRUE) #make this column shown by default
        
        n_columns = n_columns + 1
        raw_data = cbind(raw_data,new_column)
        colnames(raw_data)[n_columns] = new_column_name
      }else{ 
        #change sig NAs to no
        col_n = match("significant",correct_colnames)
        if(sum(is.na(raw_data[,col_n])) > 0){
          raw_data[(is.na(raw_data[,col_n])),col_n] = "no"
        }
      }
      
      
      #initialise column_checks
      for (i in 1:(n_columns + 1)){
        user_input[[paste0("col",i,"_checkbox")]] = starting_checks[i]
      }
      
      #work out which columns are which and do translation if possible
      data = checkNameColumns(raw_data,species,format)
      if (is.null(data)){ #species not found / Other
        data.object$file_ref = "OFFLINE MODE"
        
        #Fix gene_id column to remove "Gene:"
        raw_data[gene_id_name] = sub("Gene:","",raw_data[gene_id_name][[1]])
        raw_data[gene_id_name] = sub("gene:","",raw_data[gene_id_name][[1]])
        
        if (data_subtype == 'transcript'){
          isoform_id_position = format_details$isoform_id_position
          isoform_id_name = correct_colnames[isoform_id_position]
          raw_data[isoform_id_name] = sub("CDS:","",raw_data[isoform_id_name][[1]])      
          raw_data[isoform_id_name] = sub("Transcript:","",raw_data[isoform_id_name][[1]])     
          #raw_data[isoform_id_name] = sub("Pseudogene:","",raw_data[isoform_id_name][[1]])      
        }
        
        #reset = resetValues() #this was only working the first time a switch occured so I moved all into this function
        #Bits related to data
        #data_cache$new_column = NULL
        data_cache$lookup = NULL
        data.object$lookup = NULL
        
        #Bits related to connection
        data_cache$connectionError = NULL
        data_cache$connectionErrorTime = NULL
        data_cache$connectionErrorSS = NULL
        data_cache$connectionRetryTime = NULL
        
        #Bits related to species
        data_cache$built_in_checkboxes_values = NULL
        data_cache$built_in_checkboxes = NULL
        data_cache$built_in_data = NULL
        
        #UI bits
        user_input$pergene_text = NULL
        
        #Other bits
        
        
        data.object$modified_data = raw_data
        #data.object$linked_modified_data = raw_data
        #data_cache$new_column = data$new_column
        data_cache$lookup = NULL
        data.object$lookup = NULL
        
        if (input$species != "Other"){ 
          data_cache$geneIDsUnknown = TRUE
        }
        
        #TESTING ERROR MESSAGES
        #data_cache$connectionError = TRUE
        #data_cache$noTransTable = TRUE
        #data_cache$connectionErrorSS = TRUE
        
        #make name lists
        names_list_list = makeNameLists(data.object,data.object$lookup)
        names_list = names_list_list$names_list
        names_list2 = names_list_list$names_list2
        
        data_cache$names_list = names_list          #all_gene_names variants or transcript_name variants
        data_cache$trans_names_list = names_list2   #the associated gene_ids or transcript_ids
        
        data.object$raw_data = NULL
        return(data.object)
      }else if (is.null(data$lookup)){ #offline and no translation table
        if (is.null(data$file_ref)){
          data.object$file_ref = "OFFLINE MODE"
        }else{
          data.object$file_ref = data$file_ref
        }
          
        #Fix gene_id column to remove "Gene:"
        raw_data[gene_id_name] = sub("Gene:","",raw_data[gene_id_name][[1]])
        raw_data[gene_id_name] = sub("gene:","",raw_data[gene_id_name][[1]])
        
        if (data_subtype == 'transcript'){
          isoform_id_position = format_details$isoform_id_position
          isoform_id_name = correct_colnames[isoform_id_position]
          raw_data[isoform_id_name] = sub("CDS:","",raw_data[isoform_id_name][[1]])      
          raw_data[isoform_id_name] = sub("Transcript:","",raw_data[isoform_id_name][[1]])     
          #raw_data[isoform_id_name] = sub("Pseudogene:","",raw_data[isoform_id_name][[1]])      
        }
        
        reset = resetValues()
        
        data.object$modified_data = raw_data
        #data.object$linked_modified_data = raw_data
        #data_cache$new_column = data$new_column
        data_cache$lookup = NULL
        data.object$lookup = NULL
        
        data_cache$noTransTable = TRUE
        
        #make name lists
        names_list_list = makeNameLists(data.object,data.object$lookup)
        names_list = names_list_list$names_list
        names_list2 = names_list_list$names_list2
        
        data_cache$names_list = names_list          #all_gene_names variants or transcript_name variants
        data_cache$trans_names_list = names_list2   #the associated gene_ids or transcript_ids
        
        data.object$raw_data = NULL
        return(data.object)       
      }else{
        data.object$file_ref = data$file_ref
        
        new_raw = data$table
        connection_error = data$connection_error
        no_link_rows = data$no_link_rows
        new_column = data$new_column
        lookup = data$lookup
        test_time = data$time
        
        data.object$modified_data = new_raw
        #data_cache$new_column = new_column
        data_cache$lookup = lookup
        data.object$lookup = lookup
        data_cache$connectionError = connection_error
        data_cache$connectionErrorTime = test_time
        data_cache$connectionErrorSS = FALSE
        
        data_cache$connectionRetryTime = 5 #mins since connection failure before allow retries
        
        # #make version of modified data with links
        # if (is.null(lookup)){
        #   table = new_raw
        # }else{
        #   data = createTableLinks(species,new_raw,no_link_rows,lookup,new_column,format)
        #   table = data$table
        # }
        # data.object$linked_modified_data = table
        data.object$no_link_rows = no_link_rows

        #make name lists
        names_list_list = makeNameLists(data.object,data.object$lookup)
        names_list = names_list_list$names_list
        names_list2 = names_list_list$names_list2
        
        data_cache$names_list = names_list          #all_gene_names variants or transcript_name variants
        data_cache$trans_names_list = names_list2   #the associated gene_ids or transcript_ids
        
        data.object$raw_data = NULL
        return(data.object)
      }
    }
  })
    #Reset variables when species changes
    resetValues <- reactive({
      #Bits related to data
      #data_cache$new_column = NULL
      data_cache$lookup = NULL
      
      #Bits related to connection
      data_cache$connectionError = NULL
      data_cache$connectionErrorTime = NULL
      data_cache$connectionErrorSS = NULL
      data_cache$connectionRetryTime = NULL
      
      #Bits related to species
      data_cache$built_in_checkboxes_values = NULL
      data_cache$built_in_checkboxes = NULL
      data_cache$built_in_data = NULL
      
      #UI bits
      user_input$pergene_text = NULL
      
      #Other bits
      
    })

  #Generate data for table, summary (same as table without annotations and column filtering)
  #and download (without html links)  
  #Add annotations and filter data
  #Given: data.object from modifyData, details of filtering from collectUser, annotations from builtIn and 
  #Returns: modified data.object
  #Requires reactive functions: modifyData, collectUser
  #Requires reactive variables: user_input$format
  #Requires functions: formatLookup, addAnnotationToRaw
  filteredData <- reactive({
    if (is.null(modifyData())){ #Prior to file loading
      return(NULL)
    }else{
      data.object = modifyData()
      format = data.object$data_format #user_input$format
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      
      ui = collectUser()
      species = data.object$species
      link_root_gene = (speciesNames(species))$link_root_gene
      
      lookup = data.object$lookup #data_cache$lookup
      
      
      #Filter based on UI
      #Which comparison and multiple comparisons filter
      if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
        return("Invalid comparison")
      }else{
        data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
      }
      
      data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)
      
      #Apply Filter_file if present and selected
      if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
        data.object$filtered_data = data.object$filtered_data
        data.object$unused_genes = NULL
      }else{
        which = match(user_input$select_gene_list,geneListNames())
        gene_list = data_cache$gene_lists[[which]]
        new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
        if (ui$invertfilter_file_checkbox){
          data.object$filtered_data = new_data$nonmatches_have
          #used_genes = new_data$nonmatches_want
          data.object$unused_genes = gene_list
        }else{
          data.object$filtered_data = new_data$matches
          data.object$unused_genes = new_data$nonmatches_want
        }
      }
      
      
      #update no_link_rows to their new positions in the filtered data
      #mapping = row.match(data.object$modified_data,data.object$filtered_data)
      #data.object$no_link_rows = (mapping[data.object$no_link_rows])[!is.na(mapping[data.object$no_link_rows])]
      data.object$no_link_rows = NULL
      
      #user_input$nrows2 = dim(data.object$filtered_data)[1]
      
      #Add annotations
      #Built-in
      if (is.null(lookup)){ #no annotation
      }else{
        for (i in 1:length((data_cache$built_in_checkboxes_values))){
          value = data_cache$built_in_checkboxes_values[i]
          #value = TRUE #FALSE #NULL
          if (value){ #Want data
            update = getBuiltIn() #this caches the data if it hasn't already been done
            
            new_data = data_cache$built_in_data[[i]]
            #new_data = NULL
            if (is.null(new_data)){
            }else{
              data.object$filtered_data = addAnnotationToRaw(new_data,data.object$filtered_data,lookup=lookup,format = format, species = species)
            }
          }
        }
      }
      #Other
      #Add annotation
      update = queryAnnotationSets() #loads annotation files
      if (!is.null(global_annotation$annotation_sets)){
        if (length(global_annotation$annotation_sets) > 0){ #necessary because when removing all annot sets does not go back to null
          for (i in 1:length(global_annotation$annotation_sets)){
            if (dim(global_annotation$annotation_sets[[i]])[2] > 1){
              data.object$filtered_data = addAnnotationToRaw(global_annotation$annotation_sets[[i]],data.object$filtered_data,lookup=lookup,format = format, species = species)
            }
          }
        }
      }
      
      return(data.object)
    }
  }) 
  
  
  
  ##Generate annotation  
  #Initialise 'built-in' annotations  
  #Given: user_input$species
  #Returns: boolean of FALSE if no built-in annotations for species, else TRUE
  #Returns reactive variables: data_cache$built_in_checkboxes, data_cache$built_in_data
  #Requires reactive variables: user_input$species
  #Requires functions: speciesNames()
  initialiseBuiltIn <- reactive({
    species = user_input$species
    species_data = speciesNames(species)
    if (is.null(species_data$built_in_list)){
      return(FALSE)
    }else{
      data_cache$built_in_checkboxes = species_data$built_in_list
      data_cache$built_in_data[[length(data_cache$built_in_checkboxes_values) + 1]] = c(TRUE,FALSE) #sets up list so data is NULL so all values in built_in_checkboxes, there's probably a more sensible way of doing this
      return(TRUE)
    }
  })
    
  #Cache 'built-in' annotations for all genes, saves having to do it again
  #Given: Booleans of checkboxes of whether to load data (data_cache$built_in_checkboxes_values);
  #  info on whether data can be looked up (data_cache$lookup);
  #  currently loaded data to see if the relavant data has already been cached (data_cache$built_in_data)
  #Returns: TRUE if data_cache$built_in_checkboxes_values is NULL (no built-in data for this species) 
  #  or data_cache$lookup is NULL (no web connection or unknown species)
  #  or couldn't connect to server,
  #  else FALSE (if found web data and loaded it - n.b. may not actually be any annotation for the given genes)
  #Returns reactive variables: data_cache$built_in_data, data_cache$connectionError, data_cache$connectionErrorTime
  #Requires reactive functions: initialiseBuiltIn()
  #Requires reactive variables: data_cache$built_in_checkboxes_values, data_cache$lookup, data_cache$built_in_data
  #Requires functions: getAnnot()
  getBuiltIn <- reactive({
    init = initialiseBuiltIn() #initialises the variables
    if (is.null(data_cache$built_in_checkboxes_values)){
      return(TRUE)  
    }else{
      for (i in 1:length((data_cache$built_in_checkboxes_values))){
        value = data_cache$built_in_checkboxes_values[i]
        if (value){ #Want data
          if (is.null(data_cache$built_in_data[[i]])){ #Want data for the first time 
            species = user_input$species
            lookup = data_cache$lookup
            
            if (is.null(lookup)){ #this shouldn't ever be reached as should be handled before this function gets called
              #data_cache$connectionError = TRUE
              return(TRUE)            
            }else{
              name = (data_cache$built_in_checkboxes)[i]
              
              new_genes_data = getAnnot(name,species) #annotateBuiltIn(genes,species,name)

              if(new_genes_data$connection_error){ #couldn't connect to server
                data_cache$connectionError = TRUE
                data_cache$connectionErrorTime = new_genes_data$time
                return(TRUE)
              }else{
                data_cache$connectionError = FALSE
                data_cache$built_in_data[[i]] =  new_genes_data$table
              }
              return(FALSE)
            }
          }
        }   
      }
    }
  })
  
  #Load annotation file  
  # Check column names against existing columns and change duplicates
  #Given: vector of filename of annotations (input$annot_file),
  # data.object from modifyData()
  #Returns: NULL if input$annot_file is NULL or list of
  #  data.frame of lookup_identifier gene identifiers and annotation(s) (annotation) filtered by data.object, 
  #  vector of filename (name)
  #Requires reactive functions: modifyData()
  #Requires reactive variables: global_annotation$annotation_sets
  annotationFile <- reactive({
    inFile <- input$annot_file
    if (is.null(inFile)){ #No file selected
      return(NULL)
    }else{
      annotation = read.csv(inFile$datapath,as.is=1,header=T,sep="\t") 
      
      data.object = modifyData()
      data = data.object$modified_data
      format = data.object$data_format
      
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      extra_position = format_details$extra_position
      
      lookup = data.object$lookup #data_cache$lookup
      
      #Filter for annotation matching raw_data
      if (is.null(lookup)){ #don't know which column to merge with, guess column 1
        by.y = 1
      }else{
        if (!lookup){
          by.y=extra_position
        }else{
          by.y=gene_id_position
        }
      }
      want = data[,by.y]
      split_annotation = validateFilter(want,annotation)
      annotation = split_annotation$matches
      
      #colnames(annotation) = c(colnames(data)[2],rep("annotation",(dim(annotation)[2]-1))) #Kept incase ever want to allow no header loading this causes warnings due to duplicated column names, no errors in my code, but DateTables code may have issues
      if (is.null(lookup)){ #don't know which column to rename
      }else{
        if (!lookup){
          colnames(annotation)[1] = colnames(data)[extra_position]
        }else{
          colnames(annotation)[1] = colnames(data)[gene_id_position] #should this be data_cache$lookup dependant?
        }
      }
            
      if(!is.na(annotation)){
        #Make unique column names
        #check if other columns have allowable names
        existing_col_names = colnames(data)
        #add from existing annotation
        annotations = isolate(global_annotation$annotation_sets)
        if (is.null(annotations)){}else{ #No annotations loaded yet
          if (length(annotations) > 0){ #No annotations, due to them all being deleted
            for (i in 1:length(annotations)){ 
              new_col_names = colnames(annotations[[i]])
              new_col_names = new_col_names[-1]
              existing_col_names = c(existing_col_names,new_col_names)
            }
          }
        }
        to_change = intersect(colnames(annotations),existing_col_names) #n.b. loading doesn't allow annotation to have identical colnames so don't worry about that. It itterates with .n so don't use that format here
        if (length(to_change) > 0){
          j = 1
          for (i in 1:length(to_change)){
            new_colname = paste0(to_change[i],'_',j) #really should do this properly and check for existance of new colname but should be sufficiently uncommon and duplicated colnames doesn't actually break any of my code
            j = j + 1
            
            to_replace_index = which(colnames(annotation) == to_change[i])
            colnames(annotation)[to_replace_index] = new_colname
          }
        }
      }
      
      output = list(annotation = annotation, name = (input$annot_file)$name)
      return(output)
    }
  })  

  #Update/Retrieve global annotations  
  #Given: list of
  #  data.frame of lookup_identifier gene identifiers and annotation(s), 
  #  vector of filename (name) from annotationFile()
  #Returns: NULL if annotationFile is NULL or list of
  #  list of list of data.frames of as above with new one added (annotation_sets);
  #  list of vectors of names  with new one added (annotation_sets_names)
  #Returns reactive variables: global_annotation$annotation_sets, global_annotation$annotation_sets_names
  #Requires reactive functions: annotationFile()
  #Requires reactive variables: global_annotation$annotation_sets, global_annotation$annotation_sets_names
  queryAnnotationSets <- reactive({
    if (is.null(annotationFile())){ #No file selected
      return(NULL) 
    }else{
      input = annotationFile()
      annotation_sets = isolate(global_annotation$annotation_sets) #isolate necessary to prevent recursion
      annotation_sets_names = isolate(global_annotation$annotation_sets_names)
      
      #former combineAnnotationSets
      if(is.null(annotation_sets)){ #nothing to add
        annotation_sets = list(input$annotation)
        annotation_sets_names = list(input$name)
      }else{
        annotation_sets = append(annotation_sets,list(input$annotation))
        annotation_sets_names = append(annotation_sets_names,list(input$name))        
      }

      global_annotation$annotation_sets = annotation_sets
      global_annotation$annotation_sets_names = annotation_sets_names
      
      output = list('annotation_sets' = annotation_sets, 'annotation_sets_names' = annotation_sets_names)
      return(output)
    }
  })
  
  #Watch for annotation_removal_buttons activation - only works for the number of annotations you have below, couldn't get it to do these automatically without running into an index error when removing last annotation  
  #if you add more search for max_annotation_remove_buttons and update it with the right number
  #Given: boolean of button pressed/not (e.g. input$button1)
  #Returns reactive variables: global_annotation$annotation_sets, global_annotation$annotation_sets_names
  #  without the list element whose number corresponds to the button
  #Requires reactive variables: global_annotation$annotation_sets, global_annotation$annotation_sets_names
  observeEvent(input$button1,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 1
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()
  })
  observeEvent(input$button2,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 2
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button3,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 3
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button4,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 4
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button5,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 5
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button6,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 6
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button7,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 7
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button8,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 8
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button9,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 9
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })
  observeEvent(input$button10,{
    annotation_sets_names = (global_annotation$annotation_sets_names)
    annotation_sets = (global_annotation$annotation_sets)
    
    i = 10
    
    annotation_sets_names[[i]] = NULL
    annotation_sets[[i]] = NULL
    
    global_annotation$annotation_sets = annotation_sets
    global_annotation$annotation_sets_names = annotation_sets_names
    #anotes()()
  })

  
  ##Plots
  #Get FPKM for all comparisons featuring specific gene  
  #Given: data.object from modifyData() and specific gene name (user_input$pergene_text)
  #Returns: graph of the expression in all samples for that gene, 
  #  with asterixes to show if expression is DE relative to a given sample 
  #  and with the tallest bar removed if user_input$remove_max = TRUE
  #  with y-axis log10 scale if user_input$gene_exp_log_scale = TRUE
  #  and y values < 0 removed if user_input$gene_exp_log_scale_lt_zero = TRUE
  #Requires reactive functions: modifyData()
  #Requires reactive variables: user_input$pergene_text, user_input$remove_max, user_input$gene_exp_log_scale, user_input$gene_exp_log_scale_lt_zero
  #  (input$select_DE)
  #Requires functions: formatLookup()
  plot2_plot <- function(){
    if (is.null(modifyData())){ #Prior to file loading
      return(NULL)
    }else{
      data.object = modifyData()
      format = data.object$data_format
      lookup = data.object$lookup
      
      format_details = formatLookup(format)
      exp_unit = format_details$exp_unit
      gene_id_position = format_details$gene_id_position
      value_1_position = format_details$value_1_position
      value_2_position = format_details$value_2_position
      FC_position = format_details$FC_position
      locus_position = format_details$locus_position
      condition_1_position = format_details$condition_1_position
      condition_2_position = format_details$condition_2_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      value_1_name = correct_colnames[value_1_position]
      value_2_name = correct_colnames[value_2_position]
      
      adjusted_p_position = format_details$adjusted_p_position
      if (is.na(adjusted_p_position)){
        sig_name = "significant"
      }else{
        sig_position = format_details$sig_position
        sig_name = correct_colnames[sig_position]
      }
      gene_symbol_position = format_details$gene_symbol_position
      if (!is.na(gene_symbol_position)){
        gene_symbol_name = correct_colnames[gene_symbol_position]
      }else{
        gene_symbol_name = "gene_symbol"
      }
      
      raw_gene_names = data.object$modified_data
      gene = user_input$pergene_text
      versus = input$select_DE
      
      #Check gene in table - allow different versions of name
      isoform_id_position = format_details$isoform_id_position
      if (is.null(isoform_id_position)){
        geneLookupID = findGene(gene,raw_gene_names,format = format,lookup)
      }else{
        isoform_id_name = correct_colnames[isoform_id_position]
        geneLookupID = findTranscript(gene,raw_gene_names,format = format,lookup)
        transcript = gene
        gene = geneLookupID
        if(is.null(gene)){
          geneLookupID = NULL
        }else{
          geneLookupID = lookup
        }
      }
      
      if(is.null(geneLookupID)){ #not in table
        return(NULL)
      }else{
        rows = which(gene == raw_gene_names[gene_id_name])
      }
      
      if (length(rows) < 1){
        return(NULL)
      }else{
        data_cache$empty_graph = FALSE
        if (is.null(lookup)){
          if (!is.na(gene_symbol_position)){
            gene_symbol = raw_gene_names[(rows[1]),gene_symbol_name]
            gene_symbol = as.character(gene_symbol)
            
            #choose what name to display on plot title
            if (nchar(gene) <= nchar(gene_symbol)){
              print_gene = gene
            }else{
              print_gene = gene_symbol
            }
          }else{
            print_gene = gene
          }
        }else{
          gene_symbol = raw_gene_names[(rows[1]),gene_symbol_name]
          gene_symbol = as.character(gene_symbol)
          
          #choose what name to display on plot title
          if (nchar(gene) <= nchar(gene_symbol)){
            print_gene = gene
          }else{
            print_gene = gene_symbol
          }
        }
        
        if (!is.na(condition_1_position)){
          condition_1_name = correct_colnames[condition_1_position]
          condition_2_name = correct_colnames[condition_2_position]
          
          table = raw_gene_names[rows,c(condition_1_name,condition_2_name,value_1_name,value_2_name,sig_name)]
        }else{
          condition_1_name = "sample_1"
          condition_2_name = "sample_2"
          conditions = c(condition_1_name,condition_2_name)
          c1_col = rep(conditions[1],rows)
          c2_col = rep(conditions[2],rows)
          
          table = cbind(c1_col,c2_col,raw_gene_names[rows,c(value_1_name,value_2_name,sig_name)])
          colnames(table)[1] = condition_1_name
          colnames(table)[2] = condition_2_name
        }
        
        #set y-axis scale
        if(!is.null(user_input$gene_exp_log_scale)){
          if (user_input$gene_exp_log_scale){
            table[,value_1_name] = log10(table[,value_1_name])
            table[,value_2_name] = log10(table[,value_2_name])
            exp_unit = paste0("log10(",exp_unit,")")
            
            #reset infinites
            table[is.infinite(table[,value_1_name]),value_1_name] = 0
            table[is.infinite(table[,value_2_name]),value_2_name] = 0
          }
        }
        if(!is.null(user_input$gene_exp_log_scale_lt_zero)){
          if (user_input$gene_exp_log_scale_lt_zero){
            table[(table[,value_1_name] < 0),value_1_name] = 0
            table[(table[,value_2_name] < 0),value_2_name] = 0
          }
        }
        
        #set conditions from UI
        conditions = input$select2InputSampleOrder
        if(is.null(conditions)){
          return(NULL)
        }
        
        
        if (is.null(isoform_id_position)){
          y_vals = rep(NA,length(conditions))
          DE = rep(FALSE,length(conditions))
          cols = rep('black',length(conditions)) #black will be the one that DE is calculated versus
          if(length(conditions) > 1){
            pal = palette(rainbow(length(conditions)))
          }else{
            pal = "black"
          }
          j = 1 #which colour up to
          
          for (i in 1:length(conditions)){ #this assumes all samples normalised together
            if (is.na(match(conditions[i],table[condition_1_name][[1]]))){
              if (is.na(match(conditions[i],table[condition_2_name][[1]]))){
                y_vals[i] = NA
              }else{
                y_vals[i] = table[match(conditions[i],table[condition_2_name][[1]]),value_2_name]
              }
            }else{
              y_vals[i] = table[match(conditions[i],table[condition_1_name][[1]]),value_1_name]
            }
            #is DE vs 'versus'?
            if (conditions[i] != versus){
              row = which(table[condition_1_name] == versus &(table[condition_2_name] == conditions[i]))
              if (length(row) > 0){
                if (table[row,sig_name] == 'yes'){
                  DE[i] = TRUE
                }
              }else{
                row = which(table[condition_2_name] == versus &(table[condition_1_name] == conditions[i]))
                if (length(row) > 0){
                  if (table[row,sig_name] == 'yes'){
                    DE[i] = TRUE
                  }
                }
              }
              cols[i] = pal[j]
              j = j + 1
            }
          }
          #pair up colours for mutant/OE pairs
          withOE = grep("([A-Za-z0-9]+)OE$",conditions,perl=T) #indices of conditions ending in OE
          possibleMutantPairs = regmatches(conditions[withOE],regexpr("[^OE]+",conditions[withOE])) #those condition names minus the 'OE'
          if (length(possibleMutantPairs) > 0){
            for (i in 1:length(possibleMutantPairs)){
              withoutOE = which(conditions == possibleMutantPairs[i])
              if (length(withoutOE) > 0){ #then we have a pair
                if (((which(conditions == versus)) != withoutOE)){ #don't change if matching to the black DE marker
                  cols[withoutOE] = cols[withOE[i]]
                  if (withOE[i] == (which(conditions == versus))){
                    cols[withoutOE] = "white"
                  }
                }else{
                  cols[withOE[i]] = "white"
                }
              }
            }
          }
          
          #make those without comparisons grey
          
          for (i in 1:length(conditions)){ #this assumes all samples normalised together
            #is DE vs 'versus'?
            if (conditions[i] != versus){
              row = which(table[condition_1_name] == versus &(table[condition_2_name] == conditions[i]))
              row2 = which(table[condition_2_name] == versus &(table[condition_1_name] == conditions[i]))
              if ((length(row) == 0)&(length(row2) == 0)){
                cols[i] = "grey"
              }
            } 
          }
          
          
          
          max_y = max(y_vals,na.rm=T)
          
          if (length(user_input$remove_max) > 0){
            if (user_input$remove_max){
              old_max_y = match(max_y,y_vals)
              y_vals[old_max_y] = NA
              max_y = max(y_vals,na.rm=T)
            }
          }
          
          plot = barplot(y_vals,names.arg = conditions, ylab = exp_unit,main=paste0("Expression of ",print_gene," under conditions in dataset"),ylim=c(0,(ceiling(max_y)*1.13)),col=cols,xaxt="n")
          
          #add DE *'s
          x_vals = plot[DE,1]
          points(x = x_vals, y = rep((ceiling(max_y)*1.10),length(x_vals)),pch='*',cex=2)
          
          #add values
          value_pos = max((y_vals + (ceiling(max_y)*0.05)),0.1)
          text(x = plot[,1], y = value_pos,labels=round(y_vals,2))
          
          #add text when a bar has been removed
          if (length(user_input$remove_max) > 0){
            if (user_input$remove_max){
              text(x = plot[old_max_y,1], y = max_y / 2, labels = "Omitted") 
            }
          }
          
          if (user_input$rotate_sample){
            #set angle for names based on length of longest condition name
            degrees = 45
            if (max(nchar(conditions)) > 8){
              degrees = degrees - ((max(nchar(conditions)) - 8)*2) #empirically based
            }
            if (degrees < 10){
              degrees = 10
            }
          }else{
            degrees = 0
          }
          
          #determine position of names below graph
          name_pos = min(-(max_y/10),-0.15)
          
          text(x= plot[,1], y=name_pos, conditions,srt=degrees, xpd=TRUE, cex = user_input$pergene_cexx) #rotates the names so that they all fit on in low resolutions
          
        }else{ #isoform view
          
          if (!is.na(condition_1_position)){
            condition_1_name = correct_colnames[condition_1_position]
            condition_2_name = correct_colnames[condition_2_position]
            
            table = raw_gene_names[rows,c(condition_1_name,condition_2_name,value_1_name,value_2_name,sig_name,isoform_id_name)]
          }else{
            condition_1_name = "sample_1"
            condition_2_name = "sample_2"
            c1_col = rep(conditions[1],rows)
            c2_col = rep(conditions[2],rows)
            
            table = cbind(c1_col,c2_col,raw_gene_names[rows,c(value_1_name,value_2_name,sig_name,isoform_id_name)])
            colnames(table)[1] = condition_1_name
            colnames(table)[2] = condition_2_name
          }
          
          #set y-axis scale
          if (user_input$gene_exp_log_scale){
            table[,value_1_name] = log10(table[,value_1_name])
            table[,value_2_name] = log10(table[,value_2_name])
            #exp_unit = paste0("log10(",exp_unit,")")
          }
          if (user_input$gene_exp_log_scale_lt_zero){
            table[(table[,value_1_name] < 0),value_1_name] = 0
            table[(table[,value_2_name] < 0),value_2_name] = 0
          }
          
          alt_splice = unique(table[isoform_id_name][[1]])
          variants = length(alt_splice) * length(conditions)
          
          #convert table into useful 
          #condition  alt_splice  y-val   DE
          #e.g. WT    14a         2       TRUE
          table2 = data.frame(matrix(NA,nrow = variants,ncol = 4), stringsAsFactors=FALSE)

          for (i in 1:variants){ #this assumes all samples normalised together
            var_n = (i %% length(alt_splice))
            con_n = (i %/% length(alt_splice))
            con_n = con_n + 1
            if (var_n == 0){var_n = length(alt_splice)
            con_n = con_n - 1}
            alt = alt_splice[var_n]
            con = conditions[con_n]
            
            table2[i,1] = con
            table2[i,2] = alt
            
            #subset by alt
            rows = which(alt == table[,6])
            temp_table = table[rows,]
            
            
            #subset by con
            if(is.na(match(con,temp_table[,1]))){
              if(is.na(match(con,temp_table[,2]))){
                table2[i,3] = 0
              }else{
                table2[i,3] = temp_table[match(con,temp_table[,2]),4]
              }
            }else{
              table2[i,3] = temp_table[match(con,temp_table[,1]),3]
            }
            
            if (con != versus){
              row = which(table[condition_1_name] == versus &(table[condition_2_name] == con))
              if (length(row) > 0){
                row = (row[which(table[row,6] == alt)])
                if (length(row > 0)){
                  if (table[row,sig_name] == 'yes'){
                    table2[i,4] = TRUE
                  }
                }
              }else{
                row = which(table[condition_2_name] == versus &(table[condition_1_name] == con))
                if (length(row) > 0){
                  row = (row[which(table[row,6] == alt)])
                  if (length(row > 0)){
                    if (table[row,sig_name] == 'yes'){
                      table2[i,4] = TRUE
                    }
                  }
                }
              }
            }
          }
          
          colnames(table2)=c("condition","isoforms","yval","DE")
          
          
          max_y = max(table2[,3],na.rm = T)*1.10
          
          for (i in 1:variants){ #this assumes all samples normalised together
            if (!is.na(table2[i,4])){
              table2[i,4] = (max_y) #table2[i,3] + 1
            }
          }
          table2[,1] = factor(table2[,1], levels = conditions)

          
          if (length(user_input$remove_max) > 0){
            if (user_input$remove_max){
              row = which(table2[,3] == max(table2[,3],na.rm = T))
              row = row[1] #must be singular otherwise plotting issues
              rows = which(table2[,1] == table2[row,1]) #removes all isoforms for that condition
              table2[rows,3] = NA
              max_y = max(table2[,3],na.rm=T)*1.10 #recalculate max value
            }
          }
          
          
          
          if (length(alt_splice) < 9){ #palette Set1 only has 8 colours (+ grey which blends into background so I skip)
            plot = ggplot(table2, aes(x = condition, y = yval, fill = isoforms)) + 
              geom_bar(stat="identity", position = "dodge") + 
              scale_fill_brewer(palette = "Set1") +
              theme(axis.text.x = element_text(size=user_input$pergene_cexx)) +
              geom_point(aes(y = as.numeric(DE)),colour = "black", size = 10, position = position_dodge(width = .9),shape ="*") + #, y = DE
              geom_text(aes(label = round(yval,digits =1)),colour = "black", size = 5, position = position_dodge(width = .9), vjust=-0.5,check_overlap = TRUE) + 
              labs(y = exp_unit, x = "")
          }else{ #allows for more than 8 isoforms to be coloured
            plot = ggplot(table2, aes(x = condition, y = yval, fill = isoforms)) + 
              geom_bar(stat="identity", position = "dodge") + 
              theme(axis.text.x = element_text(size=user_input$pergene_cexx)) +
              geom_point(aes(y = as.numeric(DE)),colour = "black", size = 10, position = position_dodge(width = .9),shape ="*") + #, y = DE
              geom_text(aes(label = round(yval,digits =1)),colour = "black", size = 5, position = position_dodge(width = .9), vjust=-0.5,check_overlap = TRUE) + 
              labs(y = exp_unit, x = "")
          }
          
          
          
          if (length(user_input$remove_max) > 0){
            if (user_input$remove_max){
              plot = plot + geom_text(data = table2, aes(x = table2[row,1], y = max_y/4, label = "Omitted"),colour = "black", size = 5, vjust=-0.5,check_overlap = TRUE, angle = 90) 
            }
          }
          
        }
        
        return(plot)
        
      }
    }
  }
  output$plot2 <- renderPlot({
    plot2_plot()
  })
  
  #Volcano plot for selected data
  #Given: data.object from modifyData(), format and ui
  #Returns: graph of the log2 fold change against -log10 pvalue
  #Requires reactive functions: modifyData(), collectUser(), pAdjThreshold(), geneListNames()
  #Requires reactive variables: user_input$select_gene_list, data_cache$gene_lists
  #Requires functions: formatLookup(), combineComparisons(), applyFilters()
  plot_plot <- function(){
    data.object = modifyData()
    ui = collectUser()
    format = data.object$data_format
    p_threshold = pAdjThreshold()
    
    format_details = formatLookup(format)
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    adjusted_p_position = format_details$adjusted_p_position
    FC_position = format_details$FC_position
    gene_id_position = format_details$gene_id_position
    
    correct_colnames = format_details$correct_colnames
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    adjusted_p_name = correct_colnames[adjusted_p_position]
    FC_name = correct_colnames[FC_position]
    gene_id_name = correct_colnames[gene_id_position]
    
    #Filter based on UI
    #Which comparison and multiple comparisons filter
    if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
      return(NULL)
    }else{
      data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
    }
    
    if(user_input$volcano_plot_grey_data){ #ploting original data in background
      original_data = combineComparisons(data.object$modified_data,ui = ui,format = format,main_only=TRUE)
    }
    if(user_input$volcano_plot_no_rearrange | user_input$volcano_plot_grey_data){ #using min max from original data
      data = combineComparisons(data.object$modified_data,ui = ui,format = format,main_only=TRUE)
      
      #update the positions
      value_1_position = match(value_1_name,colnames(data))
      value_2_position = match(value_2_name,colnames(data))
      adjusted_p_position = match(adjusted_p_name,colnames(data))
      FC_position = match(FC_name,colnames(data))
      
      
      #ignore Inf where calculating max/min x-axis
      data2 = data[(as.numeric(data[,value_1_position]) > 1),]
      data2 = data2[(as.numeric(data2[,value_2_position]) > 1),]
      if (length(data2[,FC_position]) > 0){
        x_min = min(data2[,FC_position],na.rm=T)
        x_max = max(data2[,FC_position],na.rm=T)
      }else{
        x_min = min(data[,FC_position],na.rm=T)
        x_max = max(data[,FC_position],na.rm=T)
      }
      y_min = 0
      y_max = max(as.numeric(-log10(data[,adjusted_p_position])),na.rm = TRUE)
    }
    
    data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)

    if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
    }else{
      which = match(user_input$select_gene_list,geneListNames())
      gene_list = data_cache$gene_lists[[which]]
      new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
      if (ui$invertfilter_file_checkbox){
        data.object$filtered_data = new_data$nonmatches_have
        #used_genes = new_data$nonmatches_want
        data.object$unused_genes = gene_list
      }else{
        data.object$filtered_data = new_data$matches
        data.object$unused_genes = new_data$nonmatches_want
      }
    }
    
    data = data.object$filtered_data
    
    #update the positions
    value_1_position = match(value_1_name,colnames(data))
    value_2_position = match(value_2_name,colnames(data))
    adjusted_p_position = match(adjusted_p_name,colnames(data))
    FC_position = match(FC_name,colnames(data))
    
    
    if (!user_input$volcano_plot_no_rearrange){ #only using subset
      if (!user_input$volcano_plot_grey_data){
        #ignore Inf where calculating max/min x-axis
        data2 = data[(as.numeric(data[,value_1_position]) > 1),]
        data2 = data2[(as.numeric(data2[,value_2_position]) > 1),]
        if (length(data2[,FC_position]) > 0){
          x_min = min(data2[,FC_position],na.rm=T)
          x_max = max(data2[,FC_position],na.rm=T)
        }else{
          x_min = min(data[,FC_position],na.rm=T)
          x_max = max(data[,FC_position],na.rm=T)
        }
        y_min = min(as.numeric(-log10(data[,adjusted_p_position])),na.rm = TRUE)
        y_max = max(as.numeric(-log10(data[,adjusted_p_position])),na.rm = TRUE)
      }
    }
    
    main_name = paste0("Volcano plot of ",user_input$select_comp)
    
    if (user_input$volcano_plot_grey_data){
      #remove data from original data
      to_remove = which(data[gene_id_name][[1]] %in% original_data[gene_id_name][[1]])
      original_data = original_data[-to_remove,]
      
      plot(original_data[,FC_position],-log10(original_data[,adjusted_p_position]),xlim=c(x_min,x_max),ylim=c(y_min,y_max),xlab="log2FoldChange",ylab="-log10(pvalue)",main=main_name,col="grey",pch='.',cex=2)
      points(data[,FC_position],-log10(data[,adjusted_p_position]))
    }else{
      plot(data[,FC_position],-log10(data[,adjusted_p_position]),xlim=c(x_min,x_max),ylim=c(y_min,y_max),xlab="log2FoldChange",ylab="-log10(pvalue)",main=main_name)
    }
    points(data[(data[,adjusted_p_position] < p_threshold),FC_position],-log10(data[(data[,adjusted_p_position] < p_threshold),adjusted_p_position]),col="red")
    points(data[(abs(data[,FC_position]) > 1),FC_position],-log10(data[(abs(data[,FC_position]) > 1),adjusted_p_position]),col="orange")
    points(data[((data[,adjusted_p_position] < p_threshold)&(abs(data[,FC_position]) > 1)),FC_position],-log10(data[((data[,adjusted_p_position] < p_threshold)&(abs(data[,FC_position]) > 1)),adjusted_p_position]),col="green")
    
    
  }
  output$plot <- renderPlot({
    plot_plot()
  })
  
  #Hierarchical clustering for significantly expressed genes in selected data
  #Given: matrix of FC and expressions for selected genes from genes_for_h_clust()
  #Returns: Hierarchical clustering of the significantly DE genes
  #Requires reactive functions: genes_for_h_clust()
  plot3_plot <- function(){
    data2 = genes_for_h_clust()
    if (!is.null(data2)){
      clusters = hclust(dist(data2[,1]))
      tryCatch(plot(clusters,xlab=""),error=function(err){return(NULL)}) #catches a weird error returning "invalid dendrogram input" that I can't figure out - occurs with low n of genes
    }else{
      return(NULL)
    }
  }
  output$plot3 <- renderPlot({
    plot3_plot()
  })

  #Heatmap for selected genes
  #Given: matrix from genes_for_heatmap
  #Returns: Heat map of the genes
  #Requires reactive functions: genes_for_heatmap()
  #Requires reactive variables: user_input$heatmap_no_rearrange
  plot4_plot <- function(){
    data = genes_for_heatmap()

    if (is.null(data)){
      return(NULL)
    }
    data4 = data$table
    data5 = data$unmod
    
    if (user_input$heatmap_no_rearrange){
      #heatmap.2(data4,cexCol=1,margins=c(12,8),tracecol=NA,na.color='black')
      #there's a bug with giving the function NAs that I can't track down. For now tryCatch deals with this issue
      tryCatch(heatmap.2(data4,cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8),na.color='black',tracecol=NA),error=function(e){heatmap.2(data5,cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8),na.color='black',tracecol=NA)})
    }else{
      #heatmap.2(data4,cexCol=1,margins=c(12,8),tracecol=NA,na.color='black',Rowv = FALSE, Colv=FALSE,dendrogram='none')  
      tryCatch(heatmap.2(data4,cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8),tracecol=NA,na.color='black',Rowv = FALSE, Colv=FALSE,dendrogram='none'),error=function(e){heatmap.2(data4,cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8),na.color='black',tracecol=NA,Rowv = FALSE, Colv=FALSE,dendrogram='none')})
    }
    
  }
  output$plot4 <- renderPlot({
    plot4_plot()
  })
  
  #Scatter plot for selected data
  #Given: data.object from modifyData(), format and ui
  #Returns: graph of the sample 1 expression vs the sample 2 expression
  #Requires reactive functions: modifyData(), collectUser(), pAdjThreshold(), geneListNames()
  #Requires reactive variables: user_input$select_gene_list, data_cache$gene_lists
  #Requires functions: formatLookup(), combineComparisons(), applyFilters()
  plot5_plot <- function(){
    data.object = modifyData()
    ui = collectUser()
    format = data.object$data_format
    p_threshold = pAdjThreshold()
    
    format_details = formatLookup(format)
    exp_unit = format_details$exp_unit
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    adjusted_p_position = format_details$adjusted_p_position
    FC_position = format_details$FC_position
    gene_id_position = format_details$gene_id_position
    
    correct_colnames = format_details$correct_colnames
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    adjusted_p_name = correct_colnames[adjusted_p_position]
    FC_name = correct_colnames[FC_position]
    gene_id_name = correct_colnames[gene_id_position]
    
    #Filter based on UI
    #Which comparison and multiple comparisons filter
    if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
      return(NULL)
    }else{
      data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
    }
    
    original_data = combineComparisons(data.object$modified_data,ui = ui,format = format,main_only=TRUE)

    data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)
    
    if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
    }else{
      which = match(user_input$select_gene_list,geneListNames())
      gene_list = data_cache$gene_lists[[which]]
      new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
      if (ui$invertfilter_file_checkbox){
        data.object$filtered_data = new_data$nonmatches_have
        #used_genes = new_data$nonmatches_want
        data.object$unused_genes = gene_list
      }else{
        data.object$filtered_data = new_data$matches
        data.object$unused_genes = new_data$nonmatches_want
      }
    }
    
    data = data.object$filtered_data
    
    comparison = user_input$select_comp
    main_name = paste0("Scatter plot of ",comparison)
    s1_name = strsplit(as.character(comparison)," vs ")[[1]][1] #Sample 1 name
    s2_name = strsplit(as.character(comparison)," vs ")[[1]][2] #Sample 2 name
    xlab_name=paste0(s1_name," expression ","(",exp_unit ,")")
    ylab_name=paste0(s2_name," expression ","(",exp_unit ,")")
    
    if (user_input$scatter_plot_log){
      original_data[,value_1_name] = log10(original_data[,value_1_name])
      original_data[,value_2_name] = log10(original_data[,value_2_name])
      data[,value_1_name] = log10(data[,value_1_name])
      data[,value_2_name] = log10(data[,value_2_name])
      xlab_name=paste0(s1_name," expression ","(log10 ",exp_unit ,")")
      ylab_name=paste0(s2_name," expression ","(log10 ",exp_unit ,")")
    }
    
    
    x_min = 0
    x_max = max(original_data[,value_1_name],na.rm = TRUE)
    y_min = 0
    y_max = max(original_data[,value_2_name],na.rm = TRUE)
    
    #remove data from original data
    to_remove = which(data[gene_id_name][[1]] %in% original_data[gene_id_name][[1]])
    original_data = original_data[-to_remove,]
    
    plot(original_data[,value_1_name],(original_data[,value_2_name]),xlim=c(x_min,x_max),ylim=c(y_min,y_max),xlab=xlab_name,ylab=ylab_name,main=main_name,col="grey",pch='.',cex=2)
    points(data[,value_1_name],(data[,value_2_name]),col="red")
    
    if(user_input$scatter_plot_line){
      points(c(0,x_max),c(0,y_max),type="l",lty=2,col="blue")
    }
    
  }
  output$plot5 <- renderPlot({
    plot5_plot()
  })
  
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
    
    ComplexHeatmap::Heatmap(dot_plot_data,col=c("white","black"),cluster_columns = F, show_column_dend = F, row_names_side="left",column_names_side="top",show_heatmap_legend=F,
            show_row_dend = F,cluster_rows = T,
            column_names_gp = gpar(fontsize = 10),rect_gp = gpar(col= "white"),column_names_max_height= unit(1000,"mm"),
            row_names_gp = gpar(fontsize = 10))
    # heatmap.2(dot_plot_data,col=c("white","black"),dendrogram = "none",trace="none",key = F,
    #           sepwidth=c(0.1,0.1),
    #           sepcolor="white",
    #           colsep=1:ncol(dot_plot_data),
    #           rowsep=1:nrow(dot_plot_data),cexCol=user_input$heatmap_cexCol,cexRow = user_input$heatmap_cexRow,margins=c(12,8))
  }
  output$plot6 <- renderPlot({
    draw(plot6_plot())
  })
  
  #Pairwise comparison for selected data
  #Given: data.object from modifyData(), format and ui
  #Returns: graph of the sample 1 expression vs the sample 2 expression
  #Requires reactive functions: modifyData(), collectUser(), pAdjThreshold(), geneListNames()
  #Requires reactive variables: user_input$select_gene_list, data_cache$gene_lists
  #Requires functions: formatLookup(), combineComparisons(), applyFilters()
  plot7_plot <- function(){
    req(input$select_comp_pairwise)
    
    data.object = modifyData()
    ui = collectUser()
    format = data.object$data_format
    p_threshold = pAdjThreshold()
    
    format_details = formatLookup(format)
    exp_unit = format_details$exp_unit
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    adjusted_p_position = format_details$adjusted_p_position
    FC_position = format_details$FC_position
    gene_id_position = format_details$gene_id_position
    condition_1_position = format_details$condition_1_position
    condition_2_position = format_details$condition_2_position
    
    correct_colnames = format_details$correct_colnames
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    adjusted_p_name = correct_colnames[adjusted_p_position]
    FC_name = correct_colnames[FC_position]
    gene_id_name = correct_colnames[gene_id_position]
    condition_1_name = correct_colnames[condition_1_position]
    condition_2_name = correct_colnames[condition_2_position]
    
    #Filter based on UI
    #First comparison
    #Which comparison and multiple comparisons filter
    if((is.null(combineComparisons(data.object$modified_data,ui = ui,format = format)))){
      return(NULL)
    }else{
      data.object$filtered_data = combineComparisons(data.object$modified_data,ui = ui,format = format)
    }
    
    original_data = combineComparisons(data.object$modified_data,ui = ui,format = format,main_only=TRUE)
    
    data.object$filtered_data = applyFilters(data.object$filtered_data,ui = ui,format = format)
    
    if(user_input$select_gene_list == "None"){ #if no genelist or nonselected
    }else{
      which = match(user_input$select_gene_list,geneListNames())
      gene_list = data_cache$gene_lists[[which]]
      new_data = validateFilter(gene_list,data.object$filtered_data,gene_id_position)
      if (ui$invertfilter_file_checkbox){
        data.object$filtered_data = new_data$nonmatches_have
        #used_genes = new_data$nonmatches_want
        data.object$unused_genes = gene_list
      }else{
        data.object$filtered_data = new_data$matches
        data.object$unused_genes = new_data$nonmatches_want
      }
    }
    
    data = data.object$filtered_data
    
    #remove data from original data
    # to_remove = which(data[gene_id_name][[1]] %in% original_data[gene_id_name][[1]])
    # original_data = original_data[-to_remove,]
    
    
    #Second comparison
    comparison2 = input$select_comp_pairwise #input$select_comp
    ui2 = ui
    #ui2$select_comp = comparison2 #this is updating user_input$select_comp!!!!!!! WHY!!!!!!!!!!!S
    
    #original_data2 = combineComparisons(data.object$modified_data,ui = ui2,format = format,main_only=TRUE) #subsetByComparison(comparison2,data.object$modified_data,format)
    # data_cache$test = dim(original_data)
    # data_cache$test2 = dim(original_data2)
    
    #Which comparison and multiple comparisons filter
    # if((is.null(combineComparisons(data.object$modified_data,ui = ui2,format = format)))){
    #   return(NULL)
    # }else{
    #   filtered_data2 = combineComparisons(data.object$modified_data,ui = ui2,format = format)
    # }

    original_data2 = subsetByComparison(comparison2,data.object$modified_data,format) #combineComparisons(data.object$modified_data,ui = ui2,format = format,main_only=TRUE)


    
    #get in same order / same content
    allowed_ids = unique(intersect(factor(original_data[,1]),factor(original_data2[,1])))
    
    order1 = match(allowed_ids,as.character(original_data[,1]))
    od1 = original_data[order1,]
    od1 = od1[!is.na(od1[,1]),]
    
    order2 = match(allowed_ids,as.character(original_data2[,1]))
    od2 = original_data2[order2,]
    od2 = od2[!is.na(od2[,1]),]
    
    od2 = od2[order(allowed_ids),]
    od1 = od1[order(allowed_ids),]
    
    # data_cache$test = od1[1:10,1] #levels(factor(original_data[,7]))
    # data_cache$test2 = od2[1:10,1] #levels(factor(original_data2[,7]))

    #remove infinites
    od1 = od1[(!is.infinite(od1[,FC_name])),]
    od2 = od2[(!is.infinite(od2[,FC_name])),]
    allowed_ids = unique(intersect(factor(od1[,1]),factor(od2[,1])))
    
    order1 = match(allowed_ids,as.character(od1[,1]))
    od1 = od1[order1,]
    od1 = od1[!is.na(od1[,1]),]
    
    order2 = match(allowed_ids,as.character(od2[,1]))
    od2 = od2[order2,]
    od2 = od2[!is.na(od2[,1]),]
    
    od2 = od2[order(allowed_ids),]
    od1 = od1[order(allowed_ids),]   
    
    
    #get other comparison matching filtered data
    chosen_rows = match(data[,1],od2[,1])
    data2 = od2[chosen_rows,]
    data2 = applyFilters(data2,ui = ui,format = format)
    
    
    # data_cache$test = data[1:10,1] #levels(factor(data[,7]))
    # data_cache$test2 = data2[1:10,1] #levels(factor(data2[,7]))
    
    #get in same order / same content
    allowed_ids = unique(intersect(factor(data[,1]),factor(data2[,1])))
    
    order1 = match(allowed_ids,as.character(data[,1]))
    d1 = data[order1,]
    d1 = d1[!is.na(d1[,1]),]
    
    order2 = match(allowed_ids,as.character(data2[,1]))
    d2 = data2[order2,]
    d2 = d2[!is.na(d2[,1]),]
    
    d2 = d2[order(allowed_ids),]
    d1 = d1[order(allowed_ids),]
    
    #remove data from original data
    # to_remove = which(data2[gene_id_name][[1]] %in% od2[gene_id_name][[1]])
    # od2 = od2[-to_remove,]
    
    
    # data_cache$test = d1[1:10,1] #levels(factor(data[,7]))
    # data_cache$test2 = d2[1:10,1] #levels(factor(data2[,7]))
    
    
    #Perform correlation
    corp = cor.test(d1[,FC_name],d2[,FC_name],method="spearman")
    corp_r = corp$estimate
    corp_p = corp$p.value
    
    
    #draw graph
    comparison = user_input$select_comp
    comparison2 = input$select_comp_pairwise
    
    main_name = paste0("Scatter plot of ",comparison, " log2FC vs ", comparison2, " log2FC")
    xlab_name=paste0(comparison," log2FC")
    ylab_name=paste0(comparison2," log2FC")
    x_min = -6
    x_max = 6
    y_min = -6
    y_max = 6
    
    if (user_input$pairwise_plot_original){
      plot(od1[,FC_name],od2[,FC_name],xlim=c(x_min,x_max),ylim=c(y_min,y_max),xlab=xlab_name,ylab=ylab_name,main=main_name,col="grey",pch='.',cex=2)
      points(d1[,FC_name],(d2[,FC_name]),col="black",pch='.',cex=2)
    }else{
      plot(d1[,FC_name],d2[,FC_name],xlim=c(x_min,x_max),ylim=c(y_min,y_max),xlab=xlab_name,ylab=ylab_name,main=main_name,col="black",pch='.',cex=2)
    }
    text(-5,5.5,paste0("R = ",round(corp_r,2)))
    text(-5,4.5,paste0("p < ",round(corp_p,2)))
    
    if(user_input$pairwise_plot_line){
      points(c(x_min,x_max),c(y_min,y_max),type="l",lty=1,col="red")
    }
    
  }
  output$plot7 <- renderPlot({
    (plot7_plot())
  })
  
  
  
  ##ERROR HANDLING & OTHER POPUPS
  #Notifications
  #Web related
  #Deal with inability to identify gene IDs
  observeEvent(data_cache$geneIDsUnknown, {
    alert("Could not recognize the gene identifiers, the program will act as though the selected species was Other, i.e. online data will not be used. This may happen if you have chosen the wrong species.")
  })
  
  #Deal with connection error to InterMine
  observeEvent(data_cache$connectionError, {
    if (data_cache$connectionError){
      species = user_input$species
      data = speciesNames(species)
      mine = data$mine
      if (is.null(mine)){
        mine = "InterMine"
      }
      
      time = data_cache$connectionRetryTime
      
      if (!is.null(data_cache$noTransTable)){
        out = alert(paste0("Could not connect to ",mine," server and no local translation table found, the program will act as though the selected species was Other, i.e. online data will not be used"))
        return(out)
      }
      
      out = alert(paste0("Could not connect to ",mine," server - falling back to local translation table where possible"))
      
      return(out)
    }
  })
  
  #Deal with connection error to species specific method
  observeEvent(data_cache$connectionErrorSS, {
    species = user_input$species
    data = speciesNames(species)
    database = data$database
    if (is.null(database)){
      database = "species-specific"
    }
    mine = data$mine
    if (is.null(mine)){
      mine = "InterMine"
    }
    
    if (data_cache$connectionErrorSS){
      out = alert(paste0("Could not connect to ",database," server - falling back to ",mine," database where possible"))
      
      return(out)
    }
    
    
  })
  
  #Coding related
  #Can only add max_annotation_remove_buttons annotations
  observeEvent(global_annotation$annotation_sets,{
    i = length(global_annotation$annotation_sets)
    max_annotation_remove_buttons = 10
    if (i > max_annotation_remove_buttons){
      out = alert(paste0("Can only display removal buttons for ",max_annotation_remove_buttons," annotations at a time and you have just added another, to remove this annotation you must remove another annotation first"))
      return(out)
    }
  })
        
        
        
  
  ##Functions
  #Given: data.frame of cufflinks data (table);
  # vector of species name (species), as defined in speciesNames()
  # (local copy of transtable if getIdentifiers()$connection_error = TRUE);
  #  vector of the data format (format)
  #Returns: NULL if speciesNames()$lookup_use = NULL, else list of
  #  data.frame with gene name columns analysed and potentially changed, where lookup = FALSE
  #  includes new column of the lookup identifier as column 4 (table);
  #  boolean of whether was a problem connecting in getIdentifiers() (connection_error);
  #  date of the time connection was tested (time);
  #  vector of rows not found in getIdentifiers() translation table (no_link_rows);
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    couldn't get a translation table so have no idea - NULL (lookup);
  #  where lookup = FALSE a vector of the lookup_identifier in table row order (new_column);
  #  date of the translation_table used (file_ref)
  #Requires: speciesNames(), getAnnot() formatLookup()
  checkNameColumns = function(table,species,format){
    #format = user_input$format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    extra_position = format_details$extra_position
    extra_colname = format_details$extra_colname
    
    isoform_id_position = format_details$isoform_id_position #if this exists we are interested in transcripts not genes
    
    gene_symbol_position = format_details$gene_symbol_position
    if (!is.na(gene_symbol_position)){
      gene_symbol_name = correct_colnames[gene_symbol_position]
    }else{
      gene_symbol_name = "gene_symbol"
    }
    
    species_data = speciesNames(species)
    lookup_use = species_data$lookup_use
    if (is.null(lookup_use)){
      return(NULL)
    }else{
      #Fix gene_id column to remove "Gene:" in Cufflinks data
      table[gene_id_name] = sub("Gene:","",table[gene_id_name][[1]])
      table[gene_id_name] = sub("gene:","",table[gene_id_name][[1]])
      
      if (!is.null(isoform_id_position)){
        isoform_id_name = correct_colnames[isoform_id_position]
        table[isoform_id_name] = sub("CDS:","",table[isoform_id_name][[1]])      
        table[isoform_id_name] = sub("Transcript:","",table[isoform_id_name][[1]])     
        #table[isoform_id_name] = sub("Pseudogene:","",table[isoform_id_name][[1]])      
      }
      
      if (species == "Human"){ #humanmine missing a lot of annotation
        data = getAnnot("trans_table",species,offline = TRUE) #retrieveTranslation(species,offline = TRUE)
        connection_error = FALSE #otherwise thinks we're offline due to above lookup
      }else{
        data = getAnnot("trans_table",species) #retrieveTranslation(species)
        connection_error = data$connection_error
      }
      file_ref = data$file_ref
      new_table = data$table
      trans_table = new_table
      test_time = data$time
      if (is.null(data)){
        table = NULL
        new_column = NULL
        lookup = NULL
        rows_not = NULL
        
        if (species == 'Other'){
          return(NULL)
        }else{
          table = NULL
          new_column = NULL
          lookup = NULL
          rows_not = NULL
          return(list(table = table, connection_error = connection_error, no_link_rows = rows_not, lookup = lookup, new_column = new_column, time = test_time))
        }
      }
      
      
      #Check if gene_id column is (~all) lookup Identifier
      n_not = sum(is.na(match(table[gene_id_name][[1]],trans_table[,lookup_use]))) #in worm this is cases where WormMine has not be updated/gene is now obsolete, links still exist in WormBase
      proportion = 0.75 #require at least 90% to be matchable in InterMine lookup
      if ((1 - n_not/dim(table)[1]) > proportion){
        lookup = TRUE
        new_column = NULL
        rows_not = which(is.na(match(table[gene_id_name][[1]],trans_table[,lookup_use]))) #note these and exclude from link making
      }else{ 
        if (lookup_use == "primary"){
          n_not = sum(is.na(match(table[gene_id_name][[1]],trans_table[,"secondary"]))) 
          if ((1 - n_not/dim(table)[1]) > proportion){
            #generate a column that is lookup ("secondary")
            new_column = trans_table[(match(table[gene_id_name][[1]],trans_table[,"secondary"])),"primary"]
            lookup = FALSE
            rows_not = which(is.na(match(new_column,trans_table[,lookup_use]))) #note these and exclude from link making
          }else{ #no idea how to do lookup - suggest user set species to "Other"
            table = NULL
            new_column = NULL
            lookup = NULL
            rows_not = NULL
            return(NULL)
          }
        }else if (lookup_use == "secondary"){
          n_not = sum(is.na(match(table[gene_id_name][[1]],trans_table[,"primary"]))) 
          if ((1 - n_not/dim(table)[1]) > proportion){
            #generate a column that is lookup ("primary")
            new_column = trans_table[(match(table[gene_id_name][[1]],trans_table[,"primary"])),"secondary"]
            lookup = FALSE
            rows_not = which(is.na(match(new_column,trans_table[,lookup_use]))) #note these and exclude from link making
          }else{ #no idea how to do lookup - suggest user set species to "Other"
            table = NULL
            new_column = NULL
            lookup = NULL
            rows_not = NULL
            return(NULL)
          }       
        }else{ #no idea how to do lookup - suggest user set species to "Other"
          table = NULL
          new_column = NULL
          lookup = NULL
          rows_not = NULL
          return(NULL)
        }
      }
      
      if (!is.na(gene_symbol_position)){
        #Check if gene column is (all) geneSymbol
        n_not = sum(is.na(match(table[gene_symbol_name][[1]],trans_table[,"symbol"]))) 
        proportion2 = 1 #always add in gene symbols
        if ((1 - n_not/dim(table)[1]) > proportion2){
          symbol = "gene"
        }else{ #make it so
          if (is.null(lookup)){ #suggest user set species to "Other"
          }else{
            if (lookup){
              #add gene symbol where exists, else add otherIdentifer, else existing
              if (lookup_use == "primary"){
                new_genes = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"symbol"]
                new_genes2 = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"secondary"]
                new_genes[is.na(new_genes)] = (table[gene_id_name][[1]])[is.na(new_genes)]
                new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
                
                table[gene_symbol_name] = new_genes
              }else if (lookup_use == "secondary"){
                new_genes = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"symbol"]
                new_genes2 = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"primary"]
                new_genes[is.na(new_genes)] = (table[gene_id_name][[1]])[is.na(new_genes)]
                new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
                
                table[gene_symbol_name] = new_genes
              }
            }else{ #use new_column
              #add gene symbol where exists, else existing
              if (lookup_use == "primary"){
                new_genes = trans_table[(match(new_column,trans_table[,lookup_use])),"symbol"]
                #new_genes2 = trans_table[(match(new_column,trans_table[,lookup_use])),"secondary"]
                new_genes[is.na(new_genes)] = (new_column)[is.na(new_genes)]
                #new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
                
                table[gene_symbol_name] = new_genes
              }else if (lookup_use == "secondary"){
                new_genes = trans_table[(match(new_column,trans_table[,lookup_use])),"symbol"]
                #new_genes2 = trans_table[(match(new_column,trans_table[,lookup_use])),"primary"]
                new_genes[is.na(new_genes)] = (new_column)[is.na(new_genes)]
                #new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
                
                table[gene_symbol_name] = new_genes
              }          
            }
          }
        }
      }else{ #make gene column
        if (lookup){
          #add gene symbol where exists, else add otherIdentifer, else existing
          if (lookup_use == "primary"){
            new_genes = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"symbol"]
            new_genes2 = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"secondary"]
            new_genes[is.na(new_genes)] = (table[gene_id_name][[1]])[is.na(new_genes)]
            new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
            
            #table[gene_symbol_name] = new_genes
          }else if (lookup_use == "secondary"){
            new_genes = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"symbol"]
            new_genes2 = trans_table[(match(table[gene_id_name][[1]],trans_table[,lookup_use])),"primary"]
            new_genes[is.na(new_genes)] = (table[gene_id_name][[1]])[is.na(new_genes)]
            new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
            
            #table[gene_symbol_name] = new_genes
          }
        }else{ #use new_column
          #add gene symbol where exists, else existing
          if (lookup_use == "primary"){
            new_genes = trans_table[(match(new_column,trans_table[,lookup_use])),"symbol"]
            #new_genes2 = trans_table[(match(new_column,trans_table[,lookup_use])),"secondary"]
            new_genes[is.na(new_genes)] = (new_column)[is.na(new_genes)]
            #new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
            
            #table[gene_symbol_name] = new_genes
          }else if (lookup_use == "secondary"){
            new_genes = trans_table[(match(new_column,trans_table[,lookup_use])),"symbol"]
            #new_genes2 = trans_table[(match(new_column,trans_table[,lookup_use])),"primary"]
            new_genes[is.na(new_genes)] = (new_column)[is.na(new_genes)]
            #new_genes[(new_genes == "")] = new_genes2[(new_genes == "")]
            
            #table[gene_symbol_name] = new_genes
          }          
        }
        #add new column
        table = cbind(table,new_genes)
        colnames(table)[(dim(table)[2])] = gene_symbol_name #specified earlier in this function
      }
      
      #add new column
      if (is.null(lookup)){
      }else{
        if (!lookup){
          if (extra_position > 1){ #not first
            before = 1:(extra_position - 1)
            if (extra_position < dim(table)[2]){ #not last
              after = extra_position:dim(table)[2]
              table = cbind(table[,before],new_column,table[,after],stringsAsFactors=F)
            }else{ #only need to add last
              table = cbind(table,new_column,stringsAsFactors=F)
            }
          }else{
            table = cbind(new_column,table,stringsAsFactors=F)
          }
          
          #table = cbind(table[,1:3],new_column,table[,4:dim(table)[2]],stringsAsFactors=F)
          colnames(table)[gene_id_position] = extra_colname
          colnames(table)[extra_position] = gene_id_name
        }else{
          new_column = trans_table[(match(table[gene_id_name][[1]],trans_table[,"primary"])),"secondary"]
          
          if (extra_position > 1){ #not first
            before = 1:(extra_position - 1)
            if (extra_position < dim(table)[2]){ #not last
              after = extra_position:dim(table)[2]
              table = cbind(table[,before],new_column,table[,after],stringsAsFactors=F)
            }else{ #only need to add last
              table = cbind(table,new_column,stringsAsFactors=F)
            }
          }else{
            table = cbind(new_column,table,stringsAsFactors=F)
          }
          colnames(table)[extra_position] = extra_colname
        }
      }
      
      return(list(table = table, connection_error = connection_error, no_link_rows = rows_not, lookup = lookup, new_column = new_column, 
                  time = test_time, file_ref = file_ref))
      
    }
  }
  
  #Given: data object
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    couldn't get a translation table so have no idea - NULL (lookup);  
  #Returns: list of all_gene_names variants or transcript_name variants (names_list),
  #  the associated gene_ids or transcript_ids (names_list2)
  #Requires: formatLookup()
  makeNameLists <- function(data.object,lookup){
    if(is.null(data.object$modified_data)){
      data = data.object$raw_data
    }else{
      data = data.object$modified_data
    }
    format = data.object$data_format
    data_subtype = data.object$data_subtype
    
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    extra_colname = format_details$extra_colname
    gene_symbol_position = format_details$gene_symbol_position
    if (!is.na(gene_symbol_position)){
      gene_symbol_name = correct_colnames[gene_symbol_position]
    }else{
      gene_symbol_name = "gene_symbol"
    }
    
    if (data_subtype == 'transcript'){
      #Generate all possible transcript name variants for graph drawing
      isoform_id_position = format_details$isoform_id_position
      isoform_id_name = correct_colnames[isoform_id_position]
      if (is.null(lookup)){
        if (!is.na(gene_symbol_position)){
          names_list = unique(c("",as.character(data[gene_id_name][[1]]),as.character(data[isoform_id_name][[1]]),as.character(data[gene_symbol_name][[1]])))       
          
          names_list2 = names_list
          names_list2[match(data[,gene_symbol_name],names_list)] = data[match(data[,gene_symbol_name],data[,gene_symbol_name]),isoform_id_name]
          names_list2[match(data[,gene_id_name],names_list)] = data[match(data[,gene_id_name],data[,gene_id_name]),isoform_id_name]
        }else{
          names_list = unique(c("",as.character(data[gene_id_name][[1]])))
          names_list = unique(c("",as.character(data[gene_id_name][[1]]),as.character(data[isoform_id_name][[1]])))       
          
          names_list2 = names_list
          names_list2[match(data[,gene_id_name],names_list)] = data[match(data[,gene_id_name],data[,gene_id_name]),isoform_id_name]
        }
      }else{
        names_list = unique(c("",as.character(data[gene_id_name][[1]]),as.character(data[isoform_id_name][[1]]),as.character(data[gene_symbol_name][[1]])))       
        
        names_list2 = names_list
        #names_list2[match(data[,gene_symbol_name],names_list)] = data[match(data[,gene_symbol_name],data[,gene_symbol_name]),gene_id_name]
        names_list2[match(data[,gene_symbol_name],names_list)] = data[match(data[,gene_symbol_name],data[,gene_symbol_name]),isoform_id_name]
        names_list2[match(data[,gene_id_name],names_list)] = data[match(data[,gene_id_name],data[,gene_id_name]),isoform_id_name]
      }
    }else{
      #Generate all possible gene name variants for graph drawing
      if (is.null(lookup)){
        if (!is.na(gene_symbol_position)){
          names_list = unique(c("",as.character(data[gene_id_name][[1]]),as.character(data[gene_symbol_name][[1]])))
          
          names_list2 = names_list
          names_list2[match(data[,gene_symbol_name],names_list)] = data[match(data[,gene_symbol_name],data[,gene_symbol_name]),gene_id_name]
        }else{
          names_list = unique(c("",as.character(data[gene_id_name][[1]])))
          names_list2 = names_list
        }
      }else{
        names_list = unique(c("",as.character(data[gene_id_name][[1]]),as.character(data[extra_colname][[1]]),as.character(data[gene_symbol_name][[1]])))       
        
        names_list2 = names_list
        names_list2[match(data[,gene_symbol_name],names_list)] = data[match(data[,gene_symbol_name],data[,gene_symbol_name]),gene_id_name]
        names_list2[match(data[,extra_colname],names_list)] = data[match(data[,extra_colname],data[,extra_colname]),gene_id_name]
      }
    }
    
    #Remove NAs
    if (sum(is.na(names_list)) > 0){
      names_list2 = names_list2[-((1:length(names_list))[(is.na(names_list))])]
      names_list = names_list[-((1:length(names_list))[(is.na(names_list))])]
    }
    
    return(list('names_list' = names_list, 'names_list2' = names_list2))
  }
  
  #Given: vector of species name (species), as defined in speciesNames();
  #  data.frame of cufflinks data (table_names);
  #  vector of rows not found in getIdentifiers() translation table (not_rows);
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    was neither so couldn't translate to the lookup identifier - NULL (lookup);
  #  where lookup = FALSE a vector of the lookup_identifier in table row order (new_column);
  #  vector of the data format (format)
  #Returns: list of data.frame of the table with links in the "gene_id" and locus columns (table)
  #  or NULL if can't recognize species
  #Requires: createLink(), formatLookup()
  createTableLinks <- function(species,table_names,not_rows,lookup,new_column,format){
    #format = user_input$format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    extra_colname = format_details$extra_colname
    locus_position = format_details$locus_position
    
    species_data = speciesNames(species)
    table_names2 = table_names
    if (species == 'Worm'){ #some of this may be shared with other species/functionalised - particularly to join with bit in tableData
      #Wormbase links
      link_root_gene = species_data$link_root_gene
      link_root_locus = species_data$link_root_locus 
      link_root_locus_jb = species_data$link_root_locus_jb
      link_tail_locus_jb = species_data$link_tail_locus_jb
      #Convert to links
      if (lookup){
        #Gene ID links
        table_names2[gene_id_name] = createLink(paste(link_root_gene,table_names[gene_id_name][[1]],sep=""),table_names[gene_id_name][[1]])
      }else{
        table_names2[extra_colname] = createLink(paste(link_root_gene,new_column,sep=""),table_names[extra_colname][[1]])
      }
      #Locus links
      if (!is.null(link_root_locus)){
        if (!is.na(locus_position)){
          locus_name = correct_colnames[locus_position]
          
          #Gbrowse
          #table_names2[locus_name] = createLink(paste(link_root_locus,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]),sep=""),gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]))
          #Jbrowse
          link_locus_jb = paste(link_root_locus_jb,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]),link_tail_locus_jb,sep="")
          table_names2[locus_name] = createLink(link_locus_jb,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]))  
        }
      }
      #reset not_rows
      if (length(not_rows) > 0){
        table_names2[not_rows,] = table_names[not_rows,]
        if (!is.na(locus_position)){
          locus_name = correct_colnames[locus_position]
          table_names2[not_rows,locus_name] = levels(table_names[locus_name][[1]])[as.numeric(table_names2[not_rows,locus_name])]
        }
      }
      return(list(table = table_names2))
    }else{
      #Species-specific links
      link_root_gene = species_data$link_root_gene
      link_root_locus = species_data$link_root_locus 
      link_root_locus_jb = species_data$link_root_locus_jb
      link_tail_locus_jb = species_data$link_tail_locus_jb
      #Convert to links
      if (lookup){
        #Gene ID links
        table_names2[gene_id_name] = createLink(paste(link_root_gene,table_names[gene_id_name][[1]],sep=""),table_names[gene_id_name][[1]])
      }else{
        table_names2[extra_colname] = createLink(paste(link_root_gene,new_column,sep=""),table_names[extra_colname][[1]])
      }
      #Locus links
      if (!is.null(link_root_locus)){ #want to create link
        if (!is.na(locus_position)){
          locus_name = correct_colnames[locus_position]
          
          #Gbrowse
          #table_names2[locus_name] = createLink(paste(link_root_locus,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]),sep=""),gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]))
          #Jbrowse
          link_locus_jb = paste(link_root_locus_jb,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]),link_tail_locus_jb,sep="")
          table_names2[locus_name] = createLink(link_locus_jb,gsub(pattern='-',replacement = '..',x=table_names[locus_name][[1]]))  
        }
      }
      #reset not_rows
      if (length(not_rows) > 0){
        table_names2[not_rows,] = table_names[not_rows,]
        if (!is.na(locus_position)){
          locus_name = correct_colnames[locus_position]
          table_names2[not_rows,locus_name] = levels(table_names[locus_name][[1]])[as.numeric(table_names2[not_rows,locus_name])]
        }
      }
      return(list(table = table_names2))
      #return(NULL)
    }
  }
  
  #Given: data.frame of lookup_identifier gene identifiers and annotation(s) (annotation);
  #  data.frame of cufflinks data (table);
  #  [data.frame of cufflinks data without links (linkless_table)];
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    was neither so couldn't translate to the lookup identifier - NULL (lookup);
  #  vector of the data format (format);
  #  vector of species name (species), as defined in speciesNames()
  #Returns: data.frame of cufflinks data with new columns of annotation (table)
  #Requires: formatLookup()
  addAnnotationToRaw = function(annotation,table,linkless_table = NULL,lookup,format,species){
    #data = merge(table,annotation,by.x=1,by.y=1,all.x=T)
    if (is.null(lookup)){
    }else{
      #format = user_input$format
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      extra_position = format_details$extra_position
      locus_position = format_details$locus_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      extra_colname = format_details$extra_colname
      
      if (!lookup){
        if (is.null(linkless_table)){
          data = merge(table,annotation,by.x=extra_position,by.y=1,all.x=T)
          
          #fix weird bug where some things found in x but not y are ending up with multiple copies
          data = unique(data)
          
          if (extra_position > 1){ #not first
            before = 1:(extra_position - 1)
            before = before + 1
            if (extra_position < dim(data)[2]){ #not last
              after = (extra_position + 1):dim(data)[2]
              data = cbind(data[,before],data[,1],data[,after]) #needed when x!=1
            }else{ #only need to move first to last
              data = cbind(data[,before],data[,1]) #needed when x!=1
            }
          } #else don't need to rearrange
          #data = cbind(data[,2:4],data[,1],data[,5:(dim(data)[2])]) #needed when x=4
          colnames(data)[extra_position] = colnames(annotation)[1]
        }else{
          if (dim(table)[2] > dim(linkless_table)[2]){ 
            #add extra columns that from table
            extra_start = dim(linkless_table)[2] + 1
            temp_linkless = cbind(linkless_table,table[,extra_start:(dim(table)[2])])
            #fix column names
            colnames(temp_linkless)[extra_start:dim(temp_linkless)[2]] = colnames(table)[extra_start:(dim(table)[2])]
          }else{
            temp_linkless = linkless_table
          }
          
          data = merge(temp_linkless,annotation,by.x=extra_position,by.y=1,all.x=T)
          
          #fix weird bug where some things found in x but not y are ending up with multiple copies
          data = unique(data)
          
          if (extra_position > 1){ #not first
            before = 1:(extra_position - 1)
            before = before + 1
            if (extra_position < dim(data)[2]){ #not last
              after = (extra_position + 1):dim(data)[2]
              data = cbind(data[,before],data[,1],data[,after]) #needed when x!=1
            }else{ #only need to move first to last
              data = cbind(data[,before],data[,1]) #needed when x!=1
            }
          } #else don't need to rearrange
          #data = cbind(data[,2:4],data[,1],data[,5:(dim(data)[2])]) #needed when x=4
          colnames(data)[extra_position] = gene_id_name #colnames(annotation)[1]
          #restore links
          no_link_rows = NULL #data_cache$no_link_rows #because of reordering this is no longer accurate
          
          data2 = createTableLinks(species,data,no_link_rows,lookup,data[extra_position][,1],format)
          data = data2$table
          
          # data[extra_colname] = table[extra_colname]
          # if (!is.na(locus_position)){
          #   locus_name = correct_colnames[locus_position]
          #   data[locus_name] = table[locus_name]
          # }
        }
        return(data)
      }
    }
    if (is.null(linkless_table)){
      data = merge(table,annotation,by.x=gene_id_position,by.y=1,all.x=T)
      
      #fix weird bug where some things found in x but not y are ending up with multiple copies
      data = unique(data)
      if (gene_id_position > 1){ #not first
        before = 1:(gene_id_position - 1)
        before = before + 1
        if (gene_id_position < dim(data)[2]){ #not last
          after = (gene_id_position + 1):dim(data)[2]
          data = cbind(data[,before],data[,1],data[,after]) #needed when x!=1
        }else{ #only need to move first to last
          data = cbind(data[,before],data[,1]) #needed when x!=1
        }
      } #else don't need to rearrange
      #data = cbind(data[,2],data[,1],data[,3:(dim(data)[2])]) #needed when x=2
      colnames(data)[gene_id_position] = colnames(annotation)[1]
    }else{
      if (dim(table)[2] > dim(linkless_table)[2]){ 
        #add extra columns that from table
        extra_start = dim(linkless_table)[2] + 1
        temp_linkless = cbind(linkless_table,table[,extra_start:(dim(table)[2])])
        #fix column names
        colnames(temp_linkless)[extra_start:dim(temp_linkless)[2]] = colnames(table)[extra_start:(dim(table)[2])]
      }else{
        temp_linkless = linkless_table
      }
      
      data = merge(temp_linkless,annotation,by.x=gene_id_position,by.y=1,all.x=T)
      
      #fix weird bug where some things found in x but not y are ending up with multiple copies
      data = unique(data)
      if (gene_id_position > 1){ #not first
        before = 1:(gene_id_position - 1)
        before = before + 1
        if (gene_id_position < dim(data)[2]){ #not last
          after = (gene_id_position + 1):dim(data)[2]
          data = cbind(data[,before],data[,1],data[,after]) #needed when x!=1
        }else{ #only need to move first to last
          data = cbind(data[,before],data[,1]) #needed when x!=1
        }
      } #else don't need to rearrange
      #data = cbind(data[,2],data[,1],data[,3:(dim(data)[2])]) #needed when x=2
      colnames(data)[gene_id_position] = gene_id_name #colnames(annotation)[1]
      #restore links
      no_link_rows = NULL #data_cache$no_link_rows #because of reordering this is no longer accurate
      
      data2 = createTableLinks(species,data,no_link_rows,lookup,data[gene_id_position][,1],format)
      data = data2$table
      
      # data[gene_id_name] = table[gene_id_name]
      # if (!is.na(locus_position)){
      #   locus_name = correct_colnames[locus_position]
      #   data[locus_name] = table[locus_name]
      # }
    }
    return(data)
  }  

  #Given: vector of a comparison (comparison);
  #  data.frame of cufflinks data (table);
  #  vector of the data format (format)
  #Returns: NULL if comparison not in table or no comparison given,
  #  otherwise data.frame of cufflinks data only for that comparison
  #Requires: formatLookup()
  subsetByComparison = function(comparison,data,format){
    if (length(strsplit(as.character(comparison)," vs ")) > 0){ #Comparison string of correct format
      s1_name = strsplit(as.character(comparison)," vs ")[[1]][1] #Sample 1 name
      s2_name = strsplit(as.character(comparison)," vs ")[[1]][2] #Sample 2 name
      if(is.na(s2_name)){ #Invalid comparison
        return(NULL)
      }else{
        #format = user_input$format
        format_details = formatLookup(format)
        condition_1_position = format_details$condition_1_position
        condition_2_position = format_details$condition_2_position
        correct_colnames = format_details$correct_colnames
        if (!is.na(condition_1_position)){
          condition_1_name = correct_colnames[condition_1_position]
          condition_2_name = correct_colnames[condition_2_position]
          
          data = data[(data[condition_1_name] == s1_name)&(data[condition_2_name] == s2_name),]
        } #else want all data
        
        return(data)
      }
    }else{ #No comparison
      return(NULL)
    }
  }
  
  #Given: data.frame of cufflinks data (table);
  #  a vector of which columns to return, by default all (wanted_columns);
  #  list of the collected user inputs specifying how to filter the table (ui);
  #  vector of the data format (format)
  #Returns: data.frame of the filtered table
  #Requires: formatLookup()
  applyFilters = function(data,wanted_columns=NULL,ui,format){
    #format = user_input$format
    format_details = formatLookup(format)
    gene_id_position = format_details$gene_id_position
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    FC_position = format_details$FC_position
    
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    FC_name = correct_colnames[FC_position]
    
    adjusted_p_position = format_details$adjusted_p_position
    if (is.na(adjusted_p_position)){
      sig_name = "significant"
    }else{
      sig_position = format_details$sig_position
      sig_name = correct_colnames[sig_position]
    }
    
    if (is.null(wanted_columns)){
      if(ui$sig_checkbox){
        if ((ui$log2FC_up_checkbox)&(ui$log2FC_down_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)|(data[FC_name][[1]] <= -ui$log2FC_num)),]
        }else if ((ui$log2FC_up_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)),]
        }else if ((ui$log2FC_down_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] <= -ui$log2FC_num)),]
        }else{
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num),]
        }
      }else{
        if ((ui$log2FC_up_checkbox)&(ui$log2FC_down_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)|(data[FC_name][[1]] <= -ui$log2FC_num)),]
        }else if ((ui$log2FC_up_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)),]
        }else if ((ui$log2FC_down_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] <= -ui$log2FC_num)),]
        }else{
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num),]
        }
      }      
    }else{
      if(ui$sig_checkbox){
        if ((ui$log2FC_up_checkbox)&(ui$log2FC_down_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)|(data[FC_name][[1]] <= -ui$log2FC_num)),wanted_columns]
        }else if ((ui$log2FC_up_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)),wanted_columns]
        }else if ((ui$log2FC_down_checkbox)){
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] <= -ui$log2FC_num)),wanted_columns]
        }else{
          data3 = data[(data[sig_name][[1]] == "yes")&(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num),wanted_columns]
        }
      }else{
        if ((ui$log2FC_up_checkbox)&(ui$log2FC_down_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)|(data[FC_name][[1]] <= -ui$log2FC_num)),wanted_columns]
        }else if ((ui$log2FC_up_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] >= ui$log2FC_num)),wanted_columns]
        }else if ((ui$log2FC_down_checkbox)){
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num)&((data[FC_name][[1]] <= -ui$log2FC_num)),wanted_columns]
        }else{
          data3 = data[(data[value_1_name][[1]] >= ui$s1_fpkm_num)&(data[value_2_name][[1]] >= ui$s2_fpkm_num),wanted_columns]
        }
      }      
    }
    
    #altered splice patterns
    isoform_id_position = format_details$isoform_id_position
    if (!is.null(isoform_id_position)){
      if (ui$altsplice_checkbox){
        rows = duplicated(data3[gene_id_name][[1]]) + duplicated(data3[gene_id_name][[1]],fromLast=TRUE) #done in both directions as otherwise miss first/last variant of gene
        rows[rows > 1] = 1
        data3 = data3[as.logical(rows),]
      }
    }
    
    return(data3)
  }
  
  #Given: data.frame of cufflinks data (table);
  #  list of the collected user inputs specifying how to filter the table (ui);
  #  vector of the data format (format);
  #  boolean of whether to only filter by the main comparison (main_only)
  #Returns: data.frame of cufflinks data filtered by comparisons
  #Requires: subsetByComparison(); applyFilters(), formatLookup()
  combineComparisons = function(data,ui,format,main_only=FALSE){
    #format = user_input$format
    format_details = formatLookup(format)
    correct_colnames = format_details$correct_colnames
    gene_id_position = format_details$gene_id_position
    FC_position = format_details$FC_position
    gene_id_name = correct_colnames[gene_id_position]
    FC_name = correct_colnames[FC_position]
    
    main_comparison = ui$select_comp   
    other_comparisons = ui$select2InputFoundIn
    anti_other_comparisons = ui$select2InputNotFoundIn
    opposite = ui$invertFCcomp_checkbox
    
    if(main_only){
      return(subsetByComparison(main_comparison,data,format))
    }
    
    if ((length(other_comparisons) < 1)&(length(anti_other_comparisons) < 1)){ #No other comparisons, only use main comparison
      return(subsetByComparison(main_comparison,data,format))
    }else{
      data4 = subsetByComparison(main_comparison,data,format) #updating this is a pain and unnecessary
      current_genes = data4[gene_id_name][[1]]
      if ((opposite)){ 
        data4_up = data4[(data4[FC_name][[1]] > 0),] #subset by FC
        data4_down = data4[(data4[FC_name][[1]] < 0),]
      }
      if (length(other_comparisons) > 0){
        for (i in 1:length(other_comparisons)){
          new_data = subsetByComparison(other_comparisons[i],data,format)
          data3 = applyFilters(data = new_data,ui = ui,format = format)
          new_genes = data3[gene_id_name][[1]]
          
          #Check for genes that have opposite direction FC
          if ((opposite)){ 
            data3_up = data3[(data3[FC_name][[1]] > 0),] #subset by FC
            data3_down = data3[(data3[FC_name][[1]] < 0),]
            if (dim(data3_up)[1] > 0){
              new_genes_up = which(data3_up[gene_id_name][[1]] %in% data4_down[gene_id_name][[1]])
            }else{
              new_genes_up = NA
            }
            if (dim(data3_down)[1] > 0){
              new_genes_down = which(data3_down[gene_id_name][[1]] %in% data4_up[gene_id_name][[1]])
            }else{
              new_genes_down = NA
            }
            new_genes = (rbind(data3_up[new_genes_up,],data3_down[new_genes_down,]))[gene_id_name][[1]]
          }
          
          current_genes = intersect(current_genes,new_genes)  
        }
      }
      if (length(anti_other_comparisons) > 0){
        for (i in 1:length(anti_other_comparisons)){
          new_data = subsetByComparison(anti_other_comparisons[i],data,format)
          data3 = applyFilters(data = new_data,ui = ui,format = format)
          new_genes = data3[gene_id_name][[1]]
          
          #Check for genes that have opposite direction FC
          if ((opposite)){ 
            data3_up = data3[(data3[FC_name][[1]] > 0),] #subset by FC
            data3_down = data3[(data3[FC_name][[1]] < 0),]
            if (dim(data3_up)[1] > 0){
              new_genes_up = which(data3_up[gene_id_name][[1]] %in% data4_down[gene_id_name][[1]])
            }else{
              new_genes_up = NA
            }
            if (dim(data3_down)[1] > 0){
              new_genes_down = which(data3_down[gene_id_name][[1]] %in% data4_up[gene_id_name][[1]])
            }else{
              new_genes_down = NA
            }
            new_genes = (rbind(data3_up[new_genes_up,],data3_down[new_genes_down,]))[gene_id_name][[1]]
          }
          
          current_genes = setdiff(current_genes,new_genes)
        }
      }
      
      data2 = subsetByComparison(main_comparison,data,format)
      current_data = data2[which(data2[gene_id_name][[1]] %in% current_genes),]
      return(current_data)
    }
  }
  
  #Given: data.frame of cufflinks data (table);
  #  vector of sample 1 name (s1_name);
  #  vector of sample 2 name (s2_name);
  #  vector of the data format (format)
  #Returns: string warning if table empty, else matrix of summary table of means etc.
  #Requires: formatLookup()
  mySummaryTable = function(table,s1_name,s2_name,format){ 
    #format = user_input$format
    format_details = formatLookup(format)
    exp_unit = format_details$exp_unit
    gene_id_position = format_details$gene_id_position
    value_1_position = format_details$value_1_position
    value_2_position = format_details$value_2_position
    FC_position = format_details$FC_position
    
    correct_colnames = format_details$correct_colnames
    gene_id_name = correct_colnames[gene_id_position]
    value_1_name = correct_colnames[value_1_position]
    value_2_name = correct_colnames[value_2_position]
    FC_name = correct_colnames[FC_position]
    
    adjusted_p_position = format_details$adjusted_p_position
    if (is.na(adjusted_p_position)){
      sig_name = "significant"
    }else{
      sig_position = format_details$sig_position
      sig_name = correct_colnames[sig_position]
    }
    
    s_1_name = paste0("Sample 1 ",exp_unit)
    s_2_name = paste0("Sample 2 ",exp_unit)
    
    
    isoform_id_position = format_details$isoform_id_position
    if (is.null(isoform_id_position)){
      row_names = c("n genes",
                    s_1_name,
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    s_2_name,
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    "log2(fold change)",
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    "Significantly DE",
                    "yes",
                    "no")
      
      out = matrix(c(row_names,rep(NA,length(row_names))),ncol=2,nrow=length(row_names))
      if (dim(table)[1] > 0){
        out[1,2] = length(unique(table[gene_id_name][[1]],incomparables=c(NA,NULL,"-"))) #incomparables deals with cases where gene lacks primary ID
        out[2,1] = paste0("Sample 1 (",s1_name,") ",exp_unit)
        out[2,2] = ""
        out[3,2] = round(mean(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[4,2] = round(median(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[5,2] = round(max(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2) #there are warnings when this (+other max/min) have no non-missing arguments, should really handle them
        out[6,2] = round(min(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[7,1] = paste0("Sample 2 (",s2_name,") ",exp_unit)
        out[7,2] = ""
        out[8,2] = round(mean(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[9,2] = round(median(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[10,2] = round(max(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[11,2] = round(min(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[12,2] = ""
        out[13,2] = round(mean(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[14,2] = round(median(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[15,2] = round(max(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[16,2] = round(min(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[17,2] = ""
        out[18,2] = sum(table[sig_name] == "yes",na.rm=T)
        out[19,2] = sum(table[sig_name] == "no",na.rm=T)
        
        #Add spacers to improve display
        out = rbind(out[1:16,],c("",""),out[17:dim(out)[1],])
        out = rbind(out[1:11,],c("",""),out[12:dim(out)[1],])
        out = rbind(out[1:6,],c("",""),out[7:dim(out)[1],])
        out = rbind(out[1,],c("",""),out[2:dim(out)[1],])
        out = rbind(c("",""),out)
        out = rbind(out,c("",""))
      }else{
        out = "No results in table"
      }
      
    }else{
      isoform_id_name = correct_colnames[isoform_id_position]
      row_names = c("n genes",
                    "n_isoforms",
                    s_1_name,
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    s_2_name,
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    "log2(fold change)",
                    "mean value",
                    "median value",
                    "max value",
                    "min value",
                    "Significantly DE",
                    "yes",
                    "no")
      
      out = matrix(c(row_names,rep(NA,length(row_names))),ncol=2,nrow=length(row_names))
      if (dim(table)[1] > 0){
        out[1,2] = length(unique(table[gene_id_name][[1]],incomparables=c(NA,NULL,"-"))) #incomparables deals with cases where gene lacks primary ID
        out[2,2] = length((table[isoform_id_name][[1]])) 
        out[3,1] = paste0("Sample 1 (",s1_name,") ",exp_unit)
        out[3,2] = ""
        out[4,2] = round(mean(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[5,2] = round(median(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[6,2] = round(max(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2) #there are warnings when this (+other max/min) have no non-missing arguments, should really handle them
        out[7,2] = round(min(as.numeric(table[value_1_name][[1]]),na.rm=T),digits = 2)
        out[8,1] = paste0("Sample 2 (",s2_name,") ",exp_unit)
        out[8,2] = ""
        out[9,2] = round(mean(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[10,2] = round(median(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[11,2] = round(max(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[12,2] = round(min(as.numeric(table[value_2_name][[1]]),na.rm=T),digits = 2)
        out[13,2] = ""
        out[14,2] = round(mean(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[15,2] = round(median(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[16,2] = round(max(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[17,2] = round(min(as.numeric(table[FC_name][[1]]),na.rm=T),digits = 2)
        out[18,2] = ""
        out[19,2] = sum(table[sig_name] == "yes",na.rm=T)
        out[20,2] = sum(table[sig_name] == "no",na.rm=T)
        
        #Add spacers to improve display
        out = rbind(out[1:17,],c("",""),out[18:dim(out)[1],])
        out = rbind(out[1:12,],c("",""),out[13:dim(out)[1],])
        out = rbind(out[1:7,],c("",""),out[8:dim(out)[1],])
        out = rbind(out[1:2,],c("",""),out[3:dim(out)[1],])
        out = rbind(c("",""),out)
        out = rbind(out,c("",""))
      }else{
        out = "No results in table"
      }
    }
    return(out)
  }
  
  #Given: vector of a single gene;
  #  data.frame of cufflinks data (table);
  #  vector of the data format (format);
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    was neither so couldn't translate to the lookup identifier - NULL (lookup)
  #Returns: TRUE if gene in column "gene_id", FALSE, if gene in column "gene" and 
  #  NULL if in neither or if gene = NULL
  #Requires: formatLookup()
  findGene = function(gene,table,format,lookup){
    if(is.null(gene)){
      return(NULL)
    }else{
      #format = user_input$format
      format_details = formatLookup(format)
      gene_id_position = format_details$gene_id_position
      correct_colnames = format_details$correct_colnames
      gene_id_name = correct_colnames[gene_id_position]
      gene_symbol_position = format_details$gene_symbol_position
      if (!is.na(gene_symbol_position)){
        gene_symbol_name = correct_colnames[gene_symbol_position]
      }else{
        gene_symbol_name = "gene_symbol"
      }
      if (is.null(lookup)){
        if (!is.na(gene_symbol_position)){
          if (is.na(match(gene,table[gene_id_name][[1]]))){
            if (is.na(match(gene,table[gene_symbol_name][[1]]))){
              return(NULL)
            }else{
              return(FALSE)
            }
          }else{
            return(TRUE)
          }
        }else{
          if (is.na(match(gene,table[gene_id_name][[1]]))){
            # if (is.na(match(gene,table[gene_symbol_name][[1]]))){
            return(NULL)
            # }else{
            #   return(FALSE)
            # }
          }else{
            return(TRUE)
          }
        }
      }else{
        if (is.na(match(gene,table[gene_id_name][[1]]))){
          if (is.na(match(gene,table[gene_symbol_name][[1]]))){
            return(NULL)
          }else{
            return(FALSE)
          }
        }else{
          return(TRUE)
        }
      }
    }
  }
  
  #Given: vector of a single transcript;
  #  data.frame of cufflinks data (table);
  #  vector of the data format (format);
  #  boolean of whether the "gene_id" column in table was the lookup_identifier - TRUE 
  #    was the other of primary/secondary identifier as defined by InterMine - FALSE
  #    was neither so couldn't translate to the lookup identifier - NULL (lookup)
  #Returns: gene lookup_identifier if transcript in column isoform_id_position, else NULL
  #Requires: formatLookup()
  findTranscript = function(transcript,table,format,lookup){
    if(is.null(transcript)){
      return(NULL)
    }else{
      #format = user_input$format
      format_details = formatLookup(format)
      correct_colnames = format_details$correct_colnames
      
      
      isoform_id_position = format_details$isoform_id_position
      isoform_id_name = correct_colnames[isoform_id_position]
      
      if (is.na(match(transcript,table[isoform_id_name][[1]]))){
        return(NULL)
      }else{
        #return(TRUE)
        gene_id_position = format_details$gene_id_position
        gene_id_name = correct_colnames[gene_id_position]
        gene_symbol_position = format_details$gene_symbol_position
        if (!is.na(gene_symbol_position)){
          gene_symbol_name = correct_colnames[gene_symbol_position]
        }else{
          gene_symbol_name = "gene_symbol"
        }
        
        row = ((match(transcript,table[isoform_id_name][[1]])))
        
        if (is.null(lookup)){
          return(table[row,gene_id_name])
        }else{
          if(lookup){
            return(table[row,gene_id_name])
          }else{
            return(table[row,gene_symbol_name])
          }
        }
      }
    }
  }
  

  #DATA FORMATS MUST CONTAIN THESE 4 COLUMNS:
  # - gene identifier
  # - expression value in sample 1
  # - expression value in sample 2
  # - fold change (possible altered, e.g. log2) between the samples - if no p-value of this FC exists then all fold changes are assumed to be significant
  # Up to three columns will be added to them:
  # - a binary "Is DE significant?" column, will be generated if doesn't exist and adjusted p-value of DE does
  # - gene symbol, if it doesn't exist and can get translation table (local or online)
  # - second gene identifier, if can get translation table (local or online)
  #N.B. CURRENTLY ALLOWS FOR A MAXIMUM OF 14 COLUMNS IN THE BASE (without extra column) FORMAT, ANYMORE WILL CAUSE PROBLEMS
  # - extra column is not counted in this
  # - only 13 if a binary "Is DE significant?" column does not exist, or a gene symbol column doesn not exist (12 if neither exist)
  # To add more columns only need two changes - bit in collectUser() to allow the update of the selection of columns 
  # and outside of function the initialisation of the checkbox state (search for "Initialise column checkboxes")
  #N.B. gene_id_position MUST BE < extra_position
  #Given: vector of file format name (format)
  #Returns: list of values for that format used in app
  formatLookup <- function(format){
    if (is.null(format)){
      return(NULL)
    }else{
      if (format == 'cufflinks_gene'){
        n_columns = 14
        
        #Required columns positions
        gene_id_position = 2
        value_1_position = 8
        value_2_position = 9
        FC_position = 10
        
        #Filter + other column positions (NA if don't exist)
        adjusted_p_position = 13 #needed to generate a binary "Is DE significant?" column if it doesn't exist, if is.na assume all FCs are significant
        sig_position = 14 #binary "Is DE significant?" column, will be generated if doesn't exist and adjusted_p_position does
        gene_symbol_position = 3 #column that either does or should contain the gene symbol, will be generated if doesn't exist
        
        locus_position = 4
        condition_1_position = 5 #name of the first condition
        condition_2_position = 6 #name of the second condition
        
        #Column names
        correct_colnames =  c("test_id","gene_id","gene","locus","sample_1","sample_2","status",
                              "value_1","value_2","log2.fold_change.","test_stat","p_value","q_value",
                              "significant")
        long_colnames =  c("Test ID","Gene ID","Gene","Locus","Sample 1 name","Sample 2 name","Status",
                           "Sample 1 FPKM","Sample 2 FPKM","log2(fold change)","test stat","p value","adjusted p value",
                           "DE significant")
        
        #Web lookup column
        extra_position = 4 #position to add in extra column for multiple gene ids when using web lookup
        extra_colname = "gene_id2"
        extra_long_colname = "Gene ID 2"
        
        #Which columns to show when starting
        starting_checks = c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE) #includes the check for extra_column in extra_position place - will be ignored if column doesn't exist
        
        #Numeric columns
        numeric_columns = 8:13
        
        #Expression unit
        exp_unit = "FPKM"
        
        out = list(n_columns = n_columns, gene_id_position = gene_id_position, value_1_position = value_1_position, 
                   value_2_position = value_2_position, FC_position = FC_position,
                   adjusted_p_position = adjusted_p_position, sig_position = sig_position, 
                   gene_symbol_position = gene_symbol_position, locus_position = locus_position,
                   condition_1_position = condition_1_position, condition_2_position = condition_2_position,
                   correct_colnames = correct_colnames, long_colnames = long_colnames,
                   extra_position = extra_position, extra_colname = extra_colname, extra_long_colname = extra_long_colname,
                   starting_checks = starting_checks, numeric_columns = numeric_columns, exp_unit = exp_unit)
        return(out)
      }else if (format == 'cufflinks_transcript'){
        n_columns = 14
        
        #Required columns positions
        gene_id_position = 2 
        value_1_position = 8
        value_2_position = 9
        FC_position = 10
        
        #Filter + other column positions (NA if don't exist)
        adjusted_p_position = 13 #needed to generate a binary "Is DE significant?" column if it doesn't exist, if is.na assume all FCs are significant
        sig_position = 14 #binary "Is DE significant?" column, will be generated if doesn't exist and adjusted_p_position does
        gene_symbol_position = 3 #column that either does or should contain the gene symbol, will be generated if doesn't exist
        
        locus_position = 4
        condition_1_position = 5 #name of the first condition
        condition_2_position = 6 #name of the second condition
        
        isoform_id_position = 1 #needed to display multiple transcript info
        
        #Column names
        correct_colnames =  c("test_id","gene_id","gene","locus","sample_1","sample_2","status",
                              "value_1","value_2","log2.fold_change.","test_stat","p_value","q_value",
                              "significant")
        long_colnames =  c("Test ID","Gene ID","Gene","Locus","Sample 1 name","Sample 2 name","Status",
                           "Sample 1 FPKM","Sample 2 FPKM","log2(fold change)","test stat","p value","adjusted p value",
                           "DE significant")
        
        #Web lookup column
        extra_position = 4 #position to add in extra column for multiple gene ids when using web lookup
        extra_colname = "gene_id2"
        extra_long_colname = "Gene ID 2"
        
        #Which columns to show when starting
        starting_checks = c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE) #includes the check for extra_column in extra_position place - will be ignored if column doesn't exist
        
        #Numeric columns
        numeric_columns = 8:13
        
        #Expression unit
        exp_unit = "FPKM"
        
        out = list(n_columns = n_columns, gene_id_position = gene_id_position, value_1_position = value_1_position, 
                   value_2_position = value_2_position, FC_position = FC_position,
                   adjusted_p_position = adjusted_p_position, sig_position = sig_position, 
                   gene_symbol_position = gene_symbol_position, locus_position = locus_position,
                   condition_1_position = condition_1_position, condition_2_position = condition_2_position,
                   isoform_id_position = isoform_id_position,
                   correct_colnames = correct_colnames, long_colnames = long_colnames,
                   extra_position = extra_position, extra_colname = extra_colname, extra_long_colname = extra_long_colname,
                   starting_checks = starting_checks, numeric_columns = numeric_columns, exp_unit = exp_unit)
        return(out)
      }else if (format == 'pragui_gene'){
        n_columns = 14
        
        #Required columns positions
        gene_id_position = 2
        value_1_position = 8
        value_2_position = 9
        FC_position = 10
        
        #Filter + other column positions (NA if don't exist)
        adjusted_p_position = 13 #needed to generate a binary "Is DE significant?" column if it doesn't exist, if is.na assume all FCs are significant
        sig_position = 14 #binary "Is DE significant?" column, will be generated if doesn't exist and adjusted_p_position does
        gene_symbol_position = 3 #column that either does or should contain the gene symbol, will be generated if doesn't exist
        
        locus_position = 4
        condition_1_position = 5 #name of the first condition
        condition_2_position = 6 #name of the second condition
        
        #Column names
        correct_colnames =  c("test_id","gene_id","gene","locus","sample_1","sample_2","status",
                              "value_1","value_2","log2.fold_change.","test_stat","p_value","q_value",
                              "significant")
        long_colnames =  c("Test ID","Gene ID","Gene","Locus","Sample 1 name","Sample 2 name","Status",
                           "Sample 1 FPKM","Sample 2 FPKM","log2(fold change)","test stat","p value","adjusted p value",
                           "DE significant")
        
        #Web lookup column
        extra_position = 4 #position to add in extra column for multiple gene ids when using web lookup
        extra_colname = "gene_id2"
        extra_long_colname = "Gene ID 2"
        
        #Which columns to show when starting
        starting_checks = c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE) #includes the check for extra_column in extra_position place - will be ignored if column doesn't exist
        
        #Numeric columns
        numeric_columns = 8:13
        
        #Expression unit
        exp_unit = "norm read count"
        
        out = list(n_columns = n_columns, gene_id_position = gene_id_position, value_1_position = value_1_position, 
                   value_2_position = value_2_position, FC_position = FC_position,
                   adjusted_p_position = adjusted_p_position, sig_position = sig_position, 
                   gene_symbol_position = gene_symbol_position, locus_position = locus_position,
                   condition_1_position = condition_1_position, condition_2_position = condition_2_position,
                   correct_colnames = correct_colnames, long_colnames = long_colnames,
                   extra_position = extra_position, extra_colname = extra_colname, extra_long_colname = extra_long_colname,
                   starting_checks = starting_checks, numeric_columns = numeric_columns, exp_unit = exp_unit)
        return(out)
      }else if (format == 'pragui_transcript'){
        n_columns = 14
        
        #Required columns positions
        gene_id_position = 2 
        value_1_position = 8
        value_2_position = 9
        FC_position = 10
        
        #Filter + other column positions (NA if don't exist)
        adjusted_p_position = 13 #needed to generate a binary "Is DE significant?" column if it doesn't exist, if is.na assume all FCs are significant
        sig_position = 14 #binary "Is DE significant?" column, will be generated if doesn't exist and adjusted_p_position does
        gene_symbol_position = 3 #column that either does or should contain the gene symbol, will be generated if doesn't exist
        
        locus_position = 4
        condition_1_position = 5 #name of the first condition
        condition_2_position = 6 #name of the second condition
        
        isoform_id_position = 1 #needed to display multiple transcript info
        
        #Column names
        correct_colnames =  c("test_id","gene_id","gene","locus","sample_1","sample_2","status",
                              "value_1","value_2","log2.fold_change.","test_stat","p_value","q_value",
                              "significant")
        long_colnames =  c("Test ID","Gene ID","Gene","Locus","Sample 1 name","Sample 2 name","Status",
                           "Sample 1 FPKM","Sample 2 FPKM","log2(fold change)","test stat","p value","adjusted p value",
                           "DE significant")
        
        #Web lookup column
        extra_position = 4 #position to add in extra column for multiple gene ids when using web lookup
        extra_colname = "gene_id2"
        extra_long_colname = "Gene ID 2"
        
        #Which columns to show when starting
        starting_checks = c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE) #includes the check for extra_column in extra_position place - will be ignored if column doesn't exist
        
        #Numeric columns
        numeric_columns = 8:13
        
        #Expression unit
        exp_unit = "norm read count"
        
        out = list(n_columns = n_columns, gene_id_position = gene_id_position, value_1_position = value_1_position, 
                   value_2_position = value_2_position, FC_position = FC_position,
                   adjusted_p_position = adjusted_p_position, sig_position = sig_position, 
                   gene_symbol_position = gene_symbol_position, locus_position = locus_position,
                   condition_1_position = condition_1_position, condition_2_position = condition_2_position,
                   isoform_id_position = isoform_id_position,
                   correct_colnames = correct_colnames, long_colnames = long_colnames,
                   extra_position = extra_position, extra_colname = extra_colname, extra_long_colname = extra_long_colname,
                   starting_checks = starting_checks, numeric_columns = numeric_columns, exp_unit = exp_unit)
        return(out)
      }
      return(NULL)
    }
  }
  
  
  
  #Given:
  #Returns:
  #Requires:
  
  #Given:
  #Returns:
  #Returns reactive variables:
  #Requires reactive functions:
  #Requires reactive variables:
  #Requires functions: 
  
  
  ##Depreciated
}
