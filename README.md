# PEAT
PRAGUI Exploratory Analysis Tool

INSTALLATION


RUNNING


MANUAL



TERMINOLOGY

Samples - the different samples that were sequenced.<br>
Comparisons - any pairwise comparison of two samples, where sample 1 refers to the first sample in the comparison and sample 2 to the second, e.g. WT vs KO - WT is sample 1 and KO is sample 2. <br>
Sidebar - the options down the left hand side. <br>
Tabs - the options along the top. <br>
DE - differentially expressed. <br>
KO - Knockout.  <br>
OE - Overexpression.

Details of what most options do can be found in the tooltips shown when the mouse is left over them for a short while.


LOADING A FILE

PEAT starts with the Select Input sidebar option already open. To load a file select the correct species and file format. If you select an incompatible file format then PEAT will not do anything with the data. If you select the wrong species you will get an error message saying "Could not recognize the gene identifiers, the program will act as though the selected species was Other, i.e. online data will not be used. This may happen if you have chosen the wrong species." If your species is not available or you are offline you can select 'Other' for species and PEAT will function normally, but without any online annotation (specifically Enrichment analysis will be non-functional).


LOOKING AT DATA

PEAT has three ways of looking at the raw data, these comprise the first three tabs: Summary, Table, Per gene. The Summary tab starts open.

Summary - The Summary tab shows a summary of the currently chosen comparison as filtered by the currently chosen filters (see FILTERING).

Table - This shows the actual data summarised on the Summary tab. This data can be sorted by column and searched, on all columns, using the search box in the top right, or on a specific column using the boxes at the bottom. The latter can only search for exact matches, while the former allows for partial matching and for matching uppercase letters to lowercase and vice versa. By default not all columns are shown (exclusions include p-value for DE), column choice can be controlled from the Table columns sidebar option. Clicking on the blue boxes around gene IDs (when online) will take you to the species-specific webpage detailing that gene. Clicking on any row will take you to the Per gene tab for that gene.

Per gene - Unlike Summary and Table, which show a per comparison view of the data, the Per gene tab shows data for all samples for a given gene. Genes can be selected by linking from other tabs (e.g. Table) or from the Choose gene dropbox. This allows partial matching and attempts to autocomplete. Once a gene is selected the tab shows first details of the possible IDs for that gene, with a brief description if available, second a barplot of the expression of the gene in all samples, third options controlling the barplot and, optionally, fourth species-specific data on that gene (currently only in C. elegans). The barplot shows asterixes for all values significantly differently expressed with respect to the black bar. The identity of the black bar can be altered from the plot options. 


FILTERING

The Filter Input sidebar option allows the selection of a subset of data. Most options are self-explanatory or easily understood from their tooltips. (Not) In other comparisons allows simultaneous filtering of data in multiple comparisons, the rest of the filters will be applied to each comparison selected and either the intersect (for In other comparisons) or the set difference (for Not in other comparisons) of the resulting gene lists returned. This allows filtering to find genes DE in a range of comparisons, e.g. in KOs of a number of genes in a pathway in pairwise comparisons versus WT, to eliminate noise and/or DE that is not due to a shared function/pathway for these genes; or to find the effects of KO of a specific gene that are not mediated by another gene using the Not in other comparisons filter. The 'Fold-change in opposite direction in comparison' checkbox exists for case where you have positive and negative regulators or KO and OE mutants, such that (some) genes DE in one comparison, e.g. WT vs KO, should be DE in the other direction in the other comparison, WT vs OE. The final filter option allows you to only show genes in a particular list (that also pass the other filters - to see all genes in a particular list use the 'Gene lists genes' tab). Gene lists must be made using the 'Gene lists' sidebar option, see USING GENE LISTS.


ADDING ANNOTATIONS TO TABLE

Additional columns of annotation can be added to the 'Table' tab using the 'Add annotations' sidebar option. Most species have a few built-in options that add web-derived annotations, typically brief description (which is useful when you have no idea what a gene is) and GO terms (importantly this is not precisely the same GO annotation used in the 'Enrichment analysis' tab as they come form different online sources, though they should be very similar). Additionally you can add annotation of your own from any tab-separated file, where the first column consists of the gene ID. Buttons will appear in the sidebar for each file you load allowing its removal.


USING GENE LISTS

Gene lists are a list of gene IDs with associated metadata. They can be created using the 'Gene lists' sidebar option, their metadata viewed/modified using the 'Gene lists' tab and their contents and the associated expression values for all samples viewed in the 'Gene lists genes' tab. WARNING: any creation or deletion of gene lists resets the gene list filtering to 'None', as the dropbox must be regenerated with the possible options.

'Gene lists' sidebar option - 'Create gene list' makes a new gene list with the above name and metadata derived from the current input settings and filter options. 'Load gene list' can either load a saved gene list with existing metadata or just a list of gene IDs, with one per line, in which case the metadata will be set to Null.

'Compare lists' sidebar options create a new list consisting of the union, intersect or set difference of the selected lists.

'Gene lists' tab - This shows the metadata for all gene lists. Columns can be turned on/off using the 'Column visibility' option, columns can be sorted and searched as on the 'Table' tab. Most buttons are again self-explanatory or easily understood from their tooltips. 
The 'Perform enrichment test' button allows you to look if a particular gene list is enriched for genes from another gene list. Suppose you have a list of all genes that have a certain motif (list 1) and you want to see if your currently filtered data has more genes with that motif that expected (which you have created list 2 from), you would click on list 1, then list 2, then the 'Perform enrichment test' button. This will open a popup with four numbers, the first (Population size) is the total number of 'genes' in your data (this depends on how your analysis is done, but is typically all protein coding genes + all non-coding genes + a few random orphaned sequences, so may not actually be the number you want - e.g. if list 1 is actually all protein-coding genes with the motif, the number should be the count of all protein-coding genes), the second is the number of genes in list 1, the third is the number of genes in list 2 and the fourth is the number of genes in the intersect of list 1 and list 2. Clicking the 'Calculate enrichment' button will perform a hypergeometric test and return the p-value for whether you should reject the null hypothesis that list 2 is not enriched for genes in list 1. WARNING - the authors take no responsibility for the misuse of statistics, use this function at your own peril.
The 'Draw Venn diagram' button will produce a Venn diagram for the selected lists, n.b. areas are not proportional to the number of elements in them.


MULTI-GENE PLOTS

The 'Figures' tab contains a number of plots of the filtered data. When you click on 'Show X' the plot will be updated each time the filters are changed, so it is a good idea to only show one plot at a time. Apart from Heatmap the plots use per comparison data only.

Volcano plot - this plots significance versus log2 fold-change on the y and x axes, respectively. By default only filtered data is shown, but the options allow for the display of all data in grey.

Scatter plot - this plots the expression levels for each sample in the current comparison against each other.

Pairwise comparison - this plots the log2 fold-change of your current comparison (x-axis) against the log2 fold-change of another comparison (chosen here, y-axis) for each gene in which both comparisons have been made and calculates a R value for the correlation using Spearman correlation coefficient.

Clustering - this performs hierarchical clustering on the fold-changes of the currently filtered data. Clustered genes may share pathways, functions or regulation.

Heatmap - this shows per sample expression for up to 50 of the genes in the currently filtered data (an option allows display for all genes). If there are more than 50 genes in the currently filtered data the 50 most DE are shown. Options allow use of log10 expression and to turn off clustering e.g. to preserve the gene order in a list. Following the plot is a table of the data used to draw the plot.


ENRICHMENT ANALYSIS

The 'Enrichment analysis' tab provides a least one method for looking at the enrichment of certain annotations in the currently filtered data. Like on the 'Figures' tab the enrichment will be updated each time the filters are changed, so it is a good idea to only show one enrichment at a time. These calculations are slow. gProfile enrichment allows choice of what type of annotation to look at (bear in mind not all types are relavant for all species) and then shows all annotations for the currently filtered data, p-values for their enrichment and which genes have that annotation. Buttons allow the creation of gene lists that include all genes with an annotation ('Create gene list') and genes with selected annotations only ('Create gene list from selected'; selection done by clicking). Similar to the 'Table' tab clicking on the blue boxes around annotation IDs (when online & available) will take you to the annotation-specific webpage detailing that annotation. Below this table is a plot of which genes have each of the most highly enriched annotations. 

Other enrichment methods have similar functionality.


ORTHOLOGS

Choosing a species will show orthologs of the currently filtered data in that species. The gene list buttons have the same functionality as those on the 'Enrichment analysis' tab (see ENRICHMENT ANALYSIS).


DOWNLOADING DATA

Data for all tables (apart from the 'Gene lists' tab metadata table) can be downloaded from the 'Download data' sidebar option. Please note that any filtering done by searching on a given table or be reordering of columns etc. will not be reflected in the download, though any annotations added to the 'Table' tab using the 'Add annotations' sidebar option will be. Figures/plots can be downloaded by right-clicking and choosing Save as (or equivalent) or in their SVG format version from this tab. Links will appear after data has been generated. WARNING tables/plots are only updated when on the relavant tab, if you change filters without looking at tables/plots before saving them your changes will not be reflected.


OTHER

To convert between different forms of gene ID PEAT accesses web databases and caches the data locally for speed, only updating if the most recent version is older than 30 days. The Sidebar also contains the date from which this was last done in the form "Gene ID conversion table from".


BUGS/ERRORS

If you find any please report them to acrisp@mrc-lmb.cam.ac.uk including as many details as possible about what you did leading up the error and the actual error message. n.b. this does not apply to error popups, which are there deliberately.




