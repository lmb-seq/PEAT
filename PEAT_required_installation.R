#Get installed packages
packages = rownames(installed.packages())

#Required Packages
required = c("shiny","shinyjs","shinydashboard","DT","data.table","plotly","gplots","svglite","ComplexHeatmap","httr","jsonlite","RCurl","InterMineR","gProfileR","eulerr")
versions = c("1.1.0","1.0","0.7.1","0.4","1.11.8","4.8.0","3.0.1","1.2.1","1.17.1","1.3.1","1.5","1.95-4.11","1.0.0","0.6.7","5.0.0")
sources = c("CRAN","CRAN","CRAN","CRAN","CRAN","CRAN","CRAN","CRAN","Bioconductor","CRAN","CRAN","CRAN","Bioconductor","CRAN","CRAN")

#Test if installed and if not do so
for (i in 1:length(required)){
  package = required[i]
  version = versions[i]
  source = sources[i]
  if((!package %in% packages)){
    cat("WARNING: Updated package version is needed... Installing updated package version...")
    if (source == "CRAN"){
      install.packages(package,repos='http://www.stats.bris.ac.uk/R/') #http://cran.us.r-project.org
    }else if (source == "Bioconductor"){
    #  source("https://bioconductor.org/biocLite.R")
    #  biocLite(package)
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install(package)
    }
  }else if(packageDescription(package)$Version==version){
    cat("WARNING: Updated package version is needed... Installing updated package version...")
    if (source == "CRAN"){
      install.packages(package,repos='http://www.stats.bris.ac.uk/R/') #http://cran.us.r-project.org
    }else if (source == "Bioconductor"){
      #  source("https://bioconductor.org/biocLite.R")
      #  biocLite(package)
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install(package)
    }
  }
}

