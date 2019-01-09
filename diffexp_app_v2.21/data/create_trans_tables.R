library(AnnotationDbi)

##Human
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

#genes = read.csv("/home/acrisp/projects/kozik_crispr/new_version/libraries/homo_sapiens/genome/sabatini-lander/sabatini-lander-1000000095-library_revised2.csv")
#gene_symbol_list = as.vector(unique(genes$gene_symbol))

A = read.csv("~/projects/shiny/test_data_sets/ex_output_4_peat.txt",sep="\t")
gene_list = as.character(A[,2])
gene_list2 = as.character(A[,3])

geneinfo = as.data.frame(
  select(org.Hs.eg.db,
         keys = gene_list,
         columns=c("GENENAME","SYMBOL", "ENSEMBL","ENTREZID"),
         keytype = "ENSEMBL")
)


geneinfo2 = as.data.frame(
  select(org.Hs.eg.db,
         keys = gene_list2,
         columns=c("GENENAME","SYMBOL", "ENSEMBL","ENTREZID"),
         keytype = "SYMBOL")
)

#unique
gi1 = geneinfo[match(unique(geneinfo[,1]),geneinfo[,1]),]
gi2 = geneinfo2[match(unique(geneinfo2[,1]),geneinfo2[,1]),]

trans_table = gi1[,c(4,1,3)]
colnames(trans_table) = c("primary","secondary","symbol")

#geneinfo2 not in gi1
gi2 = (gi2[(is.na(match(gi2[,3],gi1[,1]))),])
gi2 = gi2[,c(4,3,1)]

#get rid of gi2 symbols,entrezID already in gi1
gi2 = (gi2[(is.na(match(gi2[,3],gi1[,3]))),])
gi2 = (gi2[(is.na(match(gi2[,1],gi1[,4]))),])

#join
colnames(gi2) = c("primary","secondary","symbol")
trans_table = rbind(trans_table,gi2)


today = Sys.Date()
filename = paste0("raw_HTT_",today,".tsv")   
write.table(trans_table,paste0("data/Human/",filename),sep="\t")

#add in old 'retired' genes
old_trans = read.csv("data/Human/Human_trans_table_2017-07-20.tsv",sep="\t")

#retired ENSG
old_trans = (old_trans[(is.na(match(old_trans[,2],trans_table[,2]))),])

trans_table2 = rbind(trans_table,old_trans)

filename = paste0("Human_trans_table_",today,".tsv")   
write.table(trans_table2,paste0("data/Human/",filename),sep="\t")

