library(AnnotationDbi)

##Human
library(org.Hs.eg.db)

A = read.csv("../Human/Human_trans_table_2018-05-01.tsv",sep="\t")
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
write.table(trans_table,paste0(filename),sep="\t")

#add in old 'retired' genes
old_trans = A

#retired ENSG
old_trans = (old_trans[(is.na(match(old_trans[,2],trans_table[,2]))),])

trans_table2 = rbind(trans_table,old_trans)


#add in intermine genes
old_trans = read.csv("../Human/intermine_Human_trans_table_2020-08-19.tsv",sep="\t")

#retired ENSG
old_trans = (old_trans[(is.na(match(old_trans[,2],trans_table2[,2]))),])

trans_table3 = rbind(trans_table2,old_trans)


filename = paste0("Human_trans_table_",today,".tsv")   
write.table(trans_table3,paste0(filename),sep="\t")


##Mouse
library(org.Mm.eg.db)
ensgns <- keys(org.Mm.eg.db, "ENSEMBL")
geneinfo = as.data.frame(
  select(org.Mm.eg.db,
         keys = ensgns,
         columns=c("GENENAME","SYMBOL", "ENSEMBL","ENTREZID"),
         keytype = "ENSEMBL")
)

trans_table = cbind()



old_table = read.csv("../Mouse/Mouse_trans_table_2019-06-27.tsv",sep="\t")
data = as.matrix(old_table)
data[,1] = as.character(old_table[,1])
old_trans_table = data

#combine
trans_table = old_trans_table
matches = match((old_trans_table[,1]),geneinfo[,1])
names(matches) = 1:nrow(trans_table)
matches = matches[!is.na(matches)]
trans_table[as.numeric(names(matches)),3] = geneinfo[matches,3]

geneinfo2 = geneinfo[-matches,]
geneinfo2 = geneinfo2[,c(1,4,3)]
colnames(geneinfo2) = colnames(trans_table)

trans_table = rbind(trans_table,geneinfo2)


##add other ENSGM
#
#A = read.csv("/data1/rnaseq_analyses/meng/meng_WT_Adh5_gene_exp.diff",sep="\t")
#A[,2] = as.character(A[,2])
#A[,2] = sub("gene:","",A[,2])
#A[,1] = as.character(A[,1])
#A[,1] = sub("gene:","",A[,1])
#A[,2] = as.character(A[,2])
#A[,2] = sub("transcript:","",A[,2])
#A[,1] = as.character(A[,1])
#A[,1] = sub("transcript:","",A[,1])
#
#A2 = A[(is.na(match(A[,2],trans_table[,1]))),]
#A2 = A2[,1:3]
#colnames(A2) = colnames(trans_table)
#A2 = A2[-1,]
#trans_table = rbind(trans_table,A2)


#today = Sys.Date()
#filename = paste0("raw_MTT_",today,".tsv")   
#write.table(trans_table,paste0(filename),sep="\t")

filename = paste0("Mouse_trans_table_",today,".tsv")   
write.table(trans_table,paste0(filename),sep="\t")


