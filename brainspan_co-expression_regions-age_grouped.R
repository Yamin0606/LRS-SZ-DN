setwd("~/Desktop/lrs-de_novo/analysis_2022")
library(dplyr)
col_meta <- read.csv("./brainspan_genes_matrix_csv/columns_metadata.csv",header = T)
row_meta <- read.csv("./brainspan_genes_matrix_csv/rows_metadata.csv",header = T)
#group into regions

annot <- xlsx::read.xlsx2("./brainspan_genes_matrix_csv/brainspan_样本例数-变量名注释.xlsx",sheetName = "brain")
a <- unique(annot$new_group)
annot[which(annot$new_group == a[1]),]
col_meta <- merge(col_meta,annot[,c("structure_id","new_group")],by="structure_id",all.x=T)

######
t1 <- strsplit(col_meta$age, " ")
age_num <- t(as.data.frame(t1))[,1]
group_ym1 <- t(as.data.frame(t1))[,2]
col_meta <- cbind(col_meta,cbind(age_num,group_ym1))
for (i in 1: nrow(col_meta)){
  if (col_meta$group_ym1[i]=="pcw"){
    col_meta$group_ym3[i] <- "pwd0-40"
  }else if (col_meta$group_ym1[i]=="mos" ){
    col_meta$group_ym3[i] <- "yrs1-12"
  }else if (col_meta$group_ym1[i]=="yrs" &&col_meta$age_num[i] <= 12){
    col_meta$group_ym3[i] <- "yrs1-12"
  }else if (col_meta$group_ym1[i]=="yrs" &&col_meta$age_num[i] <= 40){
    col_meta$group_ym3[i] <- "yrs13-40"
  }
}
table(col_meta$group_ym3,col_meta$new_group)
#################################
exp <- read.csv("./brainspan_genes_matrix_csv/expression_matrix.csv", header = FALSE)
exp <- exp[,2:ncol(exp)]
rownames(exp) <- row_meta$ensembl_gene_id
colnames(exp) <- col_meta$column_num

######select genes
# ours <- read.table("./results/gene_list_0413.txt",header = T)
# aa <- merge(ours,row_meta,by.x="gene",by.y="gene_symbol")
# nrow(ours) - nrow(aa) #47

ours <- read.table("./results/gene_list_0413_ensembl.txt",header = F)
names(ours) <- "ensembl_gene_id"
aa <- merge(ours,row_meta,by = "ensembl_gene_id")
nrow(ours) - nrow(aa) #7

#########################################
list_ours <- aa[,"ensembl_gene_id"]
exp_ours <- exp[c(list_ours),]
tt <- row_meta[!is.na(row_meta$entrez_id),"ensembl_gene_id"] #otherwise there are 50k+ genes
list_left <- setdiff(tt, list_ours)
exp_left <- exp[which(row.names(exp)  %in% list_left),]



str <- unique(col_meta$new_group)
period <- unique(col_meta$group_ym3)

dataf_ours <- list()  
dataf <- list()
for (x in str) {
  for (y in period) {
    target <- col_meta[which(col_meta$group_ym3==y & col_meta$new_group==x),"column_num"]
    dataf_ours[[paste0(y,"_",x)]] <- t(exp_ours[,target])
    dataf[[paste0(y,"_",x)]] <- t(exp_left[,target])
  }
}

#not all sample have tissue of all brain regions available.
connect_num <- function(x) {
  if (nrow(x) == 0){
    return(NA)
  }else {
    a0 <- x[,colSums(x)!=0]
    a1 <- cor(a0)
    a2 <- abs(a1) >= 0.8
    a3 <- (sum(a2,na.rm = T)-nrow(a2))/2
    return(a3)
  }
}

connect_list_ours <- list()
for (i in names(dataf)) {
  data_t <- as.data.frame(dataf_ours[[i]])
  aa <- connect_num(data_t)
  connect_list_ours[[i]] <- c(aa)
}

connect_list <- list()
n=1000
for (i in names(dataf)) {
  data_t <- as.data.frame(dataf[[i]])
  aa <- replicate(n,connect_num(data_t[,sample(row.names(exp_left),length(list_ours),FALSE)]))
  connect_list[[i]] <- c(aa)
}
write.table(connect_list,"results_brainspanco_expression_permu1000_0413.txt")

table <- list()
for (i in names(connect_list)) {
  aa <- connect_list[[i]]
  a3 <- connect_list_ours[[i]]
  if (!is.na(a3)){
    limit1 <- max(max(aa),a3)
    limit2 <- min(min(aa),a3)
    png(paste0("./results/plots/brainspan_co-exp/",i,".png"))
    plot(density(aa),main=i,xlab="# of connections",
         ylab="% of simulation", lwd=2,xlim=c(limit2,limit1))
    abline(v=a3,col="red",lwd=3)
    abline(v=quantile(aa,0.95),lty = 3,col="blue",lwd=3)
    p <- table((aa >= a3))["TRUE"]/n;p
    p1 <- ifelse(is.na(p),"< 1E-03",p)
    mid <- mean(c(range(density(aa)$y)[1],range(density(aa)$y)[2]))
    text(a3-20,mid,paste0("p: ",p1),pos=2,srt=90)
    dev.off()
    table[[i]] <- c(a3,p)
  }
} 



