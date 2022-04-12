
all_annot_flt <- read.table("~/lrs_do_novo/res_annotSV_flt3/denovo_flt3_merged.annotated0408_fltzym.tsv",sep="\t")
all_annot_flt$group <- "case"

all_annot_match <- {}
annot_all <- read.table("~/lrs_do_novo/res_annotSV/CN_sniffles.annotated.tsv",sep="\t",na.strings="",header=T,quote="")
list_var <- intersect(names(all_annot_flt),names(annot_all))
for (m in 1:1000){
  match <- read.table(paste0("~/lrs_do_novo/data/vcf/nano-denovo/population_random_match/random_match_",m),header=T)
  names(match)[1:3] = c("SV_chrom", "SV_start","SV_end")
  t2 <- merge(annot_all,match,by=c("SV_chrom", "SV_start","SV_end"))
  t2$group = "match"
  data_case_match <- rbind(all_annot_flt[,c(list_var,"group")],t2[,c(list_var,"group")])
  write.table(data_case_match,paste0("~/lrs_do_novo/data/data4test_merge/data4test_",m),quote = F,row.names = F,sep="\t")
}


quantiVar <- read.table("~/lrs_do_novo/res_annotSV/list_quantitative_variables")
qualiVar <- read.table("~/lrs_do_novo/res_annotSV/list_qualitative_variables")


res_qua <- {}
for (m in 1:1000){
  data <- read.table(paste0("~/lrs_do_novo/data/data4test_merge/data4test_",m),header = T,sep="\t",na.strings=c("",NA))
  data1 <- data[,c("group",qualiVar$V1)]
  for (var in qualiVar[c(2,4:9),1]) {
    data1[,paste0(var,"_b")] <- ifelse(is.na(data1[,var]),0,1)
  }
  data1$P_all_phen <- apply(data1[,c("P_gain_phen_b","P_loss_phen_b","P_ins_phen_b","P_snvindel_phen_b")],1,sum)
  data1$P_all_phen_b <- ifelse(data1$P_all_phen>0,1,0)
  for (var1 in c("Gene_name_b","RE_gene_b","OMIM_ID_b","P_all_phen_b")){
    t <- table(data1[,var1],data1$group)
    t1 <- chisq.test(t)
    res_qua[[paste0(var1,"_",m)]] <- c(as.vector(t),t1$p.value)
  }
}

res_qua2 <- data.frame()
for (var in c("Gene_name_b","RE_gene_b","OMIM_ID_b","P_all_phen_b")) {
  results <- data.frame()
  for (m in 1:1000){
    results <- rbind(results, res_qua[[paste0(var,"_",m)]])
  }
  names(results) <- c("case_0","case_1","match_0","match_1","p")
  attach(results)
  results$OR <- case_1*match_0/match_1/case_0
  detach(results)
  res_qua2 <- rbind(res_qua2,quantile(results$OR[results$OR!="Inf"],c(0.5,0.05,0.95)))
  write.table(results,paste0("~/lrs_do_novo/data/data4test_merge/","res_comparison_",var),quote = F,row.names = F)
}
res_qua2

res_quan <- data.frame()
for (var in quantiVar[13:15,1]) {
  res_mean <- data.frame()
  res_median <- data.frame()
  for (m in 1:1000){
    data <- read.table(paste0("~/lrs_do_novo/data/data4test_merge/data4test_",m),header = T,sep="\t")
    mean1 <- tapply(data[,var],data[,"group"],function(x) mean(na.omit(x)))
    median1 <- tapply(data[,var],data[,"group"],function(x) median(na.omit(x)))
    res_mean <- rbind(res_mean,mean1)
    res_median <- rbind(res_median,median1)
  }
  names(res_mean) <-c("case","match")
  names(res_median) <-c("case","match")
  write.table(res_mean,paste0("~/lrs_do_novo/data/data4test_merge/","res_comparison_mean_",var),quote = F,row.names = F)
  write.table(res_median,paste0("~/lrs_do_novo/data/data4test_merge/","res_comparison_median_",var),quote = F,row.names = F)
  res1 <- c(res_mean[1,"case"],quantile(res_mean$match,c(0.5,0.05,0.95)))
  res2 <- c(res_median[1,"case"],quantile(res_median$match,c(0.5,0.05,0.95)))
  res_quan <- rbind(res_quan,rbind(res1,res2))
}
  

