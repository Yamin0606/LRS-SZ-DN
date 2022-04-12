list <- read.table("~/lrs_do_novo/data/list_nano_final3.txt")
all_annot_flt <- data.frame()
for (id in list$V1){
  annot_flt <- read.table(paste0("~/lrs_do_novo/res_annotSV_flt3/",id,"_flt3.annotated.tsv"),sep="\t",na.strings="",header=T,quote="")
  keep <- read.table("~/lrs_do_novo/res_annotSV_flt3/sv_keep_be_hc_rbind_individual.txt",header=T)
  t1 <- annot_flt[annot_flt$AnnotSV_ID %in% keep[which(keep$sample_ID==id),2],]
  t1$group="case"
  if (id==list$V1[1]){
    col_name <- names(t1)
  }
  all_annot_flt <-  rbind(all_annot_flt, setNames(t1,col_name))
}


all_annot_match <- {}
for (id in list$V1){
  annot_all <- read.table(paste0("~/lrs_do_novo/res_annotSV/",id,"-0.annotated.tsv"),sep="\t",na.strings="",header=T,quote="")
  for (m in 1:1000){
    match <- read.table(paste0("~/lrs_do_novo/data/vcf/nano-denovo/inherited/random_match/",id,"_match_",m),header=T)
    names(match) = c("SV_chrom", "SV_start","SV_end")
    t2 <- merge(annot_all,match,by=c("SV_chrom", "SV_start","SV_end"))
    t2$group = "match"
    if (id==list$V1[1]){
      col_name <- names(t2)
    }
    names(t2) <- col_name
    all_annot_match[[paste0(id,"_",m)]] <- t2
  }
}

data_case_match <- {}
for (m in 1:1000){
  data_case_match[[m]] <- data.frame()
  for (id in list$V1){
    data_case_match[[m]] <- rbind(data_case_match[[m]],all_annot_match[[paste0(id,"_",m)]])
  }
  data <- rbind(all_annot_flt,data_case_match[[m]][,c(1:113,118)])
  write.table(data,paste0("~/lrs_do_novo/data/data4test/data4test_",m),quote = F,row.names = F,sep="\t")
}


quantiVar <- read.table("~/lrs_do_novo/res_annotSV/list_quantitative_variables")
qualiVar <- read.table("~/lrs_do_novo/res_annotSV/list_qualitative_variables")


res_qua <- {}
for (m in 1:1000){
  data <- read.table(paste0("~/lrs_do_novo/data/data4test/data4test_",m),header = T,sep="\t")
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


for (var in c("Gene_name_b","RE_gene_b","OMIM_ID_b","P_all_phen_b")) {
  results <- data.frame()
  for (m in 1:1000){
    results <- rbind(results, res_qua[[paste0(var,"_",m)]])
  }
  names(results) <- c("case_0","case_1","match_0","match_1","p")
  attach(results)
  results$OR <- case_1*match_0/match_1/case_0
  detach(results)
  write.table(results,paste0("~/lrs_do_novo/data/data4test/","res_comparison_",var),quote = F,row.names = F)
}

