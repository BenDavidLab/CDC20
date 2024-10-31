library(dplyr)

gene_name = 'CDC20'
drugs_sensitivity <- read.csv('data/PRISM_Repurposing_Public_23Q2_subsetted.csv')
drugs <- colnames(drugs_sensitivity)[-1]

res <- data.frame()
for(treatment_index in seq(length(drugs))){
  print(paste0('Working on ', treatment_index))
  treatment_name <- drugs[treatment_index]
  analysis <- read.csv(paste0('results/PRISM_23_correlations/',treatment_index,'.csv'))
  
  if(nrow(analysis) == 0){
    print(paste0(treatment_name, ' is empty for ', convariat))
    next
  }
  
  ranked_res <- analysis %>% dplyr::filter(ttest_cohens_d < 0) %>% arrange(ttest_pvalue)
  rank = which(ranked_res$gene == gene_name)
  if(!gene_name %in% ranked_res$gene){
    rank = -1
  }
  total = nrow(ranked_res)
  row <- data.frame(treatment = treatment_name, rank = rank, total = total, precent = 100*rank/total)
  res <- rbind(res, row)
}
  
write.csv(res, 'treatment_vs_gene_corr_ttest_ranking_neg_effect_only_PRISM_23.csv',row.names = F)