library(jtools)
library(dplyr)

###function declarations####
get_partlize_cor <- function(data, outcome, v, covariat) {
  full_model <- lm(as.formula(paste0(outcome,'~',v,'+',covariat)), data = data)
  par_model <- partialize(full_model, vars = covariat, center = TRUE, scale = "response", set.offset = 1)
  colnames(par_model)[1] <- paste0('adj_',colnames(par_model)[1])
  par_data <- merge(data, par_model[,c(1,3)], by=covariat)
  adj_cor <- cor.test(par_data[,v], par_data[,colnames(par_model)[1]], method = "spearman")
  unadj_cor <- cor.test(data[,v], data[,outcome], method = "spearman")
  
  top25_val <- quantile(par_data[,v], 0.75) %>% unname()
  bottom_25_val <- quantile(par_data[,v], 0.25) %>% unname()
  top_25 <- par_data[par_data[,v] >= top25_val,]
  bottom_25 <- par_data[par_data[,v] <= bottom_25_val,]
  
  x = top_25[,outcome]
  y = bottom_25[,outcome]
  unadj_ttest_pvalue <- t.test(x,y)$p.value
  unadj_ttest_cohens_d <- cohen.d(x,y)$estimate
  
  x = top_25[,colnames(par_model)[1]]
  y = bottom_25[,colnames(par_model)[1]]
  adj_ttest_pvalue <- t.test(x,y)$p.value
  adj_ttest_cohens_d <- cohen.d(x,y)$estimate
  
  return(list(data = par_data, cor = data.frame(adj_pval = adj_cor$p.value, adj_rho = adj_cor$estimate,
                                                unadj_pval = unadj_cor$p.value, unadj_rho = unadj_cor$estimate,
                                                adj_ttest_pval = adj_ttest_pvalue, adj_ttest_cohens_d = adj_ttest_cohens_d,
                                                unadj_ttest_pvalue = unadj_ttest_pvalue, unadj_ttest_cohens_d = unadj_ttest_cohens_d)))
}
###drug name preprocessing####
drug_info <- read.csv('data/primary-screen-replicate-treatment-info.csv')
drugs <- c('MPI-0479605','AZ3146')
drugs <- drug_info %>% dplyr::filter(tolower(name) %in% tolower(drugs)) %>% dplyr::select(name, broad_id) %>% unique()

###drug sensitivity data####
primary_screen_replicate_collapsed_logfold_change <- read.csv('data/primary-screen-replicate-collapsed-logfold-change.csv')
colnames(primary_screen_replicate_collapsed_logfold_change)[1] <- 'DepMap_ID'
for(i in 1:nrow(drugs)){
  col_indx <- which(str_detect(colnames(primary_screen_replicate_collapsed_logfold_change), gsub('-','\\.',drugs$broad_id[i])))
  colnames(primary_screen_replicate_collapsed_logfold_change)[col_indx] <- drugs$name[i]
}
drugs_sensitivity <- primary_screen_replicate_collapsed_logfold_change %>% select(DepMap_ID, any_of(drugs$name))

####expression data####
expression <- read.csv('data/CCLE_expression.csv')
colnames(expression) <- gsub('\\.\\..*','',colnames(expression))
expression_zscore <- data.frame(scale(expression[,-1]))
expression_zscore <- cbind(expression[,1], expression_zscore)
colnames(expression_zscore)[1] <- 'DepMap_ID'

###aneuploidy data####
aneuploidy_data <- read.csv('data/Depmap AS updated 1700.csv') %>% select(DepMap_ID, num_arm_events)

###merge data####
cell_line_data <- merge(expression_zscore, drugs_sensitivity, by = 'DepMap_ID', all = TRUE)
cell_line_data <- merge(cell_line_data, aneuploidy_data, by = 'DepMap_ID', all = TRUE)

###get partial correlation####
gene = 'CDC20'
for(treatment_name in drugs){
  treatments_and_exp  <- dplyr::select(cell_line_data, DepMap_ID, num_arm_events, any_of(c(gene, treatment_name))) %>% na.omit() 
  colnames(treatments_and_exp)[colnames(treatments_and_exp) == treatment_name] <- 'treatment'
  x = get_partlize_cor(treatments_and_exp, 'treatment', 'num_arm_events', gene)
  write.csv(x$data, paste0('results/',treatment_name,'_',gene,'_adjusted_drug_senstevity.csv'))
}

