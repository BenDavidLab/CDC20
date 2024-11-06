library(jtools)
library(dplyr)
library(stringr)
library(effsize)
library(xlsx)

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
drug_info <- read.csv("data/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv")
drugs <- c('MPI-0479605','AZ3146')
drugs <- drug_info %>% dplyr::filter(tolower(Drug.Name) %in% tolower(drugs)) %>% dplyr::select(name = Drug.Name, broad_id =IDs) %>% unique()

###drug sensitivity data####
primary_screen_23Q2_data <- read.csv("data/PRISM_Repurposing_Public_23Q2_subsetted.csv")
colnames(primary_screen_23Q2_data)[1] <- 'DepMap_ID'
for(i in 1:nrow(drugs)){
  col_indx <- which(str_detect(colnames(primary_screen_23Q2_data), gsub(':','\\.',gsub('-','\\.',drugs$broad_id[i]))))
  colnames(primary_screen_23Q2_data)[col_indx] <- drugs$name[i]
}
drugs_sensitivity <- primary_screen_23Q2_data %>% select(DepMap_ID, any_of(drugs$name))

####expression data####
expression <- read.csv("data/CCLE_expression.csv")
colnames(expression) <- gsub('\\.\\..*','',colnames(expression))
expression_zscore <- data.frame(scale(expression[,-1]))
expression_zscore <- cbind(expression[,1], expression_zscore)
colnames(expression_zscore)[1] <- 'DepMap_ID'

###aneuploidy data####
aneuploidy_data <- read.csv("data/Depmap AS updated 1700.csv") %>% select(DepMap_ID, num_arm_events)

###merge data####
cell_line_data <- merge(expression_zscore, drugs_sensitivity, by = 'DepMap_ID', all = TRUE)
cell_line_data <- merge(cell_line_data, aneuploidy_data, by = 'DepMap_ID', all = TRUE)

###get partial correlation####
gene = 'CDC20'

if(!dir.exists('results/drug_sensitivity_vs_as')){
  dir.create('results/drug_sensitivity_vs_as', recursive = TRUE)
}

# Loop over each drug to create Excel files
for (i in 1:nrow(drugs)) {
  # Get the drug name
  drug_name <- drugs$name[i]
  
  # Prepare the data for partial correlation adjustment
  treatments_and_exp <- dplyr::select(cell_line_data, DepMap_ID, num_arm_events, any_of(c(gene, drug_name))) %>% na.omit()
  
  # Rename the drug column to 'treatment' for compatibility with the function
  colnames(treatments_and_exp)[colnames(treatments_and_exp) == drug_name] <- 'treatment'
  
  # Perform partial correlation adjustment
  result <- get_partlize_cor(treatments_and_exp, 'treatment', 'num_arm_events', gene)
  
  # Prepare the data for saving
  output_data <- data.frame(
    DepMap_ID = result$data$DepMap_ID,
    aneuploidy_score = result$data$num_arm_events,
    drug_sensitivity_unadjusted = result$data$treatment,
    drug_sensitivity_adjusted = result$data$adj_treatment
  )
  
  # Save the results to an Excel file
  write.xlsx(output_data, file = paste0("results/drug_sensitivity_vs_as/",drug_name, "_sensitivity_adjustment.xlsx"), row.names = FALSE)
}
