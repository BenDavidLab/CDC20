profile_path <- Sys.getenv("R_PROFILE")
if (!is.na(profile_path) && file.exists(profile_path)) {
  source(profile_path)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(effsize)
  library(foreach)
  library(stringr)
})

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("expecting two arg", call.=FALSE)
}

treatment_index = as.numeric(args[1])
output_dir = args[2]

get_cor <- function(data, outcome, v) {
  top25_val <- quantile(data[,v], 0.75) %>% unname()
  bottom_25_val <- quantile(data[,v], 0.25) %>% unname()
  top_25 <- data[data[,v] >= top25_val,]
  bottom_25 <- data[data[,v] <= bottom_25_val,]
  
  x = top_25[,outcome]
  y = bottom_25[,outcome]
  ttest_pvalue <- t.test(x,y)$p.value
  ttest_cohens_d <- cohen.d(x,y)$estimate

  
  return(data.frame(ttest_pvalue = ttest_pvalue, ttest_cohens_d = ttest_cohens_d))
}


drugs_sensitivity <- read.csv('data/PRISM_Repurposing_Public_23Q2_subsetted.csv')
drugs_sensitivity <- drugs_sensitivity[,c(1,treatment_index + 1)]
colnames(drugs_sensitivity)[1] <- 'DepMap_ID'

expression <- read.csv('data/CCLE_expression.csv')
colnames(expression) <- gsub('\\.\\..*','',colnames(expression))
expression_zscore <- data.frame(scale(expression[,-1]))
expression_zscore <- cbind(expression[,1], expression_zscore)
colnames(expression_zscore)[1] <- 'DepMap_ID'
colnames(expression_zscore)[-1] <-  colnames(expression)[-1]

cell_line_data <- merge(expression_zscore, drugs_sensitivity, by = 'DepMap_ID', all = TRUE)

treatment_name = colnames(drugs_sensitivity)[2]
print(paste0('working on treatment number ',treatment_index))
treatment_gene_effects_cor <- data.frame()
i <- 0
for(gene in colnames(expression)[-1]){
  tryCatch({
    i <- i + 1
    if(i %% 100 == 0){
      print(paste0('treatment number ',treatment_index,': working on gene ',i))
    }
    gene_treatment_pair <- dplyr::select(cell_line_data, DepMap_ID, all_of(c(gene, treatment_name))) %>% na.omit()
    colnames(gene_treatment_pair)[colnames(gene_treatment_pair) == treatment_name] <- 'treatment'
    cor_res <- get_cor(gene_treatment_pair, 'treatment', gene)
    cor_res$gene <- gene
    treatment_gene_effects_cor <- rbind(treatment_gene_effects_cor, cor_res)
  }, error = function(e) {
    print(paste0('error in gene ',gene))
    print(e)
  })
}

write.csv(treatment_gene_effects_cor, paste0(output_dir,'/',treatment_index,'.csv'),
          row.names = F)

print(paste0('done with treatment number ',treatment_index))
