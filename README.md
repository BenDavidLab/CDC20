# CDC20
Code for the "CDC20 determines the sensitivity to MPS1 inhibitors" paper.

**For figure 2E:** 
  1. Donwload CCLE_expression.csv from [DepMap Public 22Q2 Files](https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap+Public+22Q2&filename=CCLE_expression_full.csv) to data/ dir.
  2. Download The "PRISM Repurposing Public" file from [DepMap custom download](https://depmap.org/portal/data_page/?tab=customDownloads) to data/ dir (we used 23Q2).
  3. Run parallel_executor.sh. You can specify .Rprofile and number of workers to use in the script. This can take a few days to complete... precomputed results for steps 3 and 4 are available at treatment_vs_gene_corr_ttest_ranking_neg_effect_only_PRISM_23.csv.
  4. Run the calculate_ranking.R script.
  5. Run the plotting_cdc20_gene_ranking script.

**For figure F2B and supplementary 3A:**
  1. Download Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv from [DepMap Public 23Q2 Files](https://depmap.org/portal/data_page/?tab=allData&releasename=PRISM+Repurposing+Public+23Q2&filename=Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv) to data/ dir.
  3. Donwload CCLE_expression.csv from [DepMap Public 22Q2 Files](https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap+Public+22Q2&filename=CCLE_expression_full.csv) to data/ dir.
  4.  Run anp_score_vs_drug_sensitivity.R script.
