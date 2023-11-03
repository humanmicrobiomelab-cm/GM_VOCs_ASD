
df_upset_xl= read_xlsx("upset_with_condition.xlsx")
df_upset_xl=df_upset_xl[1:64,]

library("UpSetR")
BiocManager::install("ComplexHeatmap")

m_upset= make_comb_mat(df_upset_xl)
order_strates= c("Condition_CTRL", "Condition_ASD","Male", "Female", ">5 years","with GI symptoms","without GI symptoms","low symptoms","Picky eaters","No picky eaters","CBCL-EXT","risk CBCL-INT","risk CBCL-EXT", "no clinical symptoms INT","no clinical symptoms EXT",
                 "Probiotics","Antibiotics","No probiotics","No antibiotics","with CI/DD","without CI/DD")
coul <- brewer.pal(11, "Spectral") 
UpSet(m_upset, set_order = order_strates, comb_order = order(comb_size(m_upset)), 
      pt_size = unit(4, "mm"), lwd = 4, comb_col = coul[comb_degree(m_upset)],
      top_annotation = upset_top_annotation(m_upset,height = unit(4, "cm"), add_numbers = TRUE),
      right_annotation = upset_right_annotation(m_upset,
                                                ylim= c(0,30),
                                                gp = gpar(fill = "darkgray"))  