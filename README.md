# VOCs_ASD

The folder named "scripts" contains the main scripts used to perform the analyses in the paper "Gut Microbiota functional profiling in autism spectrum disorders: bacterial VOCs and related metabolic pathways acting as disease biomarkers and predictors" by Vernocchi et al., i.e.:

-Fold change and Univariate analysis.R : this script was used to perform the differential analysis on VOCs concentrations and to compute their Fold change between two conditions (ASD vs CTRLs).
-Multivariate analysis.R : this script was used to perform the multivariate analysis, in particular PCA and PLS-DA on VOCs concentrations considering two conditions ASD and CTRL.
-Pheatmap.R : this script was used to plot the heatmap of VOCs distributions with annotation based on conditions ASD or CTRL, using R "pheatmap" package. The same script was used to perform the figures of all heatmaps in the               paper.
-Correlation matrix.R : this script was used to construct the correlation matrix, specifically to correlate OTUs and VOCs and to select only the statistically significant (p-value adjusted ≤ 0.05). The correlations matrix                            was visualized as correlation network with Cytoscape v3.8.2. [You need to](https://cytoscape.org/).
