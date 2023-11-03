# VOCs_ASD

The folder named "scripts" contains the main scripts used to perform the analyses in the paper "Gut Microbiota functional profiling in autism spectrum disorders: bacterial VOCs and related metabolic pathways acting as disease biomarkers and predictors" by Vernocchi et al., i.e.:

-**Fold change and Univariate analysis.R** : this script was used to perform the differential analysis on VOCs concentrations and to compute their Fold change between two conditions (ASD vs CTRLs).

-**Multivariate analysis.R** : this script was used to perform the multivariate analysis, in particular PCA and PLS-DA on VOCs concentrations considering two conditions ASD and CTRL.

-[Pheatmap.R](Pheatmap.R) : this script was used to plot the heatmap of VOCs distributions with annotation based on conditions ASD or CTRL, using R "pheatmap" package. The same script was used to perform the figures of all heatmaps in                  the paper.

-Correlation matrix.R : this script was used to construct the correlation matrix, specifically to correlate OTUs and VOCs and to select only the statistically significant (p-value adjusted ≤ 0.05) correlations. 
                        The correlations matrix was visualized as correlation network with [Cytoscape](https://cytoscape.org/) v3.8.2.
                        
-ROC_curve.R: this script was used to apply machine learning (ML) with Logistic regression model, to measure the accuracy, sensitivity and specificity of the model and to represent the Receiver operating characteristic 
              (ROC) curves with the relative area under curve (AUC) values.

All the R scripts along with a set of tutorials for performing WGCNA analysis are freely available at this [link](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/).
The [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html) R package can be downloaded from CRAN repository.      
