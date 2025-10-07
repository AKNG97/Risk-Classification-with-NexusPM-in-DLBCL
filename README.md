The scripts in this repository can be used to reproduce the results of the manuscript titled "Risk stratification models based on co-expression signatures for RCHOP and GCHOP treated DLBCL patients". The analysis has been organized into four different scripts than should be run in the following order:

1_DataProcessing.R: This script contains the code to retrieve and process data from the GDC and the GEO database. It also contains the code to construct single-sample networks using LIONESS and to perform feature selection from these networks to predict overall-survival and progression-free survival.
2_NullModels.R: This script performs the same methodological pipeline for feature selection as the previous script, but using random features.
3_EvaluationOnTestData.R: This script evaluates the selected models to predict treatment response in RCHOP and GCHOP treated patients. This evaluation is performed in the random test partitions and (for the RCHOP model) across cohorts.
4_DiffCoExpr_FEA.R: This script was used to evaluate the neighborhood of the genes in the models to predict treatment response. First, the script generates a stable networks across multiple random partitions, and uses them to perform a differential co-expession analysis between risk groups. Then, the script performs a functional enrichment analysis on communities to associate biological processes (from Gene Ontology) to the neighborhoods of relevant genes. Lastly, it calculates whether the neighborhood of these genes is biased towards intra- or inter-chromosomal associations.

The Clinical information and gene annotation data necessary for this analysis was retrieved as follows:

1. The Gene annotation file to filter the GOYA data was downloaded from https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html (Human.GRCh38.p13.annot.tsv.gz). Note that an additional gene annotation file from Biomart is used on the pipeline, but this is downloaded using the 1_DataProcessing.R script.
 
2. For the NCICCR clinical data, we used the Supplementary Appendix 2 (Tab S9 Characteristics DLBCL) from: Schmitz, R., Wright, G. W., Huang, D. W., Johnson, C. A., Phelan, J. D., Wang, J. Q., ... & Staudt, L. M. (2018). Genetics and pathogenesis of diffuse large B-cell lymphoma. New England journal of medicine, 378(15), 1396-1407.

3. For the GOYA clinical data, the Clinical information was downloaded from https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE125966



