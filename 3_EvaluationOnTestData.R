
home.dir <- file.path(path.expand("~"), "NexusPM")
setwd(home.dir)

Raw.dir <- file.path(home.dir, "RawData/")
ClinicalInfo <- file.path(home.dir, "ClinicalInfo/")
Norm.dir <- file.path(home.dir, "NormData/")
NullModels.dir <- file.path(home.dir, "NullModels/")
Validation.dir <- file.path(home.dir, "Validation/")
FEA.dir <- file.path(home.dir, "FEA/")
#NOTE: Download the files StableCommunities_RCHOP.xlsx and StableCommunities_GCHOP.xlsx 
# and save into the FEA folder before running this script

#### Get all samples SSNs real model ####

#For RCHOP

SSN_i <- read.csv("SSNs/SSN_RCHOP_EdgesInModels_Null2.csv", row.names = 1)

edgesRCHOP <- readRDS(file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))

SSN_i <- SSN_i[rownames(SSN_i) %in% edgesRCHOP, ]

write.csv(SSN_i, "SSNs/SSN_RCHOP_AllSamples.csv", quote = F, row.names = T)

#For GCHOP
library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Normal <- data.frame(Sample = NormalSamples, 
                        COO = NA,
                        PFS_time = NA,
                        PFS_status = NA,
                        OS_time = NA, 
                        OS_status = NA,
                        Gender = NA, 
                        Treatment = "Lymphoid normal",
                        Source = "GDC",
                        primary_diagnosis = "Lymphoid normal")

CI.Full <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

DLBCL.norm <- readRDS("NormData/DLBCL_NormCounts.RDS")

CI.GCHOP <- CI.Full.Normal %>% filter(Treatment == "GA101" |
                                        primary_diagnosis == "Lymphoid normal")

relevantEdges <- readRDS(file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))

GCN <- data.frame(Edge = relevantEdges, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

TARGET_NBM_ALALi_norm <- DLBCL.norm[,colnames(DLBCL.norm) %in% CI.GCHOP$Sample]
rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))

GCN_i <- GCN %>% 
  dplyr::mutate(Rho = cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,],
                          TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,],
                          method = "spearman"))

#write.csv(GCN_i, paste0(GCNs, "GCN_GOYA_Sample_", i, ".csv"), quote = F, row.names = F)
n_samples <- ncol(TARGET_NBM_ALALi_norm)

#Get Lioness
for(j in 1:nrow(GCN_i)){
  
  Source <- GCN_i$Source[j]
  Target <- GCN_i$Target[j]
  alpha <- GCN_i$Rho[j]
  SSN_j <- numeric()
  
  for(x in 1:ncol(TARGET_NBM_ALALi_norm)){
    
    j_cor <- stats::cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,-x],
                        TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,-x], 
                        method="spearman")
    SSN_cor <- round(n_samples*(alpha - j_cor) + j_cor, 4)
    #names(SSN_cor) <- paste(Source, Target, sep = "_")
    
    if(length(SSN_j) == 0){
      SSN_j <- SSN_cor
    } else {
      SSN_j <- c(SSN_j, SSN_cor)
    }
  }
  
  if(j == 1){
    SSN_i <- SSN_j
  } else {
    SSN_i <- rbind(SSN_i, SSN_j)
  }
}

SSN_i <- as.data.frame(SSN_i)
rownames(SSN_i) <- GCN_i$Edge
colnames(SSN_i) <- colnames(TARGET_NBM_ALALi_norm)

write.csv(SSN_i, "SSNs/SSN_GCHOP_AllSamples.csv", quote = F, row.names = T)

#### Test the proportionality assumption across random partitions ####

library(dplyr)
library(tidyr)
library(survival)
library(future.apply)
library(glmnet)
library(MASS)

RCHOP.edges <- readRDS(file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))
GCHOP.edges <- readRDS(file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

CI.Global <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Global$PFS_status[CI.Global$PFS_time > 1825] <- 0
CI.Global$PFS_time[CI.Global$PFS_time > 1825] <- 1825

CI.RCHOP <- CI.Global %>% filter(Treatment == "Rituximab")
CI.GCHOP <- CI.Global %>% filter(Treatment == "GA101")

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% RCHOP.edges,]
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
  Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
  
  CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
  
  CoxSurvCI <- CoxSurvCI %>% 
    inner_join(CI.RCHOP %>% 
                 dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
  
  Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                     data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
  
  Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                    data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
  
  ProportionalHazard.RCHOP.PFS <- cox.zph(Cox.PFS.1,  "identity")$table %>% as.data.frame %>%
    mutate(Partition = i, Type = "PFS")
  
  ProportionalHazard.RCHOP.PFS <- ProportionalHazard.RCHOP.PFS %>% 
    mutate(Edge = rownames(ProportionalHazard.RCHOP.PFS))
  
  ProportionalHazard.RCHOP.OS <- cox.zph(Cox.OS.1,  "identity")$table %>% as.data.frame %>%
    mutate(Partition = i, Type = "OS")
  
  ProportionalHazard.RCHOP.OS <- ProportionalHazard.RCHOP.OS %>% 
    mutate(Edge = rownames(ProportionalHazard.RCHOP.OS))
  
  if(i == 1){
    ProportionalHazard.RCHOP.PFS.summary <- ProportionalHazard.RCHOP.PFS
    ProportionalHazard.RCHOP.OS.summary <- ProportionalHazard.RCHOP.OS
  } else {
    ProportionalHazard.RCHOP.PFS.summary <- rbind(ProportionalHazard.RCHOP.PFS.summary, ProportionalHazard.RCHOP.PFS)
    ProportionalHazard.RCHOP.OS.summary <- rbind(ProportionalHazard.RCHOP.OS.summary, ProportionalHazard.RCHOP.OS)
    
  }
  
}

ProportionalHazard.RCHOP.OS.summary <- tibble(ProportionalHazard.RCHOP.OS.summary) %>% 
  group_by(Edge, Type) %>%
  mutate(adj.p = p.adjust(p, method = "fdr")) 
ProportionalHazard.RCHOP.PFS.summary <- tibble(ProportionalHazard.RCHOP.PFS.summary) %>% 
  group_by(Edge, Type) %>%
  mutate(adj.p = p.adjust(p, method = "fdr")) 

ProportionalHazard.RCHOP.OS.summary %>%  summarise(Median.adj.p = median(adj.p),
                                                   Median.p = median(p),
                                                   FreqSignificant = sum(p < 0.05) / n())
# # A tibble: 15 × 5
# # Groups:   Edge [15]
#   Edge           Type  Median.adj.p Median.p FreqSignificant
#   <chr>          <chr>        <dbl>    <dbl>           <dbl>
# 1 ANKIB1_NUP42   OS          0.313    0.190              0.1
# 2 ATRIP_GMPPB    OS          0.0948   0.0535             0.4
# 3 ATR_GFM1       OS          0.594    0.348              0  
# 4 B4GALT3_POLR3C OS          0.929    0.593              0  
# 5 CCDC85C_POLR1B OS          0.968    0.802              0  
# 6 CD4_CST3       OS          0.882    0.585              0.1
# 7 CD69_JUNB      OS          0.977    0.714              0  
# 8 CYP27A1_VNN1   OS          0.633    0.419              0  
# 9 FOLR2_SLC46A1  OS          0.767    0.561              0  
# 10 GLOBAL         OS          0.569    0.323              0.1
# 11 KIFBP_METTL14  OS          0.756    0.434              0.1
# 12 NUDCD1_OTUD6B  OS          0.939    0.706              0  
# 13 PHF20L1_POLR2K OS          0.311    0.179              0.1
# 14 SLC52A2_ZNF16  OS          0.910    0.711              0  
# 15 ZNF429_ZNF681  OS          0.845    0.530              0.1

ProportionalHazard.RCHOP.PFS.summary %>%  summarise(Median.adj.p = median(adj.p),
                                                    Median.p = median(p),
                                                    FreqSignificant = sum(p < 0.05) / n())
# A tibble: 15 × 5
# Groups:   Edge [15]
#   Edge           Type  Median.adj.p Median.p FreqSignificant
#   <chr>          <chr>        <dbl>    <dbl>           <dbl>
# 1 ANKIB1_NUP42   PFS         0.464   0.276               0.1
# 2 ATRIP_GMPPB    PFS         0.463   0.269               0  
# 3 ATR_GFM1       PFS         0.127   0.0755              0.4
# 4 B4GALT3_POLR3C PFS         0.730   0.458               0.1
# 5 CCDC85C_POLR1B PFS         0.935   0.540               0  
# 6 CD4_CST3       PFS         0.825   0.575               0  
# 7 CD69_JUNB      PFS         0.956   0.742               0  
# 8 CYP27A1_VNN1   PFS         0.267   0.184               0  
# 9 FOLR2_SLC46A1  PFS         0.0868  0.0496              0.5
# 10 GLOBAL         PFS         0.0165  0.00912             0.9
# 11 KIFBP_METTL14  PFS         0.0113  0.00675             0.8
# 12 NUDCD1_OTUD6B  PFS         0.851   0.470               0  
# 13 PHF20L1_POLR2K PFS         0.258   0.153               0.1
# 14 SLC52A2_ZNF16  PFS         0.724   0.482               0  
# 15 ZNF429_ZNF681  PFS         0.730   0.423               0.1

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% GCHOP.edges,]
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
  Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
  
  CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
  
  CoxSurvCI <- CoxSurvCI %>% 
    inner_join(CI.GCHOP %>% 
                 dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
  
  Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                     data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
  
  Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                    data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
  
  ProportionalHazard.GCHOP.PFS <- cox.zph(Cox.PFS.1,  "identity")$table %>% as.data.frame %>%
    mutate(Partition = i, Type = "PFS")
  
  ProportionalHazard.GCHOP.PFS <- ProportionalHazard.GCHOP.PFS %>% 
    mutate(Edge = rownames(ProportionalHazard.GCHOP.PFS))
  
  ProportionalHazard.GCHOP.OS <- cox.zph(Cox.OS.1,  "identity")$table %>% as.data.frame %>%
    mutate(Partition = i, Type = "OS")
  
  ProportionalHazard.GCHOP.OS <- ProportionalHazard.GCHOP.OS %>% 
    mutate(Edge = rownames(ProportionalHazard.GCHOP.OS))
  
  if(i == 1){
    ProportionalHazard.GCHOP.PFS.summary <- ProportionalHazard.GCHOP.PFS
    ProportionalHazard.GCHOP.OS.summary <- ProportionalHazard.GCHOP.OS
  } else {
    ProportionalHazard.GCHOP.PFS.summary <- rbind(ProportionalHazard.GCHOP.PFS.summary, ProportionalHazard.GCHOP.PFS)
    ProportionalHazard.GCHOP.OS.summary <- rbind(ProportionalHazard.GCHOP.OS.summary, ProportionalHazard.GCHOP.OS)
    
  }
  
}

ProportionalHazard.GCHOP.OS.summary <- tibble(ProportionalHazard.GCHOP.OS.summary) %>% 
  group_by(Edge, Type) %>%
  mutate(adj.p = p.adjust(p, method = "fdr")) 
ProportionalHazard.GCHOP.PFS.summary <- tibble(ProportionalHazard.GCHOP.PFS.summary) %>% 
  group_by(Edge, Type) %>%
  mutate(adj.p = p.adjust(p, method = "fdr")) 

ProportionalHazard.GCHOP.OS.summary %>%  summarise(Median.adj.p = median(adj.p),
                                                   Median.p = median(p),
                                                   FreqSignificant = sum(p < 0.05) / n())

# A tibble: 8 × 5
# Groups:   Edge [8]
# Edge          Type  Median.adj.p Median.p FreqSignificant
# <chr>         <chr>        <dbl>    <dbl>           <dbl>
#   1 ANKRD13D_CUL5 OS          0.884    0.558              0  
# 2 ASF1A_SEC63   OS          0.807    0.498              0  
# 3 CNOT10_RAB5A  OS          0.0295   0.0166             0.7
# 4 FRG1_PLRG1    OS          0.908    0.657              0  
# 5 GLOBAL        OS          0.118    0.0671             0.5
# 6 HMGCR_SHISA8  OS          0.833    0.475              0  
# 7 ICMT_NOL9     OS          0.830    0.672              0  
# 8 JMJD6_PSMD12  OS          0.921    0.650              0  

ProportionalHazard.GCHOP.PFS.summary %>%  summarise(Median.adj.p = median(adj.p),
                                                    Median.p = median(p),
                                                    FreqSignificant = sum(p < 0.05) / n())
# A tibble: 8 × 5
# Groups:   Edge [8]
# Edge          Type  Median.adj.p Median.p FreqSignificant
# <chr>         <chr>        <dbl>    <dbl>           <dbl>
#   1 ANKRD13D_CUL5 PFS          0.538    0.296             0  
# 2 ASF1A_SEC63   PFS          0.918    0.652             0  
# 3 CNOT10_RAB5A  PFS          0.239    0.135             0.1
# 4 FRG1_PLRG1    PFS          0.855    0.488             0  
# 5 GLOBAL        PFS          0.374    0.233             0.1
# 6 HMGCR_SHISA8  PFS          0.261    0.144             0.1
# 7 ICMT_NOL9     PFS          0.931    0.719             0  
# 8 JMJD6_PSMD12  PFS          0.634    0.363             0  

saveRDS(ProportionalHazard.GCHOP.OS.summary, file.path(Validation.dir, "ProportionalHazardGCHOP_OS.RDS"))
saveRDS(ProportionalHazard.GCHOP.PFS.summary, file.path(Validation.dir, "ProportionalHazardGCHOP_PFS.RDS"))

saveRDS(ProportionalHazard.RCHOP.OS.summary, file.path(Validation.dir, "ProportionalHazardRCHOP_OS.RDS"))
saveRDS(ProportionalHazard.RCHOP.PFS.summary, file.path(Validation.dir, "ProportionalHazardRCHOP_PFS.RDS"))

#### Norm Test data ####
library(SummarizedExperiment)

library(EDASeq)
library(dplyr)
library(NOISeq)
library(DESeq2)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(Hmisc)
library(tidyr)
library(M3C)
library(rlist)
library(jsonlite)
library(future.apply)
library(ggbiplot)

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Normal <- data.frame(Sample = NormalSamples, 
                        COO = NA,
                        PFS_time = NA,
                        PFS_status = NA,
                        OS_time = NA, 
                        OS_status = NA,
                        Gender = NA, 
                        Treatment = "Lymphoid normal",
                        Source = "GDC",
                        primary_diagnosis = "Lymphoid normal")

CI.Full <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Full <- CI.Full %>% mutate(primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

annot <- readRDS(file.path(Raw.dir, "annot_GDC.RDS"))

annot <- annot[annot$HGNC_symbol %in% rownames(Test.Data),]

coeff.R <- readRDS(file.path(ElasticNet_results.dir,"RCHOP_SelectedEdges.RDS"))
coeff.G <- readRDS(file.path(ElasticNet_results.dir,"GCHOP_SelectedEdges.RDS"))

norm.A <- function(x, y, z){
  x <- x[match(y$HGNC_symbol, rownames(x)), ]
  x <- as.matrix(x)
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factor = z)
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, factor = "primary_diagnosis",
                                  batch = F, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

NormTest <- function(x){
  Samplei <- colnames(Test.Data)[x]
  
  Expri <- cbind(Train.Data, Test.Data[,x])
  
  colnames(Expri)[ncol(Train.Data) + 1] <- Samplei
  
  CI <- CI.Full.Normal %>% filter(Sample %in% colnames(Expri))
  CI <- CI[match(colnames(Expri), CI$Sample),]
  
  CI <- CI %>%  dplyr::select(Source, primary_diagnosis, Treatment)
  
  Expri_norm <- norm.A(as.matrix(Expri), annot, CI)
  
  saveRDS(Expri_norm[rownames(Expri_norm) %in% RelevantNodes,],
          paste0(TestNormd.dir, "Expr_TestSample_", x, ".RDS"))
  
}

RelevantEdge <- unique(c(unique(coeff.R), unique(coeff.G)))

RelevantNodes <- stringr::str_split(RelevantEdge, "_") %>% unlist %>% unique

future::plan(multisession, workers = 30)
for(i in 1:10){
  
  print(i)
  
  Test.Data <- readRDS(paste0("RawData/RandomPartitionsRAW/Test_RandomPartition_", i, ".RDS"))
  Train.Data <- readRDS(paste0("RawData/RandomPartitionsRAW/Train_RandomPartition_", i, ".RDS"))
  
  TestNormd.dir <- paste0("NormData/TestSet", i, "/")
  
  dir.create(TestNormd.dir)
  
  results <- future_lapply(1:ncol(Test.Data), 
                           FUN = NormTest, future.seed = TRUE)
}

future::plan(sequential)

#### Get SSNs RCHOP ####

library(dplyr)
library(tidyr)
library(future.apply)

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Normal <- data.frame(Sample = NormalSamples, 
                        COO = NA,
                        PFS_time = NA,
                        PFS_status = NA,
                        OS_time = NA, 
                        OS_status = NA,
                        Gender = NA, 
                        Treatment = "Lymphoid normal",
                        Source = "GDC",
                        primary_diagnosis = "Lymphoid normal")

CI.Full <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Full <- CI.Full %>% mutate(primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist

relevantEdges <- readRDS(file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))

dir.create("NormData/SSNsTests")

GCN <- data.frame(Edge = relevantEdges, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

IterateGCN_Lioness <- function(file) {
  
  TARGET_NBM_ALALi_norm <- readRDS(paste0(ExprMs, file))
  
  Samplei <- colnames(TARGET_NBM_ALALi_norm)[ncol(TARGET_NBM_ALALi_norm)]
  
  if(Samplei %in% RCHOP.samples){
    TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% 
                                                     c(RCHOP.samples, NormalSamples)]
    
    rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))
    
    GCN_i <- GCN %>% 
      dplyr::mutate(Rho = cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,],
                              TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,],
                              method = "spearman"))
    
    #write.csv(GCN_i, paste0(GCNs, "GCN_GOYA_Sample_", i, ".csv"), quote = F, row.names = F)
    n_samples <- ncol(TARGET_NBM_ALALi_norm)
    
    #Get Lioness
    for(j in 1:nrow(GCN_i)){
      
      Source <- GCN_i$Source[j]
      Target <- GCN_i$Target[j]
      alpha <- GCN_i$Rho[j]
      SSN_j <- numeric()
      
      for(x in 1:ncol(TARGET_NBM_ALALi_norm)){
        
        j_cor <- stats::cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,-x],
                            TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,-x], 
                            method="spearman")
        SSN_cor <- round(n_samples*(alpha - j_cor) + j_cor, 4)
        #names(SSN_cor) <- paste(Source, Target, sep = "_")
        
        if(length(SSN_j) == 0){
          SSN_j <- SSN_cor
        } else {
          SSN_j <- c(SSN_j, SSN_cor)
        }
      }
      
      if(j == 1){
        SSN_i <- SSN_j
      } else {
        SSN_i <- rbind(SSN_i, SSN_j)
      }
    }
    
    SSN_i <- as.data.frame(SSN_i)
    rownames(SSN_i) <- GCN_i$Edge
    colnames(SSN_i) <- colnames(TARGET_NBM_ALALi_norm)
    
    names.SSN <- stringr::str_split(file, "_") %>% unlist
    names.SSN <- names.SSN[3]
    names.SSN <- gsub(".RDS", ".csv", names.SSN)
    
    write.csv(SSN_i, paste0(SSNs,"SSN_RCHOP_TestSample_Sample_", names.SSN), quote = F, row.names = T)
    
  }
  
}

future::plan(multisession, workers = 15)
for(i in 1:10){
  
  print(i)
  
  ExprMs <- paste0("NormData/TestSet", i, "/")
  SSNs <- paste0(paste0("NormData/SSNsTests/SSN_TestSet", i, "/"))
  dir.create(SSNs)
  results <- future_lapply(list.files(ExprMs),
                           FUN = IterateGCN_Lioness, future.seed = TRUE)
}
future::plan(sequential)

#### Get SSNs GCHOP ####

library(dplyr)
library(tidyr)
library(future.apply)

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Normal <- data.frame(Sample = NormalSamples, 
                        COO = NA,
                        PFS_time = NA,
                        PFS_status = NA,
                        OS_time = NA, 
                        OS_status = NA,
                        Gender = NA, 
                        Treatment = "Lymphoid normal",
                        Source = "GDC",
                        primary_diagnosis = "Lymphoid normal")

CI.Full <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Full <- CI.Full %>% mutate(primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

GCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "GA101") %>%
  dplyr::select(Sample) %>% unlist

relevantEdges <- readRDS(file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))

GCN <- data.frame(Edge = relevantEdges, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

IterateGCN_Lioness <- function(file) {
  
  TARGET_NBM_ALALi_norm <- readRDS(paste0(ExprMs, file))
  
  Samplei <- colnames(TARGET_NBM_ALALi_norm)[ncol(TARGET_NBM_ALALi_norm)]
  
  if(Samplei %in% GCHOP.samples){
    TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% 
                                                     c(GCHOP.samples, NormalSamples)]
    
    rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))
    
    GCN_i <- GCN %>% 
      dplyr::mutate(Rho = cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,],
                              TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,],
                              method = "spearman"))
    
    #write.csv(GCN_i, paste0(GCNs, "GCN_GOYA_Sample_", i, ".csv"), quote = F, row.names = F)
    n_samples <- ncol(TARGET_NBM_ALALi_norm)
    
    #Get Lioness
    for(j in 1:nrow(GCN_i)){
      
      Source <- GCN_i$Source[j]
      Target <- GCN_i$Target[j]
      alpha <- GCN_i$Rho[j]
      SSN_j <- numeric()
      
      for(x in 1:ncol(TARGET_NBM_ALALi_norm)){
        
        j_cor <- stats::cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,-x],
                            TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,-x], 
                            method="spearman")
        SSN_cor <- round(n_samples*(alpha - j_cor) + j_cor, 4)
        #names(SSN_cor) <- paste(Source, Target, sep = "_")
        
        if(length(SSN_j) == 0){
          SSN_j <- SSN_cor
        } else {
          SSN_j <- c(SSN_j, SSN_cor)
        }
      }
      
      if(j == 1){
        SSN_i <- SSN_j
      } else {
        SSN_i <- rbind(SSN_i, SSN_j)
      }
    }
    
    SSN_i <- as.data.frame(SSN_i)
    rownames(SSN_i) <- GCN_i$Edge
    colnames(SSN_i) <- colnames(TARGET_NBM_ALALi_norm)
    
    names.SSN <- stringr::str_split(file, "_") %>% unlist
    names.SSN <- names.SSN[3]
    names.SSN <- gsub(".RDS", ".csv", names.SSN)
    
    write.csv(SSN_i, paste0(SSNs,"SSN_GCHOP_TestSample_Sample_", names.SSN), quote = F, row.names = T)
    
  }
  
}

future::plan(multisession, workers = 15)
for(i in 1:10){
  
  print(i)
  
  ExprMs <- paste0("NormData/TestSet", i, "/")
  SSNs <- paste0(paste0("NormData/SSNsTests/SSN_TestSet", i, "/"))
  dir.create(SSNs)
  results <- future_lapply(list.files(ExprMs),
                           FUN = IterateGCN_Lioness, future.seed = TRUE)
}
future::plan(sequential)

#### Evaluate RCHOP model PFS ####

library(dplyr)
library(tidyr)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(caret)
library(future.apply)
library(timeROC)

SurvCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

SurvCI <- SurvCI %>% filter(Treatment == "Rituximab")

future::plan(multisession, workers = 5)
for(i in 1:10){
  
  print(i)
  
  IterationCox <- function(file){
    #i = 1
    
    SSN_i <- read.csv(paste0(SSNs,file), row.names = 1)
    SSN_i <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]
    
    Sample_i <- colnames(SSN_i)[ncol(SSN_i)]
    
    original_data_relevant <- SSN_i[, 1:(ncol(SSN_i)-1)]
    row_means <- rowMeans(original_data_relevant)
    row_sds <- apply(original_data_relevant, 1, sd)
    
    SSN_i_scaled <- (SSN_i - row_means) / row_sds
    
    e <- t(SSN_i_scaled)
    e <- e %>% as.data.frame %>% mutate(Sample = rownames(e))
    
    SurvCIi <- SurvCI %>% dplyr::select(Sample, PFS_time, PFS_status) %>%
      inner_join(e, by = c("Sample" = "Sample"))
    
    training_data <- SurvCIi %>% filter(Sample != Sample_i) %>% dplyr::select(-Sample)
    
    Cox.Ref <- coxph(Surv(PFS_time, PFS_status) ~ . +tt(KIFBP_METTL14),
                     data = training_data,
                     tt = function(x, t, ...) x * t)
    
    # Cox.Ref <- coxph(Surv(PFS_time, PFS_status) ~ . ,
    #                  data = training_data)
    
    coef.Ref  <- coef(Cox.Ref)
    se_coefs.Ref  <- sqrt(diag(vcov(Cox.Ref)))
    p_values.Ref  <- summary(Cox.Ref)$coefficients[, "Pr(>|z|)"]
    concordance.Ref  <- summary(Cox.Ref)$concordance
    logtest.Ref  <- summary(Cox.Ref)$logtest
    wald_test.Ref  <- summary(Cox.Ref)$waldtest
    score_logrank.Ref  <- summary(Cox.Ref)$sctest
    
    Summary.Ref <- list(
      coefs = coef.Ref,
      se_coefs = se_coefs.Ref,
      p_values = p_values.Ref,
      concordance = concordance.Ref,
      logtest.Ref = logtest.Ref,
      wald_test = wald_test.Ref,
      score_logrank = score_logrank.Ref
    )
    
    saveRDS(Summary.Ref, paste0(SummaryDFS.dir,"SummaryPFS_RCHOP_TestSample_", Sample_i, ".RDS"))
    
    #SurvCIi$LP <- predict(Cox.Ref, type = "lp", newdata = SurvCIi)
    SurvCIi$RiskScore <- predict(Cox.Ref, type = "risk", newdata = SurvCIi)
    
    risk_thresholds <- quantile(SurvCIi %>% filter(Sample != Sample_i) %>% 
                                  dplyr::select(RiskScore) %>% unlist,
                                probs = c(0, 0.25, 0.5, 0.75, 1))
    
    # Crea una nueva columna para clasificar el riesgo
    SurvCIi$RiskCategory <- cut(
      SurvCIi$RiskScore,
      breaks = risk_thresholds,
      labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
      include.lowest = TRUE
    )
    
    Cox.Ref.Cat <- coxph(Surv(PFS_time, PFS_status) ~ .,
                         data = SurvCIi %>% filter(Sample != Sample_i) %>% 
                           dplyr::select(RiskCategory, PFS_time, PFS_status))
    
    MP2_strat <- data.frame(Sample = Sample_i,
                            It_number = i,
                            C.Ref = summary(Cox.Ref)$concordance[1],
                            RiskScore = SurvCIi$RiskScore[SurvCIi$Sample == Sample_i],
                            RiskCat = SurvCIi$RiskCategory[SurvCIi$Sample == Sample_i],
                            C.Ref.Cat = summary(Cox.Ref.Cat)$concordance[1])
    
    if(is.na(MP2_strat$RiskCat)){
      
      MP2_strat <- MP2_strat %>% mutate(RiskCat = case_when(RiskScore > risk_thresholds[5] ~ "High",
                                                            RiskScore < risk_thresholds[1] ~ "Low",
                                                            .default = RiskCat))
      
    }
    
    
    saveRDS(SurvCIi, paste0(SummaryDFS.dir, "SurvCIi_RCHOP_TestSample_", Sample_i, ".RDS"))
    
    MP2_strat <- MP2_strat %>% inner_join(SurvCI, by = c("Sample" = "Sample"))
    
    return(MP2_strat)
    
  }
  
  SummaryDFS.dir <- paste0("NormData/SSNsTests", "/Summaries_RCHOP_Test_tt", i, "/")
  dir.create(SummaryDFS.dir)
  SSNs <- paste0("NormData/SSNsTests/SSN_TestSet", i, "/")
  nfiles <- length(list.files(SSNs)[grep("RCHOP", list.files(SSNs))])
  results <- future_lapply(list.files(SSNs)[grep("RCHOP", list.files(SSNs))],
                           FUN = IterationCox, future.seed = TRUE)
  Global_Strata <- do.call(rbind, results)
  
  write.csv(Global_Strata, paste0("NormData/SSNsTests/Global_RCHOP_Strata_tt_Val", i, ".csv"))
  
  cox_model <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore, 
                     data = Global_Strata)
  
  cox_model.Cat <- coxph(Surv(PFS_time, PFS_status) ~ RiskCat, 
                         data = Global_Strata)
  
  c.risk <- summary(cox_model)$concordance[1]
  c.cat <- summary(cox_model.Cat)$concordance[1]
  
  if(i == 1){
    C.risk <- c.risk
    C.cat <- c.cat
  } else {
    C.risk <- rbind(C.risk, c.risk)
    C.cat <- rbind(C.cat, c.cat)
  }
  
  #print(summary(cox_model)$concordance)
  
}

future::plan(sequential)

saveRDS(C.cat, file.path(Validation.dir, "C_cat_PFS_RCHOP.RDS"))
saveRDS(C.risk, file.path(Validation.dir, "C_risk_PFS_RCHOP.RDS"))


#### Evaluate RCHOP model OS ####

library(dplyr)
library(tidyr)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(caret)
library(future.apply)
library(timeROC)

SurvCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

SurvCI <- SurvCI %>% filter(Treatment == "Rituximab")

future::plan(multisession, workers = 5)
for(i in 1:10){
  
  print(i)
  
  IterationCox <- function(file){
    #i = 1
    
    SSN_i <- read.csv(paste0(SSNs,file), row.names = 1)
    SSN_i <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]
    
    Sample_i <- colnames(SSN_i)[ncol(SSN_i)]
    
    original_data_relevant <- SSN_i[, 1:(ncol(SSN_i)-1)]
    row_means <- rowMeans(original_data_relevant)
    row_sds <- apply(original_data_relevant, 1, sd)
    
    SSN_i_scaled <- (SSN_i - row_means) / row_sds
    
    e <- t(SSN_i_scaled)
    e <- e %>% as.data.frame %>% mutate(Sample = rownames(e))
    
    SurvCIi <- SurvCI %>% dplyr::select(Sample, OS_time, OS_status) %>%
      inner_join(e, by = c("Sample" = "Sample"))
    
    training_data <- SurvCIi %>% filter(Sample != Sample_i) %>% dplyr::select(-Sample)
    
    Cox.Ref <- coxph(Surv(OS_time, OS_status) ~ .,
                     data = training_data)
    
    coef.Ref  <- coef(Cox.Ref)
    se_coefs.Ref  <- sqrt(diag(vcov(Cox.Ref)))
    p_values.Ref  <- summary(Cox.Ref)$coefficients[, "Pr(>|z|)"]
    concordance.Ref  <- summary(Cox.Ref)$concordance
    logtest.Ref  <- summary(Cox.Ref)$logtest
    wald_test.Ref  <- summary(Cox.Ref)$waldtest
    score_logrank.Ref  <- summary(Cox.Ref)$sctest
    
    Summary.Ref <- list(
      coefs = coef.Ref,
      se_coefs = se_coefs.Ref,
      p_values = p_values.Ref,
      concordance = concordance.Ref,
      logtest.Ref = logtest.Ref,
      wald_test = wald_test.Ref,
      score_logrank = score_logrank.Ref
    )
    
    saveRDS(Summary.Ref, paste0(SummaryDFS.dir,"SummaryOS_RCHOP_TestSample_", Sample_i, ".RDS"))
    
    #SurvCIi$LP <- predict(Cox.Ref, type = "lp", newdata = SurvCIi)
    SurvCIi$RiskScore <- predict(Cox.Ref, type = "risk", newdata = SurvCIi)
    
    risk_thresholds <- quantile(SurvCIi %>% filter(Sample != Sample_i) %>% 
                                  dplyr::select(RiskScore) %>% unlist,
                                probs = c(0, 0.25, 0.5, 0.75, 1))
    
    # Crea una nueva columna para clasificar el riesgo
    SurvCIi$RiskCategory <- cut(
      SurvCIi$RiskScore,
      breaks = risk_thresholds,
      labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
      include.lowest = TRUE
    )
    
    Cox.Ref.Cat <- coxph(Surv(OS_time, OS_status) ~ .,
                         data = SurvCIi %>% filter(Sample != Sample_i) %>% 
                           dplyr::select(RiskCategory, OS_time, OS_status))
    
    MP2_strat <- data.frame(Sample = Sample_i,
                            It_number = i,
                            C.Ref = summary(Cox.Ref)$concordance[1],
                            RiskScore = SurvCIi$RiskScore[SurvCIi$Sample == Sample_i],
                            RiskCat = SurvCIi$RiskCategory[SurvCIi$Sample == Sample_i],
                            C.Ref.Cat = summary(Cox.Ref.Cat)$concordance[1])
    
    if(is.na(MP2_strat$RiskCat)){
      
      MP2_strat <- MP2_strat %>% mutate(RiskCat = case_when(RiskScore > risk_thresholds[5] ~ "High",
                                                            RiskScore < risk_thresholds[1] ~ "Low",
                                                            .default = RiskCat))
      
    }
    
    
    saveRDS(SurvCIi, paste0(SummaryDFS.dir, "SurvCIi_RCHOP_TestSample_", Sample_i, ".RDS"))
    
    MP2_strat <- MP2_strat %>% inner_join(SurvCI, by = c("Sample" = "Sample"))
    
    return(MP2_strat)
    
  }
  
  SummaryDFS.dir <- paste0("NormData/SSNsTests", "/OS_Summaries_RCHOP_Test", i, "/")
  dir.create(SummaryDFS.dir)
  SSNs <- paste0("SSNsTests/SSN_TestSet", i, "/")
  nfiles <- length(list.files(SSNs)[grep("RCHOP", list.files(SSNs))])
  results <- future_lapply(list.files(SSNs)[grep("RCHOP", list.files(SSNs))],
                           FUN = IterationCox, future.seed = TRUE)
  Global_Strata <- do.call(rbind, results)
  
  write.csv(Global_Strata, paste0("NormData/SSNsTests/OS_Global_RCHOP_Strata_Val", i, ".csv"))
  
  cox_model <- coxph(Surv(OS_time, OS_status) ~ RiskScore, 
                     data = Global_Strata)
  
  cox_model.Cat <- coxph(Surv(OS_time, OS_status) ~ RiskCat, 
                         data = Global_Strata)
  
  c.risk <- summary(cox_model)$concordance[1]
  c.cat <- summary(cox_model.Cat)$concordance[1]
  
  if(i == 1){
    C.risk <- c.risk
    C.cat <- c.cat
  } else {
    C.risk <- rbind(C.risk, c.risk)
    C.cat <- rbind(C.cat, c.cat)
  }
  
  #print(summary(cox_model)$concordance)
  
}

future::plan(sequential)

saveRDS(C.cat, file.path(Validation.dir, "C_cat_OS_RCHOP.RDS"))
saveRDS(C.risk, file.path(Validation.dir, "C_risk_OS_RCHOP.RDS"))


#### Evaluate GCHOP model  PFS ####

library(dplyr)
library(tidyr)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(caret)
library(future.apply)
library(timeROC)

SurvCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

SurvCI <- SurvCI %>% filter(Treatment == "GA101")

future::plan(multisession, workers = 5)
for(i in 1:10){
  
  print(i)
  
  IterationCox <- function(file){
    #i = 1
    
    SSN_i <- read.csv(paste0(SSNs,file), row.names = 1)
    SSN_i <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]
    
    Sample_i <- colnames(SSN_i)[ncol(SSN_i)]
    
    original_data_relevant <- SSN_i[, 1:(ncol(SSN_i)-1)]
    row_means <- rowMeans(original_data_relevant)
    row_sds <- apply(original_data_relevant, 1, sd)
    
    SSN_i_scaled <- (SSN_i - row_means) / row_sds
    
    e <- t(SSN_i_scaled)
    e <- e %>% as.data.frame %>% mutate(Sample = rownames(e))
    
    SurvCIi <- SurvCI %>% dplyr::select(Sample, PFS_time, PFS_status) %>%
      inner_join(e, by = c("Sample" = "Sample"))
    
    training_data <- SurvCIi %>% filter(Sample != Sample_i) %>% dplyr::select(-Sample)
    
    Cox.Ref <- coxph(Surv(PFS_time, PFS_status) ~ .,
                     data = training_data)
    
    coef.Ref  <- coef(Cox.Ref)
    se_coefs.Ref  <- sqrt(diag(vcov(Cox.Ref)))
    p_values.Ref  <- summary(Cox.Ref)$coefficients[, "Pr(>|z|)"]
    concordance.Ref  <- summary(Cox.Ref)$concordance
    logtest.Ref  <- summary(Cox.Ref)$logtest
    wald_test.Ref  <- summary(Cox.Ref)$waldtest
    score_logrank.Ref  <- summary(Cox.Ref)$sctest
    
    Summary.Ref <- list(
      coefs = coef.Ref,
      se_coefs = se_coefs.Ref,
      p_values = p_values.Ref,
      concordance = concordance.Ref,
      logtest.Ref = logtest.Ref,
      wald_test = wald_test.Ref,
      score_logrank = score_logrank.Ref
    )
    
    saveRDS(Summary.Ref, paste0(SummaryDFS.dir,"SummaryPFS_GCHOP_TestSample_", Sample_i, ".RDS"))
    
    #SurvCIi$LP <- predict(Cox.Ref, type = "lp", newdata = SurvCIi)
    SurvCIi$RiskScore <- predict(Cox.Ref, type = "risk", newdata = SurvCIi)
    
    risk_thresholds <- quantile(SurvCIi %>% filter(Sample != Sample_i) %>% 
                                  dplyr::select(RiskScore) %>% unlist,
                                probs = c(0, 0.25, 0.5, 0.75, 1))
    
    # Crea una nueva columna para clasificar el riesgo
    SurvCIi$RiskCategory <- cut(
      SurvCIi$RiskScore,
      breaks = risk_thresholds,
      labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
      include.lowest = TRUE
    )
    
    Cox.Ref.Cat <- coxph(Surv(PFS_time, PFS_status) ~ .,
                         data = SurvCIi %>% filter(Sample != Sample_i) %>% 
                           dplyr::select(RiskCategory, PFS_time, PFS_status))
    
    MP2_strat <- data.frame(Sample = Sample_i,
                            It_number = i,
                            C.Ref = summary(Cox.Ref)$concordance[1],
                            RiskScore = SurvCIi$RiskScore[SurvCIi$Sample == Sample_i],
                            RiskCat = SurvCIi$RiskCategory[SurvCIi$Sample == Sample_i],
                            C.Ref.Cat = summary(Cox.Ref.Cat)$concordance[1])
    
    if(is.na(MP2_strat$RiskCat)){
      
      MP2_strat <- MP2_strat %>% mutate(RiskCat = case_when(RiskScore > risk_thresholds[5] ~ "High",
                                                            RiskScore < risk_thresholds[1] ~ "Low",
                                                            .default = RiskCat))
      
    }
    
    
    saveRDS(SurvCIi, paste0(SummaryDFS.dir, "SurvCIi_GCHOP_TestSample_", Sample_i, ".RDS"))
    
    MP2_strat <- MP2_strat %>% inner_join(SurvCI, by = c("Sample" = "Sample"))
    
    return(MP2_strat)
    
  }
  
  SummaryDFS.dir <- paste0("NormData/SSNsTests", "/Summaries_GCHOP_Test", i, "/")
  dir.create(SummaryDFS.dir)
  SSNs <- paste0("NormData/SSNsTests/SSN_TestSet", i, "/")
  nfiles <- length(list.files(SSNs)[grep("GCHOP", list.files(SSNs))])
  results <- future_lapply(list.files(SSNs)[grep("GCHOP", list.files(SSNs))],
                           FUN = IterationCox, future.seed = TRUE)
  Global_Strata <- do.call(rbind, results)
  
  write.csv(Global_Strata, paste0("NormData/SSNsTests/Global_GCHOP_Strata_Val", i, ".csv"))
  
  cox_model <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore, 
                     data = Global_Strata)
  
  cox_model.Cat <- coxph(Surv(PFS_time, PFS_status) ~ RiskCat, 
                         data = Global_Strata)
  
  c.risk <- summary(cox_model)$concordance[1]
  c.cat <- summary(cox_model.Cat)$concordance[1]
  
  if(i == 1){
    C.risk <- c.risk
    C.cat <- c.cat
  } else {
    C.risk <- rbind(C.risk, c.risk)
    C.cat <- rbind(C.cat, c.cat)
  }
  
  #print(summary(cox_model)$concordance)
  
}

future::plan(sequential)

saveRDS(C.cat, file.path(Validation.dir, "C_cat_PFS_GCHOP.RDS"))
saveRDS(C.risk, file.path(Validation.dir, "C_risk_PFS_GCHOP.RDS"))

#### Evaluate GCHOP model  OS ####
library(dplyr)
library(tidyr)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(caret)
library(future.apply)
library(timeROC)

SurvCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

SurvCI <- SurvCI %>% filter(Treatment == "GA101")

future::plan(multisession, workers = 5)
for(i in 1:10){
  
  print(i)
  
  IterationCox <- function(file){
    #i = 1
    
    SSN_i <- read.csv(paste0(SSNs,file), row.names = 1)
    SSN_i <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]
    
    Sample_i <- colnames(SSN_i)[ncol(SSN_i)]
    
    original_data_relevant <- SSN_i[, 1:(ncol(SSN_i)-1)]
    row_means <- rowMeans(original_data_relevant)
    row_sds <- apply(original_data_relevant, 1, sd)
    
    SSN_i_scaled <- (SSN_i - row_means) / row_sds
    
    e <- t(SSN_i_scaled)
    e <- e %>% as.data.frame %>% mutate(Sample = rownames(e))
    
    SurvCIi <- SurvCI %>% dplyr::select(Sample, OS_time, OS_status) %>%
      inner_join(e, by = c("Sample" = "Sample"))
    
    training_data <- SurvCIi %>% filter(Sample != Sample_i) %>% dplyr::select(-Sample)
    
    Cox.Ref <- coxph(Surv(OS_time, OS_status) ~ . ,
                     data = training_data)
    coef.Ref  <- coef(Cox.Ref)
    se_coefs.Ref  <- sqrt(diag(vcov(Cox.Ref)))
    p_values.Ref  <- summary(Cox.Ref)$coefficients[, "Pr(>|z|)"]
    concordance.Ref  <- summary(Cox.Ref)$concordance
    logtest.Ref  <- summary(Cox.Ref)$logtest
    wald_test.Ref  <- summary(Cox.Ref)$waldtest
    score_logrank.Ref  <- summary(Cox.Ref)$sctest
    
    Summary.Ref <- list(
      coefs = coef.Ref,
      se_coefs = se_coefs.Ref,
      p_values = p_values.Ref,
      concordance = concordance.Ref,
      logtest.Ref = logtest.Ref,
      wald_test = wald_test.Ref,
      score_logrank = score_logrank.Ref
    )
    
    saveRDS(Summary.Ref, paste0(SummaryDFS.dir,"SummaryOS_GCHOP_TestSample_", Sample_i, ".RDS"))
    
    #SurvCIi$LP <- predict(Cox.Ref, type = "lp", newdata = SurvCIi)
    SurvCIi$RiskScore <- predict(Cox.Ref, type = "risk", newdata = SurvCIi)
    
    risk_thresholds <- quantile(SurvCIi %>% filter(Sample != Sample_i) %>% 
                                  dplyr::select(RiskScore) %>% unlist,
                                probs = c(0, 0.25, 0.5, 0.75, 1))
    
    # Crea una nueva columna para clasificar el riesgo
    SurvCIi$RiskCategory <- cut(
      SurvCIi$RiskScore,
      breaks = risk_thresholds,
      labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
      include.lowest = TRUE
    )
    
    Cox.Ref.Cat <- coxph(Surv(OS_time, OS_status) ~ .,
                         data = SurvCIi %>% filter(Sample != Sample_i) %>% 
                           dplyr::select(RiskCategory, OS_time, OS_status))
    
    MP2_strat <- data.frame(Sample = Sample_i,
                            It_number = i,
                            C.Ref = summary(Cox.Ref)$concordance[1],
                            RiskScore = SurvCIi$RiskScore[SurvCIi$Sample == Sample_i],
                            RiskCat = SurvCIi$RiskCategory[SurvCIi$Sample == Sample_i],
                            C.Ref.Cat = summary(Cox.Ref.Cat)$concordance[1])
    
    if(is.na(MP2_strat$RiskCat)){
      
      MP2_strat <- MP2_strat %>% mutate(RiskCat = case_when(RiskScore > risk_thresholds[5] ~ "High",
                                                            RiskScore < risk_thresholds[1] ~ "Low",
                                                            .default = RiskCat))
      
    }
    
    saveRDS(SurvCIi, paste0(SummaryDFS.dir, "SurvCIi_GCHOP_TestSample_", Sample_i, ".RDS"))
    
    MP2_strat <- MP2_strat %>% inner_join(SurvCI, by = c("Sample" = "Sample"))
    
    return(MP2_strat)
    
  }
  
  SummaryDFS.dir <- paste0("NormData/SSNsTests", "/Summaries_GCHOP_Test", i, "/")
  dir.create(SummaryDFS.dir)
  SSNs <- paste0("NormData/SSNsTests/SSN_TestSet", i, "/")
  nfiles <- length(list.files(SSNs)[grep("GCHOP", list.files(SSNs))])
  results <- future_lapply(list.files(SSNs)[grep("GCHOP", list.files(SSNs))],
                           FUN = IterationCox, future.seed = TRUE)
  Global_Strata <- do.call(rbind, results)
  
  write.csv(Global_Strata, paste0("NormData/SSNsTests/OS_Global_GCHOP_Strata_Val", i, ".csv"))
  
  cox_model <- coxph(Surv(OS_time, OS_status) ~ RiskScore, 
                     data = Global_Strata)
  
  cox_model.Cat <- coxph(Surv(OS_time, OS_status) ~ RiskCat, 
                         data = Global_Strata)
  
  c.risk <- summary(cox_model)$concordance[1]
  c.cat <- summary(cox_model.Cat)$concordance[1]
  
  if(i == 1){
    C.risk <- c.risk
    C.cat <- c.cat
  } else {
    C.risk <- rbind(C.risk, c.risk)
    C.cat <- rbind(C.cat, c.cat)
  }
  
  #print(summary(cox_model)$concordance)
  
}

future::plan(sequential)

saveRDS(C.cat, file.path(Validation.dir, "C_cat_OS_GCHOP.RDS"))
saveRDS(C.risk, file.path(Validation.dir, "C_risk_OS_GCHOP.RDS"))

#### Get C.index for Train and Test sets (from full data networks) and for IPI and COO ####

#Train C viene del analisis de AIC
#Test viene del analisis de prueba del model con tt
#Vamos a probar el C del conjunto Train y del conjuto Test para IPI y COO

#El plot tendra 6 columnas al final, dos para cada modelo

#Para RCHOP

library(dplyr)
library(tidyr)
library(survival)
library(future.apply)
library(glmnet)
library(MASS)

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

CI.Global <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Global$PFS_status[CI.Global$PFS_time > 1825] <- 0
CI.Global$PFS_time[CI.Global$PFS_time > 1825] <- 1825

CI.RCHOP <- CI.Global %>% filter(Treatment == "Rituximab")
CI.GCHOP <- CI.Global %>% filter(Treatment == "GA101")

#Train
#RCHOP
CheckFinalEdges2 <- readRDS(file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))
for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges2,]
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
  Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
  
  CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
  
  CoxSurvCI <- CoxSurvCI %>% 
    inner_join(CI.RCHOP %>% 
                 dplyr::select(Sample, PFS_status, PFS_time,
                               OS_time, OS_status, COO, IPI), by = "Sample")
  
  Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ . +tt(KIFBP_METTL14), 
                     data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time,
                                                        -IPI, -COO),
                     tt = function(x, t, ...) x * t)
  
  Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                    data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time,
                                                       -IPI, -COO))
  
  #Score for PFS
  CoxSurvCI$RiskScorePFS <- predict(Cox.PFS.1, type = "risk", newdata = CoxSurvCI)
  
  risk_thresholds <- quantile(CoxSurvCI %>% 
                                pull(RiskScorePFS),
                              probs = c(0, 0.25, 0.5, 0.75, 1))
  
  # Crea una nueva columna para clasificar el riesgo
  CoxSurvCI$RiskCatsPFS <- cut(
    CoxSurvCI$RiskScorePFS,
    breaks = risk_thresholds,
    labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
    include.lowest = TRUE
  )
  
  #Score for OS
  CoxSurvCI$RiskScoreOS <- predict(Cox.OS.1, type = "risk", newdata = CoxSurvCI)
  
  risk_thresholds <- quantile(CoxSurvCI %>% 
                                pull(RiskScoreOS),
                              probs = c(0, 0.25, 0.5, 0.75, 1))
  
  # Crea una nueva columna para clasificar el riesgo
  CoxSurvCI$RiskCatsOS <- cut(
    CoxSurvCI$RiskScoreOS,
    breaks = risk_thresholds,
    labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
    include.lowest = TRUE
  )
  
  Cox.PFS.Cat <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, RiskScorePFS))
  
  Cox.OS.Cat <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, RiskScoreOS))
  
  #IPI
  Cox.PFS.IPI <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, IPI))
  
  Cox.OS.IPI <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, IPI))
  
  #COO
  Cox.PFS.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, COO) %>% drop_na(COO)) 
  
  Cox.OS.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, COO) %>% drop_na(COO)) 
  
  if(i == 1){
    summC <- data.frame(C.PFS = summary(Cox.PFS.1)$concordance[1],
                        C.OS = concordance(Cox.OS.1)$concordance, 
                        C.PFS.Cat = concordance(Cox.PFS.Cat)$concordance,
                        C.OS.Cat = concordance(Cox.OS.Cat)$concordance, 
                        C.IPI.OS = concordance(Cox.OS.IPI)$concordance,
                        C.IPI.PFS = concordance(Cox.PFS.IPI)$concordance,
                        C.COO.OS = concordance(Cox.OS.COO)$concordance,
                        C.COO.PFS = concordance(Cox.PFS.COO)$concordance,
                        RP = i)
  } else {
    summC[nrow(summC) + 1, ] <- c(summary(Cox.PFS.1)$concordance[1],
                                  concordance(Cox.OS.1)$concordance, 
                                  concordance(Cox.PFS.Cat)$concordance,
                                  concordance(Cox.OS.Cat)$concordance, 
                                  concordance(Cox.OS.IPI)$concordance,
                                  concordance(Cox.PFS.IPI)$concordance,
                                  concordance(Cox.OS.COO)$concordance,
                                  concordance(Cox.PFS.COO)$concordance,
                                  i)
  }
  
  
}
saveRDS(summC, file.path(Validation.dir, "summC_RCHOP_Train_Cats.RDS"))

#Now GCHOP
CheckFinalEdges2 <- readRDS(file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))
for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges2,]
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
  Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
  
  CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
  
  CoxSurvCI <- CoxSurvCI %>% 
    inner_join(CI.GCHOP %>% 
                 dplyr::select(Sample, PFS_status, PFS_time,
                               OS_time, OS_status, COO, IPI), by = "Sample")
  
  Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ .,
                     data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time,
                                                        -IPI, -COO))
  
  Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                    data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time,
                                                       -IPI, -COO))
  
  #Score for PFS
  CoxSurvCI$RiskScorePFS <- predict(Cox.PFS.1, type = "risk", newdata = CoxSurvCI)
  
  risk_thresholds <- quantile(CoxSurvCI %>% 
                                pull(RiskScorePFS),
                              probs = c(0, 0.25, 0.5, 0.75, 1))
  
  # Crea una nueva columna para clasificar el riesgo
  CoxSurvCI$RiskCatsPFS <- cut(
    CoxSurvCI$RiskScorePFS,
    breaks = risk_thresholds,
    labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
    include.lowest = TRUE
  )
  
  #Score for OS
  CoxSurvCI$RiskScoreOS <- predict(Cox.OS.1, type = "risk", newdata = CoxSurvCI)
  
  risk_thresholds <- quantile(CoxSurvCI %>% 
                                pull(RiskScoreOS),
                              probs = c(0, 0.25, 0.5, 0.75, 1))
  
  # Crea una nueva columna para clasificar el riesgo
  CoxSurvCI$RiskCatsOS <- cut(
    CoxSurvCI$RiskScoreOS,
    breaks = risk_thresholds,
    labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
    include.lowest = TRUE
  )
  
  Cox.PFS.Cat <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, RiskScorePFS))
  
  Cox.OS.Cat <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, RiskScoreOS))
  
  #IPI
  Cox.PFS.IPI <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, IPI))
  
  Cox.OS.IPI <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, IPI))
  
  #IPI + PACO
  Cox.PFS.IPI.PACO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                            data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI, RiskScorePFS))
  
  Cox.OS.IPI.PACO <- coxph(Surv(OS_time, OS_status) ~ ., 
                           data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI, RiskScoreOS))
  
  #COO
  Cox.PFS.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI %>% dplyr::select(PFS_time, PFS_status, COO) %>% drop_na(COO)) 
  
  Cox.OS.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI %>% dplyr::select(OS_time, OS_status, COO) %>% drop_na(COO)) 
  
  if(i == 1){
    summC <- data.frame(C.PFS = summary(Cox.PFS.1)$concordance[1],
                        C.OS = concordance(Cox.OS.1)$concordance, 
                        C.PFS.Cat = concordance(Cox.PFS.Cat)$concordance,
                        C.OS.Cat = concordance(Cox.OS.Cat)$concordance, 
                        C.IPI.OS = concordance(Cox.OS.IPI)$concordance,
                        C.IPI.PFS = concordance(Cox.PFS.IPI)$concordance,
                        C.COO.OS = concordance(Cox.OS.COO)$concordance,
                        C.COO.PFS = concordance(Cox.PFS.COO)$concordance,
                        RP = i)
  } else {
    summC[nrow(summC) + 1, ] <- c(summary(Cox.PFS.1)$concordance[1],
                                  concordance(Cox.OS.1)$concordance, 
                                  concordance(Cox.PFS.Cat)$concordance,
                                  concordance(Cox.OS.Cat)$concordance, 
                                  concordance(Cox.OS.IPI)$concordance,
                                  concordance(Cox.PFS.IPI)$concordance,
                                  concordance(Cox.OS.COO)$concordance,
                                  concordance(Cox.PFS.COO)$concordance,
                                  i)
  }
  
  
}
saveRDS(summC, file.path(Validation.dir, "summC_GCHOP_Train_Cats.RDS"))

#Test
#RCHOP
for(i in 1:10){
  CoxSurvCI.PFS <-  read.csv(paste0("NormData/SSNsTests/Global_RCHOP_Strata_tt_Val", i, ".csv"))
  CoxSurvCI.OS <-  read.csv(paste0("NormData/SSNsTests/OS_Global_RCHOP_Strata_tt_Val", i, ".csv"))
  
  CoxSurvCI.OS <- CoxSurvCI.OS %>% drop_na(COO, IPI)
  CoxSurvCI.PFS <- CoxSurvCI.PFS %>% drop_na(COO, IPI)
  
  #PACO
  Cox.PFS.PACO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                        data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, RiskCat))
  
  Cox.OS.PACO <- coxph(Surv(OS_time, OS_status) ~ ., 
                       data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, RiskCat))
  
  #IPI
  Cox.PFS.IPI <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI))
  
  Cox.OS.IPI <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI))
  
  #COO
  Cox.PFS.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, COO))
  
  Cox.OS.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, COO))
  
  #IPI + PACO
  Cox.PFS.IPI.PACO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                            data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI, RiskScore))
  
  Cox.OS.IPI.PACO <- coxph(Surv(OS_time, OS_status) ~ ., 
                           data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI, RiskScore))
  
  #IPI + COO
  Cox.PFS.IPI.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                           data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI, COO))
  
  Cox.OS.IPI.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                          data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI, COO))
  
  if(i == 1){
    summC.Test <- data.frame(PACO.PFS.Test = concordance(Cox.PFS.PACO)$concordance,
                             PACO.OS.Test = concordance(Cox.OS.PACO)$concordance, 
                             IPI.OS.Test = concordance(Cox.OS.IPI)$concordance,
                             IPI.PFS.Test = concordance(Cox.PFS.IPI)$concordance,
                             COO.OS.Test = concordance(Cox.OS.COO)$concordance,
                             COO.PFS.Test = concordance(Cox.PFS.COO)$concordance,
                             IPI.PACO.OS.Test = concordance(Cox.OS.IPI.PACO)$concordance,
                             IPI.PACO.PFS.Test = concordance(Cox.PFS.IPI.PACO)$concordance,
                             IPI.COO.OS.Test = concordance(Cox.OS.IPI.COO)$concordance,
                             IPI.COO.PFS.Test = concordance(Cox.PFS.IPI.COO)$concordance,
                             AIC.IPI = AIC(Cox.OS.IPI),
                             AIC.PACO = AIC(Cox.OS.PACO),
                             AIC.COO = AIC(Cox.OS.COO),
                             AIC.IPI.PACO = AIC(Cox.OS.IPI.PACO),
                             AIC.IPI.COO = AIC(Cox.OS.IPI.COO),
                             Anova.LTR.IPI_to_IPI.PACO = anova(Cox.OS.IPI, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                             Anova.LTR.PACO_to_IPI.PACO = anova(Cox.OS.PACO, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                             Anova.LTR.IPI_to_IPI.COO = anova(Cox.OS.IPI, Cox.OS.IPI.COO, test = "LRT")$P[2],
                             Anova.PFS.LTR.IPI_to_IPI.PACO = anova(Cox.PFS.IPI, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                             Anova.PFS.LTR.PACO_to_IPI.PACO = anova(Cox.PFS.PACO, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                             Anova.PFS.LTR.IPI_to_IPI.COO = anova(Cox.PFS.IPI, Cox.PFS.IPI.COO, test = "LRT")$P[2],
                             RP = i)
  } else {
    summC.Test[nrow(summC.Test) + 1, ] <- c(concordance(Cox.PFS.PACO)$concordance,
                                            concordance(Cox.OS.PACO)$concordance, 
                                            concordance(Cox.OS.IPI)$concordance,
                                            concordance(Cox.PFS.IPI)$concordance,
                                            concordance(Cox.OS.COO)$concordance,
                                            concordance(Cox.PFS.COO)$concordance,
                                            concordance(Cox.OS.IPI.PACO)$concordance,
                                            concordance(Cox.PFS.IPI.PACO)$concordance,
                                            concordance(Cox.OS.IPI.COO)$concordance,
                                            concordance(Cox.PFS.IPI.COO)$concordance,
                                            AIC(Cox.OS.IPI),
                                            AIC(Cox.OS.PACO),
                                            AIC(Cox.OS.COO),
                                            AIC(Cox.OS.IPI.PACO),
                                            AIC(Cox.OS.IPI.COO),
                                            anova(Cox.OS.IPI, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.OS.PACO, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.OS.IPI, Cox.OS.IPI.COO, test = "LRT")$P[2],
                                            anova(Cox.PFS.IPI, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.PFS.PACO, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.PFS.IPI, Cox.PFS.IPI.COO, test = "LRT")$P[2],
                                            i)
  }
  
}
saveRDS(summC.Test, file.path(Validation.dir, "summC_RCHOP_Test_Cats.RDS"))

#Now GCHOP
CheckFinalEdges2 <- readRDS(file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))
for(i in 1:10){
  CoxSurvCI.PFS <-  read.csv(paste0("NormData/SSNsTests/Global_GCHOP_Strata_tt_Val", i, ".csv"))
  CoxSurvCI.OS <-  read.csv(paste0("NormData/SSNsTests/OS_Global_GCHOP_Strata_tt_Val", i, ".csv"))
  
  CoxSurvCI.OS <- CoxSurvCI.OS %>% drop_na(COO, IPI)
  CoxSurvCI.PFS <- CoxSurvCI.PFS %>% drop_na(COO, IPI)
  
  #PACO
  Cox.PFS.PACO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                        data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, RiskCat))
  
  Cox.OS.PACO <- coxph(Surv(OS_time, OS_status) ~ ., 
                       data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, RiskCat))
  
  #IPI
  Cox.PFS.IPI <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI))
  
  Cox.OS.IPI <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI))
  
  #COO
  Cox.PFS.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, COO))
  
  Cox.OS.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, COO))
  
  #IPI + PACO
  Cox.PFS.IPI.PACO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                            data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI, RiskScore))
  
  Cox.OS.IPI.PACO <- coxph(Surv(OS_time, OS_status) ~ ., 
                           data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI, RiskScore))
  
  #IPI + COO
  Cox.PFS.IPI.COO <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                           data = CoxSurvCI.PFS %>% dplyr::select(PFS_time, PFS_status, IPI, COO))
  
  Cox.OS.IPI.COO <- coxph(Surv(OS_time, OS_status) ~ ., 
                          data = CoxSurvCI.OS %>% dplyr::select(OS_time, OS_status, IPI, COO))
  
  if(i == 1){
    summC.Test <- data.frame(PACO.PFS.Test = concordance(Cox.PFS.PACO)$concordance,
                             PACO.OS.Test = concordance(Cox.OS.PACO)$concordance, 
                             IPI.OS.Test = concordance(Cox.OS.IPI)$concordance,
                             IPI.PFS.Test = concordance(Cox.PFS.IPI)$concordance,
                             COO.OS.Test = concordance(Cox.OS.COO)$concordance,
                             COO.PFS.Test = concordance(Cox.PFS.COO)$concordance,
                             IPI.PACO.OS.Test = concordance(Cox.OS.IPI.PACO)$concordance,
                             IPI.PACO.PFS.Test = concordance(Cox.PFS.IPI.PACO)$concordance,
                             IPI.COO.OS.Test = concordance(Cox.OS.IPI.COO)$concordance,
                             IPI.COO.PFS.Test = concordance(Cox.PFS.IPI.COO)$concordance,
                             AIC.IPI = AIC(Cox.OS.IPI),
                             AIC.PACO = AIC(Cox.OS.PACO),
                             AIC.COO = AIC(Cox.OS.COO),
                             AIC.IPI.PACO = AIC(Cox.OS.IPI.PACO),
                             AIC.IPI.COO = AIC(Cox.OS.IPI.COO),
                             Anova.LTR.IPI_to_IPI.PACO = anova(Cox.OS.IPI, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                             Anova.LTR.PACO_to_IPI.PACO = anova(Cox.OS.PACO, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                             Anova.LTR.IPI_to_IPI.COO = anova(Cox.OS.IPI, Cox.OS.IPI.COO, test = "LRT")$P[2],
                             Anova.PFS.LTR.IPI_to_IPI.PACO = anova(Cox.PFS.IPI, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                             Anova.PFS.LTR.PACO_to_IPI.PACO = anova(Cox.PFS.PACO, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                             Anova.PFS.LTR.IPI_to_IPI.COO = anova(Cox.PFS.IPI, Cox.PFS.IPI.COO, test = "LRT")$P[2],
                             RP = i)
  } else {
    summC.Test[nrow(summC.Test) + 1, ] <- c(concordance(Cox.PFS.PACO)$concordance,
                                            concordance(Cox.OS.PACO)$concordance, 
                                            concordance(Cox.OS.IPI)$concordance,
                                            concordance(Cox.PFS.IPI)$concordance,
                                            concordance(Cox.OS.COO)$concordance,
                                            concordance(Cox.PFS.COO)$concordance,
                                            concordance(Cox.OS.IPI.PACO)$concordance,
                                            concordance(Cox.PFS.IPI.PACO)$concordance,
                                            concordance(Cox.OS.IPI.COO)$concordance,
                                            concordance(Cox.PFS.IPI.COO)$concordance,
                                            AIC(Cox.OS.IPI),
                                            AIC(Cox.OS.PACO),
                                            AIC(Cox.OS.COO),
                                            AIC(Cox.OS.IPI.PACO),
                                            AIC(Cox.OS.IPI.COO),
                                            anova(Cox.OS.IPI, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.OS.PACO, Cox.OS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.OS.IPI, Cox.OS.IPI.COO, test = "LRT")$P[2],
                                            anova(Cox.PFS.IPI, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.PFS.PACO, Cox.PFS.IPI.PACO, test = "LRT")$P[2],
                                            anova(Cox.PFS.IPI, Cox.PFS.IPI.COO, test = "LRT")$P[2],
                                            i)
  }
  
}
saveRDS(summC.Test, file.path(Validation.dir, "summC_GCHOP_Test_Cats.RDS"))

#### RCHOP Cross-cohort evaluation ####

library(tidyr)
library(dplyr)
library(survival)

CI.FULL <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

SSN.RCHOP <- read.csv("SSNs/SSN_RCHOP_AllSamples.csv", row.names = 1)

SSN.RCHOP <- as.data.frame(scale(t(SSN.RCHOP)))

SSN.RCHOP.GEO <- SSN.RCHOP[CI.FULL %>% filter(Source == "GEO",
                                                Treatment == "Rituximab") %>% pull(Sample), ]
SSN.RCHOP.GDC <- SSN.RCHOP[CI.FULL %>% filter(Source == "TCGA") %>% pull(Sample), ]

SSN.RCHOP.GEO <- SSN.RCHOP.GEO %>% mutate(Sample = rownames(SSN.RCHOP.GEO))
SSN.RCHOP.GDC <- SSN.RCHOP.GDC %>% mutate(Sample = rownames(SSN.RCHOP.GDC))

SSN.RCHOP.GEO.OS <- SSN.RCHOP.GEO %>% inner_join(CI.FULL %>% dplyr::select(Sample, OS_time, OS_status),
                                               by = "Sample")
SSN.RCHOP.GDC.OS <- SSN.RCHOP.GDC %>% inner_join(CI.FULL %>% dplyr::select(Sample, OS_time, OS_status),
                                                 by = "Sample")

SSN.RCHOP.GEO.PFS <- SSN.RCHOP.GEO %>% inner_join(CI.FULL %>% dplyr::select(Sample, PFS_time, PFS_status),
                                                 by = "Sample")
SSN.RCHOP.GDC.PFS <- SSN.RCHOP.GDC %>% inner_join(CI.FULL %>% dplyr::select(Sample, PFS_time, PFS_status),
                                                 by = "Sample")

#Train GDC in OS
Cox.RCHOP.GDC.OS.Train <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = SSN.RCHOP.GDC.OS %>% dplyr::select(-Sample))

C.Train.OS.GDC <- summary(Cox.RCHOP.GDC.OS.Train)$concordance[1]

SSN.RCHOP.GEO.OS$RiskScore <- predict(Cox.RCHOP.GDC.OS.Train, 
                                        newdata = SSN.RCHOP.GEO.OS %>% dplyr::select(-Sample),
                                        type = "risk")
Cox.RCHOP.GEO.OS.Test <- coxph(Surv(OS_time, OS_status) ~ RiskScore, 
                          data = SSN.RCHOP.GEO.OS %>% dplyr::select(-Sample))

C.Test.OS.GEO <- summary(Cox.RCHOP.GEO.OS.Test)$concordance[1]

SSN.RCHOP.GEO.OS <- SSN.RCHOP.GEO.OS %>% dplyr::select(-RiskScore)

#Train GEO in OS
Cox.RCHOP.GEO.OS.Train <- coxph(Surv(OS_time, OS_status) ~ ., 
                                data = SSN.RCHOP.GEO.OS %>% dplyr::select(-Sample))

C.Train.OS.GEO <- summary(Cox.RCHOP.GEO.OS.Train)$concordance[1]

SSN.RCHOP.GDC.OS$RiskScore <- predict(Cox.RCHOP.GEO.OS.Train,
                                      newdata = SSN.RCHOP.GDC.OS %>% dplyr::select(-Sample),
                                      type = "risk")
Cox.RCHOP.GDC.OS.Test <- coxph(Surv(OS_time, OS_status) ~ RiskScore, 
                               data = SSN.RCHOP.GDC.OS %>% dplyr::select(-Sample))
C.Test.OS.GDC <- summary(Cox.RCHOP.GDC.OS.Test)$concordance[1]
SSN.RCHOP.GDC.OS <- SSN.RCHOP.GDC.OS %>%  dplyr::select(-RiskScore)

#Train GDC in PFS
Cox.RCHOP.GDC.PFS.Train <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                                data = SSN.RCHOP.GDC.PFS %>% dplyr::select(-Sample))

C.Train.PFS.GDC <- summary(Cox.RCHOP.GDC.PFS.Train)$concordance[1]

SSN.RCHOP.GEO.PFS$RiskScore <- predict(Cox.RCHOP.GDC.PFS.Train, 
                                      newdata = SSN.RCHOP.GEO.PFS %>% dplyr::select(-Sample),
                                      type = "risk")
Cox.RCHOP.GEO.PFS.Test <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore, 
                               data = SSN.RCHOP.GEO.PFS %>% dplyr::select(-Sample))

C.Test.PFS.GEO <- summary(Cox.RCHOP.GEO.PFS.Test)$concordance[1]

SSN.RCHOP.GEO.PFS <- SSN.RCHOP.GEO.PFS %>%  dplyr::select(-RiskScore)

#Train GEO in PFS
Cox.RCHOP.GEO.PFS.Train <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                                data = SSN.RCHOP.GEO.PFS %>% dplyr::select(-Sample))

C.Train.PFS.GEO <- summary(Cox.RCHOP.GEO.PFS.Train)$concordance[1]

SSN.RCHOP.GDC.PFS$RiskScore <- predict(Cox.RCHOP.GEO.PFS.Train,
                                      newdata = SSN.RCHOP.GDC.PFS %>% dplyr::select(-Sample),
                                      type = "risk")
Cox.RCHOP.GDC.PFS.Test <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore, 
                               data = SSN.RCHOP.GDC.PFS %>% dplyr::select(-Sample))
C.Test.PFS.GDC <- summary(Cox.RCHOP.GDC.PFS.Test)$concordance[1]
SSN.RCHOP.GDC.PFS <- SSN.RCHOP.GDC.PFS %>%  dplyr::select(-RiskScore)

CrossValidation <- data.frame(C = c(C.Train.PFS.GDC, C.Test.PFS.GEO, C.Train.PFS.GEO, C.Test.PFS.GDC, 
                                    C.Train.OS.GDC, C.Test.OS.GEO, C.Train.OS.GEO, C.Test.OS.GDC),
                              EndPoint = c("PFS", "PFS", "PFS", "PFS",
                                           "OS", "OS", "OS", "OS"),
                              Type = c("Train", "Test", "Train", "Test",
                                       "Train", "Test", "Train", "Test"),
                              DataSet = c("TCGA", "GEO", "GEO", "TCGA",
                                          "TCGA", "GEO", "GEO", "TCGA"))

write.csv(CrossValidation, file.path(Validation.dir, "RCHOP_CrossCohortEvaluation.csv"), quote = F, row.names = F)

