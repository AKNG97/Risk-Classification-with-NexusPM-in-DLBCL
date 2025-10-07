
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

#### Get Stable Edges ####
library(dplyr)
library(tidyr)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(survminer)
library(MASS)
library(ggplot2)
library(ggvenn)
library(scales)
library(viridis) 
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(cowplot)
library(grid)
library(igraph)
library(clusterProfiler)
library(gprofiler2)
library(exactRankTests)

for(i in 1:10){
  
  print(i)
  
  neti <- read.csv(paste0("GCNs/GCN_RCHOP_TrainSet_NormTrainSet_", i ,".RDS.csv"))
  
  if(i == 1){
    NetSumm <- neti %>% mutate(RP = i) %>%
      dplyr::select(RP, ID, Sp_corr, Source, Target)
  } else {
    NetSumm <- rbind(NetSumm, neti %>% mutate(RP = i) %>% dplyr::select(RP, ID, Sp_corr, Source, Target))
  }
  
}

StableEdges <- NetSumm %>% group_by(ID) %>% summarise(count = n()) %>% 
  filter(count >= 6) %>% arrange(desc(count))

MedianWeights <- NetSumm %>% group_by(ID) %>% summarise(MedianWeigth = median(Sp_corr))

StableEdges <- StableEdges %>% inner_join(MedianWeights)

StableNetwork <- StableEdges %>% group_by(ID) %>%
  mutate(Source = unlist(stringr::str_split(ID, "_"))[1],
         Target = unlist(stringr::str_split(ID, "_"))[2],
         AbsWeight = abs(MedianWeigth))

write.csv(StableNetwork, file.path(FEA.dir, "StableNetwork_RCHOP.csv"), row.names = F, quote = F)

for(i in 1:10){
  
  print(i)
  
  neti <- read.csv(paste0("GCNs/GCN_GCHOP_TrainSet_NormTrainSet_", i ,".RDS.csv"))
  
  if(i == 1){
    NetSumm <- neti %>% mutate(RP = i) %>%
      dplyr::select(RP, ID, Sp_corr, Source, Target)
  } else {
    NetSumm <- rbind(NetSumm, neti %>% mutate(RP = i) %>% dplyr::select(RP, ID, Sp_corr, Source, Target))
  }
  
}

StableEdges <- NetSumm %>% group_by(ID) %>% summarise(count = n()) %>% 
  filter(count >= 6) %>% arrange(desc(count))

MedianWeights <- NetSumm %>% group_by(ID) %>% summarise(MedianWeigth = median(Sp_corr))

StableEdges <- StableEdges %>% inner_join(MedianWeights)

StableNetwork <- StableEdges %>% group_by(ID) %>%
  mutate(Source = unlist(stringr::str_split(ID, "_"))[1],
         Target = unlist(stringr::str_split(ID, "_"))[2],
         AbsWeight = abs(MedianWeigth))

write.csv(StableNetwork, file.path(FEA.dir, "StableNetwork_GCHOP.csv"), row.names = F, quote = F)

#### Get Stable Networks LIONESS weights RCHOP ####

library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

stableEdges <- read.csv(file.path(FEA.dir, "StableNetwork_RCHOP.csv"))

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
CI.Full <- CI.Full %>% mutate(Source = case_when(Source == "GOYA" ~ "GOYA",
                                                 .default = "GDC"),
                              primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist

GCN <- data.frame(Edge = stableEdges$ID, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

TARGET_NBM_ALALi_norm <- readRDS("NormData/DLBCL_NormCounts.RDS")

TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% c(RCHOP.samples, NormalSamples)]

GCN_i <- GCN %>% group_by(Edge) %>%
  dplyr::mutate(Rho = cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,],
                          TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,],
                          method = "spearman"))

write.csv(GCN_i, file.path(FEA.dir, "GCN_RCHOP_EstableEdges.csv"), quote = F, row.names = F)

RCHOP.Cancer_GCN <- GCN_i

RCHOP.Cancer_ExprM <- TARGET_NBM_ALALi_norm

Lioness_PerSource.R <- function(j) {
  
  Source <- RCHOP.Cancer_GCN[j,"Source"] %>% pull
  Target <- RCHOP.Cancer_GCN[j,"Target"] %>% pull
  
  ID <- RCHOP.Cancer_GCN[j,"Edge"] %>% pull
  
  alpha <- RCHOP.Cancer_GCN[j,"Rho"] %>% pull
  
  n_samples <- ncol(RCHOP.Cancer_ExprM)
  
  SSN <- numeric(n_samples)
  
  for(x in 1:n_samples){
    
    #print(x)
    
    j_cor <- Hmisc::rcorr(RCHOP.Cancer_ExprM[rownames(RCHOP.Cancer_ExprM) == Source,-x],
                          RCHOP.Cancer_ExprM[rownames(RCHOP.Cancer_ExprM) == Target,-x], 
                          type = "spearman")$r[1,2]
    
    SSN[x] <- round(n_samples*(alpha - j_cor) + j_cor, 4)
    
  }
  
  names(SSN) <- colnames(RCHOP.Cancer_ExprM)
  SSN <- as.data.frame(t(SSN))
  rownames(SSN) <- ID
  SSN <- SSN %>% mutate(Edge_ID = rownames(SSN))
  
  return(SSN)
}

future::plan(multisession, workers = 60)
SSN.R.Res <- future_lapply(1:nrow(RCHOP.Cancer_GCN), FUN = Lioness_PerSource.R,
                           future.seed = TRUE)
SSN.R <- bind_rows(SSN.R.Res)
write.csv(SSN.R %>% dplyr::select(-Edge_ID), file.path(FEA.dir, "SSN_StableEdges_RCHOP.csv"), 
          quote = FALSE, row.names = T)
plan(sequential)

#### Get Stable Networks LIONESS weights GCHOP ####

library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

stableEdges <- read.csv(file.path(FEA.dir, "StableNetwork_GCHOP.csv"))

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
CI.Full <- CI.Full %>% mutate(Source = case_when(Source == "GOYA" ~ "GOYA",
                                                 .default = "GDC"),
                              primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

GCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "GA101") %>%
  dplyr::select(Sample) %>% unlist

GCN <- data.frame(Edge = stableEdges$ID, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

TARGET_NBM_ALALi_norm <- readRDS("NormData/DLBCL_NormCounts.RDS")

TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% c(GCHOP.samples, NormalSamples)]

GCN_i <- GCN %>% group_by(Edge) %>%
  dplyr::mutate(Rho = cor(TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Source,],
                          TARGET_NBM_ALALi_norm[rownames(TARGET_NBM_ALALi_norm) == Target,],
                          method = "spearman"))

write.csv(GCN_i, file.path(FEA.dir,"GCN_GCHOP_EstableEdges.csv"), quote = F, row.names = F)

GCHOP.Cancer_GCN <- GCN_i

GCHOP.Cancer_ExprM <- TARGET_NBM_ALALi_norm

Lioness_PerSource.R <- function(j) {
  
  Source <- GCHOP.Cancer_GCN[j,"Source"] %>% pull
  Target <- GCHOP.Cancer_GCN[j,"Target"] %>% pull
  
  ID <- GCHOP.Cancer_GCN[j,"Edge"] %>% pull
  
  alpha <- GCHOP.Cancer_GCN[j,"Rho"] %>% pull
  
  n_samples <- ncol(GCHOP.Cancer_ExprM)
  
  SSN <- numeric(n_samples)
  
  for(x in 1:n_samples){
    
    #print(x)
    
    j_cor <- Hmisc::rcorr(GCHOP.Cancer_ExprM[rownames(GCHOP.Cancer_ExprM) == Source,-x],
                          GCHOP.Cancer_ExprM[rownames(GCHOP.Cancer_ExprM) == Target,-x], 
                          type = "spearman")$r[1,2]
    
    SSN[x] <- round(n_samples*(alpha - j_cor) + j_cor, 4)
    
  }
  
  names(SSN) <- colnames(GCHOP.Cancer_ExprM)
  SSN <- as.data.frame(t(SSN))
  rownames(SSN) <- ID
  SSN <- SSN %>% mutate(Edge_ID = rownames(SSN))
  
  return(SSN)
}

future::plan(multisession, workers = 60)
SSN.R.Res <- future_lapply(1:nrow(GCHOP.Cancer_GCN), FUN = Lioness_PerSource.R,
                           future.seed = TRUE)
SSN.R <- bind_rows(SSN.R.Res)
write.csv(SSN.R %>% dplyr::select(-Edge_ID), file.path(FEA.dir,"SSN_StableEdges_GCHOP.csv"), 
          quote = FALSE, row.names = T)
plan(sequential)

##### Run DiffCoExpr ####

library(dplyr)
library(tidyr)
library(future.apply)
library(Hmisc)
library(rstatix)
library(survival)
library(glmnet)
library(MASS)
library(exactRankTests)

CI.FULL <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

SSN.RCHOP <- read.csv("SSNs/SSN_RCHOP_AllSamples.csv", row.names = 1)
SSN.GCHOP <- read.csv("SSNs/SSN_GCHOP_AllSamples.csv", row.names = 1)

SSN.RCHOP <- SSN.RCHOP[,colnames(SSN.RCHOP) %in% CI.FULL$Sample]
SSN.GCHOP <- SSN.GCHOP[,colnames(SSN.GCHOP) %in% CI.FULL$Sample]

SSN.RCHOP <- as.data.frame(scale(t(SSN.RCHOP)))
SSN.GCHOP <- as.data.frame(scale(t(SSN.GCHOP)))

SSN.RCHOP <- SSN.RCHOP %>% mutate(Sample = rownames(SSN.RCHOP))
SSN.GCHOP <- SSN.GCHOP %>% mutate(Sample = rownames(SSN.GCHOP))

CoxSurvCI.RCHOP.OS <- SSN.RCHOP %>% inner_join(CI.FULL %>% dplyr::select(Sample, OS_time, OS_status),
                                               by = "Sample")
CoxSurvCI.GCHOP.OS <- SSN.GCHOP %>% inner_join(CI.FULL %>% dplyr::select(Sample, OS_time, OS_status),
                                               by = "Sample")

Cox.RCHOP.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.RCHOP.OS %>% dplyr::select(-Sample))

Cox.GCHOP.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = CoxSurvCI.GCHOP.OS %>% dplyr::select(-Sample))

CoxSurvCI.RCHOP.OS$RiskScore <- predict(Cox.RCHOP.OS, 
                                        newdata = CoxSurvCI.RCHOP.OS %>% dplyr::select(-Sample),
                                        type = "risk")

CoxSurvCI.GCHOP.OS$RiskScore <- predict(Cox.GCHOP.OS, 
                                        newdata = CoxSurvCI.GCHOP.OS %>% dplyr::select(-Sample),
                                        type = "risk")

RCHOP_thresholds <- quantile(CoxSurvCI.RCHOP.OS %>% 
                               pull(RiskScore),
                             probs = c(0, 0.25, 0.5, 0.75, 1))

GCHOP_thresholds <- quantile(CoxSurvCI.GCHOP.OS %>% 
                               pull(RiskScore),
                             probs = c(0, 0.25, 0.5, 0.75, 1))

CoxSurvCI.RCHOP.OS$RiskCategory <- cut(
  CoxSurvCI.RCHOP.OS$RiskScore,
  breaks = RCHOP_thresholds,
  labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
  include.lowest = TRUE
)

CoxSurvCI.GCHOP.OS$RiskCategory <- cut(
  CoxSurvCI.GCHOP.OS$RiskScore,
  breaks = GCHOP_thresholds,
  labels = c("Low", "Low-Intermediate", "High-Intermediate", "High"),
  include.lowest = TRUE
)

CoxSurvCI.RCHOP.OS <- CoxSurvCI.RCHOP.OS %>% inner_join(CI.FULL %>% dplyr::select(-OS_status, -OS_time))
CoxSurvCI.GCHOP.OS <- CoxSurvCI.GCHOP.OS %>% inner_join(CI.FULL %>% dplyr::select(-OS_status, -OS_time))

#### DiffCoExpr GCHOP

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

GCHOP.Stable.SSN <- read.csv(file.path(FEA.dir,"SSN_StableEdges_GCHOP.csv"), row.names = 1)
GCHOP.Stable.SSN <- GCHOP.Stable.SSN[,!(colnames(GCHOP.Stable.SSN) %in% NormalSamples)]
GCHOP.Stable.SSN <- as.data.frame(scale(t(GCHOP.Stable.SSN)))
GCHOP.Stable.SSN <- GCHOP.Stable.SSN %>% mutate(Sample = rownames(GCHOP.Stable.SSN))

GCHOP.Stable.SSN.pivot <- GCHOP.Stable.SSN %>% pivot_longer(!Sample, names_to = "Edge", values_to = "ScaleWeight")

GCHOP.Stable.SSN.pivot <- GCHOP.Stable.SSN.pivot %>% inner_join(CoxSurvCI.GCHOP.OS %>% 
                                                                  dplyr::select(Sample, RiskCategory, IPI, Gender, COO, Treatment, Source),
                                                                by = "Sample")

WeightsDiff.GCHOP <- GCHOP.Stable.SSN.pivot %>% group_by(Edge) %>%
  summarise(LvsH = median(ScaleWeight[RiskCategory == "Low"]) - median(ScaleWeight[RiskCategory == "High"]))

Wilcox.GCHOP <- GCHOP.Stable.SSN.pivot %>% 
  group_by(Edge) %>% 
  summarise(pval.W = wilcox.exact(ScaleWeight[RiskCategory == "High"],
                                  ScaleWeight[RiskCategory == "Low"])$p.value) %>%
  arrange(pval.W) %>%
  mutate(padj = p.adjust(pval.W, method = "fdr", n = ncol(GCHOP.Stable.SSN))) %>% filter()

GCHOP.DiffCoExprSummary <- Wilcox.GCHOP %>% inner_join(WeightsDiff.GCHOP, by = "Edge")

saveRDS(GCHOP.DiffCoExprSummary, file.path(FEA.dir, "GCHOP.DiffCoExprSummary.RDS"))

#### DiffCoExpr RCHOP

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

RCHOP.Stable.SSN <- read.csv(file.path(FEA.dir, "SSN_StableEdges_RCHOP.csv"), row.names = 1)
RCHOP.Stable.SSN <- RCHOP.Stable.SSN[,!(colnames(RCHOP.Stable.SSN) %in% NormalSamples)]
RCHOP.Stable.SSN <- as.data.frame(scale(t(RCHOP.Stable.SSN)))
RCHOP.Stable.SSN <- RCHOP.Stable.SSN %>% mutate(Sample = rownames(RCHOP.Stable.SSN))

RCHOP.Stable.SSN.pivot <- RCHOP.Stable.SSN %>% pivot_longer(!Sample, names_to = "Edge", values_to = "ScaleWeight")

RCHOP.Stable.SSN.pivot.GDC <- RCHOP.Stable.SSN.pivot %>% inner_join(CoxSurvCI.RCHOP.OS %>% filter(Source == "GDC") %>% 
                                                                       dplyr::select(Sample, RiskCategory, IPI, Gender, COO, Treatment, Source),
                                                                     by = "Sample")
RCHOP.Stable.SSN.pivot.GEO <- RCHOP.Stable.SSN.pivot %>% inner_join(CoxSurvCI.RCHOP.OS %>% filter(Source == "GEO") %>% 
                                                                      dplyr::select(Sample, RiskCategory, IPI, Gender, COO, Treatment, Source),
                                                                    by = "Sample")
RCHOP.Stable.SSN.pivot <- RCHOP.Stable.SSN.pivot %>% inner_join(CoxSurvCI.RCHOP.OS %>% 
                                                                  dplyr::select(Sample, RiskCategory, IPI, Gender, COO, Treatment, Source),
                                                                by = "Sample")

WeightsDiff.RCHOP <- RCHOP.Stable.SSN.pivot %>% group_by(Edge) %>%
  summarise(LvsH = median(ScaleWeight[RiskCategory == "Low"]) - median(ScaleWeight[RiskCategory == "High"]))
Wilcox.RCHOP <- RCHOP.Stable.SSN.pivot %>% 
  group_by(Edge) %>% 
  summarise(pval.W = wilcox.exact(ScaleWeight[RiskCategory == "High"],
                                  ScaleWeight[RiskCategory == "Low"])$p.value) %>%
  arrange(pval.W) %>%
  mutate(padj = p.adjust(pval.W, method = "fdr", n = ncol(RCHOP.Stable.SSN))) 

RCHOP.DiffCoExprSummary <- Wilcox.RCHOP %>% inner_join(WeightsDiff.RCHOP, by = "Edge")

saveRDS(RCHOP.DiffCoExprSummary, file.path(FEA.dir, "RCHOP_DiffCoExprSummary.RDS"))

#Calculate Difference
WeightsDiff.RCHOP.GDC <- RCHOP.Stable.SSN.pivot.GDC %>% group_by(Edge) %>%
  summarise(LvsH = median(ScaleWeight[RiskCategory == "Low"]) - median(ScaleWeight[RiskCategory == "High"]))
WeightsDiff.RCHOP.GEO <- RCHOP.Stable.SSN.pivot.GEO %>% group_by(Edge) %>%
  summarise(LvsH = median(ScaleWeight[RiskCategory == "Low"]) - median(ScaleWeight[RiskCategory == "High"]))

#Wilcoxon test
Wilcox.RCHOP.GDC <- RCHOP.Stable.SSN.pivot.GDC %>% 
  group_by(Edge) %>% 
  summarise(pval.W = wilcox.exact(ScaleWeight[RiskCategory == "High"],
                                  ScaleWeight[RiskCategory == "Low"])$p.value) %>%
  arrange(pval.W) %>%
  mutate(padj = p.adjust(pval.W, method = "fdr", n = ncol(RCHOP.Stable.SSN))) 

RCHOP.DiffCoExprSummary.GDC <- Wilcox.RCHOP.GDC %>% 
  inner_join(WeightsDiff.RCHOP.GDC, by = "Edge")

Wilcox.RCHOP.GEO <- RCHOP.Stable.SSN.pivot.GEO %>% 
  group_by(Edge) %>% 
  summarise(pval.W = wilcox.exact(ScaleWeight[RiskCategory == "High"],
                                  ScaleWeight[RiskCategory == "Low"])$p.value) %>%
  arrange(pval.W) %>%
  mutate(padj = p.adjust(pval.W, method = "fdr", n = ncol(RCHOP.Stable.SSN))) 

RCHOP.DiffCoExprSummary.GEO <- Wilcox.RCHOP.GEO %>% inner_join(WeightsDiff.RCHOP.GEO, by = "Edge")

#Save
saveRDS(RCHOP.DiffCoExprSummary.GDC, file.path(FEA.dir, "RCHOP_DiffCoExprSummary_GDC.RDS"))
saveRDS(RCHOP.DiffCoExprSummary.GEO, file.path(FEA.dir, "RCHOP_DiffCoExprSummary_GEO.RDS"))

#### RCHOP communities FEA ####

library(dplyr)
library(tidyr)
library(gprofiler2)
library(exactRankTests)

SSN.RCHOP <- read.csv("SSNs/SSN_RCHOP_AllSamples.csv", row.names = 1)

RCHOPgenes <- rownames(SSN.RCHOP)
RCHOPgenes <- stringr::str_split(RCHOPgenes, "_") %>% unlist

bg_genes <- readRDS(file.path(Raw.dir, "DLBCL_RawCounts.RDS"))
bg_genes <- rownames(bg_genes)

RCHOP.communities <- readxl::read_xlsx(file.path(FEA.dir, "StableCommunities_RCHOP.xlsx"))

RelevantCommunitiesRCHOP <- c()

for(gene in RCHOPgenes){
  
  indx <- grep(gene, RCHOP.communities$Members_list)
  
  RelevantCommunitiesRCHOP <- c(RelevantCommunitiesRCHOP, indx)
  
}

RelevantCommunitiesRCHOP <- unique(RelevantCommunitiesRCHOP)

RCHOP.communities <- RCHOP.communities[RelevantCommunitiesRCHOP, ]

for(i in 1:nrow(RCHOP.communities)){
  
  print(i)
  
  goi <- NULL
  
  goi <- gost(query = RCHOP.communities$Members_list[i] %>% stringr::str_split(" ") %>% unlist, 
              organism = "hsapiens", ordered_query = F, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              custom_bg = bg_genes, domain_scope = "custom", 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)
  
  if(!is.null(goi)){
    if(!exists("SummFEA")){
      SummFEA <-  tibble(Community = RCHOP.communities$Community[i],
                         GOBP = goi$result$term_name,
                         pval = goi$result$p_value)
    } else {
      SummFEAi <- tibble(Community = RCHOP.communities$Community[i],
                         GOBP = goi$result$term_name,
                         pval = goi$result$p_value)
      
      SummFEA <- rbind(SummFEA, SummFEAi)
    }
    
  }
  
}

write.csv(SummFEA, file.path(FEA.dir,"SummFEA_RCHOP.csv"), quote = F)

for(gene in RCHOPgenes){
  
  indx <- NULL
  
  indx <- grep(gene, RCHOP.communities$Members_list)
  
  if(!is.null(indx)){
    if(!exists("SummIndx")){
      SummIndx <-  tibble(Community = RCHOP.communities$Community[indx],
                          Gene = gene) 
    } else {
      SummIndxi <-  tibble(Community = RCHOP.communities$Community[indx],
                           Gene = gene)
      SummIndx <-  rbind(SummIndx, SummIndxi)
    }
  }
  
}

SummIndx <- SummIndx %>% inner_join(SummFEA, by = "Community")

biomarkers <- stringr::str_split(rownames(SSN.RCHOP), "_") %>% unlist
sources <- biomarkers[c(T,F)]
targets <- biomarkers[c(F,T)]

#Prep for Cytoscape viz

RCHOP_FEA_Network <- tibble(Gene = sources,
                            GOBP = targets,
                            pval = 1, 
                            Type = "CoExpr")

RCHOP_FEA_Network <- rbind(RCHOP_FEA_Network, SummIndx %>% filter(pval < 1e-8) %>%
                             dplyr::select(Gene, GOBP, pval) %>% mutate(Type = "Enrichment"))

write.table(RCHOP_FEA_Network, file.path(FEA.dir, "SummFEA_RCHOP_net_pvalFilter.tsv"), quote = F, sep = "\t", row.names = F)

#### GCHOP communities FEA ####

library(dplyr)
library(tidyr)
library(gprofiler2)
library(exactRankTests)

SSN.GCHOP <- read.csv("SSNs/SSN_GCHOP_AllSamples.csv", row.names = 1)

GCHOPgenes <- rownames(SSN.GCHOP)
GCHOPgenes <- stringr::str_split(GCHOPgenes, "_") %>% unlist

bg_genes <- readRDS(file.path(Raw.dir, "DLBCL_RawCounts.RDS"))
bg_genes <- rownames(bg_genes)

GCHOP.communities <- readxl::read_xlsx(file.path(FEA.dir, "StableCommunities_GCHOP.xlsx"))

RelevantCommunitiesGCHOP <- c()

for(gene in GCHOPgenes){
  
  indx <- grep(gene, GCHOP.communities$Members_list)
  
  RelevantCommunitiesGCHOP <- c(RelevantCommunitiesGCHOP, indx)
  
}

RelevantCommunitiesGCHOP <- unique(RelevantCommunitiesGCHOP)

GCHOP.communities <- GCHOP.communities[RelevantCommunitiesGCHOP, ]

for(i in 1:nrow(GCHOP.communities)){
  
  print(i)
  
  goi <- NULL
  
  goi <- gost(query = GCHOP.communities$Members_list[i] %>% stringr::str_split(" ") %>% unlist, 
              organism = "hsapiens", ordered_query = F, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              custom_bg = bg_genes, domain_scope = "custom", 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)
  
  if(!is.null(goi)){
    if(!exists("SummFEA")){
      SummFEA <-  tibble(Community = GCHOP.communities$Community[i],
                         GOBP = goi$result$term_name,
                         pval = goi$result$p_value)
    } else {
      SummFEAi <- tibble(Community = GCHOP.communities$Community[i],
                         GOBP = goi$result$term_name,
                         pval = goi$result$p_value)
      
      SummFEA <- rbind(SummFEA, SummFEAi)
    }
    
  }
  
}

write.csv(SummFEA, file.path(FEA.dir, "SummFEA_GCHOP.csv"), quote = F)

for(gene in GCHOPgenes){
  
  indx <- NULL
  
  indx <- grep(gene, GCHOP.communities$Members_list)
  
  if(!is.null(indx)){
    if(!exists("SummIndx")){
      SummIndx <-  tibble(Community = GCHOP.communities$Community[indx],
                          Gene = gene) 
    } else {
      SummIndxi <-  tibble(Community = GCHOP.communities$Community[indx],
                           Gene = gene)
      SummIndx <-  rbind(SummIndx, SummIndxi)
    }
  }
  
}

head(SummIndx)

SummIndx <- SummIndx %>% inner_join(SummFEA, by = "Community")

biomarkers <- stringr::str_split(rownames(SSN.GCHOP), "_") %>% unlist
sources <- biomarkers[c(T,F)]
targets <- biomarkers[c(F,T)]

#Prep for Cytoscape viz

GCHOP_FEA_Network <- tibble(Gene = sources,
                            GOBP = targets,
                            pval = 1, 
                            Type = "CoExpr")

GCHOP_FEA_Network <- rbind(GCHOP_FEA_Network, SummIndx %>% filter(pval < 1e-8) %>%
                             dplyr::select(Gene, GOBP, pval) %>% mutate(Type = "Enrichment"))

write.table(GCHOP_FEA_Network, file.path(FEA.dir, "SummFEA_GCHOP_net_pvalFilter.tsv"), quote = F, sep = "\t", row.names = F)

#### Let's check Cis proportion in Stable networks ####

library(tidyr)
library(dplyr)

annot <- readRDS(file.path(Raw.dir, "annot_GDC.RDS"))

SSN.RCHOP <- read.csv("SSNs/SSN_RCHOP_AllSamples.csv", row.names = 1)

RCHOPgenes <- rownames(SSN.RCHOP)
RCHOPgenes <- stringr::str_split(RCHOPgenes, "_") %>% unlist

RCHOP_Network <- read.csv(file.path(FEA.dir, "StableNetwork_RCHOP.csv")) %>% tibble

RCHOP_Network <- RCHOP_Network %>% 
  inner_join(annot %>% dplyr::select(HGNC_symbol, Chr), by = c("Source" = "HGNC_symbol")) %>%
  inner_join(annot %>% dplyr::select(HGNC_symbol, Chr), by = c("Target" = "HGNC_symbol"), suffix = c("_Source", "_Target")) %>% 
  mutate(EdgeType = case_when(Chr_Source == Chr_Target ~ "Cis", .default = "Trans"))

RCHOP.communities <- readxl::read_xlsx(file.path(FEA.dir, "StableCommunities_RCHOP.xlsx"))

for(i in 1:nrow(RCHOP.communities)){
  
  print(i)
  
  comm_members <- RCHOP.communities$Members_list[i] %>% stringr::str_split(" ") %>% unlist
  
  RCHOP_Networki <- RCHOP_Network %>% filter(Source %in% comm_members, Target %in% comm_members)
  
  Proportionsi <- RCHOP_Networki %>% group_by(EdgeType) %>% summarise(Proportion = n()/nrow(RCHOP_Networki))
  
  cis.p <- Proportionsi %>% filter(EdgeType == "Cis") %>% pull(Proportion)
  trans.p <- Proportionsi %>% filter(EdgeType == "Trans") %>% pull(Proportion)
  
  DominantChr <- RCHOP_Networki %>% filter(EdgeType == "Cis") %>% group_by(Chr_Target) %>% 
    summarise(Proportion = n()/nrow(RCHOP_Networki)) %>% slice_max(Proportion) %>% pull(Chr_Target)
  
  DominantChr.p <- RCHOP_Networki %>% filter(EdgeType == "Cis") %>% group_by(Chr_Target) %>% 
    summarise(Proportion = n()/nrow(RCHOP_Networki)) %>% slice_max(Proportion) %>% pull(Proportion)
  
  if(length(cis.p) == 0){
    cis.p <- 0
    
    DominantChr <- NA
    
    DominantChr.p <- 0
    
  }
  
  if(length(trans.p) == 0){
    trans.p <- 0
  }
  
  if(trans.p >= unique(DominantChr.p)){
    SummEdgeTypei <- tibble(Community = RCHOP.communities$Community[i],
                            Cis.p = cis.p,
                            Trans.p = trans.p,
                            CommunityType = "Trans",
                            DominantChr = NA,
                            DominantChr.p = NA)
  } else {
    
    SummEdgeTypei <- tibble(Community = RCHOP.communities$Community[i],
                            Cis.p = cis.p,
                            Trans.p = trans.p,
                            CommunityType = "Cis",
                            DominantChr = DominantChr,
                            DominantChr.p = DominantChr.p)
    
  }
  
  if(i == 1){
    SummEdgeType.RCHOP <- SummEdgeTypei
  } else {
    SummEdgeType.RCHOP <- rbind(SummEdgeType.RCHOP, SummEdgeTypei)
  }
  
}

Gene.annot <- tibble(Community = RCHOPgenes,
                     Cis.p = NA,
                     Trans.p = NA,
                     CommunityType = "Gene",
                     DominantChr = NA,
                     DominantChr.p = NA)

Gene.annot <- Gene.annot %>% inner_join(annot %>% filter(HGNC_symbol %in% RCHOPgenes) %>% dplyr::select(Chr, HGNC_symbol),
                                        by = c("Community" = "HGNC_symbol"))
Gene.annot$DominantChr <- Gene.annot$Chr

Gene.annot <- Gene.annot %>% dplyr::select(-Chr)

NodeAnnot.RCHOP <- rbind(Gene.annot, SummEdgeType.RCHOP)

for(gene in RCHOPgenes){
  
  indx <- NULL
  
  indx <- grep(gene, RCHOP.communities$Members_list)
  
  if(!is.null(indx)){
    if(!exists("SummIndx")){
      SummIndx <-  tibble(Community = RCHOP.communities$Community[indx],
                          Gene = gene) 
    } else {
      SummIndxi <-  tibble(Community = RCHOP.communities$Community[indx],
                           Gene = gene)
      SummIndx <-  rbind(SummIndx, SummIndxi)
    }
  }
  
}

RCHOP_EdgeType_Network <- SummIndx %>% mutate(Weight = 1)

biomarkers <- stringr::str_split(rownames(SSN.RCHOP), "_") %>% unlist
sources <- biomarkers[c(T,F)]
targets <- biomarkers[c(F,T)]

RCHOP_Gene_Network <- tibble(Community = sources,
                             Gene = targets,
                             Weight = 1)

RCHOP_EdgeType_Network <- rbind(RCHOP_EdgeType_Network, RCHOP_Gene_Network)

write.table(RCHOP_EdgeType_Network, file.path(FEA.dir, "RCHOP_EdgeType_Network.tsv"), quote = F, sep = "\t", row.names = F)
write.table(NodeAnnot.RCHOP, file.path(FEA.dir, "EdgeType_NodeAnnotRCHOP.tsv"), quote = F, sep = "\t", row.names = F)

#### Let's check Cis proportion in Stable networks GCHOP ####

library(tidyr)
library(dplyr)

annot <- readRDS(file.path(Raw.dir, "annot_GDC.RDS"))

SSN.GCHOP <- read.csv("SSNs/SSN_GCHOP_AllSamples.csv", row.names = 1)

GCHOPgenes <- rownames(SSN.GCHOP)
GCHOPgenes <- stringr::str_split(GCHOPgenes, "_") %>% unlist

GCHOP_Network <- read.csv(file.path(FEA.dir, "StableNetwork_GCHOP.csv")) %>% tibble

GCHOP_Network <- GCHOP_Network %>% 
  inner_join(annot %>% dplyr::select(HGNC_symbol, Chr), by = c("Source" = "HGNC_symbol")) %>%
  inner_join(annot %>% dplyr::select(HGNC_symbol, Chr), by = c("Target" = "HGNC_symbol"), suffix = c("_Source", "_Target")) %>% 
  mutate(EdgeType = case_when(Chr_Source == Chr_Target ~ "Cis", .default = "Trans"))

GCHOP.communities <- readxl::read_xlsx(file.path(FEA.dir, "StableCommunities_GCHOP.xlsx"))

for(i in 1:nrow(GCHOP.communities)){
  
  print(i)
  
  comm_members <- GCHOP.communities$Members_list[i] %>% stringr::str_split(" ") %>% unlist
  
  GCHOP_Networki <- GCHOP_Network %>% filter(Source %in% comm_members, Target %in% comm_members)
  
  Proportionsi <- GCHOP_Networki %>% group_by(EdgeType) %>% summarise(Proportion = n()/nrow(GCHOP_Networki))
  
  cis.p <- Proportionsi %>% filter(EdgeType == "Cis") %>% pull(Proportion)
  trans.p <- Proportionsi %>% filter(EdgeType == "Trans") %>% pull(Proportion)
  
  DominantChr <- GCHOP_Networki %>% filter(EdgeType == "Cis") %>% group_by(Chr_Target) %>% 
    summarise(Proportion = n()/nrow(GCHOP_Networki)) %>% slice_max(Proportion) %>% pull(Chr_Target)
  
  DominantChr.p <- GCHOP_Networki %>% filter(EdgeType == "Cis") %>% group_by(Chr_Target) %>% 
    summarise(Proportion = n()/nrow(GCHOP_Networki)) %>% slice_max(Proportion) %>% pull(Proportion)
  
  if(length(cis.p) == 0){
    cis.p <- 0
    
    DominantChr <- NA
    
    DominantChr.p <- 0
    
  }
  
  if(length(trans.p) == 0){
    trans.p <- 0
  }
  
  if(trans.p >= unique(DominantChr.p)){
    SummEdgeTypei <- tibble(Community = GCHOP.communities$Community[i],
                            Cis.p = cis.p,
                            Trans.p = trans.p,
                            CommunityType = "Trans",
                            DominantChr = NA,
                            DominantChr.p = NA)
  } else {
    
    SummEdgeTypei <- tibble(Community = GCHOP.communities$Community[i],
                            Cis.p = cis.p,
                            Trans.p = trans.p,
                            CommunityType = "Cis",
                            DominantChr = DominantChr,
                            DominantChr.p = DominantChr.p)
    
  }
  
  if(i == 1){
    SummEdgeType.GCHOP <- SummEdgeTypei
  } else {
    SummEdgeType.GCHOP <- rbind(SummEdgeType.GCHOP, SummEdgeTypei)
  }
  
}

Gene.annot <- tibble(Community = GCHOPgenes,
                     Cis.p = NA,
                     Trans.p = NA,
                     CommunityType = "Gene",
                     DominantChr = NA,
                     DominantChr.p = NA)

Gene.annot <- Gene.annot %>% inner_join(annot %>% filter(HGNC_symbol %in% GCHOPgenes) %>% dplyr::select(Chr, HGNC_symbol),
                                        by = c("Community" = "HGNC_symbol"))
Gene.annot$DominantChr <- Gene.annot$Chr

Gene.annot <- Gene.annot %>% dplyr::select(-Chr)

NodeAnnot.GCHOP <- rbind(Gene.annot, SummEdgeType.GCHOP)

for(gene in GCHOPgenes){
  
  indx <- NULL
  
  indx <- grep(gene, GCHOP.communities$Members_list)
  
  if(!is.null(indx)){
    if(!exists("SummIndx")){
      SummIndx <-  tibble(Community = GCHOP.communities$Community[indx],
                          Gene = gene) 
    } else {
      SummIndxi <-  tibble(Community = GCHOP.communities$Community[indx],
                           Gene = gene)
      SummIndx <-  rbind(SummIndx, SummIndxi)
    }
  }
  
}

GCHOP_EdgeType_Network <- SummIndx %>% mutate(Weight = 1)

biomarkers <- stringr::str_split(rownames(SSN.GCHOP), "_") %>% unlist
sources <- biomarkers[c(T,F)]
targets <- biomarkers[c(F,T)]

GCHOP_Gene_Network <- tibble(Community = sources,
                             Gene = targets,
                             Weight = 1)

GCHOP_EdgeType_Network <- rbind(GCHOP_EdgeType_Network, GCHOP_Gene_Network)

write.table(GCHOP_EdgeType_Network, file.path(FEA.dir, "GCHOP_EdgeType_Network.tsv"), quote = F, sep = "\t", row.names = F)
write.table(NodeAnnot.GCHOP, file.path(FEA.dir, "EdgeType_NodeAnnotGCHOP.tsv"), quote = F, sep = "\t", row.names = F)


