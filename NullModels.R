
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

#### Get Pool of edges from univariate selection, select 1326 random edges x1000 times ####

library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

for(i in 1:10){
  
  Summ.i <- readRDS(paste0("SummaryKMs/SummaryKM_RCHOP_RP", i,".RDS"))
  
  Summ.i <- Summ.i %>% filter(PFS.KM.pval < 0.05 | Concordance.PFS >= 0.6 |
                                OS.KM.pval < 0.05 | Concordance.OS >= 0.6)
  
  if(i == 1){
    Summ.R <- Summ.i
  } else {
    Summ.R <- rbind(Summ.R, Summ.i)
  }
  
}

EdgesPool <- unique(Summ.R$Edge)

RCHOP.Edges <- Summ.R %>% group_by(Edge) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>%
  filter(Count >= 6) %>% dplyr::select(Edge) %>% unlist

length(RCHOP.Edges)
#[1] 1326

#Get 1000 random sets of 1326 edges

set.seed(123)
NullEdges <- lapply(1:1000, function(i) sample(EdgesPool, size = length(RCHOP.Edges)))

saveRDS(NullEdges, file.path(NullModels.dir, "NullSets.RDS"))

relevantEdges <- unique(unlist(NullEdges))

#Get SSNs

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

CI.Full <- readRDS(file.path(ClinicalInfo,"CompleteCI_5yrCensored.RDS"))
CI.Full <- CI.Full %>% mutate(Source = case_when(Source == "GOYA" ~ "GOYA",
                                                 .default = "GDC"),
                              primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist

ExprMs <- "NormData/"
SSNs <- "SSNs/NullModel_RCHOP/"
dir.create(SSNs)

GCN <- data.frame(Edge = relevantEdges, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

IterateGCN_Lioness <- function(i) {
  
  TARGET_NBM_ALALi_norm <- readRDS(paste0(ExprMs, "NormTrainSet_", i, ".RDS"))
  rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))
  
  TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% c(RCHOP.samples, NormalSamples)]
  
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
  
  write.csv(SSN_i, paste0(SSNs,"SSN_NullModel_RCHOP_TrainSet_", i, ".csv"), quote = F, row.names = T)
  
}

library(future.apply)
future::plan(multisession, workers = 10)
future_lapply(1:10, FUN = IterateGCN_Lioness, future.seed = TRUE)

future::plan(sequential)

#### Evaluate null models ####

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

NullEdges <- readRDS(file.path(NullModels.dir, "NullSets.RDS"))

ElasticNet.Null <- function(x){
  nullSet <- unlist(NullEdges[x])
  
  for(i in 1:10){
    
    Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/NullModel_RCHOP/NullModel_RCHOPSSN_NullModel_RCHOP_TrainSet_", i, ".csv"), 
                                   row.names = 1)
    
    Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
    
    Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% nullSet,]
    
    Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
    
    CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
    
    CoxSurvCI <- CoxSurvCI %>% 
      inner_join(CI.RCHOP %>% 
                   dplyr::select(Sample, PFS_status, PFS_time, Source), by = "Sample")
    
    X <- as.matrix(CoxSurvCI %>% dplyr::select(-PFS_time,  -PFS_status,
                                               -Sample, -Source))
    y <- Surv(CoxSurvCI$PFS_time, CoxSurvCI$PFS_status)
    
    set.seed(123)
    cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 0.5)
    
    selected_edges <- coef(cv_fit, s = "lambda.min") 
    
    selected_edges_df <- as.data.frame(as.matrix(selected_edges))
    
    selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
    colnames(selected_edges_df)[1] <- "lambda.min"
    
    selected_edges_df <- selected_edges_df %>% filter(lambda.min != 0) %>% arrange(abs(lambda.min))
    
    selected_edges_df.PFS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "PFS")
    
    #Now OS
    CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
    
    CoxSurvCI <- CoxSurvCI %>% 
      inner_join(CI.RCHOP %>% 
                   dplyr::select(Sample, OS_status, OS_time, Source), by = "Sample")
    
    X <- as.matrix(CoxSurvCI %>% dplyr::select(-OS_time,  -OS_status,
                                               -Sample, -Source))
    y <- Surv(CoxSurvCI$OS_time, CoxSurvCI$OS_status)
    
    set.seed(123)
    cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 0.5)
    
    selected_edges <- coef(cv_fit, s = "lambda.min") 
    
    selected_edges_df <- as.data.frame(as.matrix(selected_edges))
    
    selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
    colnames(selected_edges_df)[1] <- "lambda.min"
    
    selected_edges_df <- selected_edges_df %>% filter(lambda.min != 0) %>% arrange(abs(lambda.min))
    
    selected_edges_df.OS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "OS")
    
    selected_edges_df <- rbind(selected_edges_df.OS, selected_edges_df.PFS)
    
    if(i == 1){
      GlobalEdges <- selected_edges_df
    } else {
      GlobalEdges <- rbind(GlobalEdges, selected_edges_df)
    }
    
  }
  
  Edges.R <- GlobalEdges %>% group_by(Edge) %>% distinct(RandomSet) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count))  %>% filter(Count >= 7) %>% dplyr::select(Edge) %>% unlist %>% as.vector
  
  return(Edges.R)
  
}

future::plan(multisession, workers = 60)
ElasticNet.NullResults <- future_lapply(1:length(NullEdges), 
                                        ElasticNet.Null,
                                        future.seed = TRUE)

plan(sequential)
saveRDS(ElasticNet.NullResults, file.path(NullModels.dir, "ElasticNetNullResults.RDS"))

length(unique(unlist(ElasticNet.NullResults)))

ElasticNet.NullResults <- readRDS(file.path(NullModels.dir, "ElasticNetNullResults.RDS"))

ElasticNetNullResults.Edges <- unique(unlist(ElasticNet.NullResults))

#Subset the complete network files to improve processing time during AIC
subsetNetworks <- function(i){
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/NullModel_RCHOP/NullModel_RCHOPSSN_NullModel_RCHOP_TrainSet_", i, ".csv"), 
                                 row.names = 1)
  
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% ElasticNetNullResults.Edges, ]
  
  saveRDS(Lioness_DLBCL_NCCI,
          paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
  
}

future::plan(multisession, workers = 10)
future_lapply(1:10, 
              subsetNetworks,
              future.seed = TRUE)
plan(sequential)

#Evaluate AIC 
EvaluateAIC <- function(x){
  
  Edges.R <- ElasticNet.NullResults[[x]]
  
  if(length(Edges.R) > 1){
    
    for(i in 1:10){
      
      print(i)
      
      Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
      Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.R,]
      Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
      
      Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
      
      CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
      
      CoxSurvCI <- CoxSurvCI %>% 
        inner_join(CI.RCHOP %>% 
                     dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
      
      Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                         data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
      Cox.PFS.2 <- stepAIC(Cox.PFS.1)
      
      Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                        data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
      Cox.OS.2 <- stepAIC(Cox.OS.1)
      
      if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
        
        coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
        coeff2.PFS <- coeff2.PFS[,c(1,5)]
        colnames(coeff2.PFS) <- c("Coefficient", "p_value")
        coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                            RP = i, Type = "PFS")
      } else {
        coeff2.PFS <- NA
      }
      if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
        
        coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
        coeff2.OS <- coeff2.OS[,c(1,5)]
        colnames(coeff2.OS) <- c("Coefficient", "p_value")
        coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                          RP = i, Type = "OS")
      } else {
        coeff2.OS <- NA
      }
      
      if(i == 1){
        coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
      } else {
        coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
      }
      
    }
    
    Edges.Filt1 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
      distinct(RP) %>%
      summarise(Count = n()) %>%
      arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
    
    if(length(Edges.Filt1) > 1){
      
      for(i in 1:10){
        
        print(i)
        
        Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
        Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt1,]
        Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
        
        Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
        
        CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
        
        CoxSurvCI <- CoxSurvCI %>% 
          inner_join(CI.RCHOP %>% 
                       dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
        
        Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                           data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
        Cox.PFS.2 <- stepAIC(Cox.PFS.1)
        
        Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                          data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
        Cox.OS.2 <- stepAIC(Cox.OS.1)
        
        if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
        } else {
          coeff2.PFS <- NA
        }
        
        if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
          
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
        } else {
          coeff2.OS <- NA
        }
        
        if(i == 1){
          coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
        } else {
          coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
        }
        
      }
      
      Edges.Filt2 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
        distinct(RP) %>%
        summarise(Count = n()) %>%
        arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
      
      if(length(Edges.Filt2) > 1){
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt2,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.RCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
            
            coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
            coeff2.PFS <- coeff2.PFS[,c(1,5)]
            colnames(coeff2.PFS) <- c("Coefficient", "p_value")
            coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                                RP = i, Type = "PFS")
          } else {
            coeff2.PFS <- NA
          }
          
          if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
            
            coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
            coeff2.OS <- coeff2.OS[,c(1,5)]
            colnames(coeff2.OS) <- c("Coefficient", "p_value")
            coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                              RP = i, Type = "OS")
          } else {
            coeff2.OS <- NA
          }
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
        }
        
        Edges.Filt3 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>%
          group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        if(length(Edges.Filt3) > 1){
          
          for(i in 1:10){
            
            print(i)
            
            Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
            Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt3,]
            Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
            
            Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
            
            CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
            
            CoxSurvCI <- CoxSurvCI %>% 
              inner_join(CI.RCHOP %>% 
                           dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
            
            Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                               data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
            Cox.PFS.2 <- stepAIC(Cox.PFS.1)
            
            Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                              data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
            Cox.OS.2 <- stepAIC(Cox.OS.1)
            
            if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
              
              coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
              coeff2.PFS <- coeff2.PFS[,c(1,5)]
              colnames(coeff2.PFS) <- c("Coefficient", "p_value")
              coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                                  RP = i, Type = "PFS")
            } else {
              coeff2.PFS <- NA
            }
            
            if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
              
              coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
              coeff2.OS <- coeff2.OS[,c(1,5)]
              colnames(coeff2.OS) <- c("Coefficient", "p_value")
              coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                                RP = i, Type = "OS")
            } else {
              coeff2.OS <- NA
            }
            
            if(i == 1){
              coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
            } else {
              coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
            }
            
          }
          
          Edges.Filt4 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>%
            group_by(Edge) %>% 
            distinct(RP) %>%
            summarise(Count = n()) %>%
            arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
          
          if(length(Edges.Filt4) > 1){
            
            FinalEdges <-  coeff.summary %>%
              filter(p_value < 0.05) %>%
              dplyr::select(Edge, Type) %>%
              group_by(Edge, Type) %>%
              summarise(Count = n()) %>%
              group_by(Edge)  %>% filter(Count >= 8) %>%
              filter(n_distinct(Type) == 2) %>% # Solo conservar los Edge que están en OS y PFS
              ungroup() %>% pull(Edge) %>% unique %>% as.vector
            
            if(length(FinalEdges) > 1){
              
              for(i in 1:10){
                
                print(i)
                
                Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
                Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% FinalEdges,]
                Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
                
                Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
                
                CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
                
                CoxSurvCI <- CoxSurvCI %>% 
                  inner_join(CI.RCHOP %>% 
                               dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
                
                Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                                   data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
                Cox.PFS.2 <- stepAIC(Cox.PFS.1)
                
                Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                                  data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
                Cox.OS.2 <- stepAIC(Cox.OS.1)
                
                if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
                  
                  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
                  coeff2.PFS <- coeff2.PFS[,c(1,5)]
                  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
                  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                                      RP = i, Type = "PFS")
                } else {
                  coeff2.PFS <- NA
                }
                
                if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
                  
                  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
                  coeff2.OS <- coeff2.OS[,c(1,5)]
                  colnames(coeff2.OS) <- c("Coefficient", "p_value")
                  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                                    RP = i, Type = "OS")
                } else {
                  coeff2.OS <- NA
                }
                
                if(i == 1){
                  coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
                } else {
                  coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
                }
                
                
              }
              
              CheckFinalEdges <- coeff.summary %>%
                filter(p_value < 0.05) %>%
                dplyr::select(Edge, Type) %>%
                group_by(Edge, Type) %>%
                summarise(Count = n()) %>%
                group_by(Edge)  %>% filter(Count >= 8) %>%
                filter(n_distinct(Type) == 2) %>% # Solo conservar los Edge que están en OS y PFS
                ungroup() %>% pull(Edge) %>% unique %>% as.vector
              
              if(length(CheckFinalEdges) > 1){
                
                for(i in 1:10){
                  
                  print(i)
                  
                  Lioness_DLBCL_NCCI <- readRDS(paste0("SSNs/NullModel_RCHOP/SubsetAIC_SSN_NullModel_RCHOP_TrainSet_", i, ".RDS"))
                  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges,]
                  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
                  
                  Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
                  
                  CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
                  
                  CoxSurvCI <- CoxSurvCI %>% 
                    inner_join(CI.RCHOP %>% 
                                 dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
                  
                  Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                                     data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
                  Cox.PFS.2 <- stepAIC(Cox.PFS.1)
                  
                  Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                                    data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
                  Cox.OS.2 <- stepAIC(Cox.OS.1)
                  
                  if(nrow(as.data.frame(summary(Cox.PFS.2)$coefficient)) != 0){
                    
                    coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
                    coeff2.PFS <- coeff2.PFS[,c(1,5)]
                    colnames(coeff2.PFS) <- c("Coefficient", "p_value")
                    coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                                        RP = i, Type = "PFS")
                  } else {
                    coeff2.PFS <- NA
                  }
                  
                  if(nrow(as.data.frame(summary(Cox.OS.2)$coefficient)) != 0){
                    
                    coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
                    coeff2.OS <- coeff2.OS[,c(1,5)]
                    colnames(coeff2.OS) <- c("Coefficient", "p_value")
                    coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                                      RP = i, Type = "OS")
                  } else {
                    coeff2.OS <- NA
                  }
                  
                  if(i == 1){
                    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
                  } else {
                    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
                  }
                  
                }
                
                CheckFinalEdges2 <- coeff.summary %>%
                  filter(p_value < 0.05) %>%
                  dplyr::select(Edge, Type) %>%
                  group_by(Edge, Type) %>%
                  summarise(Count = n()) %>%
                  group_by(Edge)  %>% filter(Count >= 8) %>%
                  filter(n_distinct(Type) == 2) %>% # Solo conservar los Edge que están en OS y PFS
                  ungroup() %>% pull(Edge) %>% unique %>% as.vector
                
                return(CheckFinalEdges2)
                
              } else {
                return(CheckFinalEdges)
              }
            } else {
              return(FinalEdges)
            }
          } else {
            return(Edges.Filt4)
          }
        } else {
          return(Edges.Filt3)
        }
      } else {
        return(Edges.Filt2)
      }
    } else {
      return(Edges.Filt1)
    }
  } else {
    return(Edges.R)
  }
}

future::plan(multisession, workers = 50)
AIC.NullResults <- future_lapply(1:1000, 
                                 EvaluateAIC,
                                 future.seed = TRUE)
plan(sequential)

saveRDS(AIC.NullResults, file.path(NullModels.dir, "AICNullResults.RDS"))

#### Perform cross-validation between cohorts (inside partitions) ####

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
library(caret)
library(ggbiplot)
library(survival)
library(glmnet)
library(MASS)

DLBCL <- readRDS(file.path(Raw.dir, "DLBCL_RawCounts.RDS"))
annot <- readRDS(file.path(Raw.dir, "annot_GDC.RDS"))
annot <- annot[annot$HGNC_symbol %in% rownames(DLBCL), ]

CompleteCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Normal <- data.frame(Sample = NormalSamples, 
                        IPI = NA,
                        Gender = NA,
                        OS_status = NA,
                        OS_time = NA, 
                        PFS_status = NA,
                        PFS_time = NA, 
                        COO = NA, 
                        Genetic_subtype = NA,
                        Treatment = "Lymphoid normal",
                        Source = "GDC",
                        primary_diagnosis = "Lymphoid normal")

CompleteCI <- CompleteCI %>% mutate(primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CompleteCI, CI.Normal)

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

DLBCL.norm <- norm.A(DLBCL, annot, CI.Full.Normal)

saveRDS(DLBCL.norm, "NormData/DLBCL_NormCounts.RDS")

AIC.NullResults <- readRDS(file.path(NullModels.dir, "AICNullResults.RDS"))
AIC.NullResults <- unique(unlist(AIC.NullResults))

PACO14 <- readRDS(file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))

EdgesInModels <- c(PACO14, AIC.NullResults)

EdgesInModels <- unique(EdgesInModels)

GCN <- data.frame(Edge = EdgesInModels, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

CI.RCHOP <- CI.Full.Normal %>% filter(Treatment == "Rituximab" |
                                        primary_diagnosis == "Lymphoid normal")

TARGET_NBM_ALALi_norm <- DLBCL.norm[,colnames(DLBCL.norm) %in% CI.RCHOP$Sample]
rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))

TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% c(RCHOP.samples, NormalSamples)]

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

write.csv(SSN_i, "SSNs/SSN_RCHOP_EdgesInModels.csv", quote = F, row.names = T)

AIC.NullResults <- readRDS(file.path(NullModels.dir, "AICNullResults.RDS"))

CI.RCHOP.GOYA <- CI.RCHOP %>% filter(Source == "GEO")

CI.RCHOP.GDC <- CI.RCHOP %>% filter(Source == "GDC")

Lioness_DLBCL_NCCI <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]

Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))

CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))

CoxSurvCI <- CoxSurvCI %>% 
  inner_join(CI.RCHOP %>% 
               dplyr::select(Sample, PFS_status, PFS_time,
                             OS_status, OS_time, Source), by = "Sample")

CoxSurvCI.GDC <- CoxSurvCI %>% filter(Source == "GDC")
CoxSurvCI.GEO <- CoxSurvCI %>% filter(Source == "GEO")

AIC.NullResults[[length(AIC.NullResults) + 1]] <- PACO14

for(i in 1:length(AIC.NullResults)){
  
  print(i)
  
  edges <- AIC.NullResults[[i]]
  
  if(length(edges) > 0){
    #Train in GDC PFS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    Cox.PFS.GDC.PFS <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                              data = CoxSurvCI.GDCi)
    
    CoxSurvCI.GEOi$RiskScore.Test.PFS <- predict(Cox.PFS.GDC.PFS, type = "risk",
                                                 newdata = CoxSurvCI.GEOi)
    
    Cox.PFS.GDC.Test.PFS <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore.Test.PFS, 
                                   data = CoxSurvCI.GEOi)
    
    #Train in GEO PFS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    Cox.PFS.GEO.PFS <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI.GEOi)
    
    CoxSurvCI.GDCi$RiskScore.Test.PFS <- predict(Cox.PFS.GEO.PFS, type = "risk",
                                                  newdata = CoxSurvCI.GDCi)
    
    Cox.PFS.GEO.Test.PFS <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore.Test.PFS, 
                                  data = CoxSurvCI.GDCi)
    
    
    #Train in GDC OS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    Cox.OS.GDC.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI.GDCi)
    
    CoxSurvCI.GEOi$RiskScore.Test.OS <- predict(Cox.OS.GDC.OS, type = "risk",
                                                newdata = CoxSurvCI.GEOi)
    
    Cox.OS.GDC.Test.OS <- coxph(Surv(OS_time, OS_status) ~ RiskScore.Test.OS, 
                                 data = CoxSurvCI.GEOi)
    
    #Train in GEO OS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    Cox.OS.GEO.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                           data = CoxSurvCI.GEOi)
    
    CoxSurvCI.GDCi$RiskScore.Test.OS <- predict(Cox.OS.GEO.OS, type = "risk",
                                                 newdata = CoxSurvCI.GDCi)
    
    Cox.OS.GEO.Test.OS <- coxph(Surv(OS_time, OS_status) ~ RiskScore.Test.OS, 
                                data = CoxSurvCI.GDCi)
    
    summary <- data.frame(Model = i,
                          EdgeNumber = length(edges),
                          C.GDC.Train.PFS = summary(Cox.PFS.GDC.PFS)$concordance[1],
                          C.GEO.Test.PFS = summary(Cox.PFS.GDC.Test.PFS)$concordance[1],
                          C.GEO.Train.PFS = summary(Cox.PFS.GEO.PFS)$concordance[1],
                          C.GDC.Test.PFS = summary(Cox.PFS.GEO.Test.PFS)$concordance[1],
                          C.GDC.Train.OS = summary(Cox.OS.GDC.OS)$concordance[1],
                          C.GEO.Test.OS = summary(Cox.OS.GDC.Test.OS)$concordance[1],
                          C.GEO.Train.OS = summary(Cox.OS.GEO.OS)$concordance[1],
                          C.GDC.Test.OS = summary(Cox.OS.GEO.Test.OS)$concordance[1])
    
  } else {
    summary <- data.frame(Model = i,
                          EdgeNumber = length(edges),
                          C.GDC.Train.PFS = NA,
                          C.GEO.Test.PFS = NA,
                          C.GEO.Train.PFS = NA,
                          C.GDC.Test.PFS = NA,
                          C.GDC.Train.OS = NA,
                          C.GEO.Test.OS = NA,
                          C.GEO.Train.OS = NA,
                          C.GDC.Test.OS = NA)
    
  }
  
  if(i == 1){
    GlobalSummary <- summary
  } else {
    GlobalSummary <- rbind(GlobalSummary, summary)
  }
  
}

saveRDS(GlobalSummary, file.path(NullModels.dir, "GlobalSummary_NullModelsAndPACO_CrossValidations.RDS"))


#### Select 14 random edges from relevant edges that passed the univariate filter ####

library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)
library(survival)
library(glmnet)
library(MASS)

for(i in 1:10){
  
  Summ.i <- readRDS(paste0("SummaryKMs/SummaryKM_RCHOP_RP", i,".RDS"))
  
  Summ.i <- Summ.i %>% filter(PFS.KM.pval < 0.05 | Concordance.PFS >= 0.6 |
                                OS.KM.pval < 0.05 | Concordance.OS >= 0.6)
  
  if(i == 1){
    Summ.R <- Summ.i
  } else {
    Summ.R <- rbind(Summ.R, Summ.i)
  }
  
}

EdgesPool <- Summ.R %>% group_by(Edge) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>%
  filter(Count >= 6) %>% dplyr::select(Edge) %>% unlist %>% as.vector

#Get 1000 random sets of 14 edges

set.seed(123)
NullEdges <- lapply(1:1000, function(i) sample(EdgesPool, size = 14))

saveRDS(NullEdges, file.path(NullModels.dir, "NullSets2.RDS"))

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

CI.RCHOP <- CI.Full.Normal %>% filter(Treatment == "Rituximab" |
                                        primary_diagnosis == "Lymphoid normal")

GCN <- data.frame(Edge = EdgesPool, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

TARGET_NBM_ALALi_norm <- DLBCL.norm[,colnames(DLBCL.norm) %in% CI.RCHOP$Sample]
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

write.csv(SSN_i, "SSNs/SSN_RCHOP_EdgesInModels_Null2.csv", quote = F, row.names = T)

#Evaluate model

CI.RCHOP.GOYA <- CI.RCHOP %>% filter(Source == "GEO")

CI.RCHOP.GDC <- CI.RCHOP %>% filter(Source == "GDC")

Lioness_DLBCL_NCCI <- SSN_i[,!(colnames(SSN_i) %in% NormalSamples)]

Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))

CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))

CoxSurvCI <- CoxSurvCI %>% 
  inner_join(CI.RCHOP %>% 
               dplyr::select(Sample, PFS_status, PFS_time,
                             OS_status, OS_time, Source), by = "Sample")

CoxSurvCI.GDC <- CoxSurvCI %>% filter(Source == "GDC")
CoxSurvCI.GEO <- CoxSurvCI %>% filter(Source == "GEO")

for(i in 1:length(NullEdges)){
  
  print(i)
  
  edges <- NullEdges[[i]]
  
  if(length(edges) > 0){
    #Train in GDC PFS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    Cox.PFS.GDC.PFS <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                              data = CoxSurvCI.GDCi)
    
    CoxSurvCI.GEOi$RiskScore.Test.PFS <- predict(Cox.PFS.GDC.PFS, type = "risk",
                                                 newdata = CoxSurvCI.GEOi)
    
    Cox.PFS.GDC.Test.PFS <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore.Test.PFS, 
                                   data = CoxSurvCI.GEOi)
    
    #Train in GEO PFS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "PFS_time", "PFS_status"))
    
    Cox.PFS.GEO.PFS <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI.GEOi)
    
    CoxSurvCI.GDCi$RiskScore.Test.PFS <- predict(Cox.PFS.GEO.PFS, type = "risk",
                                                  newdata = CoxSurvCI.GDCi)
    
    Cox.PFS.GEO.Test.PFS <- coxph(Surv(PFS_time, PFS_status) ~ RiskScore.Test.PFS, 
                                  data = CoxSurvCI.GDCi)
    
    
    #Train in GDC OS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    Cox.OS.GDC.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI.GDCi)
    
    CoxSurvCI.GEOi$RiskScore.Test.OS <- predict(Cox.OS.GDC.OS, type = "risk",
                                                newdata = CoxSurvCI.GEOi)
    
    Cox.OS.GDC.Test.OS <- coxph(Surv(OS_time, OS_status) ~ RiskScore.Test.OS, 
                                 data = CoxSurvCI.GEOi)
    
    #Train in GEO OS
    
    CoxSurvCI.GDCi <- CoxSurvCI.GDC %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    CoxSurvCI.GEOi <- CoxSurvCI.GEO %>% dplyr::select(c(edges, "OS_time", "OS_status"))
    
    Cox.OS.GEO.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                           data = CoxSurvCI.GEOi)
    
    CoxSurvCI.GDCi$RiskScore.Test.OS <- predict(Cox.OS.GEO.OS, type = "risk",
                                                 newdata = CoxSurvCI.GDCi)
    
    Cox.OS.GEO.Test.OS <- coxph(Surv(OS_time, OS_status) ~ RiskScore.Test.OS, 
                                data = CoxSurvCI.GDCi)
    
    summary <- data.frame(Model = i,
                          EdgeNumber = length(edges),
                          C.GDC.Train.PFS = summary(Cox.PFS.GDC.PFS)$concordance[1],
                          C.GEO.Test.PFS = summary(Cox.PFS.GDC.Test.PFS)$concordance[1],
                          C.GEO.Train.PFS = summary(Cox.PFS.GEO.PFS)$concordance[1],
                          C.GDC.Test.PFS = summary(Cox.PFS.GEO.Test.PFS)$concordance[1],
                          C.GDC.Train.OS = summary(Cox.OS.GDC.OS)$concordance[1],
                          C.GEO.Test.OS = summary(Cox.OS.GDC.Test.OS)$concordance[1],
                          C.GEO.Train.OS = summary(Cox.OS.GEO.OS)$concordance[1],
                          C.GDC.Test.OS = summary(Cox.OS.GEO.Test.OS)$concordance[1])
    
  } else {
    summary <- data.frame(Model = i,
                          EdgeNumber = length(edges),
                          C.GDC.Train.PFS = NA,
                          C.GEO.Test.PFS = NA,
                          C.GEO.Train.PFS = NA,
                          C.GDC.Test.PFS = NA,
                          C.GDC.Train.OS = NA,
                          C.GEO.Test.OS = NA,
                          C.GEO.Train.OS = NA,
                          C.GDC.Test.OS = NA)
    
  }
  
  if(i == 1){
    GlobalSummary <- summary
  } else {
    GlobalSummary <- rbind(GlobalSummary, summary)
  }
  
}

saveRDS(GlobalSummary, file.path(NullModels.dir, "GlobalSummary_NullModels2ndApproach_CrossValidations.RDS"))




