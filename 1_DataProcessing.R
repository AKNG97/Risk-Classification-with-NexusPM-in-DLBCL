
home.dir <- file.path(path.expand("~"), "NexusPM")
dir.create(home.dir)
setwd(home.dir)

Raw.dir <- file.path(home.dir, "RawData/")
dir.create(Raw.dir)
#NOTE: Download Human.GRCh38.p13.annot.tsv into the Raw.dir folder 
#Instruction to download the file can be found on the README file of the repository

ClinicalInfo <- file.path(home.dir, "ClinicalInfo/")
dir.create(ClinicalInfo)
#NOTE: Download the clinical files as described on the README file of the repository
#Save the file from the GDC as DLBCL_NCII_clinical_paper.csv and the file from GEO as Goya_Clinical_Data.csv 
#Yu need them on the ClinicalInfo folder before running this script

NullModels.dir <- file.path(home.dir, "NullModels/")
dir.create(NullModels.dir)

Validation.dir <- file.path(home.dir, "Validation/")
dir.create(Validation.dir)

FEA.dir <- file.path(home.dir, "FEA/")
dir.create(FEA.dir)

#### Libraries ####
library(SummarizedExperiment)
library(TCGAbiolinks)
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

#### Get Biomart annot ####

ensembl <- biomaRt::useEnsembl(biomart = "genes",
                               dataset = "hsapiens_gene_ensembl")

features <- c("ensembl_gene_id", "chromosome_name",
              "start_position", "end_position", "hgnc_symbol",
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs,
               mart = ensembl)

colnames(annot)<-c("Ensembl", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
annot <- annot[annot$HGNC_symbol != "",]
annot <- annot[!duplicated(annot$Ensembl),]
dim(annot)
#[1] 41232    10

saveRDS(annot, file.path(Raw.dir, "annot_Biomart.RDS"))
#annot <- readRDS("/STORAGE/csbig/anakamura/DLBCL_A4/A12/annot_111124.RDS")

#### Download NCICC files ####

qry.rna_NCICC <- GDCquery(project = "NCICCR-DLBCL",
                          data.category= "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
GDCdownload(qry.rna_NCICC)
DLBCL.NICC <- GDCprepare(qry.rna_NCICC, summarizedExperiment = TRUE)

saveRDS(DLBCL.NICC, file.path(Raw.dir, "DLBCL_NCICCR_Raw.RDS"))

CI.NCII <- as.data.frame(colData(DLBCL.NICC))
CI.NCII <- CI.NCII[,c("sample", "submitter_id")]

ext.ci <- read.csv(file.path(ClinicalInfo, "DLBCL_NCII_clinical_paper.csv"))

ext.ci <- ext.ci[ext.ci$dbGaP.subject.ID %in% CI.NCII$submitter_id,]
length(which(ext.ci$Follow_up.Time._yrs == 0))
#[1] 247
ext.ci <- ext.ci[ext.ci$Follow_up.Time._yrs != 0,]
write.csv(ext.ci, file.path(ClinicalInfo, "NCICCR_clinicaldata_FULL.csv"), quote = F, row.names = F)

colnames(DLBCL.NICC) <- gsub("-sample", "", colnames(DLBCL.NICC))
DLBCL.NICC <- DLBCL.NICC[,colnames(DLBCL.NICC) %in% ext.ci$dbGaP.subject.ID]

#### Download Normal Lymph ####
qry.rna_BL <- GDCquery(project = "CGCI-BLGSP",
                       data.category= "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       barcode = c("BLGSP-71-19-99998",
                                   "BLGSP-71-19-99997",
                                   "BLGSP-71-19-99988",
                                   "BLGSP-71-19-99999",
                                   "BLGSP-71-19-99996",
                                   "BLGSP-71-19-99989"))
GDCdownload(qry.rna_BL)
NormalLymph <- GDCprepare(qry.rna_BL, summarizedExperiment = TRUE)

coldata <- as.data.frame(colData(NormalLymph)) %>%
  dplyr::select(submitter_id, sample)

NormalLymph <- assay(NormalLymph, 3)

coldata <- coldata[match(colnames(NormalLymph), coldata$sample),]

dds <- DESeqDataSetFromMatrix(countData = NormalLymph,
                              colData = coldata,
                              design = ~ submitter_id)

ddsColl <- collapseReplicates(dds, dds$submitter_id, dds$submitter_id)
Normal_lymph.raw.coll <- counts(ddsColl)

saveRDS(Normal_lymph.raw.coll, paste0(Raw.dir, "Normal_Lymph_BL_Raw_collapsed_counts.RDS"))

#### Get CTSP ####

qry.rna_DLBCL_CTSP <- GDCquery(project = "CTSP-DLBCL1",
                               data.category= "Transcriptome Profiling",
                               data.type = "Gene Expression Quantification",
                               workflow.type = "STAR - Counts")
GDCdownload(qry.rna_DLBCL_CTSP)
qry.rna_DLBCL_CTSP.df <- qry.rna_DLBCL_CTSP$results[[1]]
qry.rna_DLBCL_CTSP.df <- qry.rna_DLBCL_CTSP.df %>% filter(version == "2")
qry.rna_DLBCL_CTSP$results[[1]] <- qry.rna_DLBCL_CTSP.df
DLBCL.CTSP <- GDCprepare(qry.rna_DLBCL_CTSP, summarizedExperiment = TRUE) 

CI.original <- as.data.frame(colData(DLBCL.CTSP))
CI.original$sample <- gsub("-sample", "", CI.original$sample)

#Complete the clinical data using the json files 
qry.rna_DLBCL_CTSP.clinical <- GDCquery(project = "CTSP-DLBCL1",
                                        data.category= "Clinical",
                                        data.type = "Clinical Supplement")
GDCdownload(qry.rna_DLBCL_CTSP.clinical)

#List clinical json

json.files <- list.files("GDCdata/CTSP-DLBCL1/Clinical/Clinical_Supplement/")

for(i in 1:length(json.files)){
  print(i)
  json.dir <- file.path("GDCdata/CTSP-DLBCL1/Clinical/Clinical_Supplement",
                        json.files[i])
  jsoni <- list.files(json.dir)
  clinical_data <- fromJSON(file.path(json.dir, jsoni))
  
  Regimes <- unique(clinical_data$ClinicalData$`GDC Clinical Data Entities`$regimen_or_line_of_therapy)
  Regimes <- Regimes[!is.na(Regimes)]
  
  if(Regimes == "R-CHOP"){
    #print(i)
    submitter_id <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$submitter_id
    submitter_id <- submitter_id[grep("-diagnosis", submitter_id)]
    #print(paste0(i, submitter_id))
    days_to_follow_up <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$days_to_follow_up
    days_to_follow_up <- max(days_to_follow_up[!is.na(days_to_follow_up)])
    vital_status <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$vital_status
    vital_status <- vital_status[!is.na(vital_status)]
    days_to_death <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$days_to_death
    days_to_death <- days_to_death[!is.na(days_to_death)]
    cause_of_death <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$cause_of_death
    cause_of_death <- cause_of_death[!is.na(cause_of_death)]
    lost_to_followup <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$lost_to_followup
    lost_to_followup <- lost_to_followup[!is.na(lost_to_followup)]
    progression <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$progression_or_recurrence
    progression <- progression[!is.na(progression)]
    progression <- ifelse(any(progression == "Yes"), "Yes", "No")
    PFS_time <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$days_to_progression
    PFS_time <- PFS_time[!is.na(PFS_time)]
    #PFS_time <- ifelse(lenght(PFS_time) > 1), "Yes", "No") 
    PFS_time <- PFS_time[!is.na(PFS_time)]
    PFS_time <- ifelse(length(PFS_time) == 0, 0, PFS_time)
    days_to_death <- ifelse(is.null(days_to_death), NA, days_to_death)
    ipi <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$international_prognostic_index
    ipi <- ipi[!is.na(ipi)]
    gender <- clinical_data$ClinicalData$`GDC Clinical Data Entities`$gender
    gender <- gender[!is.na(gender)]
    
    UpdatedCI <- data.frame(submitter_id = submitter_id,
                            gender = gender,
                            ipi = case_when(is.null(ipi) ~ NA,
                                            .default = ipi),
                            days_to_follow_up_UP = case_when(is.null(days_to_follow_up) ~ NA,
                                                             .default = days_to_follow_up),
                            vital_status_UP = case_when(is.null(vital_status) ~ NA,
                                                        .default = vital_status),
                            days_to_death_UP = case_when(is.null(days_to_death) ~ NA,
                                                         .default = days_to_death),
                            cause_of_death_UP = case_when(is.null(cause_of_death) ~ NA,
                                                          .default = cause_of_death),
                            lost_to_followup_UP = case_when(is.null(lost_to_followup) ~ NA,
                                                            .default = lost_to_followup),
                            progression = case_when(is.null(progression) ~ NA,
                                                    .default = progression),
                            PFS_time = case_when(is.null(PFS_time) ~ NA,
                                                 PFS_time == 0 ~ days_to_death,
                                                 .default = PFS_time))
    
    if(!exists("UpdatedCI.Global")){
      UpdatedCI.Global  <- UpdatedCI
    } else {
      UpdatedCI.Global <- rbind(UpdatedCI.Global, UpdatedCI)
    }
    
  }
  
}
UpdatedCI.Global$submitter_id <- gsub("-diagnosis[0-9]*$", "", UpdatedCI.Global$submitter_id)
UpdatedCI.Global <- UpdatedCI.Global[!duplicated(UpdatedCI.Global$submitter_id),]

saveRDS(UpdatedCI.Global, file.path(ClinicalInfo, "RCHOP_CTSPDLBCL_Update_16_12_24.RDS"))

CI.Full <- UpdatedCI.Global %>% dplyr::inner_join(CI.original %>% dplyr::select(sample, international_prognostic_index, ann_arbor_pathologic_stage,
                                                                                cause_of_death, vital_status, days_to_death),
                                                  by = c("submitter_id" = "sample"))

colnames(DLBCL.CTSP) <- gsub("-sample", "", colnames(DLBCL.CTSP))

DLBCL.CTSP <- DLBCL.CTSP[,colnames(DLBCL.CTSP) %in% CI.Full$submitter_id]
saveRDS(DLBCL.CTSP, file.path(Raw.dir, "RCHOP_CTSP_DLBCL_A17.RDS"))
CI.Full <- CI.Full[CI.Full$submitter_id %in% colnames(DLBCL.CTSP),]
write.csv(CI.Full, paste0(ClinicalInfo, "CIFull_CTSP_DLBCL_RCHOP_A17.csv"),
          quote = F, row.names = F)

#### Filter annotation ####

annot_GDC <- as.data.frame(rowData(DLBCL.NICC))
annot_GDC$Ensembl <- gsub("\\.[0-9]*", "", annot_GDC$gene_id)
dim(annot_GDC)
annot_GDC <- annot_GDC %>% inner_join(annot[,c("Ensembl","Ensembl_ID_Version","HGNC_symbol", "HGNC_ID", "Type", "Chr", "GC", "Start", "End","Length")], 
                                        by = c("Ensembl" = "Ensembl"))
dim(annot)
#[1] 41232    10
dim(annot_GDC)
#[1] 40935    20

which(duplicated(annot_GDC$gene_name))
#[1] 38098 38605 38823 39947 39948 39994 40024 40061 40066 40787
which(duplicated(annot_GDC$HGNC_symbol))
# [1] 31553 35560 36385 37794 38002 38605 38751 39947 39948 39994 40024 40061
# [13] 40066 40124 40763 40787 40791 40902

duplicated_GeneNames <- annot_GDC$HGNC_symbol[which(duplicated(annot_GDC$HGNC_symbol))]
dim(annot_GDC)
#40935    19
annot_GDC <- annot_GDC[!(annot_GDC$HGNC_symbol %in% duplicated_GeneNames),]
dim(annot_GDC)
#[1] 40899    20

annot_GDC <- annot_GDC[annot_GDC$Type == "protein_coding",]
dim(annot_GDC)
# [1] 19391    20

saveRDS(annot_GDC, file.path(Raw.dir, "annot_GDC.RDS"))
annot <- readRDS("annot_GDC.RDS")

#### Get GOYA ####

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE125966", "file=GSE125966_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GOYA_raw <- as.data.frame(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

goya_annot <- read.delim(file.path(Raw.dir, "Human.GRCh38.p13.annot.tsv"))
goya_annot <- goya_annot[goya_annot$GeneType == "protein-coding",]
goya_annot <- goya_annot[goya_annot$EnsemblGeneID != "",]

goya_annot <- goya_annot[goya_annot$Symbol %in% annot$HGNC_symbol,]
annot <- annot[annot$HGNC_symbol %in% goya_annot$Symbol,]
length(which(duplicated(goya_annot$EnsemblGeneID)))
#3
length(intersect(goya_annot$EnsemblGeneID, annot$Ensembl))
#[1] 19086
length(intersect(goya_annot$Symbol, annot$HGNC_symbol))
#[1] 19092

goya_annot <- goya_annot[goya_annot$EnsemblGeneID %in% annot$Ensembl,]
rmG <- goya_annot$EnsemblGeneID[which(duplicated(goya_annot$EnsemblGeneID))]
goya_annot <- goya_annot[!(goya_annot$EnsemblGeneID %in% rmG), ]

annot <- annot[annot$Ensembl %in% goya_annot$EnsemblGeneID,]
goya_annot <- goya_annot[goya_annot$EnsemblGeneID %in% annot$Ensembl,]
length(which(duplicated(goya_annot$EnsemblGeneID)))
# 0
length(intersect(goya_annot$EnsemblGeneID, annot$Ensembl))
#[1] 19083
length(intersect(goya_annot$Symbol, annot$HGNC_symbol))
#[1] 19083

GOYA_raw <- GOYA_raw[GOYA_raw$GeneID %in% goya_annot$GeneID,]
GOYA_raw$GeneID <- goya_annot$EnsemblGeneID[match(GOYA_raw$GeneID, goya_annot$GeneID)]

rownames(GOYA_raw) <- GOYA_raw$GeneID
GOYA_raw <- GOYA_raw[,-1]

#### Filter GDC ####

Get_raw_matrix <- function(z) {
  dataFilt <- TCGAanalyze_Filtering(tabDF = z,
                                    method = "quantile",
                                    qnt.cut = 0.25)
  threshold <- round(dim(z)[2]/2)
  ridx <- rowSums(dataFilt == 0) <= threshold
  dataFilt <- dataFilt[ridx, ]
  ridx <- rowMeans(dataFilt) >= 10
  dataFilt <- dataFilt[ridx, ]
  z <- z[rownames(z) %in% rownames(dataFilt), ]
  print(dim(z))
  return(z)
}

NICC_StrandedSecond <- assay(DLBCL.NICC, 3)
CTSP_StrandedSecond <- assay(DLBCL.CTSP, 3)

GDC_Raw.F <- Get_raw_matrix(cbind(CTSP_StrandedSecond[rownames(CTSP_StrandedSecond) %in% annot$gene_id,],
                                   NICC_StrandedSecond[rownames(NICC_StrandedSecond) %in% annot$gene_id,], 
                                   Normal_lymph.raw.coll[rownames(Normal_lymph.raw.coll) %in% annot$gene_id,]))
#[1] 14159   258

annot.GDC <- annot[annot$gene_id %in% rownames(GDC_Raw.F), ]

rownames(GDC_Raw.F) <- annot.GDC$HGNC_symbol[match(annot.GDC$gene_id, rownames(GDC_Raw.F))]

GOYA_Raw.F <- Get_raw_matrix(cbind(GOYA_raw[match(annot$Ensembl, rownames(GOYA_raw)),],
                                   Normal_lymph.raw.coll[rownames(Normal_lymph.raw.coll) %in% annot$gene_id,]))
#[1] 14270   559

annot.GOYA <- annot[annot$Ensembl %in% rownames(GOYA_Raw.F), ]

rownames(GOYA_Raw.F) <- annot.GOYA$HGNC_symbol[match(annot.GOYA$Ensembl, rownames(GOYA_Raw.F))]

length(intersect(rownames(GOYA_Raw.F), rownames(GDC_Raw.F)))
#[1] 12887

GDC_Raw.F <- GDC_Raw.F[rownames(GDC_Raw.F) %in% rownames(GOYA_Raw.F),]
GOYA_Raw.F <- GOYA_Raw.F[rownames(GOYA_Raw.F) %in% rownames(GDC_Raw.F),]

GOYA_Raw.F <- GOYA_Raw.F[match(rownames(GDC_Raw.F), 
                               rownames(GOYA_Raw.F)), -c(554:559)]

#### Get Full Clinical object ####

CI.Full <- CI.Full %>% dplyr::select(submitter_id, ipi, gender,
                                     days_to_follow_up_UP, vital_status_UP,
                                     days_to_death_UP, progression,
                                     PFS_time)

colnames(CI.Full) <- c("Sample", "IPI", "Gender", "Days_to_last_follow_up",
                       "OS_status", "OS_time", "PFS_status",
                       "PFS_time")
CI.Full$OS_time <- ifelse(is.na(CI.Full$OS_time), CI.Full$Days_to_last_follow_up, CI.Full$OS_time)
CI.Full$PFS_time <- ifelse(is.na(CI.Full$PFS_time), CI.Full$Days_to_last_follow_up, CI.Full$PFS_time)

CI.Full$PFS_status <- ifelse(CI.Full$PFS_status == "Yes", 1, 0)
CI.Full$OS_status <- ifelse(CI.Full$OS_status == "Dead", 1, 0)
CI.Full$Gender <- ifelse(CI.Full$Gender == "male", "Male", "Female")

ext.ci <- ext.ci %>% mutate(IPI.c = case_when(IPI.Range == 2 ~ "Low-Intermediate Risk",
                                              IPI.Range == 3 ~ "High-Intermediate Risk",
                                              IPI.Range >= 4 ~ "High",
                                              .default = "Low"))

ext.ci <- ext.ci %>% dplyr::select(dbGaP.subject.ID, Gene.Expression.Subgroup,
                                   Genetic.Subtype, Genetic.Subtype,
                                   Gender, Status.at.Follow_up_.0.Alive_.1.Dead,
                                   Follow_up.Time._yrs, Progression_Free.Survival._PFS_.Status_.0.No.Progressoin_.1.Progression,
                                   Progression_Free.Survival._PFS_.Time._yrs, IPI.c)
colnames(ext.ci) <- c("Sample", "COO", "Genetic_subtype", "Gender",
                      "OS_status", "OS_time", "PFS_status", "PFS_time",
                      "IPI")

ext.ci$Gender <- ifelse(ext.ci$Gender == "M", "Male", "Female")
ext.ci$OS_time <- round(365*ext.ci$OS_time)
ext.ci$PFS_time <- round(365*ext.ci$PFS_time)
ext.ci$PFS_time <- ifelse(is.na(ext.ci$PFS_time), ext.ci$OS_time, ext.ci$PFS_time)

CI.GOYA <- read.csv(file.path(ClinicalInfo, "Goya_Clinical_Data.csv"))
CI.GOYA <- CI.GOYA[,c(2,8)]
CI.GOYA <- CI.GOYA %>% dplyr::rename("Sample" = "Accession",
                                     "IPI" = "Ipi",
                                     "Gender" = "Sex",
                                     "OS_status" = "Oscs"  ,
                                     "OS_time" = "Ostt"  ,
                                     "PFS_status" = "Pfscs"  ,
                                     "PFS_time" = "Pfstt"  ,
                                     "COO" = "Cell.of.origin")

CI.GOYA <- CI.GOYA %>% mutate(Genetic_Subtype = NA)

CI.GOYA <- CI.GOYA %>% dplyr::select(Sample, IPI, Gender, OS_status,
                                     OS_time, PFS_status, PFS_time,
                                     COO, Genetic_subtype, Treatment)

ext.ci <- ext.ci[,c(1,9,4,5,6,7,8,2,3)]
ext.ci <- ext.ci %>% mutate(Treatment = "Rituximab")
ext.ci$IPI <- gsub(" Risk", "", ext.ci$IPI)

CI.Full <- CI.Full[,-4]

CI.Full <- CI.Full %>% mutate(COO = NA,
                              Genetic_subtype = NA,
                              Treatment = "Rituximab")

CI.GOYA$COO <- ifelse(CI.GOYA$COO == "UNCLASSIFIED", "Unclassified", CI.GOYA$COO)
CI.GOYA$COO <- ifelse(CI.GOYA$COO == "", NA, CI.GOYA$COO)

ext.ci$COO <- ifelse(ext.ci$COO == "Unclass", "Unclassified", ext.ci$COO)

CI.Full <- CI.Full %>% mutate(Source = "GDC")
ext.ci <- ext.ci %>% mutate(Source = "GDC")
CI.GOYA <- CI.GOYA %>% mutate(Source = "GEO")

CompleteCI <- rbind(CI.Full, ext.ci, CI.GOYA)
saveRDS(CompleteCI, file.path(ClinicalInfo, "CompleteCI.RDS"))

CompleteCI$PFS_status[CompleteCI$PFS_time > 1825] <- 0
CompleteCI$PFS_time[CompleteCI$PFS_time > 1825] <- 1825

CompleteCI$OS_status[CompleteCI$OS_time > 1825] <- 0
CompleteCI$OS_time[CompleteCI$OS_time > 1825] <- 1825

saveRDS(CompleteCI, file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))

#### Split data ####

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

DLBCL <- cbind(GOYA_Raw.F, GDC_Raw.F)
saveRDS(DLBCL, file.path(Raw.dir, "DLBCL_RawCounts.RDS"))

RandomPartitions <- list()

RP.dir <- paste0(Raw.dir, "RandomPartitionsRAW")
dir.create(RP.dir)

set.seed(123) 
for (i in 1:10) {
  
  trainIndex <- createDataPartition(
    y = interaction(CompleteCI$Treatment, CompleteCI$PFS_status, CompleteCI$Source), 
    p = 0.7, 
    list = T
  )
  
  trainSamples <- c(CompleteCI[trainIndex[[1]], "Sample"], NormalSamples)
  testSamples <- CompleteCI[-trainIndex[[1]], "Sample"]
  
  TrainSet <- DLBCL[,colnames(DLBCL) %in% trainSamples]
  TestSet <- DLBCL[,!(colnames(DLBCL) %in% trainSamples)]
  
  saveRDS(TrainSet, file.path(RP.dir, paste0("Train_RandomPartition_", i,".RDS")))
  saveRDS(TestSet, file.path(RP.dir, paste0("Test_RandomPartition_", i,".RDS")))
  
}

#### Normalization ####
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

annot.Global <- annot[annot$HGNC_symbol %in% rownames(DLBCL),]
annot.Global <- annot.Global[match(rownames(DLBCL), annot.Global$HGNC_symbol),]

dir.create("NormData/")
dir.create("PCAs/")

Norm.Partitions <- function(i){
  TrainSeti <- readRDS(paste0("RawData/RandomPartitionsRAW/Train_RandomPartition_", i,".RDS"))
  
  CI.Full.Normali <- CI.Full.Normal %>%
    filter(Sample %in% colnames(TrainSeti)) %>%
    dplyr::select(Sample, Source, primary_diagnosis, Treatment)
  
  CI.Full.Normali <- CI.Full.Normali[match(colnames(TrainSeti), 
                                           CI.Full.Normali$Sample),] %>%
    dplyr::select(-Sample)
  
  TrainSeti.Norm <- norm.A(TrainSeti, annot.Global, CI.Full.Normali)
  
  saveRDS(TrainSeti.Norm, paste0("NormData/", "NormTrainSet_", i, ".RDS"))
  
  after.pca <- prcomp(t(TrainSeti.Norm),
                      center = TRUE,scale. = TRUE)
  ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali$primary_diagnosis)
  ggsave(paste0("PCAs/", "PCA_Norm_TrainSet_", i, "_Diagnosis.pdf"))
  ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali$Source)
  ggsave(paste0("PCAs/", "PCA_Norm_TrainSet_", i, "_DataBase.pdf"))
  
  before.pca <- prcomp(t(TrainSeti),
                       center = TRUE,scale. = TRUE)
  ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali$primary_diagnosis)
  ggsave(paste0("PCAs/", "PCA_Raw_TrainSet_", i, "_Diagnosis.pdf"))
  ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali$Source)
  ggsave(paste0("PCAs/", "PCA_Raw_TrainSet_", i, "_DataBase.pdf"))
  ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali$Treatment)
  ggsave(paste0("PCAs/", "PCA_Raw_TrainSet_", i, "_Treatment.pdf"))
  
  after.pca <- prcomp(t(TrainSeti.Norm[,!(colnames(TrainSeti.Norm) %in% NormalSamples)]),
                      center = TRUE,scale. = TRUE)
  ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali[CI.Full.Normali$primary_diagnosis == "DLBCL",]$Treatment)
  ggsave(paste0("PCAs/", "PCA_Norm_CancerTrainSet_", i, "_Treatment.pdf"))
  ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali[CI.Full.Normali$primary_diagnosis == "DLBCL",]$Source)
  ggsave(paste0("PCAs/", "PCA_Norm_CancerTrainSet_", i, "_DataBase.pdf"))
  
  before.pca <- prcomp(t(TrainSeti[,!(colnames(TrainSeti) %in% NormalSamples)]),
                       center = TRUE,scale. = TRUE)
  ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali[CI.Full.Normali$primary_diagnosis == "DLBCL",]$Treatment)
  ggsave(paste0("PCAs/", "PCA_Raw_CancerTrainSet_", i, "_Treatment.pdf"))
  ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, 
           groups=CI.Full.Normali[CI.Full.Normali$primary_diagnosis == "DLBCL",]$Source)
  ggsave(paste0("PCAs/", "PCA_Raw_CancerTrainSet_", i, "_DataBase.pdf"))
  
}

future::plan(multisession, workers = 10)
results <- future_lapply(1:10,
                         Norm.Partitions,
                         future.seed = TRUE)
plan(sequential)

#### Get Networks p < 1e-06 ####
library(dplyr)
library(tidyr)
library(Hmisc)
library(future.apply)

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

CompleteCI <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CompleteCI <- CompleteCI %>% mutate(primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CompleteCI, CI.Normal)

par_Sp_Corr <- function(expr_m, x){
  
  gene1 <- rownames(expr_m)[x]
  coexpr_list <- list()  
  
  for(y in 1:nrow(expr_m)){
    
    gene2 <- rownames(expr_m)[y]
    
    if(y != x && gene1 < gene2){
      
      x_expr = expr_m[x, ]
      y_expr = expr_m[y, ]
      
      net <- rcorr(x_expr, y_expr, type = "spearman")
      
      if(net$P[1,2] < 1e-06){
        
        coexpr_list[[length(coexpr_list) + 1]] <- c(gene1, gene2, net$r[1,2], net$P[1,2], abs(net$r[1,2]),
                                                    paste(pmin(gene1,gene2),pmax(gene1,gene2),sep="_"))
      }
    }
  }
  
  if(length(coexpr_list) > 0){
    coexpr <- do.call(rbind, coexpr_list)
    colnames(coexpr) <- c("Source", "Target", "Sp_corr", "p_value", "abs", "ID")  
    return(as.data.frame(coexpr, stringsAsFactors = FALSE))
  }
}

ExprMs <- list.files("NormData")

future::plan(multisession, workers = 60)

dir.create("GCNs")

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist %>% as.vector

GCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "GA101") %>%
  dplyr::select(Sample) %>% unlist %>% as.vector

for(i in ExprMs){
  print(i)
  
  DLBCL <- readRDS(paste0("NormData/", i))
  results.RCHOP <- future_lapply(1:nrow(DLBCL), 
                                 function(x) par_Sp_Corr(DLBCL[,colnames(DLBCL) %in% c(RCHOP.samples, NormalSamples)], x),
                                 future.seed = TRUE)
  GCN.RCHOP <- do.call(rbind, results.RCHOP)
  GCN.RCHOP <- GCN.RCHOP %>% as.data.frame() %>% arrange(desc(abs))
  write.csv(GCN.RCHOP, paste0("GCNs/", "GCN_RCHOP_TrainSet_", i, ".csv"), 
            quote = FALSE, row.names = FALSE)
  
  results.GCHOP <- future_lapply(1:nrow(DLBCL), 
                                 function(x) par_Sp_Corr(DLBCL[,colnames(DLBCL) %in% c(GCHOP.samples, NormalSamples)], x),
                                 future.seed = TRUE)
  GCN.GCHOP <- do.call(rbind, results.GCHOP)
  GCN.GCHOP <- GCN.GCHOP %>% as.data.frame() %>% arrange(desc(abs))
  write.csv(GCN.GCHOP, paste0("GCNs/", "GCN_GCHOP_TrainSet_", i, ".csv"), 
            quote = FALSE, row.names = FALSE)
  
}
plan(sequential)

GCNs <- list.files("GCNs")

dir.create("SSNs")

future::plan(multisession, workers = 60)
Lioness_PerSource.R <- function(j) {
  
  Source <- RCHOP.Cancer_GCN[j,"Source"] 
  Target <- RCHOP.Cancer_GCN[j,"Target"]
  
  ID <- RCHOP.Cancer_GCN[j,"ID"] 
  
  alpha <- RCHOP.Cancer_GCN[j,"Sp_corr"] 
  
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
Lioness_PerSource.G <- function(j) {
  
  Source <- GCHOP.Cancer_GCN[j,"Source"] 
  Target <- GCHOP.Cancer_GCN[j,"Target"]
  
  ID <- GCHOP.Cancer_GCN[j,"ID"] 
  
  alpha <- GCHOP.Cancer_GCN[j,"Sp_corr"] 
  
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

for(i in 1:10){
  
  print(i)
  
  Cancer_ExprM <- readRDS(paste0("NormData/NormTrainSet_", i, ".RDS"))
  
  #rownames(Lioness_DLBCL_NCCI) <- gsub("-", ".", rownames(Lioness_DLBCL_NCCI))
  
  RCHOP.Cancer_GCN <- read.csv(paste0("GCNs/GCN_RCHOP_TrainSet_NormTrainSet_", i, ".RDS.csv"))
  GCHOP.Cancer_GCN <- read.csv(paste0("GCNs/GCN_GCHOP_TrainSet_NormTrainSet_", i, ".RDS.csv"))
  #Cancer_GCN <- Cancer_GCN %>% filter(p_value < 1e-08)
  
  RCHOP.genes <- unique(c(RCHOP.Cancer_GCN$Source, RCHOP.Cancer_GCN$Target))
  GCHOP.genes <- unique(c(GCHOP.Cancer_GCN$Source, GCHOP.Cancer_GCN$Target))
  
  RCHOP.Cancer_ExprM <- Cancer_ExprM[rownames(Cancer_ExprM) %in% RCHOP.genes,
                                     colnames(Cancer_ExprM) %in% 
                                       c(RCHOP.samples, NormalSamples)]
  GCHOP.Cancer_ExprM <- Cancer_ExprM[rownames(Cancer_ExprM) %in% GCHOP.genes,
                                     colnames(Cancer_ExprM) %in% 
                                       c(GCHOP.samples, NormalSamples)]
  
  SSN.R.Res <- future_lapply(1:nrow(RCHOP.Cancer_GCN), FUN = Lioness_PerSource.R,
                             future.seed = TRUE)
  SSN.R <- bind_rows(SSN.R.Res)
  write.csv(SSN.R %>% dplyr::select(-Edge_ID), paste0("SSNs/SSN_RCHOP_TrainSet_", i, ".csv"), 
            quote = FALSE, row.names = T)
  
  SSN.G.Res <- future_lapply(1:nrow(GCHOP.Cancer_GCN), FUN = Lioness_PerSource.G,
                             future.seed = TRUE)
  SSN.G <- bind_rows(SSN.G.Res)
  write.csv(SSN.G %>% dplyr::select(-Edge_ID), paste0("SSNs/SSN_GCHOP_TrainSet_", i, ".csv"), 
            quote = FALSE, row.names = T)
  
}
future::plan(sequential)

#### Univariate selection ####

library(randomForest)
library(dplyr)
library(tidyr)
library(M3C)
library(ggsurvfit)
library(survival)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(exactRankTests)
library(future.apply)
library(glmnet)

NormalSamples <- c("BLGSP-71-19-99988", "BLGSP-71-19-99989", "BLGSP-71-19-99996",    
                   "BLGSP-71-19-99997", "BLGSP-71-19-99998", "BLGSP-71-19-99999")

CI.Global <- readRDS(file.path(ClinicalInfo, "CompleteCI_5yrCensored.RDS"))
CI.Global$PFS_status[CI.Global$PFS_time > 1825] <- 0
CI.Global$PFS_time[CI.Global$PFS_time > 1825] <- 1825

CI.Global$OS_status[CI.Global$OS_time > 1825] <- 0
CI.Global$OS_time[CI.Global$OS_time > 1825] <- 1825

# CI.GCHOP <- CI.Global %>% filter(Treatment == "GA101")
# CI.RCHOP <- CI.Global %>% filter(Treatment == "Rituximab")

dir.create("SummaryKMs")

options(future.globals.maxSize = 12 * 1024^3)
getSummary <- function(i, matrix){
  
  edge <- colnames(matrix)[i]
  
  #!!edge := matrix[,i],
  
  SurvCI.i <- tibble(EdgeWeight = matrix[,i],
                     Sample = rownames(matrix)) %>%
    inner_join(CI.Global, by = c("Sample" = "Sample")) %>%
    dplyr::select(-Sample, -COO, -Gender, -Treatment, -Source, -IPI, -Genetic_subtype)
  
  Cox.PFS <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                   data = SurvCI.i %>% dplyr::select(-OS_time, -OS_status))
  
  Cox.OS <- coxph(Surv(OS_time, OS_status) ~ ., 
                  data = SurvCI.i %>% dplyr::select(-PFS_time, -PFS_status))
  
  SurvCI.i <- SurvCI.i %>%
    mutate(EdgeWeight = case_when(EdgeWeight > 0 ~ "Positive",
                                  EdgeWeight < 0 ~ "Negative",
                                  .default = "Neutral"))
  
  Cox.PFS.Cat <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                       data = SurvCI.i %>% dplyr::select(-OS_time, -OS_status))
  
  Cox.OS.Cat <- coxph(Surv(OS_time, OS_status) ~ ., 
                      data = SurvCI.i %>% dplyr::select(-PFS_time, -PFS_status))
  
  if(length(table(SurvCI.i$EdgeWeight)) > 1){
    
    pval.pfs <- survdiff(Surv(PFS_time, PFS_status) ~ EdgeWeight, data = SurvCI.i)$pvalue
    pval.OS <- survdiff(Surv(OS_time, OS_status) ~ EdgeWeight, data = SurvCI.i)$pvalue
    
    fit.pfs <- survfit(Surv(PFS_time, PFS_status) ~ EdgeWeight, data = SurvCI.i)
    summary_fit.pfs <- summary(fit.pfs)
    survival_by_group.pfs <- tapply(summary_fit.pfs$surv, 
                                    summary_fit.pfs$strata, function(x) tail(x, 1))
    names(survival_by_group.pfs) <- gsub("EdgeWeight=", "",  names(survival_by_group.pfs))
    
    fit.OS <- survfit(Surv(OS_time, OS_status) ~ EdgeWeight, data = SurvCI.i)
    summary_fit.OS <- summary(fit.OS)
    survival_by_group.OS <- tapply(summary_fit.OS$surv, 
                                   summary_fit.OS$strata, function(x) tail(x, 1))
    names(survival_by_group.OS) <- gsub("EdgeWeight=", "",  names(survival_by_group.OS))
    
    Summary.KM <- data.frame(Edge = edge,
                             PFS.KM.pval = pval.pfs,
                             OS.KM.pval = pval.OS,
                             Concordance.PFS = summary(Cox.PFS)$concordance[1],
                             Concordance.OS = summary(Cox.OS)$concordance[1],
                             Concordance.PFS.Cat = summary(Cox.PFS.Cat)$concordance[1],
                             Concordance.OS.Cat = summary(Cox.OS.Cat)$concordance[1],
                             GoodSurvGroup.PFS = names(survival_by_group.pfs[which.max(survival_by_group.pfs)]),
                             GoodSurvGroup.OS = names(survival_by_group.OS[which.max(survival_by_group.OS)]))
    
  } else {
    Summary.KM <- data.frame(Edge = edge,
                             PFS.KM.pval = NA,
                             OS.KM.pval = NA,
                             Concordance.PFS = summary(Cox.PFS)$concordance[1],
                             Concordance.OS = summary(Cox.OS)$concordance[1],
                             Concordance.PFS.Cat = summary(Cox.PFS.Cat)$concordance[1],
                             Concordance.OS.Cat = summary(Cox.OS.Cat)$concordance[1],
                             GoodSurvGroup.PFS = NA,
                             GoodSurvGroup.OS = NA)
  }
  return(Summary.KM)
}

future::plan(multisession, workers = 50)
for(RP in 1:10){
  
  print(RP)
  
  #Load RCHOP
  Lioness.RCHOP <- data.table::fread(paste0("SSNs/SSN_RCHOP_TrainSet_", RP, ".csv"))
  Lioness.RCHOP <- as.data.frame(Lioness.RCHOP)
  rownames(Lioness.RCHOP) <- Lioness.RCHOP[,1]
  Lioness.RCHOP <- Lioness.RCHOP[,-1]
  rownames(Lioness.RCHOP) <- gsub("-", ".", rownames(Lioness.RCHOP))
  Lioness.RCHOP <- Lioness.RCHOP[,!(colnames(Lioness.RCHOP) %in% NormalSamples)]
  Lioness.RCHOP.scale <- as.data.frame(scale(t(Lioness.RCHOP)))
  
  #Load GCHOP
  Lioness.GCHOP <- data.table::fread(paste0("SSNs/SSN_GCHOP_TrainSet_", RP, ".csv"))
  Lioness.GCHOP <- as.data.frame(Lioness.GCHOP)
  rownames(Lioness.GCHOP) <- Lioness.GCHOP[,1]
  Lioness.GCHOP <- Lioness.GCHOP[,-1]
  rownames(Lioness.GCHOP) <- gsub("-", ".", rownames(Lioness.GCHOP))
  Lioness.GCHOP <- Lioness.GCHOP[,!(colnames(Lioness.GCHOP) %in% NormalSamples)]
  Lioness.GCHOP.scale <- as.data.frame(scale(t(Lioness.GCHOP)))
  
  results.RCHOP <- future_lapply(1:ncol(Lioness.RCHOP.scale), 
                                 function(x) getSummary(matrix = Lioness.RCHOP.scale, x),
                                 future.seed = TRUE)
  Summary.KM.RCHOP <- do.call(rbind, results.RCHOP)
  saveRDS(tibble(Summary.KM.RCHOP), paste0("SummaryKMs/SummaryKM_RCHOP_RP", RP,".RDS"))
  
  results.GCHOP <- future_lapply(1:ncol(Lioness.GCHOP.scale), 
                                 function(x) getSummary(matrix = Lioness.GCHOP.scale, x),
                                 future.seed = TRUE)
  Summary.KM.GCHOP <- do.call(rbind, results.GCHOP)
  saveRDS(tibble(Summary.KM.GCHOP), paste0("SummaryKMs/SummaryKM_GCHOP_RP", RP,".RDS"))
}
future::plan(sequential)

#### Select Edges ####

library(randomForest)
library(dplyr)
library(tidyr)
library(M3C)
library(ggsurvfit)
library(survival)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(exactRankTests)
library(future.apply)
library(glmnet)

results <- "UnivariateAnalysis"
dir.create(results)                           

for(i in 1:10){
  
  Summ.i <- readRDS(paste0("SummaryKMs/SummaryKM_GCHOP_RP", i,".RDS"))
  
  Summ.i <- Summ.i %>% filter(PFS.KM.pval < 0.05 | Concordance.PFS >= 0.6|
                                OS.KM.pval < 0.05 | Concordance.OS >= 0.6)
  
  if(i == 1){
    Summ.G <- Summ.i
  } else {
    Summ.G <- rbind(Summ.G, Summ.i)
  }
  
}

GCHOP.Edges <- Summ.G %>% group_by(Edge) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>%
  filter(Count >= 6) %>% dplyr::select(Edge) %>% unlist

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

RCHOP.Edges <- Summ.R %>% group_by(Edge) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>%
  filter(Count >= 6) %>% dplyr::select(Edge) %>% unlist

write.csv( Summ.R %>% group_by(Edge) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>%
             filter(Count >= 6) %>% group_by(Edge) %>% mutate(Source = stringr::str_split(Edge, "_")[[1]][1],
                                                              Target = stringr::str_split(Edge, "_")[[1]][2]),
           file.path(results, "Univariable_Analysis_Stable_Network.csv"), quote = F, row.names = F)

relevantEdges <- list("GCHOP" = GCHOP.Edges, "RCHOP" = RCHOP.Edges)

saveRDS(relevantEdges, file.path(results,"relevantEdges_UnivariateFilter.RDS"))

#### Complete GCHOP stable network across random partitions ####

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
CI.Full <- CI.Full %>% mutate(Source = case_when(Source == "GOYA" ~ "GOYA",
                                                 .default = "GDC"),
                              primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist

GCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "GA101") %>%
  dplyr::select(Sample) %>% unlist

ExprMs <- "NormData/"
SSNs <- "SSNs/"

relevantEdges <- file.path(results,"relevantEdges_UnivariateFilter.RDS")

relevantEdges <- relevantEdges$GCHOP

GCN <- data.frame(Edge = relevantEdges, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"_")[[1]][1], 
                Target = stringr::str_split(Edge,"_")[[1]][2])

IterateGCN_Lioness <- function(i) {
  
  TARGET_NBM_ALALi_norm <- readRDS(paste0(ExprMs, "NormTrainSet_", i, ".RDS"))
  rownames(TARGET_NBM_ALALi_norm) <- gsub("-", ".", rownames(TARGET_NBM_ALALi_norm))
  
  TARGET_NBM_ALALi_norm <- TARGET_NBM_ALALi_norm[,colnames(TARGET_NBM_ALALi_norm) %in% c(GCHOP.samples, NormalSamples)]
  
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
  
  write.csv(SSN_i, paste0(SSNs,"SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), quote = F, row.names = T)
  
}

library(future.apply)
future::plan(multisession, workers = 5)
future_lapply(1:10, FUN = IterateGCN_Lioness, future.seed = TRUE)

future::plan(sequential)

#### Complete RCHOP stable network across random partitions ####

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

CI.Full <- file.path(results,"relevantEdges_UnivariateFilter.RDS")
CI.Full <- CI.Full %>% mutate(Source = case_when(Source == "GOYA" ~ "GOYA",
                                                 .default = "GDC"),
                              primary_diagnosis = "DLBCL")

CI.Full.Normal <- bind_rows(CI.Full, CI.Normal)

RCHOP.samples <- CI.Full.Normal %>% filter(Treatment == "Rituximab") %>%
  dplyr::select(Sample) %>% unlist

ExprMs <- "NormData/"
SSNs <- "SSNs/"

relevantEdges <- file.path(results,"relevantEdges_UnivariateFilter.RDS")

relevantEdges <- relevantEdges$RCHOP

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
  
  write.csv(SSN_i, paste0(SSNs,"SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), quote = F, row.names = T)
  
}

library(future.apply)
future::plan(multisession, workers = 5)
future_lapply(1:10, FUN = IterateGCN_Lioness, future.seed = TRUE)

future::plan(sequential)

#### Elastic net par RCHOP ####

library(dplyr)
library(tidyr)
library(survival)
library(future.apply)
library(glmnet)
library(MASS)

NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                   "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")

CI.Global <- file.path(results,"relevantEdges_UnivariateFilter.RDS")
CI.Global$PFS_status[CI.Global$PFS_time > 1825] <- 0
CI.Global$PFS_time[CI.Global$PFS_time > 1825] <- 1825

CI.RCHOP <- CI.Global %>% filter(Treatment == "Rituximab")
CI.GCHOP <- CI.Global %>% filter(Treatment == "GA101")

ElasticNet.PFS.OS <- function(i){
  
  #print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
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
  
  selected_edges <- coef(cv_fit, s = "lambda.1se") 
  
  selected_edges_df <- as.data.frame(as.matrix(selected_edges))
  
  selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
  colnames(selected_edges_df)[1] <- "lambda.1se"
  
  selected_edges_df <- selected_edges_df %>% filter(lambda.1se != 0) %>% arrange(abs(lambda.1se))
  
  selected_edges_df.PFS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "PFS")
  
  lambda_min.PFS <- cv_fit$lambda.1se
  
  mean_cvm.PFS <- min(cv_fit$cvm)
  
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
  
  selected_edges <- coef(cv_fit, s = "lambda.1se") 
  
  selected_edges_df <- as.data.frame(as.matrix(selected_edges))
  
  selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
  colnames(selected_edges_df)[1] <- "lambda.1se"
  
  selected_edges_df <- selected_edges_df %>% filter(lambda.1se != 0) %>% arrange(abs(lambda.1se))
  
  selected_edges_df.OS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "OS")
  
  selected_edges_df <- rbind(selected_edges_df.OS, selected_edges_df.PFS)
  
  lambda_min.OS <- cv_fit$lambda.1se
  
  mean_cvm.OS <- min(cv_fit$cvm)
  
  results.EN <- data.frame(
    Partition = paste0("RP", i),
    LambdaMin.OS = lambda_min.OS,
    MeanCVError.OS = mean_cvm.OS,
    LambdaMin.PFS = lambda_min.PFS,
    MeanCVError.PFS = mean_cvm.PFS)
  
  return(list(Edges = selected_edges_df, EN = results.EN))
  
}

ElasticNet.PFS.OS.min <- function(i){
  
  #print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
  
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
  
  lambda_min.PFS <- cv_fit$lambda.min
  
  mean_cvm.PFS <- min(cv_fit$cvm)
  
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
  
  lambda_min.OS <- cv_fit$lambda.min
  
  mean_cvm.OS <- min(cv_fit$cvm)
  
  results.EN <- data.frame(
    Partition = paste0("RP", i),
    LambdaMin.OS = lambda_min.OS,
    MeanCVError.OS = mean_cvm.OS,
    LambdaMin.PFS = lambda_min.PFS,
    MeanCVError.PFS = mean_cvm.PFS)
  
  return(list(Edges = selected_edges_df, EN = results.EN))
  
}

future::plan(multisession, workers = 10)
FirstFilter <- future_lapply(1:10, 
                             ElasticNet.PFS.OS.min,
                             future.seed = TRUE)
FirstFilter.Edges <- do.call(rbind, lapply(FirstFilter, function(x) x$Edges))
FirstFilter.Results <- do.call(rbind, lapply(FirstFilter, function(x) x$EN))

ElasticNet_results.dir <- "ElasticNet_results"
dir.create(ElasticNet_results.dir)                                    

saveRDS(FirstFilter.Edges, file.path(ElasticNet_results.dir, "ElasticNet_Edges_RCHOP.RDS"))
saveRDS(FirstFilter.Results, file.path(ElasticNet_results.dir, "ElasticNet_Results_RCHOP.RDS"))

plan(sequential)

Edges.R <- FirstFilter.Edges %>% group_by(Edge) %>% distinct(RandomSet) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  %>% filter(Count >= 7) %>% dplyr::select(Edge) %>% unlist %>% as.vector

length(Edges.R)
# [1] 68

for(i in 1:10){
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 1) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 1)
  }
  
}

aic.summ1 <- aic.summ

coeff.summary1 <- coeff.summary

Edges.Filt1 <- coeff.summary1 %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
  distinct(RP) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector

length(which(Edges.R %in% Edges.Filt1)) == length(Edges.R)
#[1] FALSE

length(Edges.Filt1)
#[1] 42

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 2) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 2)
  }
  
}

aic.summ2 <- aic.summ

coeff.summary2 <- coeff.summary

Edges.Filt2 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
  distinct(RP) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector

length(which(Edges.Filt1 %in% Edges.Filt2)) == length(Edges.Filt1)
#[1] FALSE

length(Edges.Filt2)
#[1] 38

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 3) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 3)
  }
  
}

aic.summ3 <- aic.summ

coeff.summary3 <- coeff.summary

Edges.Filt3 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>%
  group_by(Edge) %>% 
  distinct(RP) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector

length(which(Edges.Filt2 %in% Edges.Filt3)) == length(Edges.Filt2)
#[1] FALSE

length(Edges.Filt3)
#[1] 36

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 3) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 3)
  }
  
}

aic.summ4 <- aic.summ

coeff.summary4 <- coeff.summary

Edges.Filt4 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>%
  group_by(Edge) %>% 
  distinct(RP) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector

length(which(Edges.Filt3 %in% Edges.Filt4)) == length(Edges.Filt3)
#[1] TRUE

length(Edges.Filt4)
#[1] 36

#No more edges can be filtered down
#Let's keep only those relevant both for PFS and OS

FinalEdges <-  coeff.summary %>%
  filter(p_value < 0.05) %>%
  dplyr::select(Edge, Type) %>%
  group_by(Edge, Type) %>%
  summarise(Count = n()) %>%
  group_by(Edge)  %>% filter(Count >= 8) %>%
  filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
  ungroup() %>% pull(Edge) %>% unique %>% as.vector

length(FinalEdges)
#[1] 16

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 4) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 4)
  }
  
}

aic.Check1 <- aic.summ

CheckFinalEdges <- coeff.summary %>%
  filter(p_value < 0.05) %>%
  dplyr::select(Edge, Type) %>%
  group_by(Edge, Type) %>%
  summarise(Count = n()) %>%
  group_by(Edge)  %>% filter(Count >= 8) %>%
  filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
  ungroup() %>% pull(Edge) %>% unique %>% as.vector

length(which(FinalEdges %in% CheckFinalEdges)) == length(FinalEdges)
#FALSE

length(CheckFinalEdges)
#[1] 14

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
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
  
  coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
  coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
  
  coeff2.PFS <- coeff2.PFS[,c(1,5)]
  colnames(coeff2.PFS) <- c("Coefficient", "p_value")
  coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                      RP = i, Type = "PFS")
  
  coeff2.OS <- coeff2.OS[,c(1,5)]
  colnames(coeff2.OS) <- c("Coefficient", "p_value")
  coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                    RP = i, Type = "OS")
  
  if(i == 1){
    coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
  } else {
    coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
  }
  
  if(i == 1){
    aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                           AIC2.PFS = AIC(Cox.PFS.2),
                           AIC1.OS = AIC(Cox.OS.1),
                           AIC2.OS = AIC(Cox.OS.2),
                           Partition = i,
                           Step = 4) 
  } else {
    aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                      AIC(Cox.OS.1), AIC(Cox.OS.2),
                                      i, 4)
  }
  
}

aic.Check2 <- aic.summ

CheckFinalEdges2 <- coeff.summary %>%
  filter(p_value < 0.05) %>%
  dplyr::select(Edge, Type) %>%
  group_by(Edge, Type) %>%
  summarise(Count = n()) %>%
  group_by(Edge)  %>% filter(Count >= 8) %>%
  filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
  ungroup() %>% pull(Edge) %>% unique %>% as.vector

length(which(CheckFinalEdges %in% CheckFinalEdges2)) == length(CheckFinalEdges)
#TRUE

length(CheckFinalEdges2)
#[1] 14

saveRDS(CheckFinalEdges2, file.path(ElasticNet_results.dir, "RCHOP_SelectedEdges.RDS"))

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_RCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges2,]
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
  
  print(concordance(Cox.PFS.1)$concordance)
  print(concordance(Cox.OS.1)$concordance)
  
  if(i == 1){
    summC <- data.frame(C.PFS = concordance(Cox.PFS.1)$concordance,
                        C.OS = concordance(Cox.OS.1)$concordance, 
                        RP = i)
  } else {
    summC[nrow(summC) + 1, ] <- c(concordance(Cox.PFS.1)$concordance,
                                  concordance(Cox.OS.1)$concordance, 
                                  i)
  }
  
  
}

saveRDS(summC, file.path(ElasticNet_results.dir, "summC_RCHOP_Train.RDS"))

        
        #### Elastic net par GCHOP ####
        
        library(dplyr)
        library(tidyr)
        library(survival)
        library(future.apply)
        library(glmnet)
        library(MASS)
        
        NormalSamples <- c("BLGSP.71.19.99988", "BLGSP.71.19.99989", "BLGSP.71.19.99996",    
                           "BLGSP.71.19.99997", "BLGSP.71.19.99998", "BLGSP.71.19.99999")
        
        CI.Global <- file.path(results,"relevantEdges_UnivariateFilter.RDS")
        CI.Global$PFS_status[CI.Global$PFS_time > 1825] <- 0
        CI.Global$PFS_time[CI.Global$PFS_time > 1825] <- 1825
        
        CI.GCHOP <- CI.Global %>% filter(Treatment == "GA101")
        
        ElasticNet.PFS.OS <- function(i){
          
          #print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, Source), by = "Sample")
          
          X <- as.matrix(CoxSurvCI %>% dplyr::select(-PFS_time,  -PFS_status,
                                                     -Sample, -Source))
          y <- Surv(CoxSurvCI$PFS_time, CoxSurvCI$PFS_status)
          
          set.seed(123)
          cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 0.5)
          
          selected_edges <- coef(cv_fit, s = "lambda.1se") 
          
          selected_edges_df <- as.data.frame(as.matrix(selected_edges))
          
          selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
          colnames(selected_edges_df)[1] <- "lambda.1se"
          
          selected_edges_df <- selected_edges_df %>% filter(lambda.1se != 0) %>% arrange(abs(lambda.1se))
          
          selected_edges_df.PFS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "PFS")
          
          lambda_min.PFS <- cv_fit$lambda.1se
          
          mean_cvm.PFS <- min(cv_fit$cvm)
          
          #Now OS
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, OS_status, OS_time, Source), by = "Sample")
          
          X <- as.matrix(CoxSurvCI %>% dplyr::select(-OS_time,  -OS_status,
                                                     -Sample, -Source))
          y <- Surv(CoxSurvCI$OS_time, CoxSurvCI$OS_status)
          
          set.seed(123)
          cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 0.5)
          
          selected_edges <- coef(cv_fit, s = "lambda.1se") 
          
          selected_edges_df <- as.data.frame(as.matrix(selected_edges))
          
          selected_edges_df <- selected_edges_df %>% mutate(Edge = rownames(selected_edges_df))
          colnames(selected_edges_df)[1] <- "lambda.1se"
          
          selected_edges_df <- selected_edges_df %>% filter(lambda.1se != 0) %>% arrange(abs(lambda.1se))
          
          selected_edges_df.OS <- tibble(selected_edges_df) %>% mutate(RandomSet = i, Type = "OS")
          
          selected_edges_df <- rbind(selected_edges_df.OS, selected_edges_df.PFS)
          
          lambda_min.OS <- cv_fit$lambda.1se
          
          mean_cvm.OS <- min(cv_fit$cvm)
          
          results.EN <- data.frame(
            Partition = paste0("RP", i),
            LambdaMin.OS = lambda_min.OS,
            MeanCVError.OS = mean_cvm.OS,
            LambdaMin.PFS = lambda_min.PFS,
            MeanCVError.PFS = mean_cvm.PFS)
          
          return(list(Edges = selected_edges_df, EN = results.EN))
          
        }
        
        ElasticNet.PFS.OS.min <- function(i){
          
          #print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          #Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% relevantEdges$GCHOP,]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
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
          
          lambda_min.PFS <- cv_fit$lambda.min
          
          mean_cvm.PFS <- min(cv_fit$cvm)
          
          #Now OS
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
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
          
          lambda_min.OS <- cv_fit$lambda.min
          
          mean_cvm.OS <- min(cv_fit$cvm)
          
          results.EN <- data.frame(
            Partition = paste0("RP", i),
            LambdaMin.OS = lambda_min.OS,
            MeanCVError.OS = mean_cvm.OS,
            LambdaMin.PFS = lambda_min.PFS,
            MeanCVError.PFS = mean_cvm.PFS)
          
          return(list(Edges = selected_edges_df, EN = results.EN))
          
        }
        
        future::plan(multisession, workers = 10)
        FirstFilter <- future_lapply(1:10, 
                                     ElasticNet.PFS.OS.min,
                                     future.seed = TRUE)
        FirstFilter.Edges <- do.call(rbind, lapply(FirstFilter, function(x) x$Edges))
        FirstFilter.Results <- do.call(rbind, lapply(FirstFilter, function(x) x$EN))
        
        saveRDS(FirstFilter.Edges, file.path(ElasticNet_results.dir, "ElasticNet_Edges_GCHOP.RDS"))
        saveRDS(FirstFilter.Results, file.path(ElasticNet_results.dir, "ElasticNet_Results_GCHOP.RDS"))
        
        plan(sequential)
        
        Edges.R <- FirstFilter.Edges %>% group_by(Edge) %>% distinct(RandomSet) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% dplyr::select(Edge) %>% unlist %>% as.vector
        
        length(Edges.R)
        #[1] 48
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.R,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 1) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 1)
          }
          
        }
        
        aic.summ1 <- aic.summ
        
        coeff.summary1 <- coeff.summary
        
        Edges.Filt1 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        length(which(Edges.R %in% Edges.Filt1)) == length(Edges.R)
        #[1] FALSE
        
        length(Edges.Filt1)
        #[1] 33
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt1,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 2) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 2)
          }
          
        }
        
        aic.summ2 <- aic.summ
        
        Edges.Filt2 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        coeff.summary2 <- coeff.summary
        
        length(which(Edges.Filt1 %in% Edges.Filt2)) == length(Edges.Filt1)
        #[1] FALSE
        
        length(Edges.Filt2)
        #[1] 24
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt2,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 3) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 3)
          }
          
        }
        
        aic.summ3 <- aic.summ
        
        Edges.Filt3 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        coeff.summary3 <- coeff.summary
        
        length(which(Edges.Filt2 %in% Edges.Filt3)) == length(Edges.Filt2)
        #[1] FALSE
        
        length(Edges.Filt3)
        #[1] 22
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt3,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 4) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 4)
          }
          
        }
        
        aic.summ4 <- aic.summ
        
        Edges.Filt4 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        length(which(Edges.Filt3 %in% Edges.Filt4)) == length(Edges.Filt3)
        #[1] FALSE
        
        length(Edges.Filt4)
        #[1] 18
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt4,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 4) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 4)
          }
          
        }
        
        aic.summ5 <- aic.summ
        
        Edges.Filt5 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        length(which(Edges.Filt4 %in% Edges.Filt5)) == length(Edges.Filt4)
        #[1] FALSE
        
        length(Edges.Filt5)
        #[1] 17
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% Edges.Filt5,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 4) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 4)
          }
          
        }
        
        aic.summ6 <- aic.summ
        
        Edges.Filt6 <- coeff.summary %>% filter(p_value < 0.05) %>% dplyr::select(Edge, RP) %>% group_by(Edge) %>% 
          distinct(RP) %>%
          summarise(Count = n()) %>%
          arrange(desc(Count))  %>% filter(Count >= 7) %>% pull(Edge) %>% as.vector
        
        length(which(Edges.Filt5 %in% Edges.Filt6)) == length(Edges.Filt5)
        #[1] TRUE
        
        length(Edges.Filt6)
        #[1] 17
        
        #Since now, all of the edges are selected at least in seven partitions
        #either in the OS or the PFS model, let's keep now only those selected in both
        #models at least 7 times
        
        FinalEdges <-  coeff.summary %>%
          filter(p_value < 0.05) %>%
          dplyr::select(Edge, Type) %>%
          group_by(Edge, Type) %>%
          summarise(Count = n()) %>%
          group_by(Edge)  %>% filter(Count >= 8) %>%
          filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
          ungroup() %>% pull(Edge) %>% unique %>% as.vector
        
        length(FinalEdges)
        #[1] 10
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% FinalEdges,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 5) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 5)
          }
          
        }
        
        CheckFinalEdges <-  coeff.summary %>%
          filter(p_value < 0.05) %>%
          dplyr::select(Edge, Type) %>%
          group_by(Edge, Type) %>%
          summarise(Count = n()) %>%
          group_by(Edge)  %>% filter(Count >= 8) %>%
          filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
          ungroup() %>% pull(Edge) %>% unique %>% as.vector
        
        length(which(FinalEdges %in% CheckFinalEdges)) == length(FinalEdges)
        #[1] FALSE
        
        length(CheckFinalEdges)
        #[1] 7
        
        for(i in 1:10){
          
          print(i)
          
          Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges,]
          Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[,!(colnames(Lioness_DLBCL_NCCI) %in% NormalSamples)]
          
          Lioness_DLBCL_NCCI.scale <- as.data.frame(scale(t(Lioness_DLBCL_NCCI)))
          
          CoxSurvCI <- Lioness_DLBCL_NCCI.scale %>% mutate(Sample = rownames(Lioness_DLBCL_NCCI.scale))
          
          CoxSurvCI <- CoxSurvCI %>% 
            inner_join(CI.GCHOP %>% 
                         dplyr::select(Sample, PFS_status, PFS_time, OS_time, OS_status), by = "Sample")
          
          Cox.PFS.1 <- coxph(Surv(PFS_time, PFS_status) ~ ., 
                             data = CoxSurvCI %>% dplyr::select(-Sample, -OS_status,-OS_time))
          Cox.PFS.2 <- stepAIC(Cox.PFS.1)
          
          Cox.OS.1 <- coxph(Surv(OS_time, OS_status) ~ ., 
                            data = CoxSurvCI %>% dplyr::select(-Sample, -PFS_status,-PFS_time))
          Cox.OS.2 <- stepAIC(Cox.OS.1)
          
          coeff2.PFS <- as.data.frame(summary(Cox.PFS.2)$coefficient)
          coeff2.OS <- as.data.frame(summary(Cox.OS.2)$coefficient)
          
          coeff2.PFS <- coeff2.PFS[,c(1,5)]
          colnames(coeff2.PFS) <- c("Coefficient", "p_value")
          coeff2.PFS <- coeff2.PFS %>% mutate(Edge = rownames(coeff2.PFS),
                                              RP = i, Type = "PFS")
          
          coeff2.OS <- coeff2.OS[,c(1,5)]
          colnames(coeff2.OS) <- c("Coefficient", "p_value")
          coeff2.OS <- coeff2.OS %>% mutate(Edge = rownames(coeff2.OS),
                                            RP = i, Type = "OS")
          
          if(i == 1){
            coeff.summary <- rbind(coeff2.PFS, coeff2.OS)
          } else {
            coeff.summary <- rbind(coeff.summary, coeff2.PFS, coeff2.OS)
          }
          
          if(i == 1){
            aic.summ <- data.frame(AIC1.PFS = AIC(Cox.PFS.1),
                                   AIC2.PFS = AIC(Cox.PFS.2),
                                   AIC1.OS = AIC(Cox.OS.1),
                                   AIC2.OS = AIC(Cox.OS.2),
                                   Partition = i,
                                   Step = 6) 
          } else {
            aic.summ[nrow(aic.summ)+1, ] <- c(AIC(Cox.PFS.1),AIC(Cox.PFS.2),
                                              AIC(Cox.OS.1), AIC(Cox.OS.2),
                                              i, 6)
          }
          
        }
        
        CheckFinalEdges2 <-  coeff.summary %>%
          filter(p_value < 0.05) %>%
          dplyr::select(Edge, Type) %>%
          group_by(Edge, Type) %>%
          summarise(Count = n()) %>%
          group_by(Edge)  %>% filter(Count >= 8) %>%
          filter(n_distinct(Type) == 2) %>% # Keep edges relevant in OS and PFS
          ungroup() %>% pull(Edge) %>% unique %>% as.vector
        
        length(which(CheckFinalEdges %in% CheckFinalEdges2)) == length(CheckFinalEdges)
        #[1] TRUE
        
        length(CheckFinalEdges2)
        #[1] 7
        
        saveRDS(CheckFinalEdges2, file.path(ElasticNet_results.dir, "GCHOP_SelectedEdges.RDS"))

for(i in 1:10){
  
  print(i)
  
  Lioness_DLBCL_NCCI <- read.csv(paste0("SSNs/SSN_Selected_GCHOP_UnivariateFilter_TrainSet_", i, ".csv"), row.names = 1)
  Lioness_DLBCL_NCCI <- Lioness_DLBCL_NCCI[rownames(Lioness_DLBCL_NCCI) %in% CheckFinalEdges2,]
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
  
  print(concordance(Cox.PFS.1)$concordance)
  print(concordance(Cox.OS.1)$concordance)
  
  if(i == 1){
    summC <- data.frame(C.PFS = concordance(Cox.PFS.1)$concordance,
                        C.OS = concordance(Cox.OS.1)$concordance, 
                        RP = i)
  } else {
    summC[nrow(summC) + 1, ] <- c(concordance(Cox.PFS.1)$concordance,
                                  concordance(Cox.OS.1)$concordance, 
                                  i)
  }
  
  
}

saveRDS(summC, file.path(ElasticNet_results.dir, "summC_GCHOP_Train.RDS"))
                                                     
        
