rm(list = ls())

# library -----------------------------------------------------------------

library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)

today <- as.character(Sys.Date())
today <- gsub("-", x=today, "")

# Set WD ------------------------------------------------------------------

office="maison"

if (office=="CIRAD"){
  root <- "D:/Mes Donnees/"
} else {
  root <- "C:/Users/cedri/"
}

setwd(paste0(root,
             "OneDrive/A_GenomicPrediction/PredictionASReml/GS_VL/S02S03_rerun"))

# Load files --------------------------------------------------------
Pred_cmpt <- NULL
for (method in c("GBLUP")){
  for (trait in c("FL", "PH", "YLD14", "ZN")){
    for (CV in c("CV2_site","IMBran", "IMBcdm")){
      for (s in c(25, 50, 100, 200)){
        if (CV == "CV2_site"){
          s <- 0.3
        }

        S02 <- read.csv(paste0(trait, "_", method, "/Pred_S02_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S02$s <- s
        S02$Predictor <- "S02"
        S02$ID <- paste0(S02$trait, S02$LOC, S02$DNAID, S02$Predictor)
        
        S03 <- read.csv(paste0(trait, "_", method, "/Pred_S03_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S03$s <- s
        S03$Predictor <- "S03"
        S03$ID <- paste0(S03$trait, S03$LOC, S03$DNAID, S03$Predictor)
        
        Pred_cmpt <- rbind(Pred_cmpt, S02, S03)
        # Pred_list <- list(S0=noGEN, S02=S02, S03=S03)
        
      }
    }
  }
}

for (method in c("RKHS")){
  for (trait in c("FL", "PH", "YLD14", "ZN")){
    for (CV in c("CV2_site","IMBran")){#, "IMBcdm")){
      for (s in c(25, 50, 100, 200)){
        if (CV == "CV2_site"){
          s <- 0.3
        }
        
        S02 <- read.csv(paste0(trait, "_", method, "/Pred_S02_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S02$s <- s
        S02$Predictor <- "S02"
        S02$ID <- paste0(S02$trait, S02$LOC, S02$DNAID, S02$Predictor)
        
        S03 <- read.csv(paste0(trait, "_", method, "/Pred_S03_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S03$s <- s
        S03$Predictor <- "S03"
        S03$ID <- paste0(S03$trait, S03$LOC, S03$DNAID, S03$Predictor)
        
        Pred_cmpt <- rbind(Pred_cmpt, S02, S03)
        # Pred_list <- list(S0=noGEN, S02=S02, S03=S03)
        
      }
    }
  }
}

Pred_cmpt$Predictor <- as.factor(Pred_cmpt$Predictor)

write.csv(Pred_cmpt, paste0("Prediction_cmpt_", today, ".csv"),
row.names = F)

# Pred_cmpt <- read.csv("Prediction_cmpt_20200721.csv", header = TRUE)
# Pred_cmpt_wide <- reshape(Pred_cmpt,
#                      direction = "wide",
#                      v.names = c("ref", "pred"),
#                      idvar = c(1,2,3,4,5,6,9),
#                      timevar = "Predictor")



# load reference ----------------------------------------------------------

S02 <- read.csv(paste0(root,
                       "OneDrive/AnalysePheno/PCT27/Output_cluster/", 
                       "S02_REFERENCE_on334_20200722.csv"),
                   header = TRUE)
colnames(S02)[5] <- paste0("S02_", colnames(S02)[5])
S02$ID <- paste0(S02$Trait, S02$LOC, S02$DNAID, "S02")

S03 <- read.csv(paste0(root,
                       "OneDrive/AnalysePheno/PCT27/Output_cluster/", 
                       "S03_REFERENCE_20200722.csv"),
                   header = TRUE)
colnames(S03)[5] <- paste0("S03_", colnames(S03)[5])
S03$ID <- paste0(S03$Trait, S03$LOC, S03$DNAID, "S03")

# Complete prediction with reference --------------------------------------

Pred_ref_S02 <- merge(Pred_cmpt[Pred_cmpt$Predictor=="S02",], S02[,c(1,5)], by="ID")
Pred_ref_S03 <- merge(Pred_cmpt[Pred_cmpt$Predictor=="S03",], S03[,c(1,5)], by="ID")

Pred_list <- list(S02=Pred_ref_S02, S03=Pred_ref_S03)

# write.csv(Pred_ref, paste0("Prediction_cmpt_REF_", date, ".csv"),
# row.names = F)

# Compute general accuracy ------------------------------------------------
summary_gen <- NULL
summary_Acc <- NULL
for (predictor in c("S02", "S03")){
  Pred_ref <- Pred_list[[predictor]]
  Acc_general <- ddply(Pred_ref, .(trait, CV, s, method, LOC, CVround), function(x) cor(x[,9],
                                                                                        x[,12], method = "pearson"))
  Acc_general$Predictor <- predictor
  summary_Acc <- rbind(summary_Acc, Acc_general)
  
  summary_Acc_gen<- ddply(Acc_general, .(trait, CV, method, s, LOC),
                          summarize,
                          mean = mean(V1, na.rm = T), 
                          sd = sd(V1, na.rm = T))
  summary_Acc_gen$Predictor <- predictor
  
  summary_gen <- rbind(summary_gen, summary_Acc_gen)
}

write.csv(summary_Acc, paste0("Accuracies_", today, ".csv"),
          row.names = F,
          quote = F)

write.csv(summary_gen, paste0("summary_Accuracy_", today, ".csv"),
          row.names = F)

