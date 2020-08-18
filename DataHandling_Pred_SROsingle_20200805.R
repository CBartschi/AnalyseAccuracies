rm(list = ls())

# library -----------------------------------------------------------------

library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)

today <- as.character(Sys.Date())
today <- gsub("-", x=today, "")

# Set WD ------------------------------------------------------------------
office="CIRAD"

if (office=="CIRAD"){
  setwd("~/OneDrive/A_GenomicPrediction/PredictionASReml/GS_VL/SROsingle")
} else {
  setwd("C:/Users/cedri/OneDrive/A_GenomicPrediction/PredictionASReml/GS_VL/SROsingle")
}

getwd()

# Load files --------------------------------------------------------
Pred_cmpt <- NULL
for (method in c("GBLUP")){
  for (trait in c("FL", "PH", "YLD14", "ZN")){
    for (CV in c("IMBcdm", "IMBran")){
      for (s in c(25, 50, 100, 200, 234)){
        if (CV == "CV2_site"){
          s <- 0.3
        }

        S02 <- read.csv(paste0(trait, "_", method, "_SROsingle/Pred_S02_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S02$s <- s
        S02$Predictor <- "S02"
        S02$ID <- paste0(S02$trait, "SR", S02$DNAID, S02$Predictor)
        
        S03 <- read.csv(paste0(trait, "_", method, "_SROsingle/Pred_S03_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S03$s <- s
        S03$Predictor <- "S03"
        S03$ID <- paste0(S03$trait, "SR", S03$DNAID, S03$Predictor)
        
        Pred_cmpt <- rbind(Pred_cmpt, S02, S03)
        # Pred_list <- list(S0=noGEN, S02=S02, S03=S03)
        
      }
    }
  }
}

for (method in c("RKHS")){
  for (trait in c("FL", "PH", "ZN", "YLD14")){
    for (CV in c("IMBran", "IMBcdm")){
      for (s in c(25, 50, 100, 200)){
        if (CV == "CV2_site"){
          s <- 0.3
        }

        S02 <- read.csv(paste0(trait, "_", method, "_SROsingle/Pred_S02_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S02$s <- s
        S02$Predictor <- "S02"
        S02$ID <- paste0(S02$trait, "SR", S02$DNAID, S02$Predictor)

        S03 <- read.csv(paste0(trait, "_", method, "_SROsingle/Pred_S03_", trait, "_", method, "_", CV, "_", s, ".csv"))
        S03$s <- s
        S03$Predictor <- "S03"
        S03$ID <- paste0(S03$trait, "SR", S03$DNAID, S03$Predictor)

        Pred_cmpt <- rbind(Pred_cmpt, S02, S03)
        # Pred_list <- list(S0=noGEN, S02=S02, S03=S03)

      }
    }
  }
}

Pred_cmpt$Predictor <- as.factor(Pred_cmpt$Predictor) 

summary(Pred_cmpt)

# write.csv(Pred_cmpt, paste0("Prediction_SROsingle_", today, ".csv"),
# row.names = F)

# Pred_cmpt <- read.csv("Prediction_cmpt_20200721.csv", header = TRUE)
# Pred_cmpt_wide <- reshape(Pred_cmpt,
#                      direction = "wide",
#                      v.names = c("ref", "pred"),
#                      idvar = c(1,2,3,4,5,6,9),
#                      timevar = "Predictor")



# load reference ----------------------------------------------------------

S02 <- read.csv("~/OneDrive/AnalysePheno/PCT27/Output_cluster/S02_REFERENCE_on384_20200722.csv",
                   header = TRUE)
colnames(S02)[5] <- paste0("S02_", colnames(S02)[5])
S02$ID <- paste0(S02$Trait, S02$LOC, S02$DNAID, "S02")

S03 <- read.csv("~/OneDrive/AnalysePheno/PCT27/Output_cluster/S03_REFERENCE_20200722.csv",
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
summary_Acc100 <- NULL
summary_scenar100 <- NULL
for (predictor in c("S02", "S03")){
  Pred_ref <- Pred_list[[predictor]]
  summary(Pred_ref)
  
  # on 100 prediction
  Acc_100 <- ddply(Pred_ref, .(trait, CV, s, method, CVround), function(x) {smpl <- sample_n(x, 100, replace = FALSE)
                                                                           cor(smpl[,8],
                                                                               smpl[,11], method = "pearson")})
  Acc_100$Predictor <- predictor
  summary_Acc100 <- rbind(summary_Acc100, Acc_100)
  
  # on all prediction
  Acc_general <- ddply(Pred_ref, .(trait, CV, s, method, CVround), function(x) cor(x[,8],
                                                                                   x[,11], method = "pearson"))
  Acc_general$Predictor <- predictor
  summary_Acc <- rbind(summary_Acc, Acc_general)
  
  # by scenario 100
  tmp_100scnar <- ddply(Acc_100, .(trait, CV, method, s),
                          summarize,
                          mean = mean(V1, na.rm = T), 
                          sd = sd(V1, na.rm = T))
  tmp_100scnar$Predictor <- predictor
  
  summary_scenar100 <- rbind(summary_scenar100, tmp_100scnar)
  
  # by scenario all
  summary_Acc_gen<- ddply(Acc_general, .(trait, CV, method, s),
                          summarize,
                          mean = mean(V1, na.rm = T), 
                          sd = sd(V1, na.rm = T))
  summary_Acc_gen$Predictor <- predictor
  
  summary_gen <- rbind(summary_gen, summary_Acc_gen)
}
# 100
write.csv(summary_Acc100, paste0("Accuracies100_SROsingle", today, ".csv"),
          row.names = F,
          quote = F)

write.csv(summary_scenar100, paste0("summary_Accuracy100_SROsingle", today, ".csv"),
          row.names = F)

# all
write.csv(summary_Acc, paste0("Accuracies_SROsingle", today, ".csv"),
          row.names = F,
          quote = F)

write.csv(summary_gen, paste0("summary_Accuracy_SROsingle", today, ".csv"),
          row.names = F)

