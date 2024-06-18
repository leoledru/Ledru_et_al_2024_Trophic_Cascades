setwd("~/Documents/Post_Doc_CARRTEL/R")

library(NetIndices) # trophic levels and omnivory
library(latex2exp)
library(plotrix)
library(expm)
library(sparsevar)
library(RColorBrewer)
library(igraph)
library(matlib)
library(rlang)
library(ggplot2)
library(ggallin)
library(ggpmisc)
library(reshape2)
library(scales)
library(fdrtool)
library(viridis)
library(extrafont)
library(tidyverse)
library(rjson)
library(dplyr)
library(gridExtra)
library(docstring)
library(roxygen2)
source("FoodWebFormat.R")
source("FigEmpiricalWebs_b.R")


########################################################
# Figure 4 & A.9 : analysis of lakes empirical food webs 
########################################################

############################################################################
# Compute for each lake the "average web" (from the 150 "jacobians for each)
# and the trophic level of each species
# This step can be avoided by loading directly the average webs below 

load("Lakes_StabilityRangeTot.Rdata") # contains the inferred "jacobians"
LakeNames <- names(Lakes_Stability)
LakeNames[5] <- "Chasseforet"
LakeNames[11] <- "Merlet Superieur"
LakeNames[15] <- "Batailleres"
LakeNames[16] <- "Pepin"
names(Lakes_Stability) <- LakeNames
FishlessNames <- c("Roche Ferran","Chasseforet","Arpont","Noir du Carro", "Blanc du Carro", "Vallette")
SalmoNames <- c("Blanc de l'Archeboc", "Pepin", "Mont Coua", "Verdet de l'Archeboc", "Noir de l'Archeboc")
SalmoForageNames <- c("Noir de Forclaz", "Merlet Superieur", "Batailleres", "Noir de La Masse", "Lou", "Cerces")
ForageNames <- c("Rond des Rochilles")
# Get the 150 "jacobians" from each lake
JacobsAllLakes <- list()
for (i in LakeNames){
  JacobsAllLakes[[i]] <- Lakes_Stability[[i]][["JacobMatBS"]]
}
# Mean interaction matrix for each lake
LakesDirect <- list()
for (i in 1:length(JacobsAllLakes)){ 
  JacobsOneLake <- JacobsAllLakes[[i]]
  # define matrix of mean direct interactions (mean over the 150 replicas)
  A <- JacobsOneLake[[1]]
  for (j in 1:length(JacobsOneLake)){
    A <- A + JacobsOneLake[[j]]
  }
  A <- A / length(JacobsOneLake)
  LakesDirect[[i]] <- A
}
names(LakesDirect) <- names(JacobsAllLakes)
# Compute trophic levels and keep only lakes with MaxTroph >= 3
# check if max trophic level is >= 3
idx <- list()
for (i in 1:length(LakesDirect)){
  A <- LakesDirect[[i]]
  flow <- -1 * A * (A < 0)
  flow[is.na(flow)] <- 0
  Troph <- TrophInd(flow)
  Trophic <- round(Troph$TL)
  MaxTroph <- max(Trophic)
  if (MaxTroph < 3){
    idx <- append(idx, i)
  }
}
LakesDirect[unlist(idx)] <- NULL
save(LakesDirect, file = "LakesDirect.RData")

LakesDirectSalmo <- LakesDirect[c(SalmoNames,SalmoForageNames)]
save(LakesDirectSalmo, file = "LakesDirectSalmo.RData")

######################################################
# Avoid step 1 by loading the average web of each lake

load("LakesDirect.RData")
load("LakesDirectSalmo.RData")

############################
# Analysis of the lakes webs

source("FigEmpiricalWebs_b.R")
docstring(FigEmpiricalWebs_b)
#debug(FigEmpiricalWebs_b)
MaxInter <- 2
MinInter <- 0.01
NbrStep <- 50
EmpiricalPlot <- FigEmpiricalWebs_b(LakesDirectSalmo, MaxInter, MinInter, NbrStep)
EmpiricalPlot[[1]]
EmpiricalPlot[[2]]
EmpiricalPlot[[3]]
EmpiricalPlot[[4]]
EmpiricalPlot[[5]]


#####################################################
# Figure A.10 : analysis of Bascompte lab's food webs
#####################################################
source("FoodWebFormat.R")
# debug(FoodWebFormat)
base_url <- "https://www.web-of-life.es/"
json_url <- paste0(base_url,"get_networks.php?interaction_type=FoodWebs&data_type=weighted")
FoodWebs <- jsonlite::fromJSON(json_url)
# correction of the 1st food-web which is inverted in the way it is annotated compared to the others
FoodWebs_correct <- FoodWebs
FoodWebs_correct$species1[FoodWebs_correct$network_name == "FW_001"] <- FoodWebs$species2[FoodWebs$network_name == "FW_001"]
FoodWebs_correct$species2[FoodWebs_correct$network_name == "FW_001"] <- FoodWebs$species1[FoodWebs$network_name == "FW_001"]
FoodWebs <- FoodWebs_correct
# format the food webs
AllFwbs <- FoodWebFormat(FoodWebs, ConvCoeff = 0.5)
FwbsNames <- unique(FoodWebs$network_name)
FwbsInfos <- read.csv(paste0(base_url,
                             paste0("get_network_info.php?network_name=",FwbsNames[[1]])))
FwbsSpeciesInfos <- list(read.csv(paste0(base_url,
                                         paste0("get_species_info.php?network_name=",FwbsNames[[1]]))))
for (i in 2:length(FwbsNames)){
  FwbsInfos <- rbind(FwbsInfos, read.csv(paste0(base_url,
                                                paste0("get_network_info.php?network_name=",FwbsNames[[i]]))))
  FwbsSpeciesInfos[[i]] <- read.csv(paste0(base_url,
                                           paste0("get_species_info.php?network_name=",FwbsNames[[i]])))
}
AllFwbs <- Filter(Negate(is.null), AllFwbs)
MaxInter <- 2
MinInter <- 0.01
NbrStep <- 50
BascomptePlot <- FigEmpiricalWebs_b(AllFwbs, MaxInter, MinInter, NbrStep)
BascomptePlot[[1]]
BascomptePlot[[2]]
BascomptePlot[[3]]
BascomptePlot[[4]]
BascomptePlot[[5]]

