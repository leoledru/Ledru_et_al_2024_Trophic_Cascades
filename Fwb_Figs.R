setwd("~/Documents/Post_Doc_CARRTEL/R")
setwd("~/Documents/Post_Doc_CARRTEL-main/R")
#############################
library(expm)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(expm)
library(sparsevar)
library(rlang)
library(ggplot2)
library(ggallin)
library(reshape2)
library(dplyr)
library(scales)
library(fdrtool)
library(viridis)
library(extrafont)
library(ggpmisc)
library(docstring)
library(roxygen2)

####################################
# Figure 1 : Simple food chain model
####################################
source("ChainAttenuation.R")
docstring(ChainAttenuation)
ChainAttenuation()

################################
# Figures : 2 / A.4 / A.6 / A.7
###############################
source("FigEachChainDirectNet_b.R")
docstring(FigEachChainDirectNet_b)
# Loading previously generated stable networks
load("StableFwbs_TL3_b.Rdata") # trophic level max = 3 (Figure 2 & A.4)
load("StableFwbs_TL4_b.Rdata") # trophic level max = 4 (Figure A.6)
load("StableFwbs_TL3_Structured_b.Rdata") # perfectly structured food webs (Figure A.7)

StableFwbs <- StableFwbs_TL4 # set the correct variable according to the file loaded
# filter too exclude too high collectivity foodwebs
FilterList <- function(sublist){all(sublist$Collect < 2)}
StableFwbs <- Filter(FilterList, StableFwbs)
A_list <- list()
Collect_list <- list()
Troph_list <- list()
Omni_list <- list()
for (i in 1:length(StableFwbs)){
  A_list[[i]] <- StableFwbs[[i]][["A"]]
  Collect_list[[i]] <- StableFwbs[[i]][["Collect"]]
  Troph_list[[i]] <- StableFwbs[[i]][["Troph"]]
  Omni_list[[i]] <- StableFwbs[[i]][["Omni"]]
}

# For networks with maximum trophic level = 3 or 4 and trophic cascade from top to top-2
Gplots <- FigEachChainDirectNet_b(A_list, Collect_list, Troph_list, Omni_list)

#  For networks with maximum trophic level = 4 and trophic cascade from top to bottom
# source("FigEachChainDirectNet_c.R")
# Gplots <- FigEachChainDirectNet_c(A_list, Collect_list, Troph_list, Omni_list)

# Loading saved analysis to see figures without running analysis # 
# load("Gplots_TL3.Rdata")
# load("Gplots_TL4.Rdata")
# load("Gplots_Structured_TL3.Rdata")

#################
# Show the graphs
# Gplots[[1]][["layers"]][[1]][["aes_params"]][["alpha"]] <- 0.2
Gplots[[1]] # Scatter plot of species-cascade diversity (all food chains)
Gplots[[2]] # Barplot of species-cascade diversity
Gplots[[3]] # Percentage of community-cascade inversion
Gplots[[4]] # Percentage of community-cascade inversion Barplot
Gplots[[5]] # Percentage of community-cascade inversion according to omnivory
Gplots[[6]] # Percentage of community-cascade inversion according to omnivory Barplot
Gplots[[7]] # Percentage of species-cascade inversion
Gplots[[8]] # Percentage of species-cascade inversion according to omnivory
Gplots[[9]] # Scatter plot of community-cascade diversity (all food webs)
Gplots[[10]] # Scatter plot of community-cascade diversity (all food webs) according to omnivory
Gplots[[11]] # Barplot of community-cascade diversity
Gplots[[12]] # Types ~ Collect with variance
par(family = "LM Roman 10")
plot(Collect_list, Omni_list, xlab = "Collectivity", ylab = "Omnivory")


###########################################
# Figures : 3 Proxies of cascade divergence
###########################################
source("DivergenceProxy_Figs.R")
docstring(DivergenceProxy_Figs)
Out <- DivergenceProxy_Figs(A_list, Collect_list, Troph_list, Omni_list)
Out[[1]] # Divergence according to food chain integration
Out[[2]] # Divergence according to food chain omnivory
Out[[3]] # Divergence according to food chain intraguild predation


################################################
# Figures A.1 & A.2 : Cascade inversion dynamics
################################################
source("InversionAndPerturb.R")
docstring(InversionAndPerturb)
InversionAndPerturb(tmax, tstep, S, C, Omni = FALSE)
  

################################################
# Figures A.5 : Compensation between net effects
################################################
source("FigEachPredToConsoNet.R")
docstring(FigEachPredToConsoNet)
FigEachPredToConsoNet(A_list, Collect_list, Troph_list)
  
  
  

