# setwd("~/Documents/Post_Doc_CARRTEL/R")
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


####################################
# Figure 1 : Simple food chain model
####################################
source("ChainAttenuation.R")
ChainAttenuation()

################################################
# Figures : 2 / 3 / A.4 / A.5 / A.6 / A.9 / A.10
################################################
source("FigEachChainDirectNet_b.R")
# source("FigEachChainDirectNet_c.R")
# Loading previously generated stable networks
load("StableFwbs_TL3_b.Rdata") # trophic level max = 3
load("StableFwbs_TL4_b.Rdata") # trophic level max = 4
load("StableFwbs_TL3_Structured_b.Rdata") # perfectly structured food webs

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
source("FigEachChainDirectNet_b.R")
Gplots <- FigEachChainDirectNet_b(A_list, Collect_list, Troph_list, Omni_list)
#  For networks with maximum trophic level = 4 and trophic cascade from top to bottom
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
plot(Collect_list, Omni_list, xlab = "Collectivity", ylab = "Omnivory") # Fig A.6


##############################################
# Figures : A.7 Proxies of cascade divergence
##############################################
source("DivergenceProxy_Figs.R")
Out <- DivergenceProxy_Figs(A_list, Collect_list, Troph_list, Omni_list)
Out[[1]] # Divergence according to food chain integration
Out[[2]] # Divergence according to food chain omnivory
Out[[3]] # Divergence according to food chain intraguild predation






