###################### Flicker Noise Generation #####################
##### This script generates and saves flicker noise simulations #####
##### for use in FlickerNoise_noGaps.R and FlickerNoise_withGaps.R ##
#####################################################################


setwd("/home/cmb15/ClockDataAnalysis/Code/Paper1")
#  

################## libraries ###############################

library(RobPer) #flicker noise generation function TK95()

############################################################

numberOfSimulations = 300
N = 2048

####### Simulating with no gaps #########

X.t_sims_flk <- matrix(NA, nrow = numberOfSimulations, ncol = N)

for(i in 1:numberOfSimulations){
  set.seed(i)
  X.t_sims_flk[i,] <-  arfima.sim(N,model = list(dfrac = 0.49))
}


####### Simulating with gaps #########
NwithGaps <- 600 + N
  
X.t_sims_flk_gaps <- matrix(NA, nrow = numberOfSimulations, ncol = NwithGaps)

for(i in 1:numberOfSimulations){
  set.seed(i + 300)
  X.t_sims_flk_gaps[i,] <-  TK95(N = NwithGaps, alpha = 1)
}


saveRDS(X.t_sims_flk,paste("Results/FlickerSims_noGaps","_N",N,"_",numberOfSimulations,".Rds",sep=""))
saveRDS(X.t_sims_flk_gaps,paste("Results/FlickerSims_withGaps","_N",NwithGaps,"_",numberOfSimulations,".Rds",sep=""))



