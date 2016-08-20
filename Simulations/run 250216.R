########### MIXTURE SIMULATIONS EXPERIMENTS ###########################
### This script calls various simulation scripts to run sequentially
rm(list = ls())

numsims <- 50
conv <- 0.1

time_start <- proc.time()

### No measures
#source("Full simulations/No measures/Spec_1 210216.R")
#source("Full simulations/No measures/Spec_2 210216.R")
#source("Full simulations/No measures/Spec_3 260216.R")

### Realistic distribution
#source("Full simulations/No measures/Spec_4 220216.R")

### With measures
#source("Full simulations/With measures/Spec_5 220216.R")
#source("Full simulations/With measures/Spec_6 220216.R")
#source("Full simulations/With measures/Spec_7 220216.R")

### Mixture of 3 normals
#source("Full simulations/No measures/Spec_2_3mix 220216.R")
#source("Full simulations/No measures/Spec_3_3mix 250216.R")

time_tot <- proc.time() - time_start
