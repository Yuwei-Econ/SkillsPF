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

### Additional specifications
#source("Full simulations/No measures/Additional specifications/Spec_A 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_A 3mix 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_B 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_B 3mix 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_C 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_C 3mix 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_D 260216.R")
#source("Full simulations/No measures/Additional specifications/Spec_D 3mix 260216.R")

### Additional assumptions on starting Xs
#source("Full simulations/3 mix/3mix Spec_1 210316.R") # have disabled the part omitting local optima
#source("Full simulations/3 mix/3mix Spec_1 3mix 210316.R") # have disabled the part omitting local optima
#source("Full simulations/chi square/chisq Spec_1 210316.R")
#source("Full simulations/chi square/chisq Spec_1 3mix 210316.R")
#source("Full simulations/No measures/Spec_1 1mix 210316.R") # have disabled the part omitting local optima

### Assumption of a single normal
#source("Full simulations/No measures/Spec_1 1mix 210316.R") # have disabled the part omitting local optima


time_tot <- proc.time() - time_start
