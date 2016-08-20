########### MIXTURE SIMULATIONS EXPERIMENTS ###########################
### This script pulls all simulation results into tables
rm(list = ls())

library("stargazer")
#library("xtable")

######################## TABLE 1: SINGLE NORMAL ################################

load("Full simulations/Output/single normal 070416.RData")
spec_1 <- list(rhoest, alphaest, betaest)
spec_1_data <- as.data.frame(spec_1)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_1_data <- spec_1_data[,!(names(spec_1_data) %in% drops)]
row.names(spec_1_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

stargazer(spec_1_data, summary = F, title = "Simulation evidence assuming a single normal", 
          notes = c("Based on 50 runs with alpha = 0.5 and beta = 0.5"),label = "tab_single_normal", digits = 4)

######################## TABLE 2: NO MEASURES ##################################

load("Full simulations/Output/Spec_1 210216.Rdata")
spec_1 <- list(rhoest, alphaest, betaest)
spec_1_data <- as.data.frame(spec_1)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_1_data <- spec_1_data[,!(names(spec_1_data) %in% drops)]
row.names(spec_1_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_2 210216.Rdata")
spec_2 <- list(rhoest, alphaest, betaest)
spec_2_data <- as.data.frame(spec_2)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_2_data <- spec_2_data[,!(names(spec_2_data) %in% drops)]
row.names(spec_2_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_3 210216.Rdata")
spec_3 <- list(rhoest, alphaest, betaest)
spec_3_data <- as.data.frame(spec_3)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_3_data <- spec_3_data[,!(names(spec_3_data) %in% drops)]
row.names(spec_3_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_data <- rbind(spec_1_data, spec_2_data, spec_3_data)
colnames(spec_data) <- c("rho","MSE(rho)","alpha","MSE(alpha)","beta","MSE(beta)")

# Copy and paste sections of following into Latex for further editing
stargazer(spec_data, summary = F, title = "Simulation evidence abstracting from measures", 
          notes = c("Based on XXX runs with alpha = 0.5 and beta = 0.5",
          "Parameters for each specification given in Appendix X"),label = "tab_sim_1", digits = 4)

###################### TABLE 3: REALISTIC DISTRIBUTION ########################

load("Full simulations/Output/Spec_4 220216.Rdata")
spec_4 <- list(rhoest, alphaest, betaest)
spec_4_data <- as.data.frame(spec_4)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_4_data <- spec_4_data[,!(names(spec_4_data) %in% drops)]
row.names(spec_4_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1")

## Need to add in a realistic distribution from my China estimation here

spec_data <- spec_4_data
colnames(spec_data) <- c("rho","MSE(rho)","alpha","MSE(alpha)","beta","MSE(beta)")

# Copy and paste sections of following into Latex for further editing
stargazer(spec_data, summary = F, title = "Simulation evidence abstracting from measures", 
          notes = c("Based on XXX runs with alpha = 0.5 and beta = 0.5",
                    "Parameters for each specification given in Appendix X"),label = "tab_sim_2", digits = 4)

###################### TABLE 4: INCLUDING MEASURES #############################
## Specification 5 excluded for now as problems running. Will need to update to include MSE if include

#load("Full simulations/Output/Spec_5 220216.Rdata")
#spec_5 <- list(rhoest, alphaest, betaest)
#spec_5_data <- as.data.frame(spec_5)
#drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
#spec_5_data <- spec_5_data[,!(names(spec_5_data) %in% drops)]
#row.names(spec_5_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_6 240216.Rdata")
spec_6 <- list(rhoest, alphaest, betaest)
spec_6_data <- as.data.frame(spec_6)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_6_data <- spec_6_data[,!(names(spec_6_data) %in% drops)]
row.names(spec_6_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_7 240216.Rdata")
spec_7 <- list(rhoest, alphaest, betaest)
spec_7_data <- as.data.frame(spec_7)
drops <- c("rhos","rhosdest","alpha","alphasdest","beta","betasdest")
spec_7_data <- spec_7_data[,!(names(spec_7_data) %in% drops)]
row.names(spec_7_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_data <- rbind(spec_6_data, spec_7_data) #spec_5_data
colnames(spec_data) <- c("rho","MSE(rho)","alpha","MSE(alpha)","beta","MSE(beta)")

# Copy and paste sections of following into Latex for further editing
stargazer(spec_data, summary = F, title = "Simulation evidence including measures", 
          notes = c("Based on XXX runs with alpha = 0.5 and beta = 0.5",
                    "Parameters for each specification given in Appendix X"),label = "tab_sim_3", digits = 4)

####################### TABLE 5: MIXTURE OF 3 NORMALS ##########################

load("Full simulations/Output/Spec_7 270216.Rdata")
spec_7 <- list(rhoest) #, alphaest, betaest)
spec_7_data <- as.data.frame(spec_7)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_7_data <- spec_7_data[,!(names(spec_7_data) %in% drops)]
row.names(spec_7_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_7 3mix 270216.Rdata")
spec_7_3mix <- list(rhoest) #, alphaest, betaest)
spec_7_3mix_data <- as.data.frame(spec_7_3mix)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_7_3mix_data <- spec_7_3mix_data[,!(names(spec_7_3mix_data) %in% drops)]
row.names(spec_7_3mix_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_7_full_data <- cbind(spec_7_data, spec_7_3mix_data)
colnames(spec_7_full_data) <- c("rho","MSE(rho)","rho 3mix","MSE(rho) 3mix") #,"alpha","MSE(alpha)","beta","MSE(beta)")


load("Full simulations/Output/chisq Spec_1 210316.Rdata")
spec_8 <- list(rhoest) #, alphaest, betaest)
spec_8_data <- as.data.frame(spec_8)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_8_data <- spec_8_data[,!(names(spec_8_data) %in% drops)]
row.names(spec_8_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/chisq Spec_1 3mix 210316.Rdata")
spec_8_3mix <- list(rhoest) #, alphaest, betaest)
spec_8_3mix_data <- as.data.frame(spec_8_3mix)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_8_3mix_data <- spec_8_3mix_data[,!(names(spec_8_3mix_data) %in% drops)]
row.names(spec_8_3mix_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_8_full_data <- cbind(spec_8_data, spec_8_3mix_data)
colnames(spec_8_full_data) <- c("rho","MSE(rho)","rho 3mix","MSE(rho) 3mix") #,"alpha","MSE(alpha)","beta","MSE(beta)")


load("Full simulations/Output/3mix Spec_1 210316.Rdata")
spec_9 <- list(rhoest) #, alphaest, betaest)
spec_9_data <- as.data.frame(spec_9)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_9_data <- spec_9_data[,!(names(spec_9_data) %in% drops)]
row.names(spec_9_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/3mix Spec_1 3mix 210316.Rdata")
spec_9_3mix <- list(rhoest) #, alphaest, betaest)
spec_9_3mix_data <- as.data.frame(spec_9_3mix)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_9_3mix_data <- spec_9_3mix_data[,!(names(spec_9_3mix_data) %in% drops)]
row.names(spec_9_3mix_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_9_full_data <- cbind(spec_9_data, spec_9_3mix_data)
colnames(spec_9_full_data) <- c("rho","MSE(rho)","rho 3mix","MSE(rho) 3mix") #,"alpha","MSE(alpha)","beta","MSE(beta)")


load("Full simulations/Output/Spec_10 270216.Rdata")
spec_10 <- list(rhoest) #, alphaest, betaest)
spec_10_data <- as.data.frame(spec_10)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_10_data <- spec_10_data[,!(names(spec_10_data) %in% drops)]
row.names(spec_10_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

load("Full simulations/Output/Spec_10 3mix 270216.Rdata")
spec_10_3mix <- list(rhoest) #, alphaest, betaest)
spec_10_3mix_data <- as.data.frame(spec_10_3mix)
drops <- c("rhos","rhosdest") #,"alpha","alphasdest","beta","betasdest")
spec_10_3mix_data <- spec_10_3mix_data[,!(names(spec_10_3mix_data) %in% drops)]
row.names(spec_10_3mix_data) <- c("rho = -1","rho = -0.5","rho = 0","rho = 0.5","rho = 1" )

spec_10_full_data <- cbind(spec_10_data, spec_10_3mix_data)
colnames(spec_10_full_data) <- c("rho","MSE(rho)","rho 3mix","MSE(rho) 3mix") #,"alpha","MSE(alpha)","beta","MSE(beta)")

### Combining specs 7 - 10

spec_data <- rbind(spec_7_full_data, spec_8_full_data, spec_9_full_data, spec_10_full_data)


# Copy and paste sections of following into Latex for further editing
stargazer(spec_data, summary = F, title = "Simulation evidence assuming a mixture of three normals", 
          notes = c("Based on XXX runs with alpha = 0.5 and beta = 0.5",
                    "Parameters for each specification given in Appendix X"),label = "tab_sim_4", digits = 4)

