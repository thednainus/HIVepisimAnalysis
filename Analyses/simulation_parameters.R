#10th March 2022

#create dataframe to run simulations in the cluster
#each line will be a combination of parameters to run the simulations

# I have 100 replicates using the same combinations of parameter values
# and I will analyse different percentages of the population
# from 5%, 10% to 90%.

params <- c(rep(1067, 1000), rep(2348, 1000))

percentage <- c(rep(0.05, 100), rep(0.1, 100),
                rep(0.2, 100), rep(0.3, 100),
                rep(0.4, 100), rep(0.5, 100),
                rep(0.6, 100), rep(0.7, 100),
                rep(0.8, 100), rep(0.9, 100))
replicate <- rep(1:100)

simulations <- data.frame(params = params, perc = rep(percentage, 2), rep = rep(replicate, 20))

saveRDS(simulations, "inst/data/simulatins_imperial_clsuter.RDS")


#15th March 2022

#create dataframe to run simulations in the cluster (XSEDE)
#each line will be a combination of parameters to run the simulations

# I have 100 replicates using the same combinations of parameter values
# and I will analyse different percentages of the population
# from 5%, 10% to 90%.

params <- c(rep(1067, 200), rep(2348, 200))

percentage <- c(rep(0.05, 100), rep(0.3, 100))
replicate <- rep(1:100)

simulations <- data.frame(params = params, perc = rep(percentage, 2),
                          rep = rep(replicate, 4))

saveRDS(simulations, "inst/data/simulations_deepseq_cluster.RDS")


#test run in order to have very few sequences
params <- c(rep(1067, 2), rep(2348, 2))

percentage <- c(rep(0.01, 4))
replicate <- rep(1:2)

simulations_test <- data.frame(params = params, perc = percentage,
                          rep = rep(replicate, 2))

saveRDS(simulations, "inst/data/simulations_deepseq_cluster_test.RDS")


