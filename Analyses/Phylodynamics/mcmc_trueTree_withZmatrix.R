# Run MCMC on phylogenetic trees for the Ethics project

library(BayesianTools)
library(phydynR)
library(ape)
library(pangeaZA)
library(HIVepisim)
library(HIVepisimAnalysis)

seed <- sample(1:1000000, 1)
print(seed)
set.seed(seed)

path_to_gene <- "../FitTrajectories/BESTFITsmallpop_combination_params2/small_msm_pop_sim_284/rep_20/cluster_results/50perc/output/vts/results_vts_merged_trees_sampling_0.059.tre"
path_to_sts <- "../FitTrajectories/BESTFITsmallpop_combination_params2/small_msm_pop_sim_284/rep_20/cluster_results/50perc/output/vts/W/merged_trees_sampling_migrant_years_1_simple_0.059.RData"
BDTR <- prepare_dated_trees(path_to_gene = path_to_gene, path_to_sts = path_to_sts)

likelihood = 'PL2'

mod0 = generate_model0_net()
# alpha2021: rate of ART initiation
# The function alpha(t) is given by a simple 1-parameter linear function which
# increases from zero starting in 2005.
ESTNAMES <- c( 'srcGrowth', 'src1980', 'i1980', 'importRate', 'alpha2021', mod0$betaNames)
ELOWER <-  c(0.01, 1.00, 0.01, 1/100, 1/24 , rep(0.001, length(mod0$betaNames) ) )
#ELOWER <-  c(0.01, 1.00, 250, 1/100, 1/24 , rep(0.001, length(mod0$betaNames) ) )
#ELOWER <-  c(0.01, 1.00, 500, 1/100, 1/24 , rep(0.001, length(mod0$betaNames) ) )
#EUPPER <-  c(0.10,  1e5, 99.9,   1/5, 1/1. , rep(2, length(mod0$betaNames) ) )
EUPPER <-  c(0.50,  1e5, 1500,   1/2, 1/1. , rep(4, length(mod0$betaNames) ) )

priordensity <- function(theta){
  #browser()
  names(theta) <- ESTNAMES
  o = dlnorm(theta['srcGrowth'], log(.035), 1 , log = TRUE ) +
    dlnorm(theta['src1980'], log(2e4), 1/4 , log = TRUE ) +
    dexp(theta['i1980'], 1/1005 , log = TRUE ) +
    dlnorm(theta['importRate'], log(1/20), 1/4 , log = TRUE ) +
    sum( dlnorm(theta[mod0$betaNames], log(2/12), 1 , log = TRUE )  )
  unname( o )
}

sampler <- function(n = 1){

  srcGrowth <- rlnorm ( n, log(.035), 1/8 )
  src1980 <- rlnorm( n, log(2e4), 1/2)
  #i1980 <- rexp ( n, 1/10 )
  #based on start conditions of simulations
  i1980 <- rexp ( n, 1/1005 )
  importRate <- rlnorm( n , log( 1/20 ), 1/4 )
  alpha2021 <- runif( n, 1/8, 1/2)
  o = cbind( srcGrowth, src1980, i1980 , importRate, alpha2021,
             matrix( rlnorm( n*length(mod0$betaNames), log(2/12), 1/4), nrow = n) )
  for ( i in 1:n){
    o[i,] <- pmax( ELOWER, o[i,] )
    o[i,] <- pmin( EUPPER, o[i,] )
  }
  o
}

prior <- createPrior(density = priordensity,
                     sampler = sampler,
                     lower = ELOWER ,
                     upper = EUPPER
)

co0 <- function(theta)
{
  #browser()
  names(theta) <- ESTNAMES
  print(theta)
  p <- mod0$parms
  p[ names(theta)] <- theta
  x0 <- c(src = p$src1980, i = p$i1980 , z = 0)
  s <-  mod0$dm(p, x0 , t0 = 1980, t1 = BDTR$maxSampleTime, res = 100 )
  y1 <- tail( s[[5]], 1 )[1,]
  zprior <- dnorm( y1['z'] / (y1['z'] + y1['i'] ) , .5, sd = .05, log = TRUE )
  #~ browser()
  o = suppressWarnings( colik( BDTR, p , mod0$dm, x0 , t0 = 1980, res = 100, maxHeight = BDTR$maxSampleTime - 1980  , likelihood = likelihood) )
  print(o + zprior )
  o + zprior
}

# After we had a run to create a z-matrix we did the follow:
# Read a previous run for creating starting values for the Z matrix
runZ <- readRDS("output_run2_newPriors_trueTree.RDS")
#runZ <- readRDS("output_run3_forcingPriotfor_i1980.RDS")

# Get a good sample (the run above is not good, however it can provide a good Z matrix)
# For more information on this: https://github.com/florianhartig/BayesianTools/issues/79
x <- getSample(runZ, start=3000)

# Get the range for the parameter estimates for the previous run
rangePost = apply(x, 2, range)

#get unique values of x
u_x <- unique(x)

#cretae new Z matrix based on previous run
# now I am estimating 11 parameters (hence ncol=11)
newZ = matrix(runif(1650, rangePost[1,], rangePost[2,]), ncol = 11, byrow = T)

# Because I will run several analysis in parallel, and to avoid the initial values to be identical
# I will provide as argument position 1 (pos1), position 2 (pos2), and position 3 (pos3)
# from the unique values of x (u_x)
pos1=sample(1:length(u_x[,1]), 1)
pos2=sample(1:length(u_x[,1]), 1)
pos3=sample(1:length(u_x[,1]), 1)


#Implementing linear chain
load("iter.rdata")

while(i < 2001){
  if(!file.exists("output_zmatrix.RDS")){
    bayesianSetup3 <- createBayesianSetup(likelihood = co0 , prior = prior, names = ESTNAMES , parallel = FALSE )
    settings3 = list(iterations = 12, nrChains = 1, thin = 1)

    out <- runMCMC(bayesianSetup = bayesianSetup3, sampler = "DEzs", settings = settings3)
    saveRDS(out, "output_zmatrix.RDS")
    save(i, file="iter.rdata")
    i = i + 1
  }else{
    out <- readRDS("output_zmatrix.RDS")
    out1 <- out
    out <- runMCMC(bayesianSetup = out1)
    saveRDS(out, "output_zmatrix.RDS")
    save(i, file="iter.rdata")
    i = i + 1
  }
}



out <- readRDS("output_zmatrix.RDS")
out_sample_smatriz <- BayesianTools::getSample(out, start = 1, coda = TRUE)
plot(out_sample_smatriz, ask = TRUE)




