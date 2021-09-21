#Sampling importance resampling (SIR) using deviance for poisson distribution

#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))

pop1_data_agg <- readRDS(system.file("data/example_data.RDS", package = "HIVepisimAnalysis"))



compute_log_importance_weight <- function ( diag_obs, diag_sim )
{
	diag_sim <- diag_sim$newDx_pop1
	#subset to remove the NA in observed data (data not available fro 1981 and 2021)
	diag_sim <- diag_sim[3:41]
	diag_obs <- diag_obs[3:41]
	
	sum( dpois( diag_sim, lambda = diag_obs , log = TRUE ) )
}
log_weights <- plyr::ddply(pop1_data_agg, "rep_param", compute_log_importance_weight,
               diag_obs = incidenceDiag$frequency)
w =   exp( log_weights[,2]  - max( log_weights[,2] ) )
resample <- sample( as.character( log_weights[, 'rep_param']  ), prob = w, replace = TRUE ) 



# this is effectively sampling only one trajectory so let's look at taht alongside the real data 
pop1_data_agg$year <- as.numeric(  pop1_data_agg$year ) + 1980
sims = split( pop1_data_agg , pop1_data_agg$rep_param )
Y = do.call( cbind, lapply( sims , '[[', 'newDx_pop1' )  )
matplot( sims[[1]]$year, Y , type = 'l', col = 'grey', lwd=.5, ylim = c( 0, max(na.omit( c(pop1_data_agg$newDx_pop1,  incidenceDiag$frequency)) ) ) )
lines( sims[[ resample[1] ]]$year , sims[[ resample[1] ]]$newDx_pop1, col = 'red', lwd = 2 )
with( incidenceDiag, points( time, frequency, pty = 20, col = 'black', cex = 2 ) )

