#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White                                                   ##
## Institut Pasteur                                                ##
## michael.white@pasteur.fr                                        ##
## m.white08@imperial.ac.uk                                        ##
##                                                                 ##
#####################################################################
#####################################################################


rm(list=ls())

#setwd("")

## call the libraries that are needed for analysis

library(MASS)
library(compiler)



###############################################
###############################################
##          ##                               ##
##   ####   ##  ####    ####  ######  ####   ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ##  ##  ##  ## ######   ##   ######  ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##   ####   ##  ####   ##  ##   ##   ##  ##  ##
##          ##                               ##
###############################################
###############################################
##



age_SD = 18					## define the age of sexual debut


data <- read.csv("", header = TRUE)		## read the csv file in

head(data)

data <- data[]							## you may need to subset the data by provience or region 

data_2 <- as.data.frame(cbind(data$age, data$conc))

data_3 <- na.omit(data_2) 				## make sure there are no NAs

AB_data <- data_3

AB_cs <- data_3							## assign the object to a new name that will be evalauated by the function

###############################################
## 0.3 Prepare data for plotting

age_bins     = seq(from=0, to=90, by=5)														## age bins will be set from the youngest age group
age_bins_mid = 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )		## in your data to the maximum age by = user defined

N_bins = length(age_bins) - 1 



GMT_bins_cs = rep(NA, N_bins)

AB_range_bins_cs = matrix(NA, nrow=N_bins, ncol=3)											## generate an empty matrix to be filled with antibody levels for each age bin and the binomial CIs
colnames(AB_range_bins_cs) = c("med", "low", "high")

for(i in 1:N_bins)
{
	index = intersect( which(AB_cs[,1]>age_bins[i]), which(AB_cs[,1]<=age_bins[i+1]) ) 
	temp  = AB_cs[index,2]

	GMT_bins_cs[i] = exp(mean(log(temp)))

	AB_range_bins_cs[i,] = quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


###############################################
## 0.4 Plot data
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins_cs, 
pch=15, cex=2,
log="y", xlim=c(0,90), ylim=c(0.01,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Cross-sectional antibody data"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins_cs[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins_cs[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}



###################################################
###################################################
##        ##                                     ##
##   ##   ##  #     #  ####  ####   ##### ##     ##
##  ###   ##  ##   ## ##  ## ## ##  ##    ##     ##
##   ##   ##  ####### ##  ## ##  ## ####  ##     ##
##   ##   ##  ## # ## ##  ## ## ##  ##    ##     ##
##  ####  ##  ##   ##  ####  ####   ##### #####  ##
##        ##                                     ##
###################################################
###################################################


################################################### 
## 1.1 MODEL specifiy the AA model that may have
## generated the data  

rr <- 0.091055 			## define the rate of antibody decay


model_M2 = function(a, t_survey, par_M2)
{
	alpha_0   = par_M2[1]
	gamma     = par_M2[2]
	time_c    = par_M2[3]
	sigma     = par_M2[4]
	alpha_STI = par_M2[5]

	alpha_c = gamma*alpha_0

	age_xx = a - time_c + t_survey
 


	if( (age_SD <= age_xx) && (age_xx <= a) )
	{
		AB_titre = ( alpha_0/rr )*( 1 - exp(-rr*age_SD) )
		
		AB_titre = AB_titre*exp(-rr*(age_xx-age_SD)) + ( (alpha_0+alpha_STI)/rr )*( 1 - exp(-rr*(age_xx-age_SD)) )

		AB_titre = AB_titre*exp(-rr*(a-age_xx)) + ( (alpha_c+alpha_STI)/rr )*( 1 - exp(-rr*(a-age_xx)) )
	}

	if( (age_SD <= a) && (a <= age_xx) )
	{
		AB_titre = ( alpha_0/rr )*( 1 - exp(-rr*age_SD) )
		
		AB_titre = AB_titre*exp(-rr*(a-age_SD)) + ( (alpha_0+alpha_STI)/rr )*( 1 - exp(-rr*(a-age_SD)) )
	}

	if( (age_xx > 0) && (age_xx <= age_SD) && (age_SD <= a) )
	{
		AB_titre = ( alpha_0/rr )*( 1 - exp(-rr*age_xx) )
		
		AB_titre = AB_titre*exp(-rr*(age_SD-age_xx)) + ( (alpha_c)/rr )*( 1 - exp(-rr*(age_SD-age_xx)) )

		AB_titre = AB_titre*exp(-rr*(a-age_SD)) + ( (alpha_c+alpha_STI)/rr )*( 1 - exp(-rr*(a-age_SD)) )
	}

	if( (age_xx > 0) && (age_xx <= a) && (a <= age_SD) )
	{
		AB_titre = ( alpha_0/rr )*( 1 - exp(-rr*age_xx) )
		
		AB_titre = AB_titre*exp(-rr*(a-age_xx)) + ( (alpha_c)/rr )*( 1 - exp(-rr*(a-age_xx)) )
	}	

	if( (age_xx > 0) && (a <= age_xx) && (a <= age_SD) )
	{
		AB_titre = ( alpha_0/rr )*( 1 - exp(-rr*a) )
	}	

	if( (age_xx <= 0) && (a <= age_SD) )
	{
		AB_titre = ( alpha_c/rr )*( 1 - exp(-rr*a) )
	}	

	if( (age_xx <= 0) && (age_SD <= a) )
	{
		AB_titre = ( alpha_c/rr )*( 1 - exp(-rr*age_SD) )

		AB_titre = AB_titre*exp(-rr*(a-age_SD)) + ( (alpha_c+alpha_STI)/rr )*( 1 - exp(-rr*(a-age_SD)) )
	}	


	AB_titre

}

model_M2 = cmpfun(model_M2, options=list(optimize=3)) 


###################################################
## 1.2 define and compute the LIKELIHOOD 


loglike_M2_cs = function( par_M2 )
{
	alpha_0   = par_M2[1]
	gamma     = par_M2[2]
	time_c    = par_M2[3]
	sigma     = par_M2[4]
	alpha_STI = par_M2[5]

	AB_model = sapply(AB_cs[,1], model_M2, t_survey=0, par_M2=par_M2)

	mu = log(AB_model) 

	loglike = -log(AB_cs[,2]) - log(2.506628*sigma) - 0.5*( (log(AB_cs[,2])-mu)/sigma )^2

      sum( loglike )
}

loglike_M2_cs = cmpfun(loglike_M2_cs, options=list(optimize=3)) 



loglike_M2_total = function( par )
{
	loglike_M2_cs( par ) 
}




###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M2 = function( par_M2 )
{
	alpha_0   = par_M2[1]
	gamma     = par_M2[2]
	time_c    = par_M2[3]
	sigma     = par_M2[4]
	alpha_STI = par_M2[5]


	######################################
	## Uniform prior on alpha_0 ~ U(0,1000)

	if( alpha_0>0 && alpha_0<1000 )
	{
		prior_alpha_0 = log(1/1000)
	}else{
		prior_alpha_0 = -LARGE
	}

	######################################
	## Uniform prior on gamma ~ U(0,1)

	if( gamma>0 && gamma<1 )
	{
		prior_gamma = log(1/1)
	}else{
		prior_gamma = -LARGE
	}


	######################################
	## Uniform prior on time_c ~ U(0,60)

	if( time_c>0 && time_c<60 )
	{
		prior_time_c = log(1/60)
	}else{
		prior_time_c = -LARGE
	}

	######################################
	## Uniform prior on sigma ~ U(0,10)

	if( sigma>0 && sigma<10 )
	{
		prior_sigma = log(1/10)
	}else{
		prior_sigma = -LARGE
	}

	######################################
	## Uniform prior on alpha_STI ~ U(0,1000)

	if( alpha_STI>0 && alpha_STI<1000 )
	{
		prior_alpha_STI = log(1/1000)
	}else{
		prior_alpha_STI = -LARGE
	}


	prior <- prior_alpha_0 + prior_gamma + prior_time_c + prior_sigma + prior_alpha_STI

	prior
}

prior_M2 = cmpfun(prior_M2, options=list(optimize=3))


#################################################
#################################################
##          ##                                 ##
##   ####   ##  #     #  ####  #     #  ####   ##
##  ##  ##  ##  ##   ## ##  ## ##   ## ##  ##  ##
##     ##   ##  ####### ##     ####### ##      ##
##    ##    ##  ## # ## ##  ## ## # ## ##  ##  ##
##   #####  ##  ##   ##  ####  ##   ##  ####   ##
##          ##                                 ##
#################################################
#################################################


N_mcmc       = 100000      ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 20000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate


#################################################
## 2.1 Robbins-munro step scaler

rm_scale = function(step_scale, mc, log_prob)
{
	dd = exp(log_prob)
	if( dd < -30 ){ dd <- 0 }
	dd = min( dd, 1 )

	rm_temp = ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
		
	out = step_scale*exp(rm_temp)
	
	out = max( out, 0.02 )
	out = min( out, 2)
	out
}


#################################################
## 2.2 Prepare object for MCMC fitting
 
MCMC_par           = matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par) = c("alpha_0", "gamma", "time_c", "sigma", "alpha_STI", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations


par_MC = c(40, 0.5, 10, 0.6, 10)       ## (alpha_0, gamma, rr, time_c, sigma, alpha_STI)


Sigma_MC = diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

                  
prior_MC   = prior_M2( par_MC )

loglike_MC = loglike_M2_total( par_MC ) + prior_MC




for(mc in 1:N_mcmc)
{
	par_MCp1 = mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)

	prior_MCp1 = prior_M2( par_MCp1 ) 

	if( prior_MCp1 > -0.5*LARGE )
	{
 		loglike_MCp1 = loglike_M2_total( par_MCp1 ) + prior_MCp1


		log_prob = min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC = par_MCp1
			
			loglike_MC  = loglike_MCp1
			prior_MC    = prior_MCp1

			MCMC_accept = MCMC_accept + 1                       
		}

		#######################################
		## RM scaling of proposal step size

		if( mc < N_adapt )
		{
			step_scale = rm_scale( step_scale, mc, log_prob)
		}

		#######################################
		## Adaptive tuning of covariance matrix

		if( (mc > N_tune_start) && (mc < N_tune_end) )
		{
			cov_MC = cov( MCMC_par[1:(mc-1),1:5] )
		}
	}


	MCMC_par[mc,1:5] = par_MC
	MCMC_par[mc,6]   = loglike_MC
	MCMC_par[mc,7]   = prior_MC
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,4))

for(k in 1:5)
{
	#####################################
	## PANEL k: 

	plot(x=1:N_mcmc, y=MCMC_par[,k], 
	pch=19, col="grey", cex=0.25,
	xlab="MCMC iteration", ylab="", 
	main=colnames(MCMC_par)[k])
}


#####################################
## PANEL 7: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,6], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
ylim=quantile( MCMC_par[,6], prob=c(0.01,1)),
main="likelihood" )

	
#####################################
## PANEL 8: prior

plot(x=1:N_mcmc, y=MCMC_par[,7], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="prior", 
main="prior" )


	




#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn = MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,3))


#####################################
## PANEL 1: MCMC posterior

for(k in 1:5)
{
	DEN = density( MCMC_burn[,k] )
	
	QUANT = quantile( MCMC_burn[,k], prob=c(0.025, 0.5, 0.975) )

	plot(x=DEN$x, y=DEN$y, type='l',
	xlim=c(0, max(DEN$x)),
	xlab="alpha_0", ylab="", 
	main=colnames(MCMC_par)[k] )


	low_index  = which(DEN$x<QUANT[1])
	mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
	high_index = which(DEN$x>QUANT[3])

	polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
            	 col="pink")

	polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
	         col="grey")

	polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
		   col="pink")

	points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)
}





#############################################
#############################################
##          ##                             ##
##   ####   ##  ###### #####  ###  ######  ##
##  ##  ##  ##    ##   ##    ##      ##    ##
##     ##   ##    ##   ####   ###    ##    ##
##  ##  ##  ##    ##   ##       ##   ##    ##
##   ####   ##    ##   #####  ###    ##    ##
##          ##                             ##
#############################################
#############################################

#############################################
## 3.1 Extract posterior medians and 
##     calculate model prediction

par_median = apply(X=MCMC_burn[,1:5], MARGIN=2, FUN=median)


age_seq = seq(from=0, to=90, by=1)

M2_predict_cs = sapply(age_seq, model_M2, t_survey=0, par=par_median)


age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )

###############################################
## 0.2 Prepare data for plotting

#age_bins     <- seq(from=0, to=60, by=5)
#age_bins_mid <- seq(from=2.5, to=57.5, by=5) 
 
N_bins <- length(age_bins) - 1 


GMT_bins      <- rep(NA, N_bins)

AB_range_bins <- matrix(NA, nrow=N_bins, ncol=3)
colnames(AB_range_bins) <- c("med", "low", "high")

for(i in 1:N_bins)
{
	index <- intersect( which(AB_data[,1]>age_bins[i]), which(AB_data[,1]<=age_bins[i+1]) ) 
	temp  <- AB_data[index,3]

	GMT_bins[i] <- exp(mean(log(temp)))

	AB_range_bins[i,] <- quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


#############################################
## 3.2 Posterior prediction intervals 

N_sam = 500
sam_seq = round(seq(from=1, to=nrow(MCMC_burn), length=N_sam))




M2_sam_cs = matrix(NA, nrow=N_sam, ncol=length(age_seq))
for(k in 1:N_sam)
{
	M2_sam_cs[k,] = sapply(age_seq, model_M2, t_survey=0, par=MCMC_burn[sam_seq[k],1:5])
}

M2_quant_cs = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
	M2_quant_cs[,j] = quantile( M2_sam_cs[,j], prob=c(0.025, 0.5, 0.975) )
}



###############################################
## 3.1 Plot data and model prediction

par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,90), ylim=c(0.01,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Antibody acquisition Model 2 with STI fix RR - Rombo CT694"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins_cs[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins_cs[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}

points(x=age_seq, y=M2_predict_cs, 
type='l', lwd=3, col="green")



######################################
######################################
##          ##                      ##
##  ##      ##  ####   ####  ####   ##
##  ## ##   ##  ## ##   ##  ##  ##  ##
##  ######  ##  ##  ##  ##  ##      ##
##     ##   ##  ## ##   ##  ##  ##  ##
##     ##   ##  ####   ####  ####   ## 
##          ##                      ##
######################################
######################################
##
## Estimate the Deviance Information Criterion.
## Note that the deviance is -2*log(likelihood)

######################################
## Mean of the posterior

theta_bar = apply( X=MCMC_burn[,1:5], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2_total( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*mean( MCMC_burn[,6] - MCMC_burn[,7] )


######################################
## Effective number of model parameters
	
pD = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC = pD + D_bar

######################################
## Median of the posterior

theta_bar = apply( X=MCMC_burn[,1:5], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2_total( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*median( MCMC_burn[,6] - MCMC_burn[,7] )


######################################
## Effective number of model parameters
	
pD_median = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC_median = pD_median + D_bar



