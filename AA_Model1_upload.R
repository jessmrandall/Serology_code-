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

## call the libraries that are needed for analysis
pacman::p_load("here", "MASS", "compiler", 
               "binom", "coda", "readr", 
               "janitor", "tidyverse", "assertr")

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
 
###############################################
## 0.1 Read in data

# creates a list of all filepaths to data to be read into the script, being explicit about which file you're using with paths
# this script can be scaled to accomodate as many files as you have
files <- list(
  file1 = here("/model/input/file1.csv"),
  file2 = here("/model/input/file2.csv"),
  file3 = here("/model/input/file3.csv"),
  file4 = here("/model/input/file4.csv")
  )

# use a unit test to check your list has the number of files you're expecting to call and none are missing
stopifnot(length(files) == 4 & is.na(files) == FALSE)

# create a list of all files as connections
fileslist <- list(files$file1, files$file2, files$file3, files$file4)

# use a unit test to check your list has the number of files you're expecting to call and none are missing
stopifnot(length(fileslist) == 16 & is.na(fileslist) == FALSE)

# create a list called dfs, (lines 63-78)
# the list holds all 4 dataframes created from csvs and cleans up their column titles to be R friendly (lines 68-69)
# use unit tests to check the dataframes look as you expect them to look with
# the right variables and dimensions (lines 71-76)
dfs <- lapply(fileslist, function(x) {
  
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
  clean_names()
  
  x_df  %>%
    verify(ncol(x_df) == 3) %>%
    verify(is.na(x_df) == FALSE) %>%
    transmute(age = age, 
              titre = titre, 
              sero_pos = sero_pos)
})

# add names for each df in the list corresponding to appropriate names for each
# spreadheet, in this case the test and the country it was performed in
# use names without spaces as these will become the file names for your exported estimates

df_names <- c("test1_country1", "test2_country1", "test1_country2", "test2_country2")

names(dfs) <- df_names

## using dfs, we extract each df,run it through the model, 
## and save and export result to the plot task ##

### calculate observed age seroprevalence ###

# loop though 4 datasets to calculate observed age seroprevalence
# start i loop 
for (i in seq_along(dfs)){
	
  # set seed for reproducibility of results
  set.seed(22315)            
  seed = 22315
  
  #messages for the user to keep them aware of model progress
  start_time <- Sys.time()
  print(paste0("Age seroprevalence calculation for dataset ",names(dfs)[i]," has now begun..."))
  
  df <- as.data.frame(pluck(dfs, i))
  
  # age bin the data for plotting
  age_bins <- seq(from=0, to=9, by=1)
  age_bins_mid <- seq(from=0.5, to=8.5, by=1) 
  
  N_bins <- length(age_bins) - 1
      
  # initialize empty dataframe to fill with binomial confidence intervals
  # from observed data 
  sp_bins <- data.frame(med = numeric(0), 
                        low_95 = numeric(0), 
                        high_95 = numeric(0))
  
  # loop thorugh data to populate the sp_bins into a 9x3 matrix
  # containing observed age seroprevalence proportions for each age group and
  # 95% confidence intervals
  # start k loop
      for(k in 1:N_bins){
        
        index <- which(df[,1]> age_bins[k] & df[,1]<=age_bins[k+1])
          
        temp  <- df[index,3]
        
        sp_bins[k,] <- as.numeric(as.vector(binom.confint(sum(temp), 
                                                            length(temp),
                                                            method="wilson",
                                                            seed = seed)[1,4:6]))
        } # close k loop

  # create a new column called age from the row numbers to join with 
  # age bin data

  sp_bins <-as.data.frame(sp_bins) %>%
    mutate(age = as.numeric(row.names(sp_bins)))

  #check that no data are missing 
  stopifnot(not_na(sp_bins) == TRUE)
    
  # create country and test specific age bin data frame
  age_bins <- as.data.frame(age_bins_mid) %>%
    filter(age_bins_mid != 9.5) %>%
    mutate(age = as.numeric(sp_bins$age))
  
  #merge observed prevalence proprtions, confidence intervals, age bins and age 
  # can export this for plotting in tidyverse
  obs<- left_join(sp_bins, age_bins, by = "age") %>%
    mutate(age = as.numeric(sp_bins$age))
  
  #message to let the user know that each iteration has completed
  print(paste0("Age seroprevalence for dataset ",names(dfs)[i]," has completed successfully."))
  
} # close i loop

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
## 1.1 MODEL   specifiy the AA model that may have
## generated the data 

par_MC <- c(2, 0.1, 0.55)  ## (alpha, rr, sigma)
 
model_M1 <- function(a, par)
{
	alpha <- par[1]
     rr    <- par[2]
	sigma <- par[3]
 
	AB_titre <- ( alpha/rr )*( 1 - exp(-rr*a) )
 
	AB_titre
}


model_M1 <- cmpfun(model_M1, options=list(optimize=3)) 


###################################################
## 1.2 LIKELIHOOD evaluate the binomial likelihood
 
loglike_M1 <- function( par )
{
	alpha <- par[1]
    rr    <- par[2]
	sigma <- par[3]

	AB_model <- model_M1( AB_data[,1], par )

	mu <- log(AB_model) 

	loglike <- - log(AB_data[,2]) - log(2.506628*sigma) - 0.5*( (log(AB_data[,2])-mu)/sigma )^2

      sum( loglike )
}

loglike_M1 <- cmpfun(loglike_M1, options=list(optimize=3))


###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M1 <- function( par )
{
	alpha <- par[1]
      rr    <- par[2]
	sigma <- par[3]

	######################################
	## Uniform prior on alpha ~ U(0,100)

	if( alpha>0 && alpha<100 )
	{
		prior_alpha <- log(1/100)
	}else{
		prior_alpha <- -LARGE
	}

	######################################
	## Uniform prior on rr ~ U(0,10)

	if( rr>0 && rr<10 )
	{
		prior_rr <- log(1/10)
	}else{
		prior_rr <- -LARGE
	}

	######################################
	## Uniform prior on sigma ~ U(0,10)

	if( sigma>0 && sigma<10 )
	{
		prior_sigma <- log(1/10)
	}else{
		prior_sigma <- -LARGE
	}

	prior <- prior_alpha + prior_rr + prior_sigma

	prior
}

prior_M1 <- cmpfun(prior_M1, options=list(optimize=3))


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


N_mcmc       <- 10000      ## Number of MCMC iterations
N_tune_start <- 300        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   <- 3000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      <- 4000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  <- 1           ## Scaler for step size
MCMC_accept <- 0           ## Track the MCMC acceptance rate

max_corr    <- 0.75        ## Maximum degree of correlation

#################################################
## 2.1 Robbins-munro step scaler


rm_scale <- function(step_scale, mc, log_prob)
{
	dd <- exp(log_prob)
	if( dd < -30 ){ dd <- 0 }
	dd <- min( dd, 1 )

	rm_temp <- ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
		
	out <- step_scale*exp(rm_temp)
	
	out <- max( out, 0.05 )
	out <- min( out, 5)
	out
}


#################################################
## 2.2 Prepare object for MCMC fitting
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par) <- c("alpha", "rr", "sigma", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations

par_MC <- c(50, 0.1, 0.5)                 ## Initial guess: (alpha, rr, sigma)

Sigma_MC <- diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

loglike_MC <- loglike_M1( par_MC ) + prior_M1( par_MC )



for(mc in 1:N_mcmc)
{
	par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)



	if( prior_M1(par_MCp1) > -0.5*LARGE  ){
 
		loglike_MCp1 <- loglike_M1( par_MCp1 ) + prior_M1( par_MCp1 )


		log_prob <- min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC <- par_MCp1

			loglike_MC  <- loglike_MCp1
			MCMC_accept <- MCMC_accept + 1                       
		}

		#######################################
		## RM scaling of proposal step size

		if( mc < N_adapt ){
			step_scale <- rm_scale( step_scale, mc, log_prob)
		}


		#######################################
		## Adaptive tuning of covariance matrix

		if( (mc > N_tune_start) && (mc < N_tune_end) )
		{
			cov_MC <- cov( MCMC_par[1:(mc-1),1:3] )

			if( min(diag(cov_MC))>1e-6 )
			{ 
				###########################
				## Check for high degree of correlation

				sd_MC_inv <- 1/sqrt(diag(cov_MC))

				corr_MC <- t(t(cov_MC*sd_MC_inv)*sd_MC_inv)

				corr_MC[intersect( which( corr_MC > max_corr ), which(corr_MC<0.99999) )] <- max_corr
				corr_MC[  which( corr_MC < -max_corr )] <- -max_corr

				cov_MC <- t(t(corr_MC*(1/sd_MC_inv))*(1/sd_MC_inv))

				Sigma_MC <- cov_MC
			}
		}
	}

	MCMC_par[mc,1:3] <- par_MC
	MCMC_par[mc,4]   <- loglike_MC
	MCMC_par[mc,5]   <- prior_M1( par_MC )
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,2))



#####################################
## PANEL 1: alpha MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,1], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="alpha", 
main="alpha")




#####################################
## PANEL 2: rr MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,2], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="rr", 
main="rr")




#####################################
## PANEL 3: sigma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,3], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="sigma", 
main="sigma" )




#####################################
## PANEL 4: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,4], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
main="likelihood" )







#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(1,3))


#####################################
## PANEL 1: alpha MCMC posterior


DEN <- density( MCMC_burn[,1] )
	
QUANT <- quantile( MCMC_burn[,1], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="alpha", ylab="", 
main="posterior: alpha" )

	
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





#####################################
## PANEL 2: rr MCMC posterior


DEN <- density( MCMC_burn[,2] )
	
QUANT <- quantile( MCMC_burn[,2], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="rr", ylab="", 
main="posterior: rr" )

	
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




#####################################
## PANEL 3: sigma MCMC posterior


DEN <- density( MCMC_burn[,3] )
	
QUANT <- quantile( MCMC_burn[,3], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="sigma", ylab="", 
main="posterior: sigma" )

	
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


par_median <- apply(X=MCMC_burn[,1:3], MARGIN=2, FUN=median)



age_seq <- seq(from=0, to=90, by=1)

M1_predict <- model_M1(age_seq, par_median )





###############################################
## 3.1 Plot data and model prediction
 
par(mfrow=c(1,1))


plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,60), ylim=c(0.1,500),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Antibody acquisition Model 1 fit"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}


points(x=age_seq, y=M1_predict, 
type='l', lwd=3, col="blue")


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

theta_bar = apply( X=MCMC_burn[,1:3], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M1( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*mean( MCMC_burn[,4] - MCMC_burn[,5] )


######################################
## Effective number of model parameters
	
pD = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC = pD + D_bar


######################################
## Median of the posterior

theta_bar = apply( X=MCMC_burn[,1:3], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M1( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*median( MCMC_burn[,4] - MCMC_burn[,5] )


######################################
## Effective number of model parameters
	
pD_median = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC_median = pD_median + D_bar

