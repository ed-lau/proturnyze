#OU_Protein functions

#--------------dexp_log-----------------------------------------
log_dexp <- function(x, rate)   # rate is decay rate, remember that you have used
{                               # the R function to for prior distribution,
   if(x <= 0)		        # and R uses -x/rate for exponential distributions
   {
     log_p <- -10^8;
   }
   else
   {
      log_p <- -x/rate - log(rate);   
   }

   log_p
}

#--------------prime_log_dexp----------------------------------
#
# derivative of the log of the exponential distribution
# function wrt to its variable
#
prime_log_dexp <- function(x, rate){ # rate is decay rate

   if(x <= 0)
   {
     prime_log <- -10^8;
   }
   else
   {
      prime_log <- -1./rate;   
   }

   log_p
}

#-------------- OU---------------------------------------
# The elements of the Covariance Matrix of the Ornstein-Uhlenbeck process
# For different stochastic processes it could vary
# for two-compartment model of protein turnover the
# formula is OH(t,s) = exp(-abs(s-t)kdeg)*sigma_g/(2*kdeg)
# sigma_g is the model fluctuation (volatility) variance
# for meaning of the kdeg see the two-compartment model
# in the formulae below d = t - s; sigma_g = par[1]; - remember
# that sigma_g is actualy the square of the sigma_g, so you do not
# need to raise it again to the power of 2.
# kdeg = par[2]
OU <- function(par, d){

    sigma_gamma <- par[1];

    kdeg <- par[2];

    OU_element  <- exp(-kdeg*abs(d))*sigma_gamma^2/(2.*kdeg);

    OU_element
}


#---Covariance Matrix for Ornstein Uhlenbek Process + White Noise----
# computes covariance function on a single grid
# par[1] is the sigma_g - variance of the model fluctuations (volatility)
# par[2] is the kdeg in the two-compartment model of protein turnover
# par[3] is the sigma_e - experimental white noise
# N_time is the number of time points, t is the vector of the time points
# where the measurements were taken
# pr_par <- c(rate_kdeg, rate_ksyn, rate_sigma_g, rate_sigma_e);
#
Cov_OU <- function(par, t)
{ 
    par1 <- par[1:2]; 

    N_time <- length(t);

    sigma_epsilon = par[3];

    epsilon = sigma_epsilon; 

    GRAM <- mat.or.vec(N_time, N_time);

    for (i in 1:N_time)
    {
        for (j in i:N_time)
        {
            GRAM[i,j] <- OU(par1, t[i]-t[j] )
            
            if (i == j)
            {
                GRAM[i, j] <- GRAM[i, j] + epsilon^2;
            }
            
            GRAM[j, i] = GRAM[i, j];
        }
    }
  
    GRAM
}


#-------------Derivative_Of_Covariance/GRAM_Function-----------------------------
# computes derivative of covariance function on a single grid
# with respect to par[3] determines the sigma - experimental white noise
# kdeg = par[1]; ksyn = par[2]; 
# sigma_g = par[3]; sigma_e = par[4]
#
Prime_GRAM_Sigma_GAMMA <- function(par, t){ 

    kdeg <- par[1];

    sigma_g <- par[3];

    sigma_e <- par[4];

    N <- length(t)

    PRIME_GRAM <- mat.or.vec(N, N) ;

    for (i in 1:N)
    {
        for (j in 1:N)
        {
            PRIME_GRAM[i,j] <- exp(-abs(t[i] - t[j])*kdeg)/(2.*kdeg);
        }
    }
    
    PRIME_GRAM;
}

#-------------Derivative_Of_Covariance/GRAM_Function-w.r.t.-kdeg----------
# computes derivative of covariance function on a single grid
# with respect to par[3] determines the sigma - experimental white noise
# kdeg = par[1]; ksyn = par[2]; 
# sigma_g = par[3]; sigma_e = par[4]
#
Prime_GRAM_KDEG <- function(par, t){ 

    kdeg <- par[1];

    sigma_g <- par[3];

    sigma_e <- par[4];

    N <- length(t)

    PRIME_GRAM_KDEG <- mat.or.vec(N, N) ;

    for (i in 1:N)
    {
        for (j in 1:N)
        {
            PRIME_GRAM_KDEG[i,j] <- -exp(-abs(t[i] - t[j])*kdeg)*sigma_g/(2.*kdeg) *
  				    (1./kdeg + abs(t[i] - t[j]) );
        }
    }
    
    PRIME_GRAM_KDEG;
}


#-----log_likelihood(data probability)_two_compartents---------------------------------
# This function is adapted to the OU covariance matrix
# the log-likelihood function uses two-compartment model 
# plus stochasticity
# to determinine the functional fit.
# par contains 4 parameters.
# kdeg = par[1]; ksyn = par[2]; 
# sigma_g = par[3]; sigma_e = par[4]
#
# pr_par <- c(rate_kdeg, rate_ksyn, rate_sigma_g, rate_sigma_e);
#
log_like_two_compartments <- function(par, Y, t1, prior){

    if(par[1] < 0 || par[2] < 0 || par[3] < 0 || par[4] < 0)
    {
       return (10^8);
    }

    N <- length(Y)

    kdeg = par[1]; ksyn = par[2]; 
    
    sigma_g = par[3]; sigma_e = par[4];

    par2 <- c(0.,0.,0.);

    par2[1] = sigma_g; par2[2] = kdeg; par2[3] = sigma_e;

    #print("PARS: ");
    print(c("Params", signif(par, digits=3) ));
   
    GRAM <- Cov_OU(par2, t1);

    GRAM_Inv <- solve(GRAM + diag(N)*0.01001, tol=1e-100)  

    Y_temp = Y - par[5]* (1.0 - (exp(-t1*kdeg)*ksyn - exp(-t1*ksyn)*kdeg)/
		                  (ksyn - kdeg) );

    # resudial sum of squares
    RSS <- t(Y_temp) %*% GRAM_Inv %*% Y_temp      

    prior_kdeg <- log_dexp(kdeg, prior[1]);

    prior_ksyn <- log_dexp(ksyn, prior[2]);

    prior_sigma_g <- log_dexp(sigma_g^2, prior[3]);   # the prior is defined for sigma_g

    prior_sigma_e <- log_dexp(sigma_e^2, prior[4]); # the prior is defined for sigma_e

    GRAM_SVD <- svd(GRAM) 

    A <- sum(log(GRAM_SVD$d)) # instead of log(det(GRAM_SVD))

    log_lkh <- 0.5*A + 0.5*RSS + 0.5*N*log(2*pi);

    log_lkh <- log_lkh - prior_kdeg - prior_ksyn - prior_sigma_g - prior_sigma_e;

    log_lkh
}

#--------------Gradient_for_Stochastic_two_compartents-------------------------------------
# a stochastic function fit  two-compartmental model (Guan approach)
# to determinine the functional fit.
# par contains 4 parameters.
# kdeg = par[1]; ksyn = par[2]; 
# sigma_g = par[3]; sigma_e = par[4]
# In driving the gradients, remember that you are minimizing -log_liklh not
# the log_liklh itself. Therefore, the gradient signs will seem to be minus
# compared to those of the log_liklh
#
Stochastic_grad_two_comparts <- function(par, Y, t1, prior){

    N <- length(Y)

    N <- length(Y)
   
    GRAM <- Cov_OU(par, t1);

    GRAM_Inv <- solve(GRAM + diag(N)*0.0001, tol=1e-100);

    kdeg = par[1]; ksyn = par[2]; 
    
    sigma_g = par[3]; sigma_e = par[4]

    Target = Y - par[5] * (1.0 - (exp(-t1*kdeg)*ksyn - exp(-t1*ksyn)*kdeg)/
		                  (ksyn - kdeg) );


    #derivative w.r.t. k_deg = par[1]  - only part of it, will need to add the rest!!!!!
    #
    grad1 = (  (t1*(par[2] - par[1]) - 1)*exp(-t1*par[1]) + exp(-t1*par[2]) )*par[2]/
		(par[2] - par[1])^2;

    grad1 <- -t(grad1) %*% GRAM_Inv %*% Target;

    # add the derivative of the prior for kdeg, since it is derivative
    # of the log of the exponential function, it is simply the value, of the rate
    # in this case the prior
    grad1 = grad1 + prior[1];
     
    temp_grad1 = Prime_GRAM_KDEG(par, t1);

    #these are components of the kdeg derivative wrt Gram/Cov matrix
    grad1 = grad1 + 0.5* sum(diag(GRAM_Inv %*% temp_grad1)) - 
		0.5 * t(Target) %*% GRAM_Inv %*% temp_grad1 %*% GRAM_Inv %*% Target;


    #derivative w.r.t. ksyn = par[2]
    #
    grad2 <- ( exp(-t1*par[1])- (1 + (par[2] - par[1])*t1)*exp(-par[2]*t1) )*par[1]/
		(par[2] - par[1])^2;

    grad2 <- -t(grad2) %*% GRAM_Inv %*% Target;

    # add the derivative of the prior for ksyn, since it is derivative
    # of the log of the exponential function, it is simply the value, prior[2]
    grad2 = grad2 + prior[2];

    #derivative w.r.t. sigma_g = par[3]
    #
    grad3 = Prime_GRAM_Sigma_GAMMA(par, t1);  # Prime_GRAM_Sigma_GAMMA

    grad3 = 0.5* sum(diag(GRAM_Inv %*% grad3)) - 
		0.5 * t(Target) %*% GRAM_Inv %*% grad3 %*% GRAM_Inv %*% Target;

    # add the derivative of the prior for sigma_g, since it is derivative
    # of the log of the exponential function, it is simply the value, prior[3]
    grad3 = grad3 + prior[3];

    #derivative w.r.t. sigma_epsilon (white noise)= par[4]
    # derivative of the GRAM matrix w.r.t. sigma_e is an Identity matrix
    #
    grad4 =  0.5* sum(diag(GRAM_Inv)) - 
		   0.5 * t(Target) %*% GRAM_Inv %*% GRAM_Inv %*% Target;

    # add the derivative of the prior for sigma_g, since it is derivative
    # of the log of the exponential function, it is simply the value, prior[4]
    #
    grad4 = grad4 + prior[4];

    
    #derivative w.r.t. asymptote parameter = par[5]
    #

    grad5 = (1.0 - (exp(-t1*kdeg)*ksyn - exp(-t1*ksyn)*kdeg)/
		                  (ksyn - kdeg) );

    grad5 =  -t(grad5) %*% GRAM_Inv %*% Target;
   
    c(grad1, grad2, grad3, grad4, grad5);
}

#--------------nonStochastic_two_compartents-------------------------------------------------
# a non-stochastic function fit  two-compartmental model (Guan approach)
# to determinine the functional fit.
# par[1] and par[2] are the k_deg and k_syn, respectively, in the 
# two-compartmental model.
# par[3] is asymptote level in the model. if not specified
# it is set equal to 1.
nonStochastic_two_compartments <- function(par, Y, t1){

    N <- length(Y)
   
    Y_temp = Y - par[3]*(1.0 - (exp(-t1*par[1])*par[2] - exp(-t1*par[2])*par[1])/
		                  (par[2] - par[1]) );
    
    M <- t(Y_temp) %*% Y_temp;      

    log_lkh = M;             #this is an important moment to correctly
  			     #define the function that you want to minimize!

    c(log_lkh)
}

#--------------Gradient_for_nonStochastic_two_compartents-------------------------------------
# a non-stochastic function fit  two-compartmental model (Guan approach)
# to determinine the functional fit.
# par[1] and par[2] are the k_deg and k_syn, respectively, in the 
# two-compartmental model.
# par[3] - the asymptotic value
grad_two_compartments <- function(par, Y, t1){

    N <- length(Y)
   
    Y_temp = Y - par[3]*(1.0 - (exp(-t1*par[1])*par[2] - exp(-t1*par[2])*par[1])/
		                  (par[2] - par[1]) );
    
    #derivative w.r.t. par[1], k_deg
    x1 = par[3]* (  (t1*(par[2] - par[1]) - 1)*exp(-t1*par[1]) + exp(-t1*par[2]) )*par[2]/
		(par[2] - par[1])^2;

    #derivative w.r.t. par[2], k_syn
    x2 = par[3]* ( exp(-t1*par[1])- (1 + (par[2] - par[1])*t1)*exp(-par[2]*t1) )*par[1]/
		(par[2] - par[1])^2;

    #derivative of the sum of squares.

    #derivative w.r.t. par[3] - the asymptotic value

    x3 = (1.0 - (exp(-t1*par[1])*par[2] - exp(-t1*par[2])*par[1])/
		                  (par[2] - par[1]) );

    c(-2.*t(x1) %*% Y_temp , -2.*t(x2) %*% Y_temp, -2*t(x3) %*% Y_temp);
    
    }