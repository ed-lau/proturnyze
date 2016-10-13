source("Rahman_et_al_Code/OU_Protein.R")

# ========Generate_Data================
brain <- read.table("data/OU-test_brain.txt", header=F, sep="\t", fill=T, as.is=T, quote="", skip = 1)
colnames(brain) <- c("Acc", "nProt", "nProtUsed", "d0", "d0.38", "d1", "d2", "d4", "d8", "d16", "d24", "d32" )


data = brain


Nrows = nrow(data);  Ncols = ncol(data);

Proteins <- data$Acc;

R2 <- array(Nrows, dim=1);  Asympt <- array(Nrows, dim=1);  Kdeg <- array(Nrows, dim=1);

Kdeg_nonStochastic <- array(Nrows, dim=1); Ksyn_nonStochastic <- array(Nrows, dim=1);

Kdeg_Stochastic <- array(Nrows, dim=1); Ksyn_Stochastic <- array(Nrows, dim=1);

Sigma_Epsilon_nonStochastic <- array(Nrows, dim=1);

Sigma_Gamma_Stochastic  <- array(Nrows, dim=1);

Sigma_Epsilon_Stochastic  <- array(Nrows, dim=1);

LogLike_Stochastic <- array(Nrows, dim=1); R2_nonStochastic <- array(Nrows, dim=1);

t <-  c(0, 0.38, 1, 2, 4, 8, 16, 24, 32);    # time course points

t1 <- t;

N1 <- length(t1);

# this module computes the decay rates using a non-stochastic curve fitting
for(j in 1:Nrows) 
  {
	Isotopes = as.numeric(data[j, 4:12]);   


	Y <- Isotopes;

        #non-stochastic, derivative driven estimation using BFGS
	      #nonStochastic_two_compartments(c(0.1,  0.2, 1),Y,t1);
	      #grad_two_compartments(c(0.1,  0.2, 1),Y,t1);
	      
        par_optim <- optim(c(0.1,  0.2, 1), fn = nonStochastic_two_compartments,
                  gr = grad_two_compartments, method = "BFGS", Y = Y, t1 = t);

        print(c("BFGS", par_optim));

        temp = par_optim$par[1];

        if(temp >  par_optim$par[2])
	{
            par_optim$par[1] = par_optim$par[2];

	    par_optim$par[2] = temp;
	}

	Kdeg_nonStochastic[j] = par_optim$par[1];

	Ksyn_nonStochastic[j] = par_optim$par[2];

        Asympt[j]              = par_optim$par[3];

        Y_fit <- 1 - (exp(-t1*par_optim$par[1])*par_optim$par[2] - 
		exp(-t1*par_optim$par[2])*par_optim$par[1])/(par_optim$par[2] - par_optim$par[1]);

        Y_fit = Asympt[j] * Y_fit;
	
	R2_nonStochastic[j] = cor(Y, Y_fit)#; 1 - sum((Y-Y_fit)^2)/sum((Y-mean(Y))^2);

        # the par_optim$value is always positive, because the likelihood
        # is defined as a residual sum of squares
        Sigma_Epsilon_nonStochastic[j] <- par_optim$value/(N1 - 1);

  }  #for(j in 1:Nrows) finishes non-stochastic estimate

#stopifnot(-1);

i_premat = 0;   i_conv = 0;

for(k in 1:Nrows)    # for every protein 
{

    Isotopes = as.numeric(data[k, 4:12]);     #protein incorporation measure   

    Y <- Isotopes;



    #--------Hyper-parameters------------------
    # set the initial values of the hyperparameters
    # equal to the the results from non_stochastic fitting
    #
    # kdeg - par[1]; ksyn = par[2]; 
    # sigma_g = par[3]; sigma_e = par[4]
    rate_kdeg <- 1/Kdeg_nonStochastic[k];

    rate_ksyn <- 1/Ksyn_nonStochastic[k];

    rate_sigma_e <- 1./Sigma_Epsilon_nonStochastic[k];  # for sigma_epsilon square determine from RSS
 
    rate_sigma_g <- rate_sigma_e;  # for sigma_gamma assume 

    asymptote = 1.; # the fifth parameter and the asymptotic value

    # prior parameters; changes for each protein, as it should
    # to reflect the different rate_kdeg, and rate_ksyn 
    pr_par <- c(rate_kdeg, rate_ksyn, rate_sigma_g, rate_sigma_e, asymptote);

    #-----------------Optimisation--------------------------------
    
    kdeg_init <- rexp(1, rate_kdeg);  #initiate the parameters

    ksyn_init <- rexp(1, rate_ksyn);  #initiate the parameters

    sigma_g_init <- rexp(1, rate_sigma_g);  #initiate the parameters

    sigma_e_init <- rexp(1, rate_sigma_e);  #initiate the parameters    

    par <- c(kdeg_init, ksyn_init, sigma_g_init, sigma_e_init, asymptote); 
  
    par_init <- par;

        par_init[1] = Kdeg_nonStochastic[k];

        par_init[2] = Ksyn_nonStochastic[k];

        par_init[3] = Sigma_Epsilon_nonStochastic[k]; # this is sigma_g
       
	par_init[4] = Sigma_Epsilon_nonStochastic[k]/2; # this is sigma_e;

	par_init[5] = Asympt[k];

 
        par_optim <- optim(par_init, fn = log_like_two_compartments,
        gr = Stochastic_grad_two_comparts, method = "BFGS", Y = Y, 
		prior = pr_par, t1 = t );


        LogLike_Stochastic[k] = par_optim$value;

        Sigma_Gamma_Stochastic[k] = par_optim$par[3];

	Sigma_Epsilon_Stochastic[k] = par_optim$par[4];

        if(par_optim$par[1] > par_optim$par[2])
        {
 	    Kdeg_Stochastic[k] = par_optim$par[2];
            Ksyn_Stochastic[k] = par_optim$par[1];
        }
        else
        {
 	    Kdeg_Stochastic[k] = par_optim$par[1];
            Ksyn_Stochastic[k] = par_optim$par[2];
        }

     par_opt <- par_optim$par;

     print("Final Optimized par_optim");


     Y_fit <- 1 - (exp(-t1*par_optim$par[1])*par_optim$par[2] - 
	    exp(-t1*par_optim$par[2])*par_optim$par[1])/(par_optim$par[2] - par_optim$par[1]);

     Y_fit <- Y_fit * par_optim$par[5];
     

     R2[k] = cor(Y, 1 - (par_opt[2]*exp(-par_opt[1]*t1) - par_opt[1]*exp(-par_opt[2]*t1))/
                           (par_opt[2] - par_opt[1]) );

     #R2[k] = 1 - sum((Y-Y_fit)^2)/sum((Y-mean(Y))^2);

     R2[k] = sqrt(t(Y - Y_fit) %*% (Y - Y_fit)/(length(t) - 1));

     R2[k] = t(Y - Y_fit) %*% (Y - Y_fit);

     R2[k] = sqrt(sum((Y-Y_fit)^2));

     R2[k] = cor(Y, Y_fit);
     
     
     if(par_optim$convergence != 0)
     {
          i_premat = i_premat + 1;
     }
     else if(par_optim$convergence == 0)
     {
	 i_conv = i_conv + 1;
     }

     if(par_optim$par[1] < 0 || par_optim$par[2] < 0)
     {
	 print("Some problem here");
        l <-  readline();
     }
     
     
 }  # finishes all proteins here! (iteration over k)

print(c("Premat number ", i_premat, "Converge number ", i_conv));

