######
######
###### Models
######
######



# The steady-state model - supply t, a, Amax, and k; return A0
SSModel <- function(t, k, a, Amax) {
        y <- a + (Amax - a) * (1 - exp(-k*t))
        return(y)
}

# Solves dK/dA analytically
SSModel_dk <- function(t, k, a, Amax, dA) {
        k*exp(1)/(Amax-a)*dA
}

# The non-steady state equation - supply t, a, k, pss, kp, N; return A0.
NSModel <- function(t, k, N, kp, pss, a){
        z <- 0
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (k/(k-n*kp))*b
                y <- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))
                z <- y+z}
        return(z)
}

# This function calculates dk of fitting, based on the analytical solution of dA/dk
NSModel_dk <-function(x, k, N, kp, pss, a, dA){
        z <- 0
        for (n in 0:N) {
                b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
                bp <- (k/(k-n*kp))*b
                y <- a *((n*kp)/(k*(k-n*kp))*bp*(exp(-k*x)-exp(-n*kp*x))-x*(1/(N+1)-bp)*exp(-k*x))                                                              
                z <- y + z}
        return(dA/z)
}


# This is an adaptaion of the Guan et al. two-compartment model - supply t, a, Amax, ksyn, and kdeg, return A0.
CCModel <- function(t, kdeg, ksyn, a, Amax) {
        
        a + (Amax-a) * (1.0 - (exp(-t*kdeg)*ksyn - exp(-t*ksyn)*kdeg)/
                                (ksyn - kdeg))
        
}


# This function calculates dk for fitting, based on the analytical solution of dA/dkdeg
CCModel_dk <- function(t, kdeg, ksyn, a, Amax, dA) {
        
        z <- (Amax-a) * ((t*(ksyn-kdeg) - 1.0) * exp(-t*kdeg) + exp(-t*ksyn))*ksyn / (ksyn-kdeg)^2
        dkdeg <- dA/z
        
        y <- (Amax-a) * (exp(-t*kdeg) - (1.0 + (ksyn - kdeg)*t) * exp(-ksyn*t))*kdeg/ (ksyn-kdeg)^2
        dksyn <- dA/y
        
        return(dkdeg)
}




######
######
###### Helpers
######
######


# Calculate N
calcLabelingSite <- function(string){
        
        N <- 0.0
        carb <- 0
        nitr <- 0
        oxyg <- 0 
        sulf <- 0
        
        # Deal with modifications first
        if(grepl("\\(15\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr + 0; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(42\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 2 ; nitr = nitr + 0; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(0\\.984\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr - 1; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(79\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr + 0; oxyg = oxyg + 4; sulf = sulf + 0}
        if(grepl("\\(114\\..+?\\)",string) == T) {N = N + 4.12; carb = carb + 4 ; nitr = nitr + 3; oxyg = oxyg + 1; sulf = sulf + 0}
        
        string <- gsub( " *\\(.*?\\) *", "", string)
        
        # Remove one oxygen per peptide bond
        oxyg <- oxyg - nchar(string) +1
        
        # Each character - add labeling sites from Commerford et al., add atoms.
        for (i in 1:nchar(string)){
                char <- substring(string,i,i)
                if(char == "A"){N = N + 4.00; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "C"){N = N + 1.62; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 1}
                if(char == "D"){N = N + 1.89; carb = carb + 4 ; nitr = nitr + 1; oxyg = oxyg + 4; sulf = sulf + 0}
                if(char == "E"){N = N + 3.95; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 4; sulf = sulf + 0}
                if(char == "F"){N = N + 0.32; carb = carb + 9 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "G"){N = N + 2.06; carb = carb + 2 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "H"){N = N + 2.88; carb = carb + 6 ; nitr = nitr + 3; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "I"){N = N + 1.00; carb = carb + 6 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "K"){N = N + 0.54; carb = carb + 6 ; nitr = nitr + 2; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "L"){N = N + 0.69; carb = carb + 6 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "M"){N = N + 1.12; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 1}
                if(char == "N"){N = N + 1.89; carb = carb + 4 ; nitr = nitr + 2; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "P"){N = N + 2.59; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "Q"){N = N + 3.95; carb = carb + 5 ; nitr = nitr + 2; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "R"){N = N + 3.34; carb = carb + 6 ; nitr = nitr + 4; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "S"){N = N + 2.61; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "T"){N = N + 0.20; carb = carb + 4 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "V"){N = N + 0.56; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "W"){N = N + 0.08; carb = carb + 11; nitr = nitr + 2; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "Y"){N = N + 0.42; carb = carb + 9 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
        }
        
        N <- round(N)
        return(c(N,carb,nitr,oxyg,sulf))
}


# Calculate natural abundance of isotopes from NIST and http://www.iupac.org/publications/pac/83/2/0397/
calc_thr_a <- function(string){
        atoms <- calcLabelingSite(string)
        c <- atoms[2]
        n <- atoms[3]
        o <- atoms[4]
        s <- atoms[5]
        a <- 0.9893^c * 0.99636^n * 0.99757^n * 0.9499^s
        
        return(a)
}



#
# Calculate the essentiality of a peptide sequence based on mammalian requirements of essential vs. non-essential amino acid
# At present it should scale from 0 to 1, although we may want to take into account certain log length factors later on,
# or the quantify the individual amino acid's essentiality by quantifying how much nutritional requirement is needed for a particular species
# Based on Dispensable and indispensable amino acids for humans Reeds J Nutr 2000 PMID 10867060
#

calc_essentiality <- function(string) {
        
        # essentiallity counter
        ess <- 0
        
        len <- nchar(string)
        
        # Each character - add labeling sites from Commerford et al., add atoms.
        for (i in 1:nchar(string)) {
                char <- substring(string,i,i)
                
                # Essential
                if (char %in% c('L','I','V','K','T','M','W','F','H','P')){ess = ess + 1}
                
                # Non-essential with no other information
                if (char %in% c('Y', 'N', 'Q', 'C')){ess = ess + 0}
                
                # Conditionally essential, or non-essential but with synthesis only accounting partially for flux
                if (char == 'E'){ess = ess + 0.02}
                if (char == 'S'){ess = ess + 0.05}
                if (char == 'D'){ess = ess + 0.22}
                if (char == 'A'){ess = ess + 0.54}
                if (char == 'G'){ess = ess + 0.65}
                if (char == 'R'){ess = ess + 0.86}
                
        }
        return(ess/len)
        
}