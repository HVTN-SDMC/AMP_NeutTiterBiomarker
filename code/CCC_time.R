# ************************************************************************************ #
# Description:
#   Calculate concordance correlation coefficient across two timepoints
#   
# Code file:  CCC_time.R
# Directory:  /trials/vaccine/Lab_Program/CoV_2_nAb_concordance_202006/manuscript/macro
#
# Author:  Jin Kee, 2021/04/21
#
# Dependencies:
# -------------
#   Inputs Files:
#
# Output: 
#
# Code History
# ------------------------------------------------------------------------------
# 21Apr2021  J Kee        adapted from CCC.R in /trials/vaccine/Lab_Program/CoV_2_nAb_concordance_202006/macro/CCC.R
## -----------------------------------------------------------------------------

##############################################################################################################
# Input: 
#      x: target values (random or fixed)
#      y: measures of the same length of x
#      time: a vector of timepoints containing 1s and 2s (1=first timepoint, 2=secondpoint)
# Output: 
#       CCC: concordance correlation coefficient
#       CCC.l: asymptotic lower bound of CCC at alpha level
#       CCC.u: asymptotic upper bound of CCC at alpha level 
###############################################################################################################


CCC_time <- function(x,y,time,alpha=0.05){
    if(length(x)!=length(y)) stop("The length of x and y must be equal")
    if(all(is.na(x))| all(is.na(y))) {ccc=NA}
    else{
        
    # function to calculate various components of CCC: Z, Z_l, Z_u
    ccc_comp <- function(x,y) {
        ok <- complete.cases(x,y)
        x <- x[ok]
        y <- y[ok]
        n <- length(x)
        sdx <- sqrt((var(x)+0.0001)*(n-1)/n)
        sdy <- sqrt((var(y)+0.0001)*(n-1)/n)
        rho <- cor(x,y,use="complete.obs",method="pearson") # pearson correlation coefficient -- precision
        v2 <- (mean(x)-mean(y))^2/(sdx*sdy)
        wbar <- sdy/sdx
        Chi_a <- 2/(wbar+1/wbar+v2)# accuracy
        rho_c <- Chi_a * rho # concordance correlation coefficient
        
        # inverse hyperbolic tangent transformation
        Z <- 0.5*log((1+rho_c)/(1-rho_c))
        
        # Asymptotic inference for CCC
        temp1 <- (1-rho^2)*(rho_c)^2/(1-(rho_c)^2)/rho^2
        temp2 <- 2*v2*(1-rho_c)*(rho_c)^3/(1-(rho_c)^2)^2/rho
        temp3 <- v2^2*(rho_c)^4/2/(1-(rho_c)^2)^2/rho^2
        sdZ <- sqrt((temp1+temp2-temp3)/(n-2))
        Z_l <- Z - qnorm(1-alpha/2)*sdZ
        Z_u <- Z + qnorm(1-alpha/2)*sdZ
        
        return(list(Z=Z, Z_l=Z_l, Z_u=Z_u))
    }
    
    # get inverse hyperbolic tangent transformation of ccc for timepoint 1 and timepoint 2
    ccc_comp1 <- ccc_comp(x=x[time==1], y=y[time==1])
    ccc_comp2 <- ccc_comp(x=x[time==2], y=y[time==2])
    z1 <- ccc_comp1$Z
    z2 <- ccc_comp2$Z
    # calculate timepoint averaged z
    meanZ <- mean(c(z1,z2))
    # calculate timepoint integrated ccc
    CCC <- (exp(2*meanZ)-1)/(exp(2*meanZ)+1)
    
    # calculate timepoint integrated lower bound of ccc
    z1_l <- ccc_comp1$Z_l
    z2_l <- ccc_comp2$Z_l
    meanZ_l <- mean(c(z1_l,z2_l))
    CCC.l <- (exp(2*meanZ_l)-1)/(exp(2*meanZ_l)+1)
    
    # calculate timepoint integrated upper bound of ccc
    z1_u <- ccc_comp1$Z_u
    z2_u <- ccc_comp2$Z_u
    meanZ_u <- mean(c(z1_u,z2_u))
    CCC.u <- (exp(2*meanZ_u)-1)/(exp(2*meanZ_u)+1)
    
    return(list(CCC=CCC, CCC.l=CCC.l, CCC.u=CCC.u))
    }
}
