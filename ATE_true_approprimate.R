##---The Monte Carlo approximation of ATE true value and related true parameters of model---
##--OLS IE functions estimate---
IE_estfun<-function(beta_hat,#Local variables m*p
                    phi_hat,#m*p*p
                    alpha_hat#m*1
){
  
  #beta_hat=beta1; phi_hat=phi1; alpha_hat=alphas
  
  IE_iter2=0
  
  for(tau in 2:nrow(beta_hat)){
    
    IE_iter1=0
    for(k in 1:(tau-1)){
      
      library(tidyverse)
      
      if((k+1)>(tau-1)){
        
        IE_iter1=alpha_hat[k,]+IE_iter1
        
      }else{
        
        IE_iter1=Reduce("%*%", phi_hat[(k+1):(tau-1)])%*%alpha_hat[k,]+IE_iter1
        
      }
      
    }
    
    IE_iter2=as.vector(t(beta_hat[tau,])%*%IE_iter1)+IE_iter2
    
  }
  
  return(IE_iter2)
  
}



library(KernSmooth)
library(ks)
library(splines)
library(mvtnorm)
##--data generatinf process---
p=3 #the dimension of state variables 
#set the number of days
Ndays=seq(16,52,4)


#set the time span intervals
TI=c(1,3,6,12,24,48)

##--set the true coef---
set.seed(2023)## random seed it may be used to verify the consistent numerical results！

taus=c(48)  #taus=m in paper 1hour--24  30min--48

beta0=rep(0,taus)

for(i in 1:taus){
  
  beta0[i]=ifelse(rbinom(1,size = 1,prob = 0.5)>0.5,round(runif(1,min = -1, max=-0.5),digits = 3),round(runif(1,min = 0.5,max= 1) ,digits = 3))
  
}


beta1=matrix(nrow = taus,ncol = p)

for(i in 1:taus){
  for(j in 1:p){
    beta1[i,j]<-ifelse(runif(1)>0.5,round(runif(1,min = -0.3, max=-0.1),digits = 3),round(runif(1,min = 0.1,max= 0.3) ,digits = 3))
  }
}



gammas=runif(taus,min = 0.5,max=0.8)

#abs(rnorm(taus,mean=0.8, sd=0.15)) #rep(0.8,taus)# 


DE0<-sum(gammas)
DE0<-DE0/taus

phi0=matrix(nrow = taus,ncol = p)


for(j in 1:p){
  for(i in 1:taus){
    
    phi0[i,j]=ifelse(rbinom(1,size = 1,prob = 0.5)>0.5,round(runif(1,min = -1, max=-0.5),digits = 3),round(runif(1,min = 0.5,max= 1) ,digits = 3))   
    
  }
}



phi1=vector("list",length = taus)

is_or_eigen=numeric(length = taus)
for(i in 1:taus){
  
  #phi1[[i]]<-diag(x=runif(3,min = -0.5,max=0.5),nrow = p)
  
  phi1.mat=matrix(nrow = p,ncol = p)
  
  for(j in 1:p){
    for(k_l in 1:p){
      
      ## for the linear DGP
      #phi1.mat[j, k_l]=runif(1,min = -0.3, max = 0.3)
      
      ## for the non-linear DGP
      phi1.mat[j, k_l]=runif(1,min = -0.6, max = 0.6)
      
    }
  }
  
  
  phi1[[i]]<-phi1.mat
  
  phi_eigenmode=eigen(phi1[[i]])
  
  #phi_eigenmode$values
  #Mod(phi_eigenmode$values)
  
  is_or_eigen[i]=all(Mod(phi_eigenmode$values)<1)
  
  #phi1[[i]]= 0.4*phi1.mat
}

is_or_eigen # check the stability of state transition


alphas=matrix(nrow = taus,ncol = p)
for(j in 1:p){
  
  ##--different strength of carryover effects--
  alphas[,j]=rnorm(taus,mean=0, sd=0.3)#0.3
  
  # alphas[,j]=rnorm(taus,mean=0.8, sd=0.3)
  
  #alphas[,j]=rnorm(taus,mean=1.6, sd=0.3)
  
  # alphas[,j]=rnorm(taus,mean=2.4, sd=0.3)
  
}  



IE0<- IE_estfun(beta_hat = beta1, phi_hat = phi1, alpha_hat = alphas)
IE0<-IE0/taus

ATE0<-DE0+IE0 #OLS parameteric ATE true value



#--the different correlations of random effects
rho0=c(0.3,0.5,0.7,0.9)


##--outcome model ----
##--generating outcome conditional expectation
emission_distribution<-function(beta0_scal,#local variables
                                beta1_vec,states_vec,acts_scal,gamma_scal){
  # #--linear model--
  
  #emission_val=beta0_scal+t(beta1_vec)%*%states_vec+acts_scal*gamma_scal
  
  
  #--nonlinear model --
  
  emission_val=beta0_scal+t(beta1_vec)%*%(2*(sin(states_vec*acts_scal)+cos(states_vec))^2+3*states_vec*gamma_scal*acts_scal)+(acts_scal*gamma_scal+cos(acts_scal*gamma_scal))^2
  
  return(emission_val)
  
  
}



trans_distribution<-function(phi0_vec,phi1_mat,state_vec,act_scal,alpha_vec){
  
  states_next<-phi0_vec+phi1_mat%*%state_vec+act_scal*alpha_vec
  
  return(states_next)
  
}

eps_cors=diag(x=rep(1.5,p),nrow = p)


Iter_n=c(500,1000)
bias_difference=rep(0,length(Iter_n))
Policy1_true=rep(0,length(Iter_n))
Policy0_true=rep(0,length(Iter_n))

for(k in 1:length(Iter_n)){
  # k=1
  #taus,p is global 
  Y.mat1<-matrix(nrow = taus,ncol = Iter_n[k])
  
  Y.mat0<-matrix(nrow = taus,ncol = Iter_n[k])
  
  states.mat1=array(dim = c(taus+1, p, Iter_n[k]))
  
  #matrix(nrow = taus+1,ncol = Iter_n[k])
  states.mat0=array(dim = c(taus+1, p, Iter_n[k]))
  
  for(n in 1:Iter_n[k]){
    #n=1
    #--initial state variable--
    S.ini0<-as.vector(rmvnorm(1,mean = rep(0,p),sigma = diag(p)) )
    
    states.mat1[1, ,n]<-S.ini0
    states.mat0[1, ,n]<-S.ini0
    
    
    for(i in 1:taus){
      # i=1
      eps.1<-rnorm(1)
      eps.2<-rmvnorm(1, mean = rep(0,p),sigma = eps_cors)
      eps.2<-as.vector(eps.2)
      
      Y.mat1[i,n]=as.vector(emission_distribution(beta0_scal = beta0[i],
                                                  beta1_vec = beta1[i,],
                                                  states_vec = as.vector(states.mat1[i,,n]),
                                                  acts_scal = 1,
                                                  gamma_scal = gammas[i]))+eps.1
      
      
      states.mat1[i+1,,n]<-as.vector(trans_distribution(phi0_vec = phi0[i,],
                                                        phi1_mat = phi1[[i]],
                                                        state_vec = as.vector(states.mat1[i,,n]),
                                                        act_scal = 1,
                                                        alpha_vec = alphas[i,]))
      
      Y.mat0[i,n]=as.vector(emission_distribution(beta0_scal = beta0[i],
                                                  beta1_vec = beta1[i,],
                                                  states_vec = as.vector(states.mat0[i,,n]),
                                                  acts_scal = 0,
                                                  gamma_scal = gammas[i]))+eps.1
      
      
      states.mat0[i+1,,n]<-as.vector(trans_distribution(phi0_vec = phi0[i,],
                                                        phi1_mat = phi1[[i]],
                                                        state_vec = as.vector(states.mat0[i,,n]),
                                                        act_scal =0,
                                                        alpha_vec = alphas[i,]))
      
    }
    
    
  }
  
  
  
  ATE_true=mean(colSums(Y.mat1))-mean(colSums(Y.mat0))
  ATE_true=ATE_true/taus
  
  bias_difference[k]=ATE_true-ATE0
  Policy1_true[k]=mean(colSums(Y.mat1))
  Policy0_true[k]=mean(colSums(Y.mat0))
  
}


ATE0
ATE_true


