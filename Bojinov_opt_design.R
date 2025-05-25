library(KernSmooth)
library(ks)
library(splines)
library(mvtnorm)
library(KernSmooth)
library(ks)
library(splines)
library(fda)



source("ATE_true_approprimate.R") ##---parametric setting--
source("actions_dgp.R")##--treatment assignments---


rho0=c(0.3,0.5,0.7,0.9)
TI=c(1, 3, 6, 12, 24, 48)

library(KernSmooth)
library(ks)
library(splines)
library(tidyverse)

different_cases=expand.grid(rho0=rho0,TI=TI,Ndays=Ndays)

#is_or_eigen

library(parallel)
library(parallelMap)
library(foreach)
library(doParallel)
detectCores()
cl<-makeCluster(detectCores())
registerDoParallel(cl)

#parallelStartSocket(cpus = detectCores())

combin0<-function(List1,List2){
  output_ATEmse<-rbind(List1[[1]],List2[[1]])
  output_bias<-rbind(List1[[2]],List2[[2]])
  output_sd<-rbind(List1[[3]],List2[[3]])
  
  return(list(output_ATEmse,output_bias,
              output_sd
              
  ))
}




calculate_mseparameters<-function(loops, different_expcases, M, states_dim, taus ){
  
  # loops=91;  different_expcases=different_cases;M=200;states_dim=p;taus=taus
  
  library(KernSmooth)
  library(ks)
  library(splines)
  
  library(mvtnorm)
  library(KernSmooth)
  library(ks)
  library(splines)
  library(fda)
  library(tidyverse)
  
  rho0_loop=different_expcases[loops,1]
  TI_loop=different_expcases[loops,2]
  Ndays_loop=different_expcases[loops,3]
  
  eta_cors=matrix(nrow = taus,ncol = taus)
  
  
  for(i in 1:taus){
    for(j in 1:taus){
      eta_cors[i,j]=rho0_loop^(abs(i-j))
    }
  }
  
  
  eta_cors_rescale=1.5*eta_cors
  eps_cors=diag(x=rep(1.5,p),nrow = p)
  
  
  
  MSE_ATE=rep(0,M)
  
  ATE_est_vec=rep(0,M)
  
  
  
  
  for(m in 1:M){
    #m=1
    set.seed(2023+m)
    #generating actions--
    
    acts_mat=actions_dgp_Boinov(taus = taus,ndays = Ndays_loop, ti=TI_loop,seed = 2023+m)
    #acts_mat=actions_dgp2(taus = taus,ndays = Ndays_loop, ti=TI_loop )
    
    
    Y_mat=matrix(nrow = taus, ncol = Ndays_loop)#m*n
    states_array=array(dim = c(taus, states_dim, Ndays_loop))#m*p*n
    
    #matrix(nrow = taus,ncol = Ndays[ndays])
    
    
    for(n in 1:Ndays_loop){
      
      #n=1
      Y<-rep(0,taus)
      
      
      states_store=matrix(nrow = taus+1,ncol = states_dim)
      
      states_store[1,]<-as.vector(rmvnorm(1,mean = rep(0,states_dim),sigma = diag(x=states_dim)))
      
      eta_taus<-as.vector(rmvnorm(1,mean = rep(0,taus),sigma =eta_cors_rescale))
      
      
      
      for(i in 1:taus){
        #i=1
        
        eps1<-rnorm(1)
        eps2<-rmvnorm(1, mean = rep(0,p),sigma = eps_cors )
        eps2<-as.vector(eps2)
        #i=1
        Y[i]=emission_distribution(beta0_scal = beta0[i],
                                   beta1_vec = beta1[i,],
                                   states_vec = as.vector( states_store[i,] ),
                                   acts_scal = acts_mat[i,n] ,
                                   gamma_scal = gammas[i])+eta_taus[i]+eps1
        
        
        
        
        
        states_store[i+1,]<-trans_distribution(phi0_vec = phi0[i,],
                                               phi1_mat = phi1[[i]],
                                               state_vec = as.vector( states_store[i,] ),
                                               act_scal =acts_mat[i,n],
                                               alpha_vec = alphas[i,])+eps2  #+Eta_taus[,i]
        
        
        
        
      }
      
      Y_mat[,n]=Y
      
      states_array[,,n]<-states_store[1:taus,]
      
      
    }
    
    
    ##---estimating ATE with IS-HT estimator----
    
    IS_Y_mat=rep(0,Ndays_loop)
    for(i in 1:Ndays_loop){
      #i=1
      
      action_vec=acts_mat[,i]
      
      Y_vec=Y_mat[,i]
      
      #TI_loop
      
      IS_Y_vec=rep(0,taus)
      IS_Y1=0
      IS_Y0=0
      weight_act1=0
      weight_act0=0
      
      for(t in 1:taus){
        #t=1
        
        
        if(t<=TI_loop){
          
          IS_Y_vec[t] =(Y_vec[t]*all(action_vec[t]==1)/(1/2)) - all( Y_vec[t]*(action_vec[t]==0)/(1/2))
          
          
          IS_Y1=IS_Y1+(Y_vec[t]*all(action_vec[t]==1)/(1/2))
          IS_Y0=IS_Y0+(Y_vec[t]*all(action_vec[t]==0)/(1/2))
          
          weight_act1=weight_act1+all(action_vec[t]==1)/(1/2)
          weight_act0=weight_act0+all(action_vec[t]==0)/(1/2)
          
        }else{
          
          IS_Y_vec[t] =(Y_vec[t]*all(action_vec[(t-TI_loop):t]==1)/(1/4)) - (Y_vec[t]*all(action_vec[(t-TI_loop):t]==0)/(1/4))  
          
          IS_Y1=IS_Y1+(Y_vec[t]*all(action_vec[(t-TI_loop):t]==1)/(1/4))
          IS_Y0=IS_Y0+(Y_vec[t]*all(action_vec[(t-TI_loop):t]==0)/(1/4))
          
          weight_act1= weight_act1+(all(action_vec[(t-TI_loop):t]==1)/(1/4))
          weight_act0= weight_act0+(all(action_vec[(t-TI_loop):t]==0)/(1/4))
          
        }
        
        
      }
      
      
      
      # 
      #       if(TI_loop<taus){
      # 
      # 
      #        # IS_Y_mat[i] =mean(IS_Y_vec[(TI_loop+1):taus])
      # 
      #         #IS_Y_mat[i] =mean(IS_Y_vec)
      # 
      #          Dy1=sum(IS_Y1)/sum(weight_act1)
      #          Dy0=sum(IS_Y0)/sum(weight_act0)
      # 
      #          IS_Y_mat[i]=Dy1-Dy0
      # 
      #       }else{
      # 
      #        # IS_Y_mat[i] =mean(IS_Y_vec)
      # 
      #          Dy1=sum(IS_Y1)/sum(weight_act1)
      #          Dy0=sum(IS_Y0)/sum(weight_act0)
      # 
      #         if(sum(weight_act1)!=0){
      #           IS_Y_mat[i]=Dy1
      #         }else{
      #           IS_Y_mat[i]=-Dy0
      #         }
      # 
      # 
      #       }
      
      
      IS_Y_mat[i] =mean(IS_Y_vec,na.rm =T)
      
    }
    
    ATE_est=mean(IS_Y_mat,na.rm =T)
    
    
    # Y_longvec=vec(Y_mat)
    # act_longvec=vec(acts_mat)
    # 
    # ATE_est=mean(4*(act_longvec-1/2)*Y_longvec)
    # 
    
    ATE_est_vec[m]=ATE_est
    
    
    
    MSE_ATE[m]=(ATE_est-ATE_true)^2
    
    
    
  }
  
  
  mean_MSE_ATE=mean(MSE_ATE,na.rm =T)
  
  mean_bias_ATE=mean(ATE_est_vec,na.rm =T)-ATE_true
  
  mean_SD_ATE=sd(ATE_est_vec,na.rm =T)
  
  
  return(list(
    
    mean_MSE_ATE,
    mean_bias_ATE,mean_SD_ATE
    
  ))
}


nrow(different_cases)


timestart<-Sys.time()
final_results_list<-foreach(m_loop=1:nrow(different_cases), .combine = "combin0")%dopar%{
  
  calculate_mseparameters(loops =m_loop, different_expcases = different_cases,
                          M=200,states_dim = p,taus = taus)
  
}




timeend<-Sys.time()
runningtime<-timeend-timestart
runningtime
stopCluster(cl)



data_ate_mse<-cbind(different_cases, MSE_ate=final_results_list[[1]])

data_ate_bias<-cbind(different_cases, bias_ate=final_results_list[[2]])

data_ate_sd<-cbind(different_cases, SD_ate=final_results_list[[3]])


data_ate_mse=as_tibble(data_ate_mse)

data_ate_bias=as_tibble(data_ate_bias)
data_ate_sd=as_tibble(data_ate_sd)



data_ate_mse<-mutate_at( data_ate_mse,c("rho0","TI"),as_factor)
data_ate_mse<-mutate( data_ate_mse, ate_rmse=sqrt(MSE_ate) )



data_ate_bias<-mutate_at(data_ate_bias,c("rho0","TI"),as_factor)
data_ate_sd<-mutate_at(data_ate_sd,c("rho0","TI"),as_factor)





data_ate_mse=as_tibble(data_ate_mse)
data_ate_mse<-mutate_at( data_ate_mse,c("rho0","TI"),as_factor)

names(data_ate_mse)<-c("rho0","m","n","MSE","RMSE")
data_ate_mse


ggplot(data_ate_mse, aes(x=n,y=RMSE, color=m) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  theme_bw()+
  ggtitle(paste("RMSE of HT estimator")) +
  theme(plot.title = element_text(vjust = 2) ) +
  theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
                                                      axis.text=element_text(size=14,face = "bold"),
                                                      axis.title.x=element_text(size=14,face="bold"),
                                                      axis.title.y=element_text(size=14,face='bold'),
                                                      axis.line = element_line(size=1),
                                                      axis.text.x = element_text(size=14),
                                                      legend.text = element_text(size=15),
                                                      legend.title=element_text(size=15),
                                                      strip.text = element_text(size=rel(1.5)))





data_ate_bias=as_tibble(data_ate_bias)
data_ate_bias<-mutate_at( data_ate_bias,c("rho0","TI"),as_factor)

names(data_ate_bias)<-c("rho0","m","n","Bias")
data_ate_bias

ggplot(data_ate_bias, aes(x=n,y=Bias, color=m) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  theme_bw()+
  ggtitle("Bias of HT estimator") +
  theme(plot.title = element_text(vjust = 2) ) +
  theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
                                                      axis.text=element_text(size=14,face = "bold"),
                                                      axis.title.x=element_text(size=14,face="bold"),
                                                      axis.title.y=element_text(size=14,face='bold'),
                                                      axis.line = element_line(size=1),
                                                      axis.text.x = element_text(size=14),
                                                      legend.text = element_text(size=15),
                                                      legend.title=element_text(size=15),
                                                      strip.text = element_text(size=rel(1.5)))







data_ate_sd<-mutate_at( data_ate_sd,c("rho0","TI"),as_factor)

names(data_ate_sd)<-c("rho0","m","n","SD")

ggplot(data_ate_sd, aes(x=n,y=SD, color=m) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  theme_bw()+
  ggtitle("SD of HT estimator") +
  theme(plot.title = element_text(vjust = 2)) +
  theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
                                                      axis.text=element_text(size=14,face = "bold"),
                                                      axis.title.x=element_text(size=14,face="bold"),
                                                      axis.title.y=element_text(size=14,face='bold'),
                                                      axis.line = element_line(size=1),
                                                      axis.text.x = element_text(size=14),
                                                      legend.text = element_text(size=15),
                                                      legend.title=element_text(size=15),
                                                      strip.text = element_text(size=rel(1.5)))






#rmse1：all-HT
#rmse2: IS-HT
#rmse3: all-Hajek



write.csv(data_ate_mse,"result_HT_rmse_linear.csv")
write.csv(data_ate_bias,"result_HT_bias_linear.csv")
write.csv(data_ate_sd,"result_HT_SD_linear.csv")

