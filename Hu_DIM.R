library(KernSmooth)
library(ks)
library(splines)
library(mvtnorm)
library(KernSmooth)
library(ks)
library(splines)
library(fda)

compute_T_DM <- function(Y, Z, l, b) {
  
  # n=1
  # Y=Y_mat[,n]
  # Z=acts_mat[,n]
  # l=TI_loop; b= burn_in
  # 
  # 
  # 输入检查
  stopifnot(length(Y) == length(Z),
            l > 0, b >= 0, b < l,
            all(Z %in% c(0, 1)))
  
  # 总时间点数
  T_total <- length(Y)
  
  # 计算有效区块数量
  n_blocks <- T_total / l
  
  # 初始化存储结构
  block_sums_1 <- numeric(0)
  block_sums_0 <- numeric(0)
  
  count_Z<-0
  count_Z_vec<-rep(0,n_blocks)
  
  # 遍历每个区块
  for (i in 1:n_blocks) {
    
    #i=3
    # 提取当前区块的Y数据
    start_idx <- (i-1)*l + 1
    
    end_idx <- i*l
    
    current_Y <- Y[start_idx:end_idx]
    
    # 提取当前区块的Z值
    current_Z <- Z[start_idx]
    count_Z_vec[i]<-Z[start_idx]
    
    count_Z<-count_Z+current_Z
    
    # 计算b+1到l的和
    if (length(current_Y) >= b) {
      
      
      sum_val <- sum(current_Y[(b+1):l])
      
      # 按Z值分类存储
      if (current_Z == 1) {
        
        block_sums_1 <- c(block_sums_1, sum_val)
        
      } else {
        
        block_sums_0 <- c(block_sums_0, sum_val)
      }
      
    }
    
    
  }
  
  # 计算k值
  k1 <-  count_Z
  k0 <-  n_blocks- count_Z
  
  # # 计算最终结果
  # if (k1 == 0 || k0 == 0) {
  #   warning("One group has no observations. Returning NA")
  #   return(NA)
  # }
  # 
  if(k1==0){
    
    term1 <- 0
    term0 <- sum(block_sums_0) / (k0 * (l - b))
    
    T_DM <- term1 - term0
    
  }else if(k0==0){
    
    term1 <- sum(block_sums_1) / (k1 * (l - b))
    term0 <- 0
    
    T_DM <- term1 - term0
  }else{
    
    term1 <- sum(block_sums_1) / (k1 * (l - b))
    term0 <- sum(block_sums_0) / (k0 * (l - b))
    
    T_DM <- term1 - term0
    
  }
  
  
  return(T_DM)
}


#Open this folder directly under the current folder via R
source("ATE_true_approprimate.R") ##---parametric setting--
source("actions_dgp.R")##--treatment assignments---

library(mvtnorm)

rho0=c(0.3,0.5,0.7,0.9) # #---negative structure--- rho0=-c(0.3,0.5,0.7,0.9)

#TI<-c( 3, 6, 12,24,48)



library(KernSmooth)
library(ks)
library(splines)
library(tidyverse)

different_cases=expand.grid(rho0=rho0,TI=TI,Ndays=Ndays)

is_or_eigen






##parallel computation
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
  
  return(list(output_ATEmse, output_bias,
              output_sd
              
  ))
}



calculate_mseparameters<-function(loops, different_expcases, M, states_dim, taus ){
  
  # loops=197; different_expcases=different_cases;M=200;states_dim=p;taus=taus
  
  
  library(KernSmooth)
  library(ks)
  library(splines)
  library(tidyverse)
  library(mvtnorm)
  
  rho0_loop=different_expcases[loops,1]
  TI_loop=different_expcases[loops,2]
  Ndays_loop=different_expcases[loops,3]
  
  burn_in<-0
  
  eta_cors=matrix(nrow = taus,ncol = taus)
  
  #--autoregressive structure---
  
  for(i in 1:taus){
    for(j in 1:taus){
      eta_cors[i,j]=rho0_loop^(abs(i-j))
    }
  }
  
  #--moving average structure---
  
  # for(i in 1:taus){
  #   for(j in 1:taus){
  # 
  #     kappa_abs=abs(i-j)
  # 
  #     eta_cors[i,j]=((taus-kappa_abs)/taus)*rho0_loop
  # 
  #   }
  # }
  # 
  
  # #---exchangeable structure---
  # 
  #   for(i in 1:taus){
  #     for(j in 1:taus){
  # 
  #       if(i!=j){
  #         eta_cors[i,j]=rho0_loop
  #       }else{
  #         eta_cors[i,j]=1
  #       }
  # 
  #     }
  #   }
  
  
  #---uncorrelated structure---
  # 
  # for(i in 1:taus){
  #   for(j in 1:taus){
  # 
  #     if(i!=j){
  #       eta_cors[i,j]=0
  #     }else{
  #       eta_cors[i,j]=rho0_loop
  #     }
  # 
  #   }
  # }
  
  
  
  
  eta_cors_rescale= 1.5*eta_cors
  
  
  eps_cors=diag(x=rep(1.5,p),nrow = p)
  
  
  
  
  
  MSE_ATE=rep(0,M)
  
  ATE_est_vec=rep(0,M)
  
  
  
  ##set the same random seed
  
  
  for(m in 1:M){
    
    #m=10
    
    
    set.seed(2023+m)
    #generating the actions
    #act_init
    #act_init=rbinom(1,1,0.5)
    acts_mat=actions_dgp_Hu(taus = taus,ndays = Ndays_loop, ti=TI_loop)
    
    
    
    # m=1
    #generating data 
    Y_mat=matrix(nrow = taus, ncol = Ndays_loop)
    states_array=array(dim = c(taus, states_dim, Ndays_loop))#m*p*n
    
    
    
    
    for(n in 1:Ndays_loop){
      
      #n=1
      Y<-rep(0,taus)
      
      
      states_store=matrix(nrow = taus+1,ncol = states_dim)
      
      states_store[1,]<-as.vector(rmvnorm(1,mean = rep(0,states_dim),sigma = diag(x=states_dim)))
      
      eta_taus<-as.vector(rmvnorm(1,mean = rep(0,taus),sigma =eta_cors_rescale))
      
      
      for(i in 1:taus){
        
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
                                               alpha_vec = alphas[i,])+eps2
        
        
        
        
      }
      
      Y_mat[,n]=Y
      
      states_array[,,n]<-states_store[1:taus,]
      
      
    }
    
    
    
    ##---estimating ATE----
    #--estimating theta---
 
    
    ATE_burn_in<-rep(0, Ndays_loop)
    
    for(n in 1:Ndays_loop){
      
      
      ATE_burn_in[n]<-compute_T_DM(Y_mat[,n], acts_mat[,n], TI_loop,  burn_in  )
      
    }
    
    ATE_est<- mean(ATE_burn_in)
    
    
    
    
    
    
    ATE_est_vec[m]=ATE_est
    
    
    
    MSE_ATE[m]=(ATE_est-ATE_true)^2
    
    
    
    
  }
  
  mean_MSE_ATE=mean(MSE_ATE)
  
  
  mean_bias_ATE=mean(ATE_est_vec)-ATE_true
  mean_SD_ATE=sd(ATE_est_vec)
  
  
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
  ggtitle(paste("RMSE of Burn-in-DIM estimator")) +
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
  ggtitle("Bias of Burn-in-DIM estimator") +
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
  ggtitle("SD of Burn-in-DIM estimator") +
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









write.csv(data_ate_mse,"result_Hu_DIM_rmse_burn0_linear.csv")
write.csv(data_ate_bias,"result_Hu_DIM_bias_burn0_linear.csv")
write.csv(data_ate_sd,"result_Hu_DIM_SD_burn0_linear.csv")
