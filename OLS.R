library(KernSmooth)
library(ks)
library(splines)
library(mvtnorm)
library(KernSmooth)
library(ks)
library(splines)
library(fda)

#Open this folder directly under the current folder via R
source("ATE_true_approprimate.R") ##---parametric setting--
source("actions_dgp.R")##--treatment assignments---

library(mvtnorm)

rho0=c(0.3,0.5,0.7,0.9) # #---negative structure--- rho0=-c(0.3,0.5,0.7,0.9)


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
  output_DEmse<-rbind(List1[[2]],List2[[2]])
  output_IEmse<-rbind(List1[[3]],List2[[3]])
  output_bias<-rbind(List1[[4]],List2[[4]])
  output_sd<-rbind(List1[[5]],List2[[5]])
  
  return(list(output_ATEmse, output_DEmse, output_IEmse,output_bias,
              output_sd
              
  ))
}



calculate_mseparameters<-function(loops, different_expcases, M, states_dim, taus ){
  
 # loops=1; different_expcases=different_cases;M=200;states_dim=p;taus=taus
  
  library(KernSmooth)
  library(ks)
  library(splines)
  library(tidyverse)
  library(mvtnorm)
  
  rho0_loop=different_expcases[loops,1]
  TI_loop=different_expcases[loops,2]
  Ndays_loop=different_expcases[loops,3]
  
  
  
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
  MSE_DE=rep(0,M)
  MSE_IE=rep(0,M)
  
  ATE_est_vec=rep(0,M)
  
  
  
  ##set the same random seed
  
  
  for(m in 1:M){
    
    set.seed(2023+m)
    #generating the actions
    act_init=rbinom(1,1,0.5)
    #act_init
    acts_mat=actions_dgp(taus = taus,ndays = Ndays_loop, ti=TI_loop,act_init=act_init )

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
    est_theta=matrix(nrow = taus,ncol = states_dim+2)#m*(p+2)
    
    for(i in 1:taus){
      #i=1
      #(p+2)*n
      Z_col=rbind(rep(1,Ndays_loop ),states_array[i,,],acts_mat[i,])
      
      erro_ridges=0.
      
      ridges_mat=diag(x=erro_ridges,nrow = nrow(Z_col))
      
      est_thetai=solve((Z_col%*%t(Z_col)+ridges_mat))%*%(Z_col%*%Y_mat[i,])
      est_theta[i,]=est_thetai
      
      
    }
    
    
    
    
    
    #  theta0=cbind(beta0,beta1,gammas)
    
  
    gamma_est=est_theta[,states_dim+2]
    
    DE_est=sum(gamma_est);#DE_est
    DE_est=DE_est/taus
  
    
    
    est_theta_state=vector("list",length = taus-1)
    
    
    for(i in 1:(taus-1)){
      #i=1
      #(p+2)*n
      Z_col=rbind(rep(1, Ndays_loop),states_array[i,,], acts_mat[i,])
      
      erro_ridges=0.0
      
      ridges_mat=diag(x=erro_ridges,nrow = nrow(Z_col))
      # p*(p+2)
      est_thetai_state=solve((Z_col%*%t(Z_col)+ridges_mat))%*%(Z_col%*%t(states_array[i+1,,]))#(p+2)*p
      
      
      est_theta_state[[i]]=t(est_thetai_state)#p*(p+2)
      
      
    }
    
    #est_theta_state
    
    beta_est<-est_theta[,2:(states_dim+1)]
    
    phi1_est<-vector("list",length = taus-1)
    alphas_est<-matrix(nrow = taus-1, ncol = states_dim )
    phi0_est<-matrix(nrow = taus-1, ncol = states_dim )
    
    for(i in 1:(taus-1)){
      
      
      
      phi0_est[i,]<-as.vector(est_theta_state[[i]][,1])
      phi1_est[[i]]=est_theta_state[[i]][, 2:(states_dim+1) ]
      
      alphas_est[i,]<-as.vector(est_theta_state[[i]][, states_dim+2 ])
      
    }  
    
    
    

    
    IE_est=IE_estfun(beta_hat = beta_est, phi_hat = phi1_est, alpha_hat = alphas_est )
    IE_est=IE_est/taus
    
   
    ATE_est=DE_est+IE_est
    
    
    ATE_est_vec[m]=ATE_est
    
    
    
    MSE_ATE[m]=(ATE_est-ATE_true)^2
    
    MSE_DE[m]=(DE_est-DE0)^2
    
    MSE_IE[m]=(IE_est-IE0)^2
    
    
  }
  
  mean_MSE_ATE=mean(MSE_ATE)
  mean_MSE_DE=mean(MSE_DE)
  mean_MSE_IE=mean(MSE_IE)
  
  mean_bias_ATE=mean(ATE_est_vec)-ATE_true
  mean_SD_ATE=sd(ATE_est_vec)
  
  
  return(list(
    
    mean_MSE_ATE,mean_MSE_DE,mean_MSE_IE,
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

data_DE_mse<-cbind(different_cases, MSE_DE=final_results_list[[2]])

data_IE_mse<-cbind(different_cases, MSE_IE=final_results_list[[3]])

data_ate_bias<-cbind(different_cases, bias_ate=final_results_list[[4]])

data_ate_sd<-cbind(different_cases, SD_ate=final_results_list[[5]])


data_ate_mse=as_tibble(data_ate_mse)
data_DE_mse=as_tibble(data_DE_mse)
data_IE_mse=as_tibble(data_IE_mse)
data_ate_bias=as_tibble(data_ate_bias)
data_ate_sd=as_tibble(data_ate_sd)



data_ate_mse<-mutate_at( data_ate_mse,c("rho0","TI"),as_factor)
data_ate_mse<-mutate( data_ate_mse, ate_rmse=sqrt(MSE_ate) )



data_DE_mse<-mutate_at(data_DE_mse,c("rho0","TI"),as_factor)
data_DE_mse<-mutate(data_DE_mse,DE_rmse=sqrt(MSE_DE))

data_IE_mse<-mutate_at(data_IE_mse,c("rho0","TI"),as_factor)

data_IE_mse<-mutate(data_IE_mse,IE_rmse=sqrt(MSE_IE) )


data_ate_bias<-mutate_at(data_ate_bias,c("rho0","TI"),as_factor)
data_ate_sd<-mutate_at(data_ate_sd,c("rho0","TI"),as_factor)



library(tidyverse)

data_ate_mse=as_tibble(data_ate_mse)
data_ate_mse<-mutate_at( data_ate_mse,c("rho0","TI"),as_factor)

names(data_ate_mse)<-c("rho0","m","n","MSE","RMSE")
data_ate_mse


library(tidyverse)
library(ggplot2)
#data_ate_mse=read.csv("numerical_datasets_unit/nonlinear_model/result_OLS_rmse.csv")[,-1]

#data_ate_mse=read.csv("/Users/qianglin/Desktop/Switch_MDP_codes/corollary_figs/nonlinear_model/result_OLS_rmse_uncorrelated.csv")[,-1]



# 假设 data_ate_mse 数据已经准备好了
data_ate_mse <- mutate_at(data_ate_mse, c("rho0", "m"), as_factor)

# 为每个 rho0 取值定义数学公式字符串
# f_names <- list( ##--negative correlated---
#   "-0.9" = "rho == -0.9",
#   "-0.7" = "rho == -0.7",
#   "-0.5" = "rho == -0.5",
#   "-0.3" = "rho == -0.3"
#   
# )

f_names <- list(
  "0.3" = "rho == 0.3",
  "0.5" = "rho == 0.5",
  "0.7" = "rho == 0.7",
  "0.9" = "rho == 0.9"
)


# 定义 labeller 函数
f_labeller <- function(variable, value) {
  return(f_names[value])
}

# 绘制图形
ggplot(data_ate_mse, aes(x = n, y = RMSE, color = m)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  facet_grid(. ~ rho0, 
             labeller = labeller(rho0 = f_labeller,
                                 .default = label_parsed),  # 使用自定义的 labeller
             scales = "free") +
  ggtitle("RMSE of OLS Estimator") +
  theme_bw() +
  theme(plot.title = element_text(vjust = 2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.line = element_line(size = 1),
        axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = rel(1.5)))  # 使用 parse = TRUE 来解析公式



ggplot(data_DE_mse, aes(x=Ndays,y=DE_rmse, color=TI) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  ggtitle(paste("DE for effects (OLS without kernelsmooth): DE0=",round(DE0,digits = 3),sep="")) +
  theme_bw()+
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





ggplot(data_IE_mse, aes(x=Ndays,y=IE_rmse, color=TI) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  ggtitle(paste("IE for effects (OLS without kernelsmooth): IE0=",round(IE0,digits = 3),sep="")) +
  theme_bw()+
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





data_ate_bias=as_tibble(data_ate_bias)
data_ate_bias<-mutate_at( data_ate_bias,c("rho0","TI"),as_factor)

names(data_ate_bias)<-c("rho0","m","n","Bias")
data_ate_bias
ggplot(data_ate_bias, aes(x=n,y=Bias, color=m) )+geom_line(size=1)+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  ggtitle("Bias of OLS estimator") +
  theme_bw()+
  #theme(plot.title = element_text(vjust = 2),legend.position = "none") +
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







data_ate_sd<-mutate_at( data_ate_sd,c("rho0","TI"),as_factor)

names(data_ate_sd)<-c("rho0","m","n","SD")

ggplot(data_ate_sd, aes(x=n,y=SD, color=m) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  ggtitle("SD of OLS estimator") +
  theme_bw()+
  #theme(plot.title = element_text(vjust = 2),legend.position = "none") +
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






##--changing different output with corresponding setting--

write.csv(data_ate_mse,"result_OLS_rmse_3fold_IE.csv")
write.csv(data_ate_bias,"result_OLS_bias_3fold_IE.csv")
write.csv(data_ate_sd,"result_OLS_SD_3fold_IE.csv")

# write.csv(data_ate_mse,"result_OLS_rmse_uncorrelated.csv")
# write.csv(data_ate_bias,"result_OLS_bias_uncorrelated.csv")
# write.csv(data_ate_sd,"result_OLS_SD_uncorrelated.csv")














