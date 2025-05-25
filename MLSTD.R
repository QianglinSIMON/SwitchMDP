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


#--LSTD with some auxiliary functions--

source("self_basis_functions.R")

rho0=c(0.3,0.5,0.7,0.9)

library(KernSmooth)
library(ks)
library(splines)
library(tidyverse)

different_cases=expand.grid(rho0=rho0,TI=TI,Ndays=Ndays)

is_or_eigen

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
  
  # loops=97;  different_expcases=different_cases;M=200;states_dim=p;taus=taus
  
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
      # m=1
    
    set.seed(2023+m)
    #生成行动策略
    act_init=rbinom(1,1,0.5)
    #act_init
    acts_mat=actions_dgp(taus = taus,ndays = Ndays_loop, ti=TI_loop,act_init=act_init )
 
    
     #m=1
    

    Y_mat=matrix(nrow = taus, ncol = Ndays_loop)#m*n
    states_array=array(dim = c(taus, states_dim, Ndays_loop))#m*p*n
    
    #matrix(nrow = taus,ncol = Ndays[ndays])
    
    
    for(n in 1:Ndays_loop){
      
      #n=1
      Y<-rep(0,taus)
      
      
      states_store=matrix(nrow = taus+1,ncol = states_dim)
      
      states_store[1,]<-as.vector(rmvnorm(1,mean = rep(0,states_dim),sigma = diag(x=states_dim)))
      
      eta_taus<-as.vector(rmvnorm(1,mean = rep(0,taus),sigma =eta_cors_rescale))
      
      # Eta_taus<-rmvnorm(p,mean = rep(0,taus),sigma =Eta_cors_rescale)#状态变量中random effects 
      
      
      
      
      
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
    
    
    ##---estimating ATE----
    library(fda)
    Ln=5; Kn=3
   
    Y_mat_all_list=vector("list",taus)
    
    # states_all_list=state_polybasis_time(states_arrays = states_array,taus = taus,
    #                                      Ln=4, Kn=2)
    
    states_all_list=state_orthogonal_basis_time(states_arrays = states_array,taus = taus,
                                                Ln=Ln, Kn=Kn)
    
    
   # toel_rence=0.5
   
    # toel_rence=0.1
    
    for(i in 1:taus){
      Y_mat_all_list[[i]] =Y_mat[i,]
      
    }
    
    sbasismat_all_list=states_all_list
    
    
    
    
    row_numbers1_list=vector("list", taus)
    row_numbers0_list=vector("list", taus)
    
    
    for(i in 1:taus){
      
      # i=1
      row_numbers1_list[[i]]=which(acts_mat[i,]==1)
      row_numbers0_list[[i]]=which(acts_mat[i,]==0)
      
    }
    
    
    ##---action =1 ----
    
    Hat_sigma1_list=vector("list",length = taus)
    Hat_sigma1_corss_list=vector("list",length = taus)
    
    hat_response1_list=vector("list",length = taus)
    
    ##------ back reduction  --
    
    for( i in 1:(taus)){
      
      #i=1
      
      
      state_biasis1_mat=sbasismat_all_list[[i]][row_numbers1_list[[i]],] 
      
    
      
      hat_sigma1=t(state_biasis1_mat)%*%state_biasis1_mat
      
      Hat_sigma1_list[[i]]=hat_sigma1
      
      if(i < taus){
        
        #state_biasis1_mat_lag=sbasismat_all_list[[i+1]][row_numbers1_list[[i+1]],]
        
        state_biasis1_mat_lag=sbasismat_all_list[[i+1]][row_numbers1_list[[i]],]
        
        Hat_sigma1_corss_list[[i]]=t(state_biasis1_mat)%*%state_biasis1_mat_lag
        
        #Hat_sigma1_corss_list[[i]]=colMeans(state_biasis1_mat)%*%t(colMeans(state_biasis1_mat_lag))
        
        #Hat_sigma1_corss_list[[i]]=t(state_biasis1_mat)%*%state_biasis1_mat_lag
        
      }else{
        
        #state_biasis1_mat_lag=sbasismat_all_list[[i]][row_numbers1_list[[i]],]
        #state_biasis1_mat_lag=0
        
        Hat_sigma1_corss_list[[i]]=0
        
      }
        
     
      
      hat_response1=as.vector(t(state_biasis1_mat)%*%Y_mat_all_list[[i]][row_numbers1_list[[i]]])
      

      hat_response1_list[[i]]= hat_response1
      
      
    }
    

    
    library(purrr)
    
    Y_response1= Reduce("+",hat_response1_list)/(Ndays_loop*taus)
    
    Hat_sigma1_zero=Reduce("+", Hat_sigma1_list)/(Ndays_loop*taus)
    
    Hat_sigma1_one=Reduce("+", Hat_sigma1_corss_list )/(Ndays_loop*taus)
    
    Hat_sigma1_zero_inverse=solve(Hat_sigma1_zero)
   
    #round( det(Hat_sigma1_zero),digits = 4)
    
    Hat_beta1_zero=Hat_sigma1_zero_inverse%*%Y_response1
    
    unit_mat1= Hat_sigma1_zero_inverse%*%Hat_sigma1_one
    
    library(expm)
    
   # x%^%0
    
    unit1_mat_list=vector("list",length = taus)
    
    for(i in 1:taus){
      
      
      unit1_mat_list[[i]]=unit_mat1%^%(i-1)
      
    }
    
    unit1_mat_sum=Reduce("+",unit1_mat_list)
    
    Hat_beta1_final=unit1_mat_sum%*%Hat_beta1_zero
    
    
    ##---action =0 ----
    
    Hat_sigma0_list=vector("list",length = taus)
    Hat_sigma0_corss_list=vector("list",length = taus)
    
    hat_response0_list=vector("list",length = taus)
    
    ##------ back reduction  --
    
    for( i in 1:(taus)){
      
      #i=1
      
      
      state_biasis0_mat=sbasismat_all_list[[i]][row_numbers0_list[[i]],] #第 i 个时间点行动为0 的基函数矩阵 
      
      
      
      hat_sigma0=t(state_biasis0_mat)%*%state_biasis0_mat
      
      Hat_sigma0_list[[i]]=hat_sigma0
      
      if(i < taus){
        
        #state_biasis0_mat_lag=sbasismat_all_list[[i+1]][row_numbers0_list[[i+1]],]
        
        state_biasis0_mat_lag=sbasismat_all_list[[i+1]][row_numbers0_list[[i]],]
        Hat_sigma0_corss_list[[i]]=t(state_biasis0_mat)%*%state_biasis0_mat_lag
        
        #Hat_sigma0_corss_list[[i]]=colMeans(state_biasis0_mat)%*%t(colMeans(state_biasis0_mat_lag))
        
        #Hat_sigma0_corss_list[[i]]=t(state_biasis0_mat)%*%state_biasis0_mat_lag
        
      }else{
        
        #state_biasis0_mat_lag=sbasismat_all_list[[i]][row_numbers0_list[[i]],]
       # state_biasis0_mat_lag=0
        
        Hat_sigma0_corss_list[[i]]=0
        
      }
      
      
      
      hat_response0=as.vector(t(state_biasis0_mat)%*%Y_mat_all_list[[i]][row_numbers0_list[[i]]])
      
      
      hat_response0_list[[i]]= hat_response0
      
      
    }
    
    
    
    library(purrr)
    
    Y_response0=Reduce("+",hat_response0_list)/(Ndays_loop*taus)
    
    Hat_sigma0_zero=Reduce("+", Hat_sigma0_list)/(Ndays_loop*taus)
    
   # det( Hat_sigma0_zero)
      
    Hat_sigma0_one=Reduce("+", Hat_sigma0_corss_list )/(Ndays_loop*taus)
    
    Hat_sigma0_zero_inverse=solve(Hat_sigma0_zero)
    
    Hat_beta0_zero=Hat_sigma0_zero_inverse%*%Y_response0
    
    unit_mat0= Hat_sigma0_zero_inverse%*%Hat_sigma0_one
    
    library(expm)
    
    # x%^%0
    
    unit0_mat_list=vector("list",length = taus)
    
    for(i in 1:taus){
      
      
      unit0_mat_list[[i]]=unit_mat0%^%(i-1)
      
    }
    
    unit0_mat_sum=Reduce("+",unit0_mat_list)
    
    Hat_beta0_final=unit0_mat_sum%*%Hat_beta0_zero
    
    
    
    
      # Ndays_loop
    
    
    ate_npar_est=mean(sbasismat_all_list[[1]]%*%(Hat_beta1_final- Hat_beta0_final))
    
    

    # ate_npar_est
    
    ATE_est=ate_npar_est/taus
    
    
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



# library(tidyverse)
# library(ggplot2)
# data_ate_mse=read.csv("numerical_datasets_unit/nonlinear_model/result_MLSTD_rmse.csv")[,-1]



# 假设 data_ate_mse 数据已经准备好了
data_ate_mse <- mutate_at(data_ate_mse, c("rho0", "m"), as_factor)

# 为每个 rho0 取值定义数学公式字符串
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
  ggtitle("RMSE of LSTD Estimator") +
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


data_ate_bias=as_tibble(data_ate_bias)
data_ate_bias<-mutate_at( data_ate_bias,c("rho0","TI"),as_factor)

names(data_ate_bias)<-c("rho0","m","n","Bias")
data_ate_bias

ggplot(data_ate_bias, aes(x=n,y=Bias, color=m) )+geom_point(size=5)+geom_line(size=2)+
  facet_grid(.~rho0)+
  theme_bw()+
  ggtitle("Bias of LSTD estimator") +
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
  ggtitle("SD of LSTD estimator") +
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










write.csv(data_ate_mse,"result_MLSTD_rmse_3fold_IE.csv")
write.csv(data_ate_bias,"result_MLSTD_bias_3fold_IE.csv")
write.csv(data_ate_sd,"result_MLSTD_SD_3fold_IE.csv")

