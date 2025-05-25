##---estimating the marginal density functions---
estimate_state_density<-function(state_vec, phi1_lists, alpha_vecs, phi0_vecs, A_vecs, time_point){
  
  library(tidyverse)
  
  if(time_point>1){
    
    meanstate_iter1=0
    meanstate_iter2=0
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        meanstate_iter1=phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
      }else{
        
        meanstate_iter1=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
        
      }
      
    }
    
    meanstate=meanstate_iter1+meanstate_iter2
    
    phi1_product=as.matrix(Reduce("%*%", phi1_lists[1:(time_point-1)]))
    varstate_iter1=phi1_product%*%t(phi1_product)
    
    varstate_iter2=0
    
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        varstate_iter2=eps_cors+varstate_iter2
        
      }else{
        
        phi1_1_product=as.matrix(Reduce("%*%", phi1_lists[(k+1):(time_point-1)]))
        
        varstate_iter2=phi1_1_product%*%eps_cors%*%t(phi1_1_product)+varstate_iter2
        
      }
      
    }
    
    varstate=varstate_iter1+varstate_iter2
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =meanstate, sigma =varstate  )
  }else{
    
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =rep(0,length(state_vec)), sigma =diag(length(state_vec)) )
  }
  
  
  return(margin_density)
  
}


estimate_state_density_sample<-function(state_vec, phi1_lists, alpha_vecs, phi0_vecs, A_vecs, time_point,
                                        var_state_init,
                                        mean_state_init,
                                        eps_cors_list,
                                        mean_cors_list
                                        
                                        ){
  
  library(tidyverse)
  
  if(time_point>1){
    
    meanstate_iter1=0
    meanstate_iter2=0
    meanstate_iter4=0
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        meanstate_iter1=phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
        meanstate_iter4=mean_cors_list[[k]]+meanstate_iter4
        
      }else{
        
        meanstate_iter1=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
        meanstate_iter4=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%mean_cors_list[[k]]+meanstate_iter4
        
      }
      
    }

    
    phi1_product=as.matrix(Reduce("%*%", phi1_lists[1:(time_point-1)]))
    
    meanstate_iter3=phi1_product%*%mean_state_init
      
    meanstate=meanstate_iter1+meanstate_iter2+meanstate_iter3+meanstate_iter4
    
    
    varstate_iter1=phi1_product%*%var_state_init%*%t(phi1_product)
    
    varstate_iter2=0
    
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        varstate_iter2=eps_cors_list[[k]]+varstate_iter2
        
      }else{
        
        phi1_1_product=as.matrix(Reduce("%*%", phi1_lists[(k+1):(time_point-1)]))
        
        varstate_iter2=phi1_1_product%*%eps_cors_list[[k]]%*%t(phi1_1_product)+varstate_iter2
        
      }
      
    }
    
    varstate=varstate_iter1+varstate_iter2
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =meanstate, sigma =varstate  )
  }else{
    
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =mean_state_init, sigma =var_state_init )
  }
  
  
  return(margin_density)
  
}



DRL_estfunction<-function( taus, ndays, states_data, acts_data,Y_matdata,Kn,Ln){
  
  # k=1
  # taus = taus; ndays = Ndays_loop/2; states_data = states_block__list[[k]]
  # 
  # acts_data = acts_block_list[[k]]; 
  # Y_matdata =Y_block_list[[k]]
  # 
  
  

  
  row_numbers1_list=vector("list", taus)
  row_numbers0_list=vector("list", taus)
  
  
  for(i in 1:taus){
    
    # i=1
    row_numbers1_list[[i]]=which(acts_data[i,]==1) 
    row_numbers0_list[[i]]=which(acts_data[i,]==0)
    
  }
  
  ##--- estimating value function ----

  library(fda)
  #Ln=3; Kn=1

  Y_mat_all_list=vector("list",taus)
  
  
  
  
  states_all_list=state_orthogonal_basis_time(states_arrays = states_data,taus = taus,
                                              Ln=Ln, Kn=Kn)
  
  toel_rence=0.01
  
  # toel_rence=0.1
  
  for(i in 1:taus){
    Y_mat_all_list[[i]] =Y_matdata[i,]
    
  }
  
  sbasismat_all_list=states_all_list
  
  
  
  ##---action =1 ----
  
  state_biasis1_mat_m=sbasismat_all_list[[taus]][row_numbers1_list[[taus]],]
  
  hat_sigma1_m=(1/nrow(state_biasis1_mat_m))*t(state_biasis1_mat_m)%*%state_biasis1_mat_m 
  
  hat_response1_m=(1/nrow(state_biasis1_mat_m))*as.vector(t(state_biasis1_mat_m)%*%Y_mat_all_list[[taus]][row_numbers1_list[[taus]]])
  
  library(matrixcalc)
  if(is.singular.matrix(hat_sigma1_m)==T){
    
    hat_beta1_m=solve(hat_sigma1_m + diag(x=toel_rence,nrow =nrow(hat_sigma1_m) ))%*%hat_response1_m
  }else{
    
    hat_beta1_m=solve(hat_sigma1_m + diag(x=0.,nrow =nrow(hat_sigma1_m) ))%*%hat_response1_m
  }
  
  
  
  
  
  Hat_beta1_mat=matrix(nrow = taus,ncol = length(hat_beta1_m))
  
  Hat_beta1_mat[1,]=hat_beta1_m
  
  
  
  Hat_sigma1_list=vector("list",length = taus)
  Hat_sigma1_corss_list=vector("list",length = taus)
  
  Hat_sigma1_list[[1]]=hat_sigma1_m
  
  
  ##------ back reduction  --
  
  for( i in 1:(taus-1)){
    #i=1
    
    
    state_biasis1_mat=sbasismat_all_list[[taus-i]][row_numbers1_list[[taus-i]],]
    
    state_biasis1_mat_lag=sbasismat_all_list[[taus-i+1]][row_numbers1_list[[taus-i+1 ]],]
    
    hat_sigma1=(1/nrow(state_biasis1_mat))*t(state_biasis1_mat)%*%state_biasis1_mat
    
    Hat_sigma1_list[[i+1]]=hat_sigma1
    
    Hat_sigma1_corss_list[[i]]=(1/nrow(state_biasis1_mat))*t(state_biasis1_mat)%*%state_biasis1_mat_lag
    
    hat_response1=(1/nrow(state_biasis1_mat))*as.vector(t(state_biasis1_mat)%*%Y_mat_all_list[[taus-i]][row_numbers1_list[[taus-i]]])
    
    
    if(is.singular.matrix( hat_sigma1 )){
      
      hat_sigma1_inv=solve(hat_sigma1 + diag(x=toel_rence,nrow =nrow(hat_sigma1) ))
    }else{
      
      hat_sigma1_inv=solve(hat_sigma1 + diag(x=0.,nrow =nrow(hat_sigma1) ))
    }
    
    
    
    
    
    hat_beta1=hat_sigma1_inv%*% hat_response1+ 
      hat_sigma1_inv%*%Hat_sigma1_corss_list[[i]]%*%Hat_beta1_mat[i,]
    
    hat_beta1=as.vector(hat_beta1)
    
    Hat_beta1_mat[i+1,]=hat_beta1
    
  }
  
  
  
  
  ##---action =0 ----
  state_biasis0_mat_m=sbasismat_all_list[[taus]][row_numbers0_list[[taus]],]
  
  hat_sigma0_m=(1/nrow(state_biasis0_mat_m))*t(state_biasis0_mat_m)%*%state_biasis0_mat_m 
  
  hat_response0_m=(1/nrow(state_biasis0_mat_m))*as.vector(t(state_biasis0_mat_m)%*%Y_mat_all_list[[taus]][row_numbers0_list[[taus]]])
  
  if(is.singular.matrix(hat_sigma0_m)==T){
    
    hat_beta0_m=solve(hat_sigma0_m + diag(x=toel_rence,nrow =nrow(hat_sigma0_m) ))%*%hat_response0_m
  }else{
    
    hat_beta0_m=solve(hat_sigma0_m + diag(x=0.,nrow =nrow(hat_sigma0_m) ))%*%hat_response0_m
  }
  
  
  
  Hat_beta0_mat=matrix(nrow = taus,ncol = length(hat_beta0_m))
  
  Hat_beta0_mat[1,]=hat_beta0_m
  
  
  
  Hat_sigma0_list=vector("list",length = taus)
  Hat_sigma0_corss_list=vector("list",length = taus)
  
  Hat_sigma0_list[[1]]=hat_sigma0_m
  
  
  ##------ back reduction  --
  
  for( i in 1:(taus-1)){
    #i=1
    
    
    state_biasis0_mat=sbasismat_all_list[[taus-i]][row_numbers0_list[[taus-i]],]
    
    state_biasis0_mat_lag=sbasismat_all_list[[taus-i+1]][row_numbers0_list[[taus-i+1]],]
    
    hat_sigma0=(1/nrow(state_biasis0_mat))*t(state_biasis0_mat)%*%state_biasis0_mat
    
    Hat_sigma0_list[[i+1]]=hat_sigma0
    
    Hat_sigma0_corss_list[[i]]=(1/nrow(state_biasis0_mat))*t(state_biasis0_mat)%*%state_biasis0_mat_lag
    
    hat_response0=(1/nrow(state_biasis0_mat))*as.vector(t(state_biasis0_mat)%*%Y_mat_all_list[[taus-i]][row_numbers0_list[[taus-i]]])
    
    
    
    if(is.singular.matrix(hat_sigma0)==T){
      
      hat_sigma0_inv=solve(hat_sigma0 + diag(x=toel_rence,nrow =nrow(hat_sigma0) ))
    }else{
      
      hat_sigma0_inv=solve(hat_sigma0 + diag(x=0.,nrow =nrow(hat_sigma0) ))
    }
    
    # 
    #       if(det(hat_sigma0)<1e-05){
    # 
    #         hat_sigma0_inv=solve(hat_sigma0 + diag(x=toel_rence, nrow =nrow(hat_sigma0) ))
    #       }else{
    # 
    #         hat_sigma0_inv=solve(hat_sigma0 + diag(x=0.,nrow =nrow(hat_sigma0) ))
    #       }
    
    
    
    hat_beta0=hat_sigma0_inv%*% hat_response0+ 
      hat_sigma0_inv%*%Hat_sigma0_corss_list[[i]]%*%Hat_beta0_mat[i,]
    
    hat_beta0=as.vector(hat_beta0)
    
    Hat_beta0_mat[i+1,]=hat_beta0
    
  }
  
  
  
  results_list=list(
    
    Hat_beta1_mat=Hat_beta1_mat,Hat_beta0_mat=Hat_beta0_mat
  )
  
  
  
  
  return(results_list)
  
}


library(MASS)

actions_covmat<-function(ti_span,time_dim){
  cov_A<-matrix(nrow = time_dim,ncol = time_dim)
  for(i in 1:time_dim){
    for(j in 1:time_dim){
      v_1=i%%ti_span; u_1=i%/%ti_span
      k_1=ifelse(v_1==0,u_1,u_1+1)
      
      v_2=j%%ti_span; u_2=j%/%ti_span
      k_2=ifelse(v_2==0,u_2,u_2+1)
      
      if(abs(k_1-k_2)%%2==0){
        
        cov_A[i,j]=1/4
      }else{
        cov_A[i,j]=-1/4
      }
      
    }
  }
  return(cov_A)
}


estimate_state_density_var_target<-function(state_vec, phi1_lists, alpha_vecs, phi0_vecs, A_vecs, time_point
){
  
  # state_vec = states_input_list[[i]][2,]
  # phi1_lists = phi1_est_after
  # alpha_vecs = alphas_est_smoothed
  # phi0_vecs = phi0_est_smoothed; A_vecs= rep(1,48)
  # time_point = 48
  # ti_span=TI_loop
  
  library(tidyverse)
  # A_cov=actions_covmat(ti_span = ti_span,time_dim = taus)
  
  if(time_point>1){
    
    meanstate_iter1=0
    meanstate_iter2=0
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        meanstate_iter1=phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
      }else{
        
        meanstate_iter1=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
        
      }
      
    }
    
    meanstate=meanstate_iter1+meanstate_iter2
    
    phi1_product=as.matrix(Reduce("%*%", phi1_lists[1:(time_point-1)]))
    varstate_iter1=phi1_product%*%t(phi1_product)
    
    varstate_iter2=0
    
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        varstate_iter2=eps_cors+varstate_iter2
        
      }else{
        
        phi1_1_product=as.matrix(Reduce("%*%", phi1_lists[(k+1):(time_point-1)]))
        
        varstate_iter2=phi1_1_product%*%eps_cors%*%t(phi1_1_product)+
          varstate_iter2
        
      }
      
    }
    
    varstate_iter3=0
    
    for(i in 1:(time_point-1)){
      for(j in 1:(time_point-1)){
        # i=47;j=47
        
        if(((i+1)>(time_point-1))&((j+1)<=(time_point-1))){
          
          phi_01=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          phi_02=as.matrix(Reduce("%*%", phi1_lists[(j+1):(time_point-1)]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*(1/4)
        }else if(((j+1)>time_point-1)&((i+1)<=time_point-1)){
          phi_01=as.matrix(Reduce("%*%", phi1_lists[(i+1):(time_point-1)]))
          
          phi_02=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*(1/4)
        }else if(((j+1)>time_point-1)&((i+1)>time_point-1)){
          phi_01=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          phi_02=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*(1/4)
        }else{
          
          phi_01=as.matrix(Reduce("%*%", phi1_lists[(i+1):(time_point-1)]))
          
          phi_02=as.matrix(Reduce("%*%", phi1_lists[(j+1):(time_point-1)]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*(1/4)
          
        }
        
      }
    }
    
    varstate=varstate_iter1+varstate_iter2+varstate_iter3
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =meanstate, sigma =varstate  )
  }else{
    
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =rep(0,length(state_vec)), sigma =diag(length(state_vec)) )
  }
  
  
  return(margin_density)
  
}



estimate_state_density_var_behavior<-function(state_vec, phi1_lists, alpha_vecs, phi0_vecs, A_vecs, time_point,
                                              ti_span
){
  
  # state_vec = states_input_list[[i]][2,]
  # phi1_lists = phi1_est_after
  # alpha_vecs = alphas_est_smoothed
  # phi0_vecs = phi0_est_smoothed; A_vecs= rep(1,48)
  # time_point = 48
  # ti_span=TI_loop
  
  library(tidyverse)
  A_cov=actions_covmat(ti_span = ti_span,time_dim = taus)
  
  if(time_point>1){
    
    meanstate_iter1=0
    meanstate_iter2=0
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        meanstate_iter1=phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
      }else{
        
        meanstate_iter1=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%phi0_vecs[k,]+meanstate_iter1
        meanstate_iter2=Reduce("%*%", phi1_lists[(k+1):(time_point-1)])%*%alpha_vecs[k,]*A_vecs[k]+meanstate_iter2
        
      }
      
    }
    
    meanstate=meanstate_iter1+meanstate_iter2
    
    phi1_product=as.matrix(Reduce("%*%", phi1_lists[1:(time_point-1)]))
    varstate_iter1=phi1_product%*%t(phi1_product)
    
    varstate_iter2=0
    
    
    for(k in 1:(time_point-1)){
      
      if((k+1)>(time_point-1)){
        
        varstate_iter2=eps_cors+varstate_iter2
        
      }else{
        
        phi1_1_product=as.matrix(Reduce("%*%", phi1_lists[(k+1):(time_point-1)]))
        
        varstate_iter2=phi1_1_product%*%eps_cors%*%t(phi1_1_product)+
          varstate_iter2
        
      }
      
    }
    
    varstate_iter3=0
    
    for(i in 1:(time_point-1)){
      for(j in 1:(time_point-1)){
        # i=47;j=47
        
        if(((i+1)>(time_point-1))&((j+1)<=(time_point-1))){
          
          phi_01=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          phi_02=as.matrix(Reduce("%*%", phi1_lists[(j+1):(time_point-1)]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*A_cov[i,j]
        }else if(((j+1)>time_point-1)&((i+1)<=time_point-1)){
          phi_01=as.matrix(Reduce("%*%", phi1_lists[(i+1):(time_point-1)]))
          
          phi_02=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*A_cov[i,j]
        }else if(((j+1)>time_point-1)&((i+1)>time_point-1)){
          phi_01=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          phi_02=diag(x=1,nrow=nrow(phi1_lists[[1]]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*A_cov[i,j]
        }else{
          
          phi_01=as.matrix(Reduce("%*%", phi1_lists[(i+1):(time_point-1)]))
          
          phi_02=as.matrix(Reduce("%*%", phi1_lists[(j+1):(time_point-1)]))
          
          corss_mat=(phi_01%*%alpha_vecs[i,])%*%t(phi_02%*%alpha_vecs[j,])
          varstate_iter3=varstate_iter3+corss_mat*A_cov[i,j]
          
        }
        
      }
    }
    
    varstate=varstate_iter1+varstate_iter2+varstate_iter3
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =meanstate, sigma =varstate  )
  }else{
    
    
    library(mvtnorm)
    margin_density=dmvnorm(state_vec,mean =rep(0,length(state_vec)), sigma =diag(length(state_vec)) )
  }
  
  
  return(margin_density)
  
}



