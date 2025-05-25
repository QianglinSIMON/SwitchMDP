
##---The Legendre basis---
Legerand_basis_funcs<-function(x,degree){
  #degree>=1
  if(degree==1){
    
    y=c(x) #--linear basis--
    
  }else if(degree==2){
    
    y=c(x,0.5*(3*x^2-1) )
  }else if(degree==3){
    
    y=c(x, 0.5*(3*x^2-1), 0.5*(5*x^3-3*x))
  }else if(degree==4){
    y=c(x, 0.5*(3*x^2-1), 0.5*(5*x^3-3*x), (1/8)*(35*x^4-30*x^2+3))
  }else if(degree==5){
    
    y=c(x, 0.5*(3*x^2-1), 0.5*(5*x^3-3*x), (1/8)*(35*x^4-30*x^2+3),
        (1/8)*(63*x^5-70*x^3+15*x))
  }else if(degree==6){
    y=c(x, 0.5*(3*x^2-1), 0.5*(5*x^3-3*x), (1/8)*(35*x^4-30*x^2+3),
        (1/8)*(63*x^5-70*x^3+15*x),(1/16)*(231*x^6 -315*x^4+105*x^2-5))
    
  }else if(degree==7){
    
    y=c(x, 0.5*(3*x^2-1), 0.5*(5*x^3-3*x), (1/8)*(35*x^4-30*x^2+3),
        (1/8)*(63*x^5-70*x^3+15*x),(1/16)*(231*x^6 -315*x^4+105*x^2-5),
        (1/16)*(429*x^7-693*x^5+315*x^3-35*x))
  }else{
    
    print("Repeat input degree!")
  }
  
  return(y)
}

##---triangle basis---
Tri_orgonal_basis<-function(x,degree){
  
  w=2*pi
  
  if(degree==1){
    y=c(cos(w*x))
  }else if(degree==2){
    
    y=c(cos(w*x),sin(w*x))
  }else if(degree==3){
    
    y=c(cos(w*x), sin(w*x), cos(2*w*x))
  }else if(degree==4){
    y=c(cos(w*x), sin(w*x), cos(2*w*x), sin(2*w*x))
  }else if(degree>=5){
    
    y=c(cos(w*x), sin(w*x), cos(2*w*x), sin(2*w*x), cos(3*w*x))
  }else{
    
    print("Repeat input degree!")
  }
  
  return(y)
  
}


state_orthogonal_basis<-function(state_vec,time_point,Ln,Kn){
  
  
  L=length(state_vec)
  state_basis_vecs=NULL
  for(l in 1:L){
    add_legd=Legerand_basis_funcs(x=state_vec[l],degree =Ln )
    state_basis_vecs=c( state_basis_vecs, add_legd)
    
  }
  
  


  
  if(Kn==0){
    
    state_basis_vecs=c(1,state_basis_vecs)
    
  }else{
    time_basis=Legerand_basis_funcs(time_point,degree = Kn)
    state_basis_vecs=c(1,state_basis_vecs, time_basis )
  }
  
  
  return(state_basis_vecs)
}


state_orthogonal_basis_time<-function(states_arrays,taus,Ln,Kn){
  
  
  #states_arrays=states_array;Kn=2;Ln=4
  
  
  final_basis_list=vector("list",length = taus)
  
  for(i in 1:taus ){
    
    #i=1
    s_nvecs=NULL
    for(j in 1:dim(states_arrays)[2]){
      
      s_nvec=states_arrays[i,j,]
      
      s_nvec=(s_nvec-min(s_nvec))/(max(s_nvec)-min(s_nvec))#---[0,1] scaled 
      
      s_nvecs=(as.matrix(rbind(s_nvecs,s_nvec)))
    }
    
    s_nvecs=t(s_nvecs)
    
    final_basis=NULL
    
    for(k in 1:nrow(s_nvecs)){
      
      #k=1
      
      final_basis=rbind(final_basis,as.vector(
        
        state_orthogonal_basis(state_vec = s_nvecs[k,],time_point = i/taus, Ln=Ln,Kn=Kn)))
      
    }
    
    
    
    #final_basis=t(final_basis)
    
    final_basis_list[[i]]=final_basis
  }
  
  return(final_basis_list)
  
}


##--interaction terms--
state_corr_basis<-function(state_vec,time_point,Ln,Kn){
  
  state_basis_vecs=c( state_vec[1]^seq(1,Ln), state_vec[2]^seq(1,Ln),
                      state_vec[1]*state_vec[2]
  )
  
  
  
  
  if(Kn==0){
    
    state_basis_vecs=c(1,state_basis_vecs)
    
  }else{
    
    time_basis=Legerand_basis_funcs(time_point,degree = Kn)
    state_basis_vecs=c(1,state_basis_vecs, time_basis )
  }
  
  
  return(state_basis_vecs)
  
}


state_corr_basis_time<-function(states_arrays,taus,Ln,Kn){
  
  
  #states_arrays=states_array;Kn=2;Ln=4
  
  
  final_basis_list=vector("list",length = taus)
  
  for(i in 1:taus ){
    
    #i=1
    s_nvecs=NULL
    for(j in 1:dim(states_arrays)[2]){
      
      s_nvec=states_arrays[i,j,]
      
      s_nvec=(s_nvec-min(s_nvec))/(max(s_nvec)-min(s_nvec))#---[0,1] scaled 
      
      #  s_nvec=(s_nvec-mean(s_nvec))/(sd(s_nvec)) #---mean-sd scaled 
      
      s_nvecs=(as.matrix(rbind(s_nvecs,s_nvec)))
    }
    
    s_nvecs=t(s_nvecs)
    
    #snvecs_aug=cbind(s_nvecs,s_nvecs[,1]*s_nvecs[,2])
    
    final_basis=NULL
    
    for(k in 1:nrow(s_nvecs)){
      
      #k=1
      
      final_basis=rbind(final_basis,as.vector(
        
        #state_orthogonal_basis(state_vec = s_nvecs[k,],time_point = i/taus, Ln=Ln,Kn=Kn)
        state_corr_basis(state_vec = s_nvecs[k,],time_point = i/taus, Ln=Ln,Kn=Kn)))
      
    }
    
    
    
    #final_basis=t(final_basis)
    
    final_basis_list[[i]]=final_basis
  }
  
  return(final_basis_list)
  
}






