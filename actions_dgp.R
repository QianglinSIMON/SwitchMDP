actions_dgp<-function(taus,ndays,ti,act_init){
  
  ##--output our experimental design with switchback design---
  acts_mat=matrix(nrow = taus,ncol = ndays )
  
  #act_init =rbinom(1,1,0.5)
  
  if( ti <=taus &(ti>0)){
    
    if(act_init==1){
      
      acts1=rep(c(1,0),each=ti,length.out=taus)
      acts2=rep(c(0,1),each=ti,length.out=taus)
      
    }else{
      acts1=rep(c(0,1),each=ti,length.out=taus)
      acts2=rep(c(1,0),each=ti,length.out=taus)
    }
    
    
    
    for(i in 1:ndays ){
      if(i%%2!=0){
        
        acts_mat[,i]=acts1
      }else{
        acts_mat[,i]=acts2
      }
    }
    
    
  }else if(ti==0){
    
    ##--UR design---##
    
    acts1=rbinom(taus,1,0.5)
    acts2=1-acts1
    
    for(i in 1:ndays ){
      
      if(i%%2!=0){
        
        acts_mat[,i]=acts1
      }else{
        acts_mat[,i]=acts2
      }
      
    }
    
    
  }else{
    break
  }
  
  
  
  
  return(acts_mat)
  
  
}

##---Bojinov (MS, 2023) optimal design---
actions_dgp_Boinov<-function(taus,ndays,ti,seed){
  
  set.seed(seed)
  #seed =2023;
  
  # taus=48;ndays=16;ti=3
  
  Time_horizon=taus
  
  itera_num=Time_horizon/ti
  
  if(itera_num>=4){ 
    
    change_points=rep(0,itera_num-2)
    change_points[1]=1
    
    for(i in 2:length(change_points)){
      change_points[i]=i*ti+1
    }
    
    # change_points
    
    action_paths=NULL
    
    for(i in 1:length(change_points)){
      #i=254
      if(i<length(change_points)){
        act_num=change_points[i+1]-change_points[i]
        act_vec=rep(rbinom(1,size = 1,prob = 0.5),act_num)
        action_paths=c(action_paths,act_vec)
      }else{
        act_num=Time_horizon-change_points[i]+1
        act_vec=rep(rbinom(1,size = 1,prob = 0.5),act_num)
        action_paths=c(action_paths,act_vec)
      }
    }
    
  }else{
    
    action_paths=rep(rbinom(1,size = 1,prob = 0.5), Time_horizon)
    
  }
  
  #change_points
  
  #action_paths
  
  action_mat=matrix(nrow = taus,ncol = ndays)
  
  for(i in 1:ndays ){
    
    if(rbinom(1,1,prob = 0.5)==1){
      action_mat[,i]=action_paths
    }else{
      action_mat[,i]=1-action_paths
    }
    
    # if(i%%2!=0){
    #   
    #   action_mat[,i]=action_paths
    # }else{
    #   action_mat[,i]=1-action_paths
    # }
    
  }
  
  return(action_mat)
  
}



actions_dgp_Hu<-function(taus,ndays,ti){
  
  #taus=48; ndays=24; ti=2; 
  ##--output our experimental design with regular switchabck design---
  acts_mat=matrix(nrow = taus,ncol = ndays )
  
  
  
  if( ti <=taus &(ti>0)){
    
    act_vec0<-rep(0,taus)
    
    for(t in 1:taus){
      
      
      
      if( (t - 1) %% ti== 0 ){
        
        act_vec0[t]<-rbinom(1,1,0.5)
        
      }else{
        
        act_vec0[t]<-act_vec0[t-1]
        
      }
      
    }
    
    
    
    for(i in 1:ndays ){
      
      if(i%%2!=0){
        
        acts_mat[,i]=act_vec0
      }else{
        acts_mat[,i]=1-act_vec0
      }
      
    }
    
    
  }else{
    break
  }
  
  
  
  
  return(acts_mat)
  
  
}







actions_dgp_Xiong<-function(taus,ndays,ti){
  
  #taus=48; ndays=24; ti=2; 
  ##--output our experimental design with regular switchabck design---
  acts_mat=matrix(nrow = taus,ncol = ndays )
  
  
  if( ti <=taus &(ti>0)){
    
    act_vec0<-rep(0,taus)
    
    for(t in 1:taus){
      
      
      
      if( (t - 1) %% ti== 0 ){
        
        act_vec0[t]<-rbinom(1,1,0.5)
        
      }else{
        
        act_vec0[t]<-act_vec0[t-1]
        
      }
      
    }
    
    
    
    for(i in 1:ndays ){
      
      if(i%%2!=0){
        
        acts_mat[,i]=act_vec0
      }else{
        acts_mat[,i]=1-act_vec0
      }
      
    }
    
    
  }else{
    break
  }
  
  
  return(acts_mat)
  
  
}
