#steal some functions for the ordered beta package!
rordbeta <- function(n=100,
                     mu=0.5,
                     phi=1,
                     cutpoints=c(-1,1)) {
  
  if(!all(mu>0 & mu<1) ) {
    
    stop("Please pass a numeric value for mu that is between 0 and 1.")
    
  }
  
  if(!all(phi>0)) {
    
    stop("Please pass a numeric value for phi that is greater than 0.")
    
  }
  
  if(!(length(mu) %in% c(1,n))) {
    
    stop("Please pass a vector for mu that is either length 1 or length N.")
    
  }
  
  if(!(length(phi) %in% c(1,n))) {
    
    stop("Please pass a vector for phi that is either length 1 or length N.")
    
  }
  
  mu_ql <- qlogis(mu)
  
  if(length(mu_ql)==1) {
    mu_ql <- rep(mu_ql, n)
  }
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(mu_ql - cutpoints[1])
  middle <- plogis(mu_ql - cutpoints[1]) - plogis(mu_ql - cutpoints[2])
  high <- plogis(mu_ql - cutpoints[2])
  
  # we'll assume the same eta was used to generate outcomes
  
  out_beta <- rbeta(n = n,mu * phi, (1 - mu) * phi)
  
  # now determine which one we get for each observation
  outcomes <- sapply(1:n, function(i) {
    
    sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
    
  })
  
  # now combine binary (0/1) with proportion (beta)
  
  final_out <- sapply(1:n,function(i) {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })
  
  
  return(final_out)
  
}

# make predictions

ordered_beta = function(x,realx,beta,lapse,alpha, conf_int,conf_beta,conf_prec,conf_shift, conf_cutzero, conf_cutone,rt_int,rt_beta,rt_sd,rt_ndt,rt_shift, minRT,id = 1){
  
  
  probab = function(x,lapse,beta,alpha){
    return((inv_logit_scaled(lapse) / 2) + (1 - 2 *  (inv_logit_scaled(lapse) / 2)) * (0.5+0.5 * pracma::erf((x-alpha) / (exp(beta) * sqrt(2)))))
    
  }
  prob = probab(x,lapse,beta,alpha)
  
  conf_mu = inv_logit_scaled(conf_int + conf_beta * (probab(x-conf_shift,lapse,beta,alpha) * (1-probab(x-conf_shift,lapse,beta,alpha))))
  
  rt = exp(rt_int + rt_beta * (probab(x-rt_shift,lapse,beta,alpha) * (1-probab(x-rt_shift,lapse,beta,alpha)))) + inv_logit_scaled(rt_ndt) * minRT
  
  
  
  
  # predictions = rordbeta(n = 101, mu = conf_mu, phi = conf_prec, cutpoints = c(conf_cutzero, conf_cutone))
  # 
  # nested_predictions <- replicate(10, rordbeta(n = 101, mu = conf_mu, phi = conf_prec, cutpoints = c(conf_cutzero, conf_cutone)), simplify = FALSE)
  # 
  # df_predictions <- do.call(rbind, lapply(seq_along(nested_predictions), function(i) {
  #   data.frame(iteration = i, index = 1:101, predictions = nested_predictions[[i]])
  # }))
  # 
  # dd = df_predictions %>% group_by(index) %>% dplyr::summarize(mean = mean(predictions), sd = sd(predictions))
  
  return(data.frame(prob = prob,rt = rt, conf  = conf_mu,x = x, id = id))
}



ordered_beta_dens = function(realx,beta,lapse,alpha, conf_int,conf_beta,conf_prec,conf_shift, conf_cutzero, conf_cutone,rt_int,rt_beta,rt_sd,rt_ndt,rt_shift, minRT,id = 1){
  
  
  probab = function(x,lapse,beta,alpha){
    return((inv_logit_scaled(lapse) / 2) + (1 - 2 *  (inv_logit_scaled(lapse) / 2)) * (0.5+0.5 * pracma::erf((x-alpha) / (exp(beta) * sqrt(2)))))
    
  }
  prob = probab(realx,lapse,beta,alpha)
  
  conf_mu = inv_logit_scaled(conf_int + conf_beta * (probab(realx-conf_shift,lapse,beta,alpha) * (1-probab(realx-conf_shift,lapse,beta,alpha))))
  
  rt = exp(rt_int + rt_beta * (probab(realx-rt_shift,lapse,beta,alpha) * (1-probab(realx-rt_shift,lapse,beta,alpha)))) + inv_logit_scaled(rt_ndt) * minRT
  
  
  conf <- tryCatch(
    {
      rordbeta(n = 1, mu = conf_mu, phi = conf_prec, cutpoints = c(conf_cutzero, conf_cutone))
    },
    error = function(e) {
      NA
    }
  )
  return(data.frame(prob = prob,rt = rt, conf  = conf,x = realx, id = id))
}




inv_logit_scaled = function(x){
  return(exp(x)/(1+exp(x)))
}


prep_data = function(raw_hrd){
  
  second_ses = raw_hrd %>% filter(session == 2)
  first_ses = raw_hrd %>% filter(!participant_id %in% unique(second_ses$participant_id))
  
  raw_hrd = rbind(second_ses,first_ses)
  
  hrd<-raw_hrd%>%
    #condensate the table so that there is one line per intensity
    group_by(nTrials,Alpha,participant_id,Modality)%>%
    summarise(x=mean(Alpha),
              rt = DecisionRT,
              Confidence = Confidence,
              x_ratio=mean(responseBPM/listenBPM),
              n=sum(nTrials<1000),
              y=sum(Decision=="More"),
              s=0,
              ID=mean(as.numeric(str_replace(participant_id,'sub-',''))))%>%
    ungroup()%>%
    select(x,x_ratio,n,y,rt,ID,s,Modality,participant_id,nTrials, Confidence) %>%
    #exclude participant for which there was an error during data acquisition (Ask Nice or Leah about this IDK)
    filter((ID<263|ID>295))
  
  #change id so that it start at 1 and goes up to N_participant in increment of 1  
  old_ID<-unique(hrd$ID)
  
  for(idx in 1:length(old_ID)){
    hrd$s[hrd$ID==old_ID[idx]]<-idx
  }
  hrd<-select(hrd,x,x_ratio,n,rt,y,s,Modality,nTrials, Confidence)
  
  #save dataframe
  #write_csv(hrd,'preped_hrd_data2.csv')
  return(hrd)
  
  
}

