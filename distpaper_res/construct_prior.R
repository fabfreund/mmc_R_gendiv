#' Draw prior parameter sets for simulating 

#' Draw theta as the Watterson estimate 
#' IMPORTANT: WE USE ms TO SIMULATE, SO GROWTH RATES g ARE ON 4N SCALE, WHICH CORRESPONDS TO 0.5g ON THE SN SCALE 
#' USED IN REBiepg
#' Models. 1: K+EXP, 2: BETA, 3: Dirac, 4: BetaXi 6: K
#' IMPORTANT: for Model 1, call this function only from the directory where ext_funs.R (and CREPiepg.o) are located
prior_obs_s <- function(n_ind,models=c(1,2,6),nsimul=c(100,100,0,0,0,100),
                        ranges = list(c(0.5,1,2.5,4,7,10,25,50,75,100,500,1000),seq(1,1.9,0.1),NULL,NULL,
                                      NULL,0),
                                      s_obs = c(5,10,15,20,40,60,75)){
theta <- vector("list",6)  
sample1 <- vector("list",6)

if (1 %in% models){
theta[[1]] <- sapply(ranges[[1]],function(g){if (g> 0){th1 <- 2*s_obs/sum(REBiepg(n_ind,g*0.5))} else
                                             {th1 <- s_obs/sum((1:(n_ind-1))^(-1))}
                                            return(th1)})
}
if (2 %in% models){
theta[[2]] <- sapply(ranges[[2]],function(a){2*s_obs/eln_f_lambda(n_ind,function(n){beta_bc_rates(n,2-a,a)})})}

if (3 %in% models){
theta[[3]] <- sapply(ranges[[3]],function(p){2*s_obs/eln_f_lambda(n_ind,function(n){dirac_bc_rates(n,p)})})}   

if (4 %in% models){
  theta[[4]] <- sapply(ranges[[4]],function(a){2*s_obs/eln_f_xi(n_ind,function(n){beta4xi_rates(n,2-a,a)})})}   

if (6 %in% models){
theta[[6]] <- sapply(ranges[[6]],function(a){s_obs/sum((1:(n_ind-1))^(-1))})} #KM RULEZ

for (m in models){sample1[[m]] <- sample(seq(along=ranges[[m]]),nsimul[m],TRUE)}
sample_s <- sample(seq(along=s_obs),sum(nsimul),TRUE)
sample_theta <- rep(0,sum(nsimul))
count1 <- 0
for (m in models){
for (i in 1:nsimul[m]){sample_theta[count1+i] <- theta[[m]][sample_s[count1+i],sample1[[m]][i]]}
count1 <- count1+nsimul[m]}  
mod_vec <- NULL 
coal_p <- NULL 
for (i in models){mod_vec <- c(mod_vec,rep(i,nsimul[i]))
                  coal_p <- c(coal_p,ranges[[i]][sample1[[i]]])}
prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                s_obs[sample_s],sample_theta)
colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt")
return(prior1)
}


#prior function for fixed theta, to partially replicate Kato's study
prior_theta <- function(n_ind,nsimul=c(2500,2500,2500),g_range = seq(10,100,10),
                        beta_range = seq(1,1.9,0.1),
                        theta1 = seq(10,100,10)){
  sample_alpha <- sample(seq(along=beta_range),nsimul[2],TRUE)
  sample_g <- sample(seq(along=g_range),nsimul[1],TRUE)
  sample_theta <- sample(theta1,sum(nsimul),TRUE)
  prior1 <- cbind(rep(n_ind,sum(nsimul)),c(rep(1,nsimul[1]),rep(2,nsimul[2]),rep(6,nsimul[3])),
                  c(g_range[sample_g],beta_range[sample_alpha],rep(0,nsimul[3])),
                  sample_theta)
  colnames(prior1) <- c("n_ind","model","coal_param","theta_watt")
  return(prior1)
}

#' as prior_obs_s, just with continuous ranges (given as vector 
#' (parameter min, parameter max))
#' Argument log_growth = FALSE:  uniform draw from the growth range
#' log_growth = TRUE: Draw uniformly on a log scaled range, i.e. 
#' draw uniformly from log(range) and then use exp(drawn values)
#' if include_g0>0, replace this many growth parameters drawn
#' by 0.
#' if integer mc1 > 1: Gives the number of cores for parallel computation of the 
#' Watterson estimators (in this case, you need to preload library(parallel))  
prior_obs_s_cont <- function(n_ind,models=c(1,2,6),nsimul=c(100,100,0,0,0,100),
                        ranges = list(c(log(0.5),log(1000)),c(1,1.99),c(0.01,0.99),
                                      NULL,NULL,c(0,0)),
                        s_obs = c(5,10,15,20,40,60,75),log_growth=TRUE,
                        include_g0=5,mc1=1){
  theta_f <- vector("list",6)  #' Function to compute theta
  #sample1 <- vector("list",6)
  if (1 %in% models){
    theta_f[[1]] <- function(g,s){if (g> 0){th1 <- 2*s/sum(REBiepg(n_ind,g*0.5))} else
    {th1 <- s/sum((1:(n_ind-1))^(-1))}
      return(th1)}
  }
  if (2 %in% models){
    theta_f[[2]] <- function(a,s){2*s/eln_f_lambda(n_ind,function(n){beta_bc_rates(n,2-a,a)})}}
  
  if (3 %in% models){
    theta_f[[3]] <- function(p,s){2*s/eln_f_lambda(n_ind,function(n){dirac_bc_rates(n,p)})}}   
  
  if (4 %in% models){
    theta_f[[4]] <- function(a,s){2*s/eln_f_xi(n_ind,function(n){beta4xi_rates(n,2-a,a)})}}   
  
  if (6 %in% models){
    theta_f[[6]] <- function(a,s){s/sum((1:(n_ind-1))^(-1))}} #KM RULEZ
  mod_vec <- NULL 
  coal_p <- NULL
  s_psobs <- NULL
  sample_theta <- NULL
  for (m in models){
  if (ranges[[m]][1]<ranges[[m]][2]){
  sample_param <- runif(min = ranges[[m]][1],max = ranges[[m]][2],
                        n = nsimul[m])} else {
                          sample_param <- rep(ranges[[m]][1],nsimul[[m]])}
  if (m==1 & include_g0>0 & log_growth){sample_param[sample(nsimul[m],include_g0)] <- -Inf}
  if (m==1 & include_g0>0 & !log_growth){sample_param[sample(nsimul[m],include_g0)] <- 0}
  if (m==1 & log_growth){sample_param <- exp(sample_param)}
  sample_s <- sample(s_obs,nsimul[m],TRUE)
  thef <- function(i){theta_f[[m]](sample_param[i],sample_s[i])}
  if (mc1==1){
  theta_drawn <- sapply(1:nsimul[m],thef)} else {
  theta_drawn <- unlist(mclapply(1:nsimul[m],thef,mc.cores = mc1))}
  mod_vec <- c(mod_vec,rep(m,nsimul[m]))
  coal_p <- c(coal_p,sample_param)
  s_psobs <- c(s_psobs,sample_s)
  sample_theta <- c(sample_theta,theta_drawn)}
  prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                  s_psobs,sample_theta)
  colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt")
  return(prior1)
}

#' As above, added option to include BSZ as model
prior_obs_s_cont2 <- function(n_ind,models=c(1,2,6),nsimul=c(100,100,0,0,0,100),
                             ranges = list(c(log(0.5),log(1000)),c(1,1.99),c(0.01,0.99),
                                           NULL,NULL,c(0,0)),
                             s_obs = c(5,10,15,20,40,60,75),log_growth=TRUE,
                             include_g0=5,include_bsz=5,mc1=1){
  theta_f <- vector("list",6)  #' Function to compute theta
  #sample1 <- vector("list",6)
  if (1 %in% models){
    theta_f[[1]] <- function(g,s){if (g> 0){th1 <- 2*s/sum(REBiepg(n_ind,g*0.5))} else
    {th1 <- s/sum((1:(n_ind-1))^(-1))}
      return(th1)}
  }
  if (2 %in% models){
    theta_f[[2]] <- function(a,s){2*s/eln_f_lambda(n_ind,function(n){beta_bc_rates(n,2-a,a)})}}
  
  if (3 %in% models){
    theta_f[[3]] <- function(p,s){2*s/eln_f_lambda(n_ind,function(n){dirac_bc_rates(n,p)})}}   
  
  if (4 %in% models){
    theta_f[[4]] <- function(a,s){2*s/eln_f_xi(n_ind,function(n){beta4xi_rates(n,2-a,a)})}}   
  
  if (6 %in% models){
    theta_f[[6]] <- function(a,s){s/sum((1:(n_ind-1))^(-1))}} #KM RULEZ
  mod_vec <- NULL 
  coal_p <- NULL
  s_psobs <- NULL
  sample_theta <- NULL
  for (m in models){
    if (ranges[[m]][1]<ranges[[m]][2]){
      sample_param <- runif(min = ranges[[m]][1],max = ranges[[m]][2],
                            n = nsimul[m])} else {
                              sample_param <- rep(ranges[[m]][1],nsimul[[m]])}
    if (m==1 & include_g0>0 & log_growth){sample_param[sample(nsimul[m],include_g0)] <- -Inf}
    if (m==1 & include_g0>0 & !log_growth){sample_param[sample(nsimul[m],include_g0)] <- 0}
    if (m==1 & log_growth){sample_param <- exp(sample_param)}
    if (m==2 & include_bsz>0){sample_param[sample(nsimul[m],include_bsz)] <- 1}
    sample_s <- sample(s_obs,nsimul[m],TRUE)
    thef <- function(i){theta_f[[m]](sample_param[i],sample_s[i])}
    if (mc1==1){
      theta_drawn <- sapply(1:nsimul[m],thef)} else {
        theta_drawn <- unlist(mclapply(1:nsimul[m],thef,mc.cores = mc1))}
    mod_vec <- c(mod_vec,rep(m,nsimul[m]))
    coal_p <- c(coal_p,sample_param)
    s_psobs <- c(s_psobs,sample_s)
    sample_theta <- c(sample_theta,theta_drawn)}
  prior1 <- cbind(rep(n_ind,sum(nsimul)),mod_vec,coal_p,
                  s_psobs,sample_theta)
  colnames(prior1) <- c("n_ind","model","coal_param","s_obs","theta_watt")
  return(prior1)
}