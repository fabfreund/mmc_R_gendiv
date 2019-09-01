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



