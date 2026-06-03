#Code to run Age-length Models for blue whales
#Zoe Rand
#Last updated: June 3 2026
library(tidyverse)
library(here)
library(RTMB)
library(RTMBdist)
library(tmbstan)
library(bayesplot)
library(posterior)
library(patchwork)
library(loo)
library(parallel)


# Data --------------------------------------------------------------------
#japanese and NP blue whale data
BW_earplug_dat<-read_csv(here("data", "BW_Earplugs_Maeda_and_Mizroch.csv"))

#soviet blue whale data (only pygmy)
PyDat_S<-read_csv(here("data", "Earplug_pygmy_S.csv"))

#data wrangling
BW_earplug <- BW_earplug_dat %>% mutate(Sex = ifelse(Sex == "F", "Female", "Male"), Length = LengthFT * 0.3048)

S_dat<-PyDat_S %>% add_column(BluePopulation = "Pygmy") %>% rename(OriginalAge = AgeOrig) %>% mutate(OriginalAge = as.character(OriginalAge))


Earplug_dat<-BW_earplug %>% bind_rows(S_dat)%>% filter(is.na(AccurateAge) | AccurateAge == "Yes")

#for censored data 
min_length<-21.336
#number of data points
ndat<-Earplug_dat %>% group_by(BluePopulation, Sex) %>% summarise(n = n())

#facet labels
labs<-unique(ndat$BluePopulation)
labs_2<-paste0(labs, "\n # F = ", ndat$n[ndat$Sex == "Female"], "\n # M = ",ndat$n[ndat$Sex == "Male"] )
names(labs_2)<-labs

#plotting data
ggplot(Earplug_dat) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  facet_wrap(~BluePopulation, labeller = as_labeller(labs_2)) + 
  theme_classic()

write_csv(Earplug_dat, here("data", "Earplug_dat_all.csv"))
# Models ------------------------------------------------------------------
#generalized growth model from Schnute
gen_growth_mod_single<-function(parms, dat){
  #b = 1 is Von Bertalanffy
  #b = 0.0001 is Gompertz
  require(RTMB)
  getAll(dat, parms, warn = FALSE)
  #parameters
  L_2_F<-exp(log_L2) #vector of 3
  L_1_F<-exp(log_L1) #vector of 3
  L_1_M<-exp(log_L1_M)
  L_2_M<-exp(log_L2_M)
  
  k<-exp(log_k) #matrix (sex on row pop on column)
  #b is a matrix (sex on row, pop on column)
  
  sigma<-exp(log_sigma) #vector of 3
  
  
  
  L_2<-matrix(c(L_2_F, L_2_M), ncol = 3, byrow = TRUE)
  #print(L_2)
  
  
  L_1<-matrix(c(L_1_F, L_1_M), ncol = 3, byrow = TRUE)
  #print(L_1)
  temp_L_inf<-matrix(0, nrow = 2, ncol = 3)
  L_inf<-matrix(0, nrow = 2, ncol = 3)
  for(i in 1:3){
    for(s in 1:2){
      Lmin<-L_1[s,i]^b[s,i]
      Lmax<-L_2[s,i]^b[s,i]
      temp_L_inf[s,i]<-Lmin + ((Lmax - Lmin)/(1-exp(-k[s,i]*(t2-t1)))) #t2 and t1 are data
      L_inf[s,i]<-temp_L_inf[s,i]^(1/b[s,i])
    }
  }
  #print(L_inf)
  #priors
  
  #jacobian for uniform priors on log scale
  pri<-sum(log_L1, log_L1_M, log_k, log_L2, log_L2_M)
  
  
  #predictions
  
  mu<-rep(0, nrow(dat))
  for(i in 1:nrow(dat)){
    Lmin<-L_1[dat$sex[i], dat$pop[i]]^b[dat$sex[i], dat$pop[i]]
    Lmax<-L_2[dat$sex[i], dat$pop[i]]^b[dat$sex[i], dat$pop[i]]
    Linf<-temp_L_inf[dat$sex[i], dat$pop[i]]
    temp<-Linf + (Lmin - Linf) *exp(-k[dat$sex[i], dat$pop[i]]*(dat$Age[i]-t1))
    mu[i]<-temp^(1/b[dat$sex[i], dat$pop[i]])
  }
  #print(mu)
  #pointwise likelihood
  log_ll<-rep(0, nrow(dat))
  for(i in 1:nrow(dat)){
    if(dat$pop[i] != 1){ #non pygmy
      log_ll[i]<-dnorm(dat$Length[i],
                       mu[i], sigma[dat$pop[i]], log = TRUE)
    }else if(dat$Source[i] == "Sazhinov"){ #soviet
      log_ll[i]<-dnorm(dat$Length[i],
                       mu[i], sigma[dat$pop[i]], log = TRUE)
    } else if(dat$Censored[i] == 1){ #censored
      log_ll[i]<-log(pnorm(min_length, mu[i], sigma[1]))
    }else { #other japanese data
      log_ll[i]<-log(dnorm(dat$Length[i], 
                           mu[i], sigma[1]))
    }
    
  }
  nll<-sum(-1*log_ll)
  
  #posterior
  post<-pri
  
  #non-pymgy
  post<- post - sum(dnorm(dat$Length[dat$pop!= 1],
                          mu[dat$pop!= 1], sigma[dat$pop[dat$pop!= 1]], log = TRUE))
  #print(post)
  #soviet data
  
  post<- post - sum(dnorm(dat$Length[dat$pop == 1 & dat$Source == "Sazhinov"], 
                          mu[dat$pop == 1 & dat$Source == "Sazhinov"], sigma[1], log = TRUE))
  #print(post)
  
  #japanese data (censored likelihood)
  #for data at min length
  cum_dist<-pnorm(min_length, mu[dat$pop == 1 & dat$Source == "Maeda" & dat$Censored == 1], sigma[1])
  
  #for data above min length
  dist_norm<-dnorm(dat$Length[dat$pop == 1 & dat$Source == "Maeda" & dat$Censored == 0], 
                   mu[dat$pop == 1 & dat$Source == "Maeda" & dat$Censored == 0], sigma[1])
  
  #update posterior
  post<-post - sum(log(dist_norm))  - sum(log(cum_dist))
  
  
  
  REPORT(L_inf)
  REPORT(k)
  REPORT(sigma)
  
  REPORT(mu)
  REPORT(log_ll)
  REPORT(nll)
  REPORT(pri)
  
  return(post)
}

#function to run model with differen subsets of data
cmb <- function(f, d) function(p) f(p, d)


# Simulation function -----------------------------------------------------

sim_length<-function(pars, dat){
  L2_F<-exp(pars[1:3])
  L1_F<-exp(pars[4:6])
  L2_M<-exp(pars[7:9])
  L1_M<-exp(pars[10:12])
  k<-exp(pars[13])
  sigma<-exp(pars[14:16])
  b<-pars[17] 
  
  #combine
  L_2<-matrix(c(L2_F, L2_M), ncol = 3, byrow = TRUE)
  #print(L_2)
  
  L_1<-matrix(c(L1_F, L1_M), ncol = 3, byrow = TRUE)
  #print(L_1)
  
  #get Linf to check that it's reasonable
  temp_L_inf<-matrix(0, nrow = 2, ncol = 3)
  L_inf<-matrix(0, nrow = 2, ncol = 3)
  for(i in 1:3){
    for(s in 1:2){
      Lmin<-L_1[s,i]^b
      Lmax<-L_2[s,i]^b
      temp_L_inf[s,i]<-Lmin + ((Lmax - Lmin)/(1-exp(-k*(t2-t1)))) #t2 and t1 are data
      L_inf[s,i]<-temp_L_inf[s,i]^(1/b)
    }
  }
  #print(L_inf)
  #predictions
  mu<-rep(0, nrow(dat))
  for(i in 1:nrow(dat)){
    Lmin<-L_1[dat$sex[i], dat$pop[i]]^b
    Lmax<-L_2[dat$sex[i], dat$pop[i]]^b
    Linf<-temp_L_inf[dat$sex[i], dat$pop[i]]
    temp<-Linf + (Lmin - Linf) *exp(-k*(dat$Age[i]-t1))
    mu[i]<-temp^(1/b)
    #print(mu)
  }
  #simulate lengths
  preds<-rep(NA, nrow(dat))
  new_preds<-rep(NA, nrow(dat))
  for(i in 1:nrow(dat)){
    preds[i]<-rnorm(1, mu[i], sigma[dat$pop[i]])
    new_preds[i]<-preds[i]
    if(dat$Censored[i] == 1){
      new_preds[i]<-ifelse(round(preds[i], 3) <= min_length, min_length, preds[i])
    } 
  }

  return(new_preds)
}
# Data Wrangling for Model ----------------------------------------------------------

mod_dat<-Earplug_dat %>% filter(BluePopulation != "W N Pacific") %>%
  mutate(pop = ifelse(BluePopulation == "Pygmy", 1, 
                      ifelse(BluePopulation == "E N Pacific", 2, 3))) %>%
  mutate(sex = ifelse(Sex == "Female", 1, 2)) %>%
  mutate(Censored = ifelse(Source == "Maeda" & BluePopulation == "Pygmy" & LengthFT == 70,1, 0)) %>%
  mutate(BluePopulation2 = ifelse(BluePopulation == "Pygmy", 
                                  paste0(BluePopulation, "_", Source), BluePopulation))


# Get model objects -------------------------------------------------------

get_report<-function(draw, obj){
  rep<-obj$report(draw)
  return(rep)
  #return(rep$mu)
}


# Base Case (Richards, only Soviet data for pygmy blue whales) --------------------------------------------------
#using simplest formulation for k and b
#see model comparison below with all data (or refer to appendix of paper)
map16<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(1, 6), nrow = 2)))

t1<-1 #minimum age overall (want to estimte the early growth curve well)
t2<-46 #maximum age in ENP data

p1<-list(
  log_L2 = c(log(22), log(22), log(25)), 
  log_L1 = c(log(15), log(12), log(18)),
  log_L2_M = c(log(21), log(21), log(24)), 
  log_L1_M = c(log(13), log(12), log(16)),
  log_k = log(matrix(rep(0.17, 6), nrow = 2)),
  log_sigma = c(log(1.1), log(0.2), log(0.3)), 
  b = matrix(rep(2.5, 6), nrow = 2))

#getting just Soviet data for pygmy blue whales
soviet_mod_dat<-mod_dat %>% filter(BluePopulation2 != "Pygmy_Maeda") 

ggplot(soviet_mod_dat) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  facet_wrap(~BluePopulation2) + 
  theme_classic()

Sov_obj<-MakeADFun(cmb(gen_growth_mod_single, soviet_mod_dat), p1, map = map16)

init_fun<-function(){
  L2<- rnorm(3, 23, 2)
  L2_M<-rnorm(3, 22, 2)
  return(list(
    log_L2 = log(L2),
    log_L1 = log(0.5*L2), #so that L1 is below L2
    log_L2_M = log(L2_M),
    log_L1_M = log(0.5*L2_M),
    log_k = rep(log(rnorm(1, 0.17, 0.01)), 1),
    log_sigma = c(log(1.1), log(1.2), log(1.3)), 
    b = rep(runif(1, 1, 5), 1)
  )
  )
}

fit_Sov<-tmbstan(Sov_obj, iter = 6000, chains = 4,
                 warmup = 5000,
                 #thin = 10,
                 control = list(adapt_delta = 0.97, max_treedepth = 17), 
                 init = init_fun, 
                 seed = 123, 
                 cores = 1,
                 lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), rep(log(0.05), 1), rep(-100, 3), rep(0, 1)), 
                 upper = c(rep(log(10000), 3), log(25), log(25), log(25), rep(log(10000), 3),log(25), log(25), log(25), rep(log(0.3), 1), rep(log(10000), 3), rep(15, 1)))
fit_Sov

mcmc_trace(fit_Sov)
mcmc_dens(fit_Sov)
mcmc_acf(fit_Sov)
pairs(fit_Sov)

#draws
Sov_draws<-as_draws_matrix(fit_Sov)
Sov_draws<-Sov_draws[, -ncol(Sov_draws)]
saveRDS(Sov_draws, here("results", "Richards_m16_Soviet_draws.RDS"))

#reports
Sov_reps<-apply(Sov_draws, 1, get_report, obj = Sov_obj)

saveRDS(Sov_reps, here("results", "Soveit_m16_reports.RDS"))

#get loglikelihood
loglike<-t(sapply(Sov_reps, function(x){x$log_ll}))
elpd(loglike)
R_loo_Sov<-loo(loglike)

#get predictions
#new data
new_dat<-expand_grid(pop = 1:3, Age = 1:85, sex = 1:2)
new_dat$Length<-rep(0, nrow(new_dat))
new_dat$Censored<-rep(0, nrow(new_dat))
new_dat$Source<-rep("Soviet", nrow(new_dat)) #so no censoring

#make new object
S_R_new<-MakeADFun(cmb(gen_growth_mod_single, new_dat), p1, map = map16)

#get predictions
reps_preds_S<-apply(Sov_draws, 1, get_report, obj = S_R_new)
preds_S<-sapply(reps_preds_S, function(x){x$mu})

saveRDS(preds_S, here("results", "Richards_preds_Soviet.RDS"))



# Von Bertalanffy (Soviet data for pygmy blue whales) -------------------------------------------
t1<-1

t2<-46
map<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(NA, 6), nrow = 2))) #b = 1for VB

p<-list(
  log_L2 = c(log(22), log(22), log(25)), 
  log_L1 = c(log(15), log(12), log(18)),
  log_L2_M = c(log(21), log(21), log(24)), 
  log_L1_M = c(log(13), log(12), log(16)),
  log_k = log(matrix(rep(0.17, 6), nrow = 2)),
  log_sigma = c(log(1.1), log(0.2), log(0.3)), 
  b = matrix(rep(1, 6), nrow = 2)) #VB, b = 1

VB_obj_S<-MakeADFun(cmb(gen_growth_mod_single, soviet_mod_dat), p, map = map)
VB_obj_S$fn()



orig_preds<-VB_obj_S$report()$mu
mod_dat_test<-soviet_mod_dat %>% add_column(preds = orig_preds)
ggplot(mod_dat_test) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  geom_point(aes(x = Age, y = preds), color = "black") +
  facet_wrap(~BluePopulation, labeller = as_labeller(labs_2), ncol = 2) + 
  theme_classic()

init_funVB<-function(){
  L2<- rnorm(3, 23, 2)
  L2_M<-rnorm(3, 22, 2)
  return(list(
    log_L2 = log(L2),
    log_L1 = log(0.5*L2), #so that L1 is below L2
    log_L2_M = log(L2_M),
    log_L1_M = log(0.5*L2_M),
    log_k = log(rnorm(1, 0.17, 0.01)),
    log_sigma = c(log(1.1), log(1.2), log(1.3)), 
    b = 1
  )
  )
}



VB_fit_S<-tmbstan(VB_obj_S, iter = 2000, chains = 4,
                  warmup = 1000,
                  #thin = 10,
                  control = list(adapt_delta = 0.97, max_treedepth = 17), 
                  init = init_funVB,
                  lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), log(0.05), rep(-100, 3)),
                  upper = c(rep(log(10000), 3), log(25), log(25), log(25), rep(log(10000), 3),log(25), log(25), log(25), log(0.2), rep(log(10000), 3)))
saveRDS(VB_fit_S, here("results", "VB_fit_Sov.RDS"))

mcmc_trace(VB_fit_S)
mcmc_dens(VB_fit_S)
mcmc_acf(VB_fit_S)
pairs(VB_fit_S)

#model draws
VB_draws_S<-as_draws_matrix(VB_fit_S)

#reports
VB_reps_S<-apply(VB_draws_S, 1, get_report, obj = VB_obj_S)

saveRDS(VB_reps_S, here("results", "VB_reports_Soviet.RDS"))

#get loglikelihood
loglike<-t(sapply(VB_reps_S, function(x){x$log_ll}))
elpd(loglike)
VB_loo<-loo(loglike)


# Gompertz (Soviet data) --------------------------------------------------

t1<-1
t2<-46
map<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(NA, 6), nrow = 2))) #b = 0.0001 for Gompertz

p<-list(
  log_L2 = c(log(22), log(22), log(25)), 
  log_L1 = c(log(15), log(12), log(18)),
  log_L2_M = c(log(21), log(21), log(24)), 
  log_L1_M = c(log(13), log(12), log(16)),
  log_k = log(matrix(rep(0.17, 6), nrow = 2)),
  log_sigma = c(log(1.1), log(0.2), log(0.3)), 
  b = matrix(rep(0.0001, 6), nrow = 2)) #Gompertz, b = 0.0001


G_obj_S<-MakeADFun(cmb(gen_growth_mod_single, soviet_mod_dat), p, map = map)
G_obj_S$fn()

orig_preds<-G_obj_S$report()$mu
mod_dat_test<-soviet_mod_dat %>% add_column(preds = orig_preds)
ggplot(mod_dat_test) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  geom_point(aes(x = Age, y = preds), color = "black") +
  facet_wrap(~BluePopulation, labeller = as_labeller(labs_2), ncol = 2) + 
  theme_classic()

init_funG<-function(){
  L2<- rnorm(3, 23, 2)
  L2_M<-rnorm(3, 22, 2)
  return(list(
    log_L2 = log(L2),
    log_L1 = log(0.5*L2), #so that L1 is below L2
    log_L2_M = log(L2_M),
    log_L1_M = log(0.5*L2_M),
    log_k = log(rnorm(1, 0.17, 0.01)),
    log_sigma = c(log(1.1), log(1.2), log(1.3)), 
    b = 0.0001
  )
  )
}



G_fit_S<-tmbstan(G_obj_S, iter = 3000, chains = 4,
                 warmup = 2000,
                 #thin = 10,
                 control = list(adapt_delta = 0.97, max_treedepth = 17), 
                 init = init_funG, 
                 lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), log(0.05), rep(-100, 3)),
                 upper = c(rep(log(10000), 3), log(25), log(25), log(25), rep(log(10000), 3),log(25), log(25), log(25), log(0.2), rep(log(10000), 3)))

saveRDS(G_fit_S, here("results", "G_fit_Soviet.RDS"))

mcmc_trace(G_fit_S)
mcmc_dens(G_fit_S)
mcmc_acf(G_fit_S)
pairs(G_fit_S)

#model draws
G_draws_S<-as_draws_matrix(G_fit_S)

#reports
G_reps_S<-apply(G_draws_S, 1, get_report, obj = G_obj_S)

saveRDS(G_reps_S, here("results", "G_reports_Soviet.RDS"))

loglike_G<-t(sapply(G_reps_S, function(x){x$log_ll}))
elpd(loglike_G)
G_loo<-loo(loglike_G)


# Base Case Model Comparison ----------------------------------------------
loo_compare(R_loo_Sov, VB_loo, G_loo)

#Full results table
VB_draws_S<-VB_fit_S %>% as_draws_matrix()
VB_out_S<-apply(VB_draws_S, 2, function(x){quantile(exp(x), c(0.025, 0.5, 0.975))}) %>% 
  as_tibble() %>% select(-`lp__`)
colnames(VB_out_S)<-c("L2_P", "L2_ENP", "L2_A", "L1_P", "L1_ENP", "L1_A", 
                      "L2M_P", "L2M_ENP", "L2M_A", "L1M_P", "L1M_ENP", "L1M_A",
                      "k", "sigma_P", "sigma_ENP", "sigma_A")
VB_out_S_2 <-VB_out_S%>% add_column(b = 1) %>% t() %>% as_tibble(rownames = "Parameter")
colnames(VB_out_S_2)[2:4]<-c("lwr", "med", "upr")


G_draws_S<-G_fit_S %>% as_draws_matrix()
G_out_S<-apply(G_draws_S, 2, function(x){quantile(exp(x), c(0.025, 0.5, 0.975))}) %>% 
  as_tibble() %>% select(-`lp__`)
colnames(G_out_S)<-c("L2_P", "L2_ENP", "L2_A", "L1_P", "L1_ENP", "L1_A", 
                     "L2M_P", "L2M_ENP", "L2M_A", "L1M_P", "L1M_ENP", "L1M_A",
                     "k", "sigma_P", "sigma_ENP", "sigma_A")
G_out_S_2 <-G_out_S%>% add_column(b= 0.0001) %>% t() %>% as_tibble(rownames = "Parameter")
colnames(G_out_S_2)[2:4]<-c("lwr", "med", "upr")

R_draws_S<-fit_Sov %>% as_draws_matrix()
R_out_S<-apply(R_draws_S[,1:16], 2, function(x){quantile(exp(x), c(0.025, 0.5, 0.975))}) %>% 
  as_tibble() 
b_out<-apply(R_draws_S[,17], 2, quantile, c(0.025, 0.5, 0.975)) %>% as_tibble()
R_out_S<-bind_cols(R_out_S, b_out)
colnames(R_out_S)<-c("L2_P", "L2_ENP", "L2_A", "L1_P", "L1_ENP", "L1_A", 
                     "L2M_P", "L2M_ENP", "L2M_A", "L1M_P", "L1M_ENP", "L1M_A",
                     "k", "sigma_P", "sigma_ENP", "sigma_A", "b")
R_out_S_2 <-R_out_S %>% t() %>% as_tibble(rownames = "Parameter")
colnames(R_out_S_2)[2:4]<-c("lwr", "med", "upr")

#combine
out_all<-R_out_S_2 %>% left_join(VB_out_S_2, by = "Parameter") %>% 
  left_join(G_out_S_2, by = "Parameter") %>% 
  select(Parameter, med.x, lwr.x, upr.x, med.y, lwr.y, upr.y, med, lwr, upr)

#put rows in better order
out_all$Parameter<-factor(out_all$Parameter, 
                          levels = c("L2_P", "L2M_P", "L2_ENP", "L2M_ENP", "L2_A","L2M_A",
                                     "L1_P", "L1M_P", "L1_ENP", "L1M_ENP", "L1_A", "L1M_A",
                                     "k",
                                     "sigma_P", "sigma_ENP", "sigma_A", 
                                     "b"))
out_tab<-out_all %>% arrange(Parameter)
write_csv(out_tab, here("results","Model_result_table_Soviet_basecase.csv"))



# Japanese data for pygmy blue whales (Richards) ------------------------------------------------
map16<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(1, 6), nrow = 2)))
p1<-list(
  log_L2 = c(log(22), log(22), log(25)), 
  log_L1 = c(log(15), log(12), log(18)),
  log_L2_M = c(log(21), log(21), log(24)), 
  log_L1_M = c(log(13), log(12), log(16)),
  log_k = log(matrix(rep(0.17, 6), nrow = 2)),
  log_sigma = c(log(1.1), log(0.2), log(0.3)), 
  b = matrix(rep(2.5, 6), nrow = 2))


japanese_mod_dat<-mod_dat %>% filter(BluePopulation2 != "Pygmy_Sazhinov")

ggplot(japanese_mod_dat) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  facet_wrap(~BluePopulation2) + 
  theme_classic()

J_obj<-MakeADFun(cmb(gen_growth_mod_single, japanese_mod_dat), p1, map = map16)

init_fun<-function(){
  L2<- rnorm(3, 23, 2)
  L2_M<-rnorm(3, 22, 2)
  return(list(
    log_L2 = log(L2),
    log_L1 = log(0.5*L2), #so that L1 is below L2
    log_L2_M = log(L2_M),
    log_L1_M = log(0.5*L2_M),
    log_k = rep(log(rnorm(1, 0.17, 0.01)), 1),
    log_sigma = c(log(1.1), log(1.2), log(1.3)), 
    b = rep(runif(1, 1, 5), 1)
  )
  )
}

fit_J<-tmbstan(J_obj, iter = 9000, chains = 4,
               warmup = 8000,
               #thin = 10,
               control = list(adapt_delta = 0.97, max_treedepth = 17), 
               init = init_fun, 
               seed = 444, 
               cores = 1,
               lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), rep(log(0.05), 1), rep(-100, 3), rep(0, 1)), 
               upper = c(rep(log(10000), 3), log(25), log(25), log(25), rep(log(10000), 3),log(25), log(25), log(25), rep(log(0.3), 1), rep(log(10000), 3), rep(15, 1)))
fit_J

mcmc_trace(fit_J)
mcmc_dens(fit_J)
mcmc_acf(fit_J)



# All data (Richards) -----------------------------------------------------
#b is an estimated parameter
#16 models tested
#in all L1 is by population and sex
#Models 1-4
#1) K by pop*sex, b by sex*pop
#2) K by pop*sex, b by population
#3) K by pop*sex, b by sex
#4) K by pop*sex, one b
#Models 5-8 same but with k by population
#Models 9-12 same but with k by sex
#Models 13-16 same but with only one k

map1<-list(log_k = factor(matrix(1:6, nrow =2)), b = factor(matrix(1:6, nrow = 2)))
map2<-list(log_k = factor(matrix(1:6, nrow =2)), b = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)))
map3<-list(log_k = factor(matrix(1:6, nrow =2)), b = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)))
map4<-list(log_k = factor(matrix(1:6, nrow =2)), b = factor(matrix(rep(1, 6), nrow = 2)))
map5<-list(log_k = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)), b = factor(matrix(1:6, nrow = 2)))
map6<-list(log_k = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)), b = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)))
map7<-list(log_k = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)), b = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)))
map8<-list(log_k = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)), b = factor(matrix(rep(1, 6), nrow = 2)))
map9<-list(log_k = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)), b = factor(matrix(1:6, nrow = 2)))
map10<-list(log_k = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)), b = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)))
map11<-list(log_k = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)), b = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)))
map12<-list(log_k = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)), b = factor(matrix(rep(1, 6), nrow = 2)))
map13<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(1:6, nrow = 2)))
map14<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(1:3, 2), nrow = 2, byrow = TRUE)))
map15<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(1:2, 3), nrow = 2, byrow = FALSE)))
map16<-list(log_k = factor(matrix(rep(1, 6), nrow = 2)), b = factor(matrix(rep(1, 6), nrow = 2)))


list_of_maps<-list(map1, map2, map3, map4, map5, map6, map7, map8,
                   map9, map10, map11, map12, map13, map14, map15, map16)

t1<-1 #minimum age overall (want to estimte the early growth curve well)
t2<-46 #maximum age in ENP data



p1<-list(
  log_L2 = c(log(22), log(22), log(25)), 
  log_L1 = c(log(15), log(12), log(18)),
  log_L2_M = c(log(21), log(21), log(24)), 
  log_L1_M = c(log(13), log(12), log(16)),
  log_k = log(matrix(rep(0.17, 6), nrow = 2)),
  log_sigma = c(log(1.1), log(0.2), log(0.3)), 
  b = matrix(rep(2.5, 6), nrow = 2))

#list of objects
R_obj_list<-lapply(list_of_maps, function(x){out<-MakeADFun(cmb(gen_growth_mod_single, mod_dat), p1, map = x); return(out)})



#making sure starting values are reasonable
orig_preds<-R_obj_list[[1]]$report()$mu
mod_dat_test<-mod_dat %>% add_column(preds = orig_preds)
ggplot(mod_dat_test) + geom_point(aes(x = Age, y = Length, color = Sex)) + 
  geom_point(aes(x = Age, y = preds), color = "black") +
  facet_wrap(~BluePopulation, labeller = as_labeller(labs_2), ncol = 2) + 
  theme_classic()

#run model
run_mod<-function(obj, map){
  nk<-length(levels(map$log_k))
  nb<-length(levels(map$b))
  #print(nk)
  #print(nb)
  init_fun<-function(){
    L2<- rnorm(3, 23, 2)
    L2_M<-rnorm(3, 22, 2)
    return(list(
      log_L2 = log(L2),
      log_L1 = log(0.5*L2), #so that L1 is below L2
      log_L2_M = log(L2_M),
      log_L1_M = log(0.5*L2_M),
      log_k = rep(log(rnorm(1, 0.17, 0.01)), nk),
      log_sigma = c(log(1.1), log(1.2), log(1.3)), 
      b = rep(runif(1, 1, 5), nb)
    )
    )
  }
  fit<-tmbstan(obj, iter = 20000, chains = 4,
               warmup = 15000,
               thin = 10,
               control = list(adapt_delta = 0.97, max_treedepth = 17), 
               init = init_fun, 
               seed = 555, 
               cores = 1,
               lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), rep(log(0.05), nk), rep(-100, 3), rep(0, nb)), 
               upper = c(rep(log(10000), 3), log(25), log(25), log(25), rep(log(10000), 3),log(25), log(25), log(25), rep(log(0.3), nk), rep(log(10000), 3), rep(15, nb)))
  return(fit)
}



R_fit_list<-mcmapply(run_mod, R_obj_list, list_of_maps)


#rerun models that didn't converge with more iterations
run_mod2<-function(obj, map){
  nk<-length(levels(map$log_k))
  nb<-length(levels(map$b))
  #print(nk)
  #print(nb)
  init_fun<-function(){
    L2<- rnorm(3, 23, 1)
    L2_M<-rnorm(3, 22, 1)
    return(list(
      log_L2 = log(L2),
      log_L1 = log(0.5*L2), #so that L1 is below L2
      log_L2_M = log(L2_M),
      log_L1_M = log(0.5*L2_M),
      log_k = rep(log(rnorm(1, 0.17, 0.01)), nk),
      log_sigma = c(log(1.1), log(1.2), log(1.3)), 
      b = rep(runif(1, 1, 5), nb)
    )
    )
  }
  fit<-tmbstan(obj, iter = 45000, chains = 4,
               warmup = 40000,
               thin = 10,
               control = list(adapt_delta = 0.97, max_treedepth = 17), 
               init = init_fun, 
               seed = 651, 
               cores = 1,
               lower = c(rep(-100, 3), log(5), log(5), log(5), rep(-100, 3),log(5), log(5), log(5), rep(log(0.05), nk), rep(-100, 3), rep(0, nb)), 
               upper = c(rep(log(1000), 3), log(25), log(25), log(25), rep(log(1000), 3),log(25), log(25), log(25), rep(log(0.3), nk), rep(log(1000), 3), rep(15, nb)))
  return(fit)
}

R_fit_list2<-mcmapply(run_mod2, R_obj_list[c(2,6,9,11,13,14)], list_of_maps[c(2,6,9,11,13,14)])

R_fit_list[[2]]<-R_fit_list2[[1]]
R_fit_list[[6]]<-R_fit_list2[[2]]
R_fit_list[[9]]<-R_fit_list2[[3]]
R_fit_list[[11]]<-R_fit_list2[[4]]
R_fit_list[[13]]<-R_fit_list2[[5]]
R_fit_list[[14]]<-R_fit_list2[[6]]

saveRDS(R_fit_list, here("results", "Richards_fit_list.RDS"))
#R_fit_list<-readRDS(here("results", "Richards_fit_list.RDS"))


destination<-here("results", "R_fits_plots.pdf")
pdf(file = destination)
for(i in 1:length(R_fit_list)){
  plot.new()
  text(x=.1, y=.1, paste("Model", i))
  print(mcmc_trace(R_fit_list[[i]]))
  print(mcmc_dens(R_fit_list[[i]]))
  print(mcmc_acf(R_fit_list[[i]]))
  print(stan_rhat(R_fit_list[[i]]))
}
dev.off()

#some models still did not converge and were removed


#outputs
draws_out<-lapply(R_fit_list, function(x){y<-x %>% as_draws_matrix; out<-y[,-ncol(y)]})
saveRDS(draws_out, here("results", "Richards_fit_draws.RDS"))

#draws_out<-readRDS(here("results", "Richards_fit_draws.RDS"))



reps<-list()
loglikes<-list()
for(i in 1:length(draws_out)){
  reps[[i]]<-apply(draws_out[[i]], 1, get_report, obj = R_obj_list[[i]])
  loglikes[[i]]<-t(sapply(reps[[i]], function(x){x$log_ll}))
}

reps<-mcmapply(function(x, y){apply(x, 1, get_report, obj = y)}, draws_out, R_obj_list, SIMPLIFY = FALSE)
loglikes<-mclapply(reps, function(x){t(sapply(x, function(y){y$log_ll}))})


saveRDS(reps, here("results", "Richards_reports_list.RDS"))
saveRDS(loglikes, here("results", "Richards_loglike_list.RDS"))

loos<-lapply(loglikes, loo)
#check diagnostics
for(i in 1:length(loos)){
  print(loos[[i]])
}
loo_compare(loos)

#get pointwise elpds
pointwise_elpd<-lapply(loos, function(x){x$pointwise[,"elpd_loo"]})
lpd_point<-bind_cols(pointwise_elpd) %>% as.matrix()


loo_model_weights(loos) #default method is stacking


#plot parameter posteriors
get_posteriors_to_plot<-function(draws){
  test<-draws %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
    rename_with(~ str_remove(., "exp_log_"), 
                cols = starts_with("exp_log_")) %>% 
    select(!starts_with("log_")) %>% 
    pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
    group_by(Parameter)%>%
    summarise(median = quantile(value, 0.5), lwr = quantile(value, 0.025), upr = quantile(value, 0.975)) %>%
    mutate(Group_number = str_extract(Parameter, pattern = "(?<=\\[)\\d+(?=\\])")) %>%
    mutate(Parameter_group = str_split_i(Parameter, "\\[", 1))
  
  nk<-nrow(test %>% filter(Parameter_group == "k"))
  nb<-nrow(test %>% filter(Parameter_group == "b"))
  
  if(nb == 6){
    grp<-c("F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant",
           "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant", 
           "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant")
  }else if(nb == 3){
    grp<-c("F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant",
           "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant", 
           "Pygmy", "ENP", "Ant")
  }else if(nb == 2){
    grp<-c("F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant",
           "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant", 
           "F", "M")
  }else if(nb == 1){
    grp<-c("F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant",
           "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant", 
           "Single")
  }
  
  if(nk == 6){
    grp<-c(grp, "F*Pygmy", "F*ENP", "F*Ant","M*Pygmy", "M*ENP", "M*Ant", 
           "Pygmy", "ENP", "Ant")
  }else if(nk == 3){
    grp<-c(grp, "Pygmy", "ENP", "Ant", 
           "Pygmy", "ENP", "Ant")
  }else if(nk == 2){
    grp<-c(grp, "F", "M",  
           "Pygmy", "ENP", "Ant")
  }else if(nk == 1){
    grp<-c(grp, "Single",  
           "Pygmy", "ENP", "Ant")
  }
  
  key<-tibble(Parameter = unique(test$Parameter),
              Group = grp)
  out<-test %>% left_join(key, by = "Parameter") %>% 
    mutate(Parameter_group = str_split_i(Parameter_group, "_", 1))
  return(out)
}

#out1<-get_posteriors_to_plot(draws_out[[1]])
#out2<-get_posteriors_to_plot(draws_out[[2]])
mod_names<-1:16
out_df<-lapply(draws_out, get_posteriors_to_plot)
out_df2<-mapply(function(x, y){x %>% add_column(Model = y)}, out_df, mod_names, SIMPLIFY = FALSE) 
out_all<-bind_rows(out_df2)

out_all$Group <-factor(out_all$Group, levels = c("Single" ,"Ant", "ENP", "Pygmy", "F" , "M", 
                                                 "F*Ant","F*ENP", "F*Pygmy",  "M*Ant", "M*ENP", "M*Pygmy"))
p_all<-out_all %>% filter(!Model %in% c(2, 9, 11, 14)) %>% ggplot(aes(x = Model, y = median)) +   
  geom_linerange(aes(x = Model, ymin = lwr, ymax = upr, color = Group),position=position_dodge(width=0.5))+
  geom_point(aes(x = Model, y = median, color = Group), position=position_dodge(width=0.5)) +
  facet_wrap(~Parameter_group, ncol = 1, scales = "free_y") + 
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  theme_minimal() + 
  labs(x = "Model", y = "Estimate")
p_all
#just K and b
out_all %>% filter(Parameter_group %in% c("k", "b")) %>% filter(!Model %in% c(1, 9)) %>% ggplot(aes(x = Model, y = median)) +   
  geom_linerange(aes(x = Model, ymin = lwr, ymax = upr, color = Group),position=position_dodge(width=0.5))+
  geom_point(aes(x = Model, y = median, color = Group), position=position_dodge(width=0.5)) +
  facet_wrap(~Parameter_group, ncol = 1, scales = "free_y") + 
  theme_minimal() + 
  labs(x = "Model", y = "Estimate")


#ggsave("Figures/Posteriors_allmodels.png", p_all, width = 6, height = 7, units = "in", dpi = 600)

#since all the posteriors are very similar we are using simplest model
#model 16
pars_out<-draws_out[[16]]

#get predictions
#new data
new_dat<-expand_grid(pop = 1:3, Age = 1:85, sex = 1:2)
new_dat$Length<-rep(0, nrow(new_dat))
new_dat$Censored<-rep(0, nrow(new_dat))
new_dat$Source<-rep("Soviet", nrow(new_dat)) #so no censoring

#make new object
R_obj16_new<-MakeADFun(cmb(gen_growth_mod_single, new_dat), p1, map = map16)

#get predictions
reps_preds_R16<-apply(pars_out, 1, get_report, obj = R_obj16_new)
preds_R16<-sapply(reps_preds_R16, function(x){x$mu})

saveRDS(preds_R16, here("results", "Richards_preds_m16.RDS"))

quants_R<-apply(preds_R16, 1, quantile, c(0.025, 0.5, 0.975))

new_dat$pred_med<-quants_R[2,]
new_dat$pred_lwr<-quants_R[1,]
new_dat$pred_upr<-quants_R[3,]
new_dat$BluePopulation<-ifelse(new_dat$pop == 1, "Pygmy", ifelse(new_dat$pop == 2, "E N Pacific", "Antarctic"))
new_dat$Sex<-ifelse(new_dat$sex == 1, "Female", "Male")

ggplot() + 
  geom_point(data = Earplug_dat[Earplug_dat$BluePopulation != "W N Pacific",], aes(x = Age, y = Length, color = Sex)) + 
  geom_line(data = new_dat, aes(x = Age, y = pred_med, group = Sex), color = "black") + 
  geom_ribbon(data = new_dat, aes(x = Age, ymin = pred_lwr, ymax = pred_upr, group = Sex), alpha = 0.5) +
  facet_wrap(Sex~BluePopulation) + 
  theme_classic()




# Plot posteriors for sensitivity tests ---------------------------------------

pars_out_all<-R_draws[[16]] %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

pars_out_all$Parameter_fct<-factor(pars_out_all$Parameter, 
                                   levels = c("L1[1]","L1[2]","L1[3]", 
                                              "L1_M[1]","L1_M[2]","L1_M[3]",
                                              "L2[1]","L2[2]", "L2[3]",
                                              "L2_M[1]","L2_M[2]","L2_M[3]",
                                              "k","b","sigma[1]","sigma[2]","sigma[3]"), 
                                   labels = c(expression(atop("L1 Female", "Pygmy")), expression(atop("L1 Female", "ENP")), expression(atop("L1 Female", "Antarctic")), 
                                              expression(atop("L1 Male", "Pygmy")), expression(atop("L1 Male", "ENP")), expression(atop("L1 Male", "Antarctic")), 
                                              expression(atop("L2 Female", "Pygmy")), expression(atop("L2 Female", "ENP")), expression(atop("L2 Female", "Antarctic")), 
                                              expression(atop("L2 Male", "Pygmy")), expression(atop("L2 Male", "ENP")), expression(atop("L2 Male", "Antarctic")), 
                                              "k", "b", "sigma[P]","sigma[E]","sigma[A]"))

pars_out_all$Model<-"All data"
pars_out_S<- Sov_draws %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

pars_out_S$Parameter_fct<-factor(pars_out_S$Parameter, 
                                 levels = c("L1[1]","L1[2]","L1[3]", 
                                            "L1_M[1]","L1_M[2]","L1_M[3]",
                                            "L2[1]","L2[2]", "L2[3]",
                                            "L2_M[1]","L2_M[2]","L2_M[3]",
                                            "k","b","sigma[1]","sigma[2]","sigma[3]"), 
                                 labels = c(expression(atop("L1 Female", "Pygmy")), expression(atop("L1 Female", "ENP")), expression(atop("L1 Female", "Antarctic")), 
                                            expression(atop("L1 Male", "Pygmy")), expression(atop("L1 Male", "ENP")), expression(atop("L1 Male", "Antarctic")), 
                                            expression(atop("L2 Female", "Pygmy")), expression(atop("L2 Female", "ENP")), expression(atop("L2 Female", "Antarctic")), 
                                            expression(atop("L2 Male", "Pygmy")), expression(atop("L2 Male", "ENP")), expression(atop("L2 Male", "Antarctic")), 
                                            "k", "b", "sigma[P]","sigma[E]","sigma[A]"))

pars_out_S$Model<-"Soviet only"
pars_out_J<- J_draws %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

pars_out_J$Parameter_fct<-factor(pars_out_J$Parameter, 
                                 levels = c("L1[1]","L1[2]","L1[3]", 
                                            "L1_M[1]","L1_M[2]","L1_M[3]",
                                            "L2[1]","L2[2]", "L2[3]",
                                            "L2_M[1]","L2_M[2]","L2_M[3]",
                                            "k","b","sigma[1]","sigma[2]","sigma[3]"), 
                                 labels = c(expression(atop("L1 Female", "Pygmy")), expression(atop("L1 Female", "ENP")), expression(atop("L1 Female", "Antarctic")), 
                                            expression(atop("L1 Male", "Pygmy")), expression(atop("L1 Male", "ENP")), expression(atop("L1 Male", "Antarctic")), 
                                            expression(atop("L2 Female", "Pygmy")), expression(atop("L2 Female", "ENP")), expression(atop("L2 Female", "Antarctic")), 
                                            expression(atop("L2 Male", "Pygmy")), expression(atop("L2 Male", "ENP")), expression(atop("L2 Male", "Antarctic")), 
                                            "k", "b", "sigma[P]","sigma[E]","sigma[A]"))
pars_out_J$Model<-"Japanese only"

pars_out<-bind_rows(pars_out_all, pars_out_S, pars_out_J)

ggplot() + geom_density(data = pars_out, aes(x = value, group = Model, color = Model), fill = NA) + 
  facet_wrap(~Parameter_fct, scales = "free", labeller = label_parsed, ncol = 6) + 
  labs(x = "Estimate", y = "density") +
  scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb')) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.spacing = unit(7, "pt"))

# Posterior predictive ----------------------------------------------------
## Base Model-----
#posterior predictive
ndraws <-4000
postpred_S<-list()
for(i in 1:ndraws){
  postpred_S[[i]]<-sim_length(Sov_draws[i,], soviet_mod_dat)
}

pp_df_S<-do.call(rbind,postpred_S)
pp_df_quants_S<-apply(pp_df_S, 2, quantile, c(0.025, 0.5, 0.975))
pp_df_quants_S<-t(pp_df_quants_S) %>% as_tibble() %>% 
  add_column(Age = soviet_mod_dat$Age, BluePopulation = soviet_mod_dat$BluePopulation, Sex = soviet_mod_dat$Sex)


ndraws<-4000
props_S<-rep(NA, nrow(soviet_mod_dat))
for(i in 1:ncol(pp_df_S)){
  tot<-sum(pp_df_S[, i] < round(soviet_mod_dat$Length[i], 6))
  props_S[i]<-tot/ndraws
}

p1_S<-ggplot()+ geom_histogram(aes(x = props_S), color = "white") + 
  labs(x = "Proportion of posterior predictve < data point")
p1_S

pp_df_quants_S$p_val<-props_S
pp_df_quants_S$Source<-soviet_mod_dat$Source
pp_df_quants_S$Censored<-soviet_mod_dat$Censored
pp_df_quants_S<-pp_df_quants_S %>% mutate(BluePopulation2 = ifelse(BluePopulation == "Pygmy", 
                                                                   paste0(BluePopulation, "_", Source), BluePopulation))

saveRDS(pp_df_quants_S, here("results", "Soviet_Posterior_Predictive.RDS"))

p2_S<-ggplot(pp_df_quants_S) + geom_jitter(aes(x = Age, y = p_val, color = as.factor(Censored)), alpha = 0.3) + 
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~BluePopulation2) + 
  labs(y = "Proportion of posterior predictve < data point") + 
  theme(legend.position = "none")
p2_S                     

p3_S<-ggplot() + 
  geom_point(data = soviet_mod_dat, 
             aes(x = Age, y = Length, shape = Sex), color = "gray20") + 
  geom_errorbar(data = pp_df_quants_S, aes(x = Age, ymin = `2.5%`, ymax = `97.5%`, color = BluePopulation), alpha = 0.5) + 
  geom_point(data = pp_df_quants_S, aes(x = Age, y = `50%`, color = BluePopulation), shape = 4, alpha = 0.5) +
  #scale_color_manual(values = pal) +
  scale_y_continuous(breaks = seq(1, 27, by = 5)) +
  facet_grid(Sex~BluePopulation2) + 
  labs(y = "Length (m)") +
  theme_classic() +
  theme(legend.position = "none", 
        strip.background = element_blank())
p3_S


S_post_pred_appendix<-p1_S + p2_S & theme_classic() & theme(legend.position = "none")
S_post_pred_appendix

## Japanese data only for pygmy blue whales ----
#posterior predictive
J_draws<-as_draws_matrix(fit_J)
J_draws<-J_draws[, -ncol(J_draws)]
saveRDS(J_draws, here("results", "Richards_draws_Japanese.RDS"))

ndraws <-4000
postpred_J<-list()
for(i in 1:ndraws){
  postpred_J[[i]]<-sim_length(J_draws[i,], japanese_mod_dat)
}

pp_df_J<-do.call(rbind,postpred_J)
pp_df_quants_J<-apply(pp_df_J, 2, quantile, c(0.025, 0.5, 0.975))
pp_df_quants_J<-t(pp_df_quants_J) %>% as_tibble() %>% 
  add_column(Age = japanese_mod_dat$Age, BluePopulation = japanese_mod_dat$BluePopulation, Sex = japanese_mod_dat$Sex)


ndraws<-4000
props_J<-rep(NA, nrow(japanese_mod_dat))
for(i in 1:ncol(pp_df_J)){
  tot<-sum(pp_df_J[, i] < round(japanese_mod_dat$Length[i], 6))
  props_J[i]<-tot/ndraws
}

props_tib<-tibble(p_val=props_J, Censored=japanese_mod_dat$Censored)

p1_J<-ggplot(props_tib)+ geom_histogram(aes(x = p_val, group = as.factor(Censored), fill = as.factor(Censored)), position = "stack") + 
  labs(x = "Proportion of posterior predictve < data point") + 
  scale_fill_manual(values = c("gray", "red")) + 
  theme(legend.position = "none")
p1_J

pp_df_quants_J$p_val<-props_J
pp_df_quants_J$Source<-japanese_mod_dat$Source
pp_df_quants_J$Censored<-japanese_mod_dat$Censored
pp_df_quants_J<-pp_df_quants_J %>% mutate(BluePopulation2 = ifelse(BluePopulation == "Pygmy", 
                                                                   paste0(BluePopulation, "_", Source), BluePopulation))


p2_J<-ggplot(pp_df_quants_J) + geom_jitter(aes(x = Age, y = p_val, color = as.factor(Censored)), alpha = 0.3) + 
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~BluePopulation2) + 
  labs(y = "Proportion of posterior predictve < data point") + 
  theme(legend.position = "none")
p2_J                     

exp1<-japanese_mod_dat %>% filter(BluePopulation == "Pygmy") %>% 
  group_by(LengthFT) %>% summarise(n = n()) %>%
  ggplot() + geom_col(aes(x = LengthFT, y = n)) 

exp2<-japanese_mod_dat %>% filter(BluePopulation == "Pygmy") %>% 
  ggplot() + geom_point(aes(x = Age, y = LengthFT, shape = as.factor(Sex)))

exp1+exp2 + plot_annotation(title = "Japanese pygmy blue whale data")
p3_J<-ggplot() + 
  geom_point(data = japanese_mod_dat, 
             aes(x = Age, y = Length, shape = Sex), color = "gray20") + 
  geom_errorbar(data = pp_df_quants_J, aes(x = Age, ymin = `2.5%`, ymax = `97.5%`, color = BluePopulation), alpha = 0.5) + 
  geom_point(data = pp_df_quants_J, aes(x = Age, y = `50%`, color = BluePopulation), shape = 4, alpha = 0.5) +
  #scale_color_manual(values = pal) +
  scale_y_continuous(breaks = seq(1, 27, by = 5)) +
  facet_grid(Sex~BluePopulation2) +  #labeller = labeller(BluePopulation2 = as_labeller(lbls))) + 
  labs(y = "Length (m)") +
  theme_classic() +
  theme(legend.position = "none", 
        strip.background = element_blank())
p3_J

p2_J + p3_J + plot_annotation(title = "Model with only Japanese data for pygmy blue whales")

J_post_pred_appendix<-p1_J + p2_J &theme_classic() & theme(legend.position = "none")
J_post_pred_appendix

## All data ----
#read in fit if necessary
R_draws<-readRDS(here("results", "Richards_fit_draws.RDS"))

#just Richards because best model
ndraws<-2000
set.seed(400)
t1<-1
t2<-46

draws_mat<-R_draws[[16]] #using simplest model
postpred<-list()
for(i in 1:ndraws){
  postpred[[i]]<-sim_length(draws_mat[i,], mod_dat)
}

pp_df<-do.call(rbind,postpred)
pp_df_quants<-apply(pp_df, 2, quantile, c(0.025, 0.5, 0.975))
pp_df_quants<-t(pp_df_quants) %>% as_tibble() %>% 
  add_column(Age = mod_dat$Age, BluePopulation = mod_dat$BluePopulation, Sex = mod_dat$Sex)


#posterior predictive p-values
#plotting

ndraws<-2000
props<-rep(NA, nrow(mod_dat))
for(i in 1:ncol(pp_df)){
  tot<-sum(pp_df[, i] < round(mod_dat$Length[i], 6))
  props[i]<-tot/ndraws
}

props_tib<-tibble(p_val=props, Censored=mod_dat$Censored)

p1<-ggplot(props_tib)+ geom_histogram(aes(x = p_val, group = as.factor(Censored), fill = as.factor(Censored)), position = "stack") + 
  labs(x = "Proportion of posterior predictve < data point") + 
  scale_fill_manual(values = c("gray", "red")) + 
  theme(legend.position = "none")
p1

pp_df_quants$p_val<-props
pp_df_quants$Source<-mod_dat$Source
pp_df_quants$Censored<-mod_dat$Censored
pp_df_quants<-pp_df_quants %>% mutate(BluePopulation2 = ifelse(BluePopulation == "Pygmy", 
                                                               paste0(BluePopulation, "_", Source), BluePopulation))

pp_df_quants %>% filter(Censored == 1) %>% ggplot() + geom_histogram(aes(x = p_val))

p2<-ggplot(pp_df_quants) + geom_jitter(aes(x = Age, y = p_val, color = as.factor(Censored)), alpha = 0.3) + 
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~BluePopulation2) + 
  labs(y = "Proportion of posterior predictve < data point") + 
  theme(legend.position = "none")
p2


p1+p2 & theme_classic() & theme(legend.position = "none")







