#Ovarian Corpora Model
#Zoe Rand
#Last updated: 12/16/25
library(tidyverse)
library(here)
library(posterior)
library(bayesplot)
library(RTMB)
library(RTMBdist)
library(tmbstan)
library(loo)
library(brms)

# Data --------------------------------------------------------------------

#Ichihara 1966 data for pygmy blue whales
Py_corpora<-read_csv(here("data", "Pygmy_corpora_ichihara.csv"))

#other blue whale corpoa
BW_earplug<-read_csv(here("data", "BW_Earplugs_Maeda_and_Mizroch.csv"), col_types = list(TotalCorpora = "n"))

#convert length to meters and combine data 
#remove uncertain ages

BW_corpora<-BW_earplug %>% filter(!is.na(TotalCorpora)) %>% 
  mutate(Length = LengthFT * 0.3048) %>% 
  filter(AccurateAge == "Yes")

Py_corpora <- Py_corpora %>% mutate(Length = `Length (ft)` * 0.3048) %>% 
  mutate(Age = `Earplug laminae`) %>%
  rename(TotalCorpora = `Ovarian corpora`) %>%
  add_column(Source = "Ichihara", BluePopulation = "Pygmy")

BW_corpora_all<-bind_rows(BW_corpora, Py_corpora)

write_csv(BW_corpora_all, here("data", "BW_corpora_all.csv"))

#plot
BW_corpora_all %>% filter(BluePopulation == "Pygmy") %>% ggplot() + geom_point(aes(x = Age, y = TotalCorpora ))
BW_corpora_all %>% ggplot() + 
  geom_point(aes(x = Age, y = TotalCorpora, color = BluePopulation)) + 
  theme_classic()

BW_corpora_all %>% ggplot() + 
  geom_point(aes(x = Age, y = TotalCorpora-mean(TotalCorpora), color = BluePopulation)) + 
  theme_classic()

BW_corpora_all %>% group_by(BluePopulation, Source) %>% summarise(n = n())


# Data Wrangling ----------------------------------------------------------

#model data
mod_dat<-BW_corpora_all %>%  
  mutate(pop = ifelse(BluePopulation == "Pygmy", 1, 
                      ifelse(BluePopulation == "E N Pacific", 2, 
                             ifelse(BluePopulation == "Antarctic", 3, 4)))) %>% 
  filter(!BluePopulation %in% c("W N Pacific", "Antarctic")) %>%
  mutate(corpora = TotalCorpora) %>% select(Age, pop, corpora)

#write_csv(mod_dat, here("data", "corpora_model_data.csv"))


# Model functions ---------------------------------------------------------
#inverse gamma for prior
dinv_gamma <- function(y, a, b){b^a / gamma(a) * y^(-a-1)*exp(-b/y)}

#smooth max to prevent negative predictions
smooth_max <- function(x, HI, steepness = 100) {
  x + log(1 + exp(steepness * (HI - x))) / steepness
}


corpora_model_NB_noRE<-function(pars){
  require(RTMB)
  #require(RTMBdist)
  getAll(pars, mod_dat, warn = FALSE)

  
  #parameters
  b_Age<-exp(log_b_Age) #slope constrained to be positive--can be up to 1 per population
  theta<-exp(log_theta) #negative binomial variance
  
  
  #priors
  pri<-0
  pri<- pri - log(dgamma(theta, 0.4, 0.3)) #inverse gamma prior on overdispersion
  
  #jacobian for uniform priors on log scale
  pri<- pri + sum(log_b_Age)
  
  post<-pri
  
  
  #model
  predCorpora<-rep(0, nrow(mod_dat))
  for(i in 1:nrow(mod_dat)){
    pred<-a_Pop[mod_dat$pop[i]] + b_Age[mod_dat$pop[i]]*mod_dat$Age[i]
    predCorpora[i]<-smooth_max(pred, 0.001) #making sure response is positive
  }
  
  
  #data likelihood
  post<-post - sum(dnbinom2(corpora, predCorpora,theta, log = TRUE))
  
  #pointwise log likelihood
  log_ll<-rep(0, nrow(mod_dat))
  for(i in 1:nrow(mod_dat)){
    log_ll[i]<-dnbinom2(mod_dat$corpora[i],
                        predCorpora[i], theta, log = TRUE)
  }
  nll<-sum(-1*log_ll)
  
  #transformed parameters
  rate <- 1/b_Age #rate of corpora production
  Age_mat<- -1*(a_Pop)/b_Age
  
  REPORT(b_Age)
  REPORT(predCorpora)
  REPORT(rate)
  REPORT(Age_mat)
  REPORT(log_ll)
  REPORT(nll)
  
  return(post)
  
}

simulate_corpora_model_NB_noRE<-function(pars){
  #parameters
  b_Age<-exp(pars[1:2]) #slope constrained to be positive
  theta<-exp(pars[3]) #negative binomial variance
  a_Pop<-pars[4:5]
  
  
  
  #model
  predCorpora<-rep(0, nrow(mod_dat))
  for(i in 1:nrow(mod_dat)){
    pred<-a_Pop[mod_dat$pop[i]] + b_Age[mod_dat$pop[i]]*mod_dat$Age[i]
    predCorpora[i]<-smooth_max(pred, 0.001) #making sure response is positive
  }
  
  #data likelihood
  corpora_sim<-rnbinom2(length(predCorpora), predCorpora,theta)
  
  #transformed parameters
  rate <- 1/b_Age #rate of corpora production
  Age_mat<- -1*(a_Pop)/b_Age
  
  
  return(corpora_sim)
}

get_report<-function(draw, obj){
  rep<-obj$report(draw)
  return(rep)
  #return(rep$mu)
}


# Run Model ---------------------------------------------------------------
#model options:
#1) Intercept and slope are population specific
#2) Intercept is population specific but one slope
#3) One intercept and slope is population specific
#4) One intercept and one slope

map1<-list(log_b_Age = factor(c(1,2)), a_Pop = factor(c(1,2)))
map2<-list(log_b_Age = factor(c(1,1)), a_Pop = factor(c(1,2)))
map3<-list(log_b_Age = factor(c(1,2)), a_Pop = factor(c(1,1)))
map4<-list(log_b_Age = factor(c(1,1)), a_Pop = factor(c(1,1)))


p<-list(log_b_Age = c(log(0.3), log(0.2)),
        log_theta = log(1),
        a_Pop = c(-3.2,-3)
)

corpora_obj1<-MakeADFun(corpora_model_NB_noRE, p, map = map1)  
corpora_obj1$fn()
corpora_obj1$report()

corpora_obj2<-MakeADFun(corpora_model_NB_noRE, p, map = map2)  
corpora_obj2$fn()
corpora_obj2$report()

corpora_obj3<-MakeADFun(corpora_model_NB_noRE, p, map = map3)  
corpora_obj3$fn()
corpora_obj3$report()

corpora_obj4<-MakeADFun(corpora_model_NB_noRE, p, map = map4)  
corpora_obj4$fn()
corpora_obj4$report()

## fitting model 1 ----
init_fun1<-function(){
  list(log_b_Age = log(rnorm(2, 0.3, 0.01)),
       log_theta = log(runif(1, 1, 10)),
       a_Pop = rnorm(2, -3, 0.2)
  )
}

corpora_fit1<-tmbstan(corpora_obj1, iter = 10000, chains = 4, 
                      warmup = 5000,
                      thin = 10,
                      control = list(adapt_delta = 0.98, max_treedepth = 15), init = init_fun1, 
                      lower = c(-100, -100, -Inf, -1000, -1000), 
                      upper = c(log(10000), log(10000), Inf, 1000, 1000))

mcmc_trace(corpora_fit1, pars = c("log_b_Age[1]", "log_b_Age[2]",  'log_theta', "a_Pop[1]","a_Pop[2]"))

mcmc_dens(corpora_fit1, pars = c("log_b_Age[1]", "log_b_Age[2]",  'log_theta', "a_Pop[1]","a_Pop[2]"))

mcmc_acf(corpora_fit1, pars = c("log_b_Age[1]", "log_b_Age[2]",  'log_theta', "a_Pop[1]","a_Pop[2]"))

saveRDS(corpora_fit1, here("results", "Corpora_earplug_fit_mod1.RDS"))

post1<-corpora_fit1 %>% as_draws_matrix()
post1<-post1[,-ncol(post1)]


#get reports
ndraws<-500*4
reps1<-apply(post1, 1, get_report, obj = corpora_obj1)
saveRDS(reps1, here("results", "mod1_reports.RDS"))

## fitting model 2 ----
init_fun2<-function(){
  list(log_b_Age = log(rnorm(1, 0.3, 0.01)),
       log_theta = log(runif(1, 1, 10)),
       a_Pop = rnorm(2, -3, 0.2)
  )
}
corpora_fit2<-tmbstan(corpora_obj2, iter = 10000, chains = 4, 
                      warmup = 5000,
                      thin = 10,
                      control = list(adapt_delta = 0.98, max_treedepth = 15), init = init_fun2, 
                      lower = c(-100, -Inf, -1000, -1000), 
                      upper = c(log(10000), Inf, 1000, 1000))

mcmc_trace(corpora_fit2, pars = c("log_b_Age",'log_theta', "a_Pop[1]","a_Pop[2]"))

mcmc_dens(corpora_fit2, pars = c("log_b_Age",'log_theta', "a_Pop[1]","a_Pop[2]"))

mcmc_acf(corpora_fit2, pars = c("log_b_Age",'log_theta', "a_Pop[1]","a_Pop[2]"))

saveRDS(corpora_fit2, here("results", "Corpora_earplug_fit_mod2.RDS"))

post2<-corpora_fit2 %>% as_draws_matrix()
post2<-post2[, -ncol(post2)]

#get reports
ndraws<-500*4
reps2<-apply(post2, 1, get_report, obj = corpora_obj2)
saveRDS(reps2, here("results", "mod2_reports.RDS"))


## fitting model 3 ----
init_fun3<-function(){
  list(log_b_Age = log(rnorm(2, 0.3, 0.01)),
       log_theta = log(runif(1, 1, 10)),
       a_Pop = rnorm(1, -3, 0.2)
  )
}
corpora_fit3<-tmbstan(corpora_obj3, iter = 10000, chains = 4, 
                      warmup = 5000,
                      thin = 10,
                      control = list(adapt_delta = 0.98, max_treedepth = 15), init = init_fun3, 
                      lower = c(-100, -100, -Inf, -1000), 
                      upper = c(log(10000), log(10000), Inf, 1000))

mcmc_trace(corpora_fit3, pars = c("log_b_Age[1]","log_b_Age[2]",'log_theta', "a_Pop"))

mcmc_dens(corpora_fit3, pars = c("log_b_Age[1]","log_b_Age[2]",'log_theta', "a_Pop"))

mcmc_acf(corpora_fit3, pars = c("log_b_Age[1]","log_b_Age[2]",'log_theta', "a_Pop"))

saveRDS(corpora_fit3, here("results", "Corpora_earplug_fit_mod3.RDS"))

post3<-corpora_fit3 %>% as_draws_matrix()
post3<-post3[, -ncol(post3)]

#get reports
ndraws<-500*4
reps3<-apply(post3, 1, get_report, obj = corpora_obj3)
saveRDS(reps3, here("results", "mod3_reports.RDS"))

## fitting model 4 ----
init_fun4<-function(){
  list(log_b_Age = log(rnorm(1, 0.3, 0.01)),
       log_theta = log(runif(1, 1, 10)),
       a_Pop = rnorm(1, -3, 0.2)
  )
}
corpora_fit4<-tmbstan(corpora_obj4, iter = 10000, chains = 4, 
                      warmup = 5000,
                      thin = 10,
                      control = list(adapt_delta = 0.98, max_treedepth = 15), init = init_fun4, 
                      lower = c(-100, -Inf, -1000), 
                      upper = c(log(10000), Inf, 1000))

mcmc_trace(corpora_fit4, pars = c("log_b_Age",'log_theta', "a_Pop"))

mcmc_dens(corpora_fit4, pars = c("log_b_Age",'log_theta', "a_Pop"))

mcmc_acf(corpora_fit4, pars = c("log_b_Age",'log_theta', "a_Pop"))

saveRDS(corpora_fit4, here("results", "Corpora_earplug_fit_mod4.RDS"))

post4<-corpora_fit4 %>% as_draws_matrix()
post4<-post4[, -ncol(post4)]

#get reports
ndraws<-500*4
reps4<-apply(post4, 1, get_report, obj = corpora_obj4)
saveRDS(reps4, here("results", "mod4_reports.RDS"))

# Model comparison --------------------------------------------------------
loglike1<-t(sapply(reps1, function(x){x$log_ll}))
elpd(loglike1)
mod1_loo<-loo(loglike1)

loglike2<-t(sapply(reps2, function(x){x$log_ll}))
elpd(loglike2)
mod2_loo<-loo(loglike2)

loglike3<-t(sapply(reps3, function(x){x$log_ll}))
elpd(loglike3)
mod3_loo<-loo(loglike3)

loglike4<-t(sapply(reps4, function(x){x$log_ll}))
elpd(loglike4)
mod4_loo<-loo(loglike4)

loo_compare(mod1_loo, mod2_loo, mod3_loo, mod4_loo)

# Results -------------------------------------------------------------------
#using model 1 since it is the most flexible

## plot predictions----
corpora_fit1<-readRDS(here("results" "Corpora_earplug_fit_mod1.RDS"))
ndraws<-500*4
fit_draws<-as_draws_matrix(corpora_fit1)
# 2 is ENP and 1 is pygmy
preds_ENP<-list()
preds_P<-list()
ENP_rates<-rep(NA, nrow(fit_draws))
P_rates<-rep(NA, nrow(fit_draws))
P_Age_mats<-rep(NA, nrow(fit_draws))
ENP_Age_mats<-rep(NA, nrow(fit_draws))
Age_pred<-seq(2, 50, by = 1)
P_slope<-rep(NA, nrow(fit_draws))
ENP_slope<-rep(NA, nrow(fit_draws))


for(i in 1:nrow(fit_draws)){
  r <- corpora_obj1$report(fit_draws[i,-ncol(fit_draws)])
  preds_ENP[[i]] <- fit_draws[i, "a_Pop[2]"][1] + r$b_Age[2]*Age_pred
  ENP_rates[i]<-r$rate[2]
  ENP_slope[i]<-r$b_Age[2]
  ENP_Age_mats[i]<-r$Age_mat[2]
}

preds_mat_ENP<-matrix(unlist(preds_ENP), ncol = length(Age_pred), byrow = TRUE)
colnames(preds_mat_ENP)<-as.character(Age_pred)


for(i in 1:nrow(fit_draws)){
  r <- corpora_obj1$report(fit_draws[i,-ncol(fit_draws)])
  preds_P[[i]] <- fit_draws[i, "a_Pop[1]"][1] + r$b_Age[1]*Age_pred
  P_rates[i]<-r$rate[1]
  P_slope[i]<-r$b_Age[1]
  P_Age_mats[i]<-r$Age_mat[1]
}

preds_mat_P<-matrix(unlist(preds_P), ncol = length(Age_pred), byrow = TRUE)
colnames(preds_mat_P)<-as.character(Age_pred)


preds_ENP<-preds_mat_ENP %>% as_tibble() %>% add_column(BluePopulation = "E N Pacific", draw = 1:ndraws)
preds_P<-preds_mat_P %>% as_tibble() %>% add_column(BluePopulation = "Pygmy", draw = 1:ndraws)


preds_all<-bind_rows(preds_ENP, preds_P) %>% pivot_longer(-c(BluePopulation, draw), names_to = "Age", values_to = "pred") %>%
  group_by(BluePopulation, Age) %>% summarise(med = quantile(pred, 0.5), lwr = quantile(pred, 0.025), upr = quantile(pred, 0.975)) %>% 
  mutate(Age = as.numeric(Age))

#save for manuscript plotting
#write_csv(preds_all, "Results/Corpora_earplug_preds.csv")

ggplot() + geom_point(data = BW_corpora_all[BW_corpora_all$BluePopulation != "W N Pacific",], aes(x = Age, y = TotalCorpora, color = BluePopulation)) + 
  geom_ribbon(data = preds_all, aes(x = Age, ymin = lwr, ymax = upr, fill = BluePopulation), alpha = 0.5) + 
  geom_line(data = preds_all, aes(x = Age, y = med, color = BluePopulation)) + 
  facet_wrap(~BluePopulation) + theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 20))

## numerical results ----
#rate of corpora estimate
#ENP
quantile(ENP_rates, c(0.025, 0.5, 0.975))
sd(ENP_rates)
summary(ENP_rates)

#pygmy
quantile(P_rates, c(0.025, 0.5, 0.975))
sd(P_rates)
summary(P_rates)
#getting approximate posterior
rsnorm<-rskew_normal(2000, 2.899, 0.548, 1)
dsnorm<-dskew_normal(rsnorm, 2.899, 0.548, 1)

#posterior of rates
ggplot() + geom_histogram(aes(x = P_rates, y = after_stat(density)), bins = 50) + 
  geom_density(aes(x = rsnorm, y = dsnorm), stat = "identity", color = "blue") + 
  theme_classic() + 
  labs(x = "Pymgy blue whale corpora formation rate")


#age at maturity estimates
Age_matrix<-matrix(c(P_Age_mats,ENP_Age_mats), ncol = 2, byrow = FALSE)
colnames(Age_matrix)<-c("Pygmy", "E N Pacific")

Age_quants<-apply(Age_matrix, 2, quantile, c(0.025, 0.5, 0.975))
Age_quants

#slope estimates
quantile(P_slope, c(0.025, 0.5, 0.975))
quantile(ENP_slope, c(0.025, 0.5, 0.975))

## posterior predictions ----

post_pred<-apply(fit_draws[, -ncol(fit_draws)], 1, simulate_corpora_model_NB_noRE) %>%
  apply(., 1, quantile, c(0.025, 0.5, 0.975)) %>% t() %>% as_tibble()


post_pred_plot<-BW_corpora_all %>% filter(!BluePopulation %in% c("W N Pacific", "Antarctic")) %>% 
  bind_cols(post_pred)


ggplot(post_pred_plot) + geom_point(aes(x = Age, y = TotalCorpora)) + 
  geom_point(aes(x = Age, y = `50%`, color = BluePopulation), shape = 8) + 
  geom_errorbar(aes(x = Age, ymin = `2.5%`, ymax = `97.5%`, color = BluePopulation)) + 
  facet_wrap(~BluePopulation)


