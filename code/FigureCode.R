#Figure Code
#Zoe Rand
#Last Updated: June 3, 2026
#Used to create all figures in manuscript
#Will require running models in other files and saving results in order to run this
#Note figure 1 requires access to full whaling data from IWC in order to recreate, therefore
#is not included in this file

library(tidyverse)
library(here)
library(posterior)
library(patchwork)
library(ggridges)
library(ggpointdensity)
library(ggh4x)

# Read in data ------------------------------------------------------------
BW_corpora_all<-read_csv(here("data", "BW_corpora_all.csv"))
Earplug_dat<-read_csv(here("data", "Earplug_dat_all.csv"))

# Read in model fits ------------------------------------------------------
#base case
Agelength_draws<-read_rds(here("results", "Richards_m16_Soviet_draws.RDS"))
AgeLength_preds<-read_rds(here("results", "Richards_preds_Soviet.RDS"))
AgeLength_reps<-read_rds(here("results", "Soveit_m16_reports.RDS"))
AgeLength_postpreds<-read_rds(here("results", "Soviet_Posterior_Predictive.RDS"))
Corpora_mod_preds<-read_csv(here("results", "Corpora_earplug_preds.csv"))
Corpora_post_pred<-read_csv(here("results", "corpora_posterior_predictive.csv"))

#Age length sensitivity
R_fit_all<-read_rds(here("results", "Richards_fit_draws.RDS"))[[16]]
R_fit_J<-read_rds(here("results", "Richards_draws_Japanese.RDS"))

# Colors ------------------------------------------------------------------

#blue whale subspecies/population
#Pygmy, ENP, Antarctic, WNP
pal<-c('#702963','#1f78b4','#8bcf4e','darkgreen')
names(pal)<-c("Pygmy", "Antarctic", "E N Pacific", "W N Pacific")

#sensitivity palette
sens_pal<-c('#66c2a5','#fc8d62','#8da0cb')


# Figure 2 ----------------------------------------------------------------
#fixing names
Earplug_dat<-Earplug_dat %>% mutate(BluePopulation2 = ifelse(BluePopulation == "Pygmy", 
                                                             paste0(BluePopulation, "_", Source), BluePopulation))


Figure2<-function(){
  #adding data sample sizes
  sample_size_labs<-Earplug_dat %>% group_by(BluePopulation2) %>% summarise(n = n()) %>%
    mutate(BluePopulation_labs = ifelse(BluePopulation2 == "Pygmy_Maeda", "Pygmy \nJapanese data", ifelse(
      BluePopulation2 == "Pygmy_Sazhinov", "Pygmy \nSoviet data", BluePopulation2))) %>%
    mutate(labs = paste0(BluePopulation_labs, "\n n = ", n))
  
  lbls_names<-sample_size_labs$labs
  names(lbls_names)<-sample_size_labs$BluePopulation2
  
  lbls1<-sample_size_labs$BluePopulation_labs
  names(lbls1)<-sample_size_labs$BluePopulation2
  
  #plotting
  Earplug_plot <- ggplot(Earplug_dat) + 
    geom_pointdensity(aes(x = Age, y = Length, shape = Sex), adjust = 4) +
    facet_wrap(~BluePopulation2, labeller = as_labeller(lbls_names), ncol = 2, axes = "all" ) + 
    labs(y = "Length (m)") +
    theme_classic() + 
    scale_color_viridis_c() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = rel(1)), 
          legend.position = "right")
  
  return(print(Earplug_plot))
}

fig2<-Figure2()

#ggsave(here("figures", "Figure2.png"), fig2, dpi = 600, width = 6, height = 6, units = "in")


# Figure 3 ----------------------------------------------------------------

Figure3<-function(){
  #get sample sizes
  sample_size_labs2<-BW_corpora_all %>% group_by(BluePopulation) %>% summarise(n = n()) %>%
    mutate(labs = paste0(BluePopulation, "\n n = ", n))
  
  lbls_names2<-sample_size_labs2$labs
  names(lbls_names2)<-sample_size_labs2$BluePopulation
  
  #plot
  
  OC_plot<-ggplot(BW_corpora_all) + 
    geom_pointdensity(aes(x = Age, y = TotalCorpora), adjust = 4) +
    facet_wrap(~BluePopulation, axes = "all", labeller = as_labeller(lbls_names2)) + 
    labs(y = "# of Corpora") +
    theme_classic() + 
    scale_color_viridis_c() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = rel(1)))
  return(print(OC_plot))
}

fig3<-Figure3()

#ggsave(here("figures", "Figure3.png"),fig3, dpi = 600, width = 6, height = 5, units = "in")


# Figure 4 ----------------------------------------------------------------
#base case
#posterior distributions of parameters
pars_out<-Agelength_draws %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

Figure4<-function(){
  #making labels in the right order and pretty for plot
  pars_out$Parameter_fct<-factor(pars_out$Parameter, 
                                 levels = c("L1[3]","L1[2]","L1[1]", 
                                            "L1_M[3]","L1_M[2]","L1_M[1]",
                                            "L2[3]","L2[2]", "L2[1]",
                                            "L2_M[3]","L2_M[2]","L2_M[1]",
                                            "k","b","sigma[3]","sigma[2]","sigma[1]"), 
                                 labels = c(expression(atop("L1 Female", "Antarctic")), expression(atop("L1 Female", "ENP")), expression(atop("L1 Female", "Pygmy")), 
                                            expression(atop("L1 Male", "Antarctic")), expression(atop("L1 Male", "ENP")), expression(atop("L1 Male", "Pygmy")), 
                                            expression(atop("L2 Female", "Antarctic")), expression(atop("L2 Female", "ENP")), expression(atop("L2 Female", "Pygmy")), 
                                            expression(atop("L2 Male", "Antarctic")), expression(atop("L2 Male", "ENP")), expression(atop("L2 Male", "Pygmy")), 
                                            "k", "b", "sigma[A]","sigma[E]","sigma[P]"))
  pars_out<-pars_out %>% mutate(BluePopulation =ifelse(Parameter %in% c("L1[1]", "L1_M[1]","L2_M[1]","L2[1]","sigma[1]" ), "Pygmy", 
                                                       ifelse(Parameter %in% c("L1[2]", "L1_M[2]","L2_M[2]","L2[2]","sigma[2]"), "E N Pacific", 
                                                              ifelse(Parameter %in% c("L1[3]", "L1_M[3]","L2_M[3]","L2[3]","sigma[3]"), "Antarctic", "All"))))                               
  #add priors
  
  L1_pri_r<-runif(10000, 5, 25)
  L1_pri_d<-dunif(L1_pri_r, 5, 25)
  k_pri_r<-runif(10000, 0.05, 0.3)
  k_pri_d<-dunif(k_pri_r, 0.05, 0.3)
  b_pri_r<-runif(10000, 0, 15)
  b_pri_d<-dunif(b_pri_r, 0, 15)
  prior_df<-tibble(Parameter = rep(c("L1[1]","L1[2]","L1[3]", "L1_M[1]","L1_M[2]","L1_M[3]", "k", "b"), each = 10000), 
                   r_vals = c(rep(L1_pri_r, 6), k_pri_r, b_pri_r),
                   d_vals = c(rep(L1_pri_d, 6), k_pri_d +1, b_pri_d))
  prior_df$Parameter_fct<-factor(prior_df$Parameter, 
                                 levels = c("L1[3]","L1[2]","L1[1]", 
                                            "L1_M[3]","L1_M[2]","L1_M[1]",
                                            "L2[3]","L2[2]", "L2[1]",
                                            "L2_M[3]","L2_M[2]","L2_M[1]",
                                            "k","b","sigma[3]","sigma[2]","sigma[1]"), 
                                 labels = c(expression(atop("L1 Female", "Antarctic")), expression(atop("L1 Female", "ENP")), expression(atop("L1 Female", "Pygmy")), 
                                            expression(atop("L1 Male", "Antarctic")), expression(atop("L1 Male", "ENP")), expression(atop("L1 Male", "Pygmy")), 
                                            expression(atop("L2 Female", "Antarctic")), expression(atop("L2 Female", "ENP")), expression(atop("L2 Female", "Pygmy")), 
                                            expression(atop("L2 Male", "Antarctic")), expression(atop("L2 Male", "ENP")), expression(atop("L2 Male", "Pygmy")), 
                                            "k", "b", "sigma[A]","sigma[E]","sigma[P]"))
  
  #custom x-scales for L2
  x_scales<-list(
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(5, 25), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(20, 27), expand = c(0,0)),
    scale_x_continuous(limits=c(0, 0.4), expand = c(0,0)),
    scale_x_continuous(limits=c(0, 15), expand = c(0,0)),
    scale_x_continuous(limits=c(0, 2), expand = c(0,0)),
    scale_x_continuous(limits=c(0, 2),expand = c(0,0)),
    scale_x_continuous(limits=c(0, 2), expand = c(0,0))
  )
  
  #colors for population
  pal3<-c(pal, "All" = "black")
  
  #plot
  post_plot<-ggplot() + geom_density(data = prior_df, aes(x = r_vals, y = d_vals, group = as.factor(Parameter)), linetype =  "dashed") + 
    geom_density(data = pars_out, aes(x = value, group = Parameter_fct, color = BluePopulation, fill = BluePopulation), alpha = 0.5) + 
    facet_wrap(~Parameter_fct, scales = "free", labeller = label_parsed, ncol = 6) + 
    facetted_pos_scales(x = x_scales) +
    labs(x = "Estimate", y = "density") +
    scale_color_manual(values = pal3) +
    scale_fill_manual(values = pal3) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() + 
    theme(strip.background = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.spacing = unit(12, "pt"), 
          legend.position = "none", 
          strip.text = element_text(margin = margin(t = 0, b = 0, unit = "lines")))
  return(print(post_plot))
}

fig4<-Figure4()

#ggsave(here("figures", "Figure4.png"),fig4 width = 6, height = 5, units = "in", dpi = 600)


# Figure 5 ----------------------------------------------------------------
#data for predictions
new_dat<-expand_grid(pop = 1:3, Age = 1:85, sex = 1:2)
new_dat$Length<-rep(0, nrow(new_dat))
new_dat$Censored<-rep(0, nrow(new_dat))
new_dat$Source<-rep("Soviet", nrow(new_dat)) #so no censoring

#get prediction summary
quants_R<-apply(AgeLength_preds, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

#add to data frame
new_dat$pred_med<-quants_R[2,]
new_dat$pred_lwr<-quants_R[1,]
new_dat$pred_upr<-quants_R[3,]
new_dat$BluePopulation<-ifelse(new_dat$pop == 1, "Pygmy", ifelse(new_dat$pop == 2, "E N Pacific", "Antarctic"))
new_dat$Sex<-ifelse(new_dat$sex == 1, "Female", "Male")

#removing WNP and Japanese data from plotting because not in model
AgeLength_postpreds$BluePopulation2<-Earplug_dat$BluePopulation2[!Earplug_dat$BluePopulation2 %in% c("W N Pacific", "Pygmy_Maeda")]

Figure5<-function(){
  plot1<-
    Earplug_dat[Earplug_dat$BluePopulation != "W N Pacific",] %>% 
    filter(BluePopulation2 != "Pygmy_Maeda") %>% 
    ggplot() + 
    geom_pointdensity(aes(x = Age, y = Length, shape = Sex), alpha = 0.7) + 
    geom_line(data = new_dat, aes(x = Age, y = pred_med, group = Sex), color = "black") + 
    geom_ribbon(data = new_dat, aes(x = Age, ymin = pred_lwr, ymax = pred_upr, group = Sex), alpha = 0.5) +
    #scale_color_manual(values = pal) +
    scale_color_viridis_c() +
    scale_y_continuous(breaks = seq(1, 27, by = 5)) +
    scale_shape_discrete(guide = "none") +
    facet_grid(Sex~BluePopulation, labeller = labeller(BluePopulation = as_labeller(lbls)), axes = "all_x") + 
    labs(y = "Length (m)") +
    theme_classic() +
    theme(strip.background = element_blank(), 
          strip.text.y = element_text(angle = 0, size = rel(1.2)), 
          strip.text.x = element_text(size = rel(1.2)))
  plot1
  
  plot2<-Earplug_dat[Earplug_dat$BluePopulation != "W N Pacific",] %>% 
    filter(BluePopulation2 != "Pygmy_Maeda") %>% 
    ggplot() + 
    geom_point(aes(x = Age, y = Length, shape = Sex), color = "gray20") + 
    geom_errorbar(data = AgeLength_postpreds, aes(x = Age, ymin = `2.5%`, ymax = `97.5%`, color = BluePopulation), alpha = 0.5) + 
    geom_point(data = AgeLength_postpreds, aes(x = Age, y = `50%`, color = BluePopulation), shape = 4, alpha = 0.5) +
    scale_color_manual(values = pal) +
    scale_y_continuous(breaks = seq(1, 27, by = 5)) +
    facet_grid(Sex~BluePopulation, labeller = labeller(BluePopulation = as_labeller(lbls)), axes = "all_x") + 
    labs(y = "Length (m)") +
    theme_classic() +
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text.y = element_text(angle = 0, size = rel(1.2)), 
          strip.text.x = element_text(size = rel(1.2)))
  plot2
  
  plot_tog<-plot1 + plot2 + plot_annotation(tag_levels = "a", tag_suffix = ")") +
    plot_layout(axes = "collect")
  return(print(plot_tog))
}

fig5<-Figure5()

#ggsave(here("figures", "Figure5.png"), fig5, dpi = 600, width = 8.5, height = 5, units = "in")


# Figure 6 ----------------------------------------------------------------

Linf_R<-t(sapply(AgeLength_reps, function(x){x$L_inf}))
colnames(Linf_R)<-c("P_F", "P_M", "ENP_F", "ENP_M", "A_F", "A_M")

Figure6<-function(){
  #making a table and prettier label names for plotting
  Linf_tib<-Linf_R %>% as_tibble() %>% 
    pivot_longer(everything(), names_to = "Pop_Sex", values_to  = "val") %>% 
    mutate(Pop = str_split_i(Pop_Sex, "_", 1), 
           S = str_split_i(Pop_Sex, "_", 2), 
           BluePopulation = ifelse(Pop == "P", "Pygmy", ifelse(Pop == "ENP", "EN Pacific", "Antarctic")), 
           Sex = ifelse(S == "F", "Female", "Male"))
  
  Linf_tib$BluePopulation<-factor(Linf_tib$BluePopulation, levels = c("Pygmy", 
                                                                      "EN Pacific", "Antarctic"))
  Linf_tib$Sex <- factor(Linf_tib$Sex, levels = c("Female", "Male"))
  
  plot3<-ggplot(Linf_tib, aes(x = val, group = Pop_Sex)) + 
    geom_density_ridges(aes(y = BluePopulation, fill= Sex, color = Sex), scale = 0.99, quantile_lines = TRUE, quantiles = c(0.5), 
                        alpha = 0.8, rel_min_height=0.01) + 
    scale_y_discrete(expand = c(0, 0)) + 
    scale_x_continuous(limits = c(20, 29)) +
    labs(x = "Asymptotic length (m)") +
    theme_classic() + 
    theme(axis.ticks.y = element_blank(), 
          axis.text.y = element_text(vjust = -1), 
          axis.title.y = element_blank())
  return(print(plot3))
}

fig6<-Figure6()

#ggsave(here("figures", "Figure6.png"), fig6, dpi = 600, width = 5, height = 4, units = "in")
#
#

# Figure 7 ----------------------------------------------------------------

#posteriors for sensitivity model runs
pars_out_all<-R_fit_all %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
rename_with(~ str_remove(., "exp_log_"), 
            cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

pars_out_S<- Agelength_draws %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

pars_out_J<- R_fit_J %>% as.matrix() %>% as_tibble() %>% mutate(across(starts_with("log"), exp, .names = "exp_{col}")) %>%
  rename_with(~ str_remove(., "exp_log_"), 
              cols = starts_with("exp_log_")) %>% 
  select(!starts_with("log_")) %>% 
  pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>% 
  group_by(Parameter)

Figure7<-function(){
  #making better labels and then combining data frames
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
  
  #combine
  pars_out<-bind_rows(pars_out_all, pars_out_S, pars_out_J)
  
  #plot
  sens_plot<-ggplot() + geom_density(data = pars_out, aes(x = value, group = Model, color = Model), 
                                     fill = NA, key_glyph = draw_key_path) + 
    facet_wrap(~Parameter_fct, scales = "free", labeller = label_parsed, ncol = 6) + 
    labs(x = "Estimate", y = "density") +
    scale_color_manual(values = sens_pal) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() + 
    theme(strip.background = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.spacing = unit(12, "pt"), 
          legend.position = "bottom", 
          strip.text = element_text(margin = margin(t = 0, b = 0, unit = "pt")))
  return(print(sens_plot))
}

fig7<-Figure7()
#ggsave(here("figures", "Figure7.png"), fig7, width = 10, height = 6, units = "in", dpi = 600)


# Figure 8 ----------------------------------------------------------------

Figure8<-function(){
  #plot predictions
  
  corp_1<-ggplot() + geom_point(data = BW_corpora_all[BW_corpora_all$BluePopulation %in% c("E N Pacific","Pygmy"),], 
                                aes(x = Age, y = TotalCorpora, color = BluePopulation)) + 
    geom_ribbon(data = Corpora_mod_preds, aes(x = Age, ymin = lwr, ymax = upr, fill = BluePopulation), alpha = 0.5) + 
    geom_line(data = Corpora_mod_preds, aes(x = Age, y = med, color = BluePopulation)) + 
    facet_wrap(~BluePopulation, labeller = as_labeller(lbls2)) + 
    scale_fill_manual(values = pal) + 
    scale_color_manual(values = pal) +
    theme_classic() + 
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 25)) + 
    labs(y = "# of Corpora") +
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_text(size = rel(1.1)))
  
  corp_1
  
  #posterior predictive
  
  corp_2<-ggplot(Corpora_post_pred) + geom_point(aes(x = Age, y = TotalCorpora)) + 
    geom_point(aes(x = Age, y = `50%`, color = BluePopulation), shape = 8) + 
    geom_errorbar(aes(x = Age, ymin = `2.5%`, ymax = `97.5%`, color = BluePopulation)) + 
    facet_wrap(~BluePopulation, labeller = as_labeller(lbls2)) +
    scale_color_manual(values = pal) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "# of Corpora") +
    theme_classic() + 
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_text(size = rel(1.1)))
  
  
  corp_2
  
  #plot together
  corp_together<-corp_1 / corp_2 + plot_annotation(tag_levels = "a", tag_suffix = ")") +
    plot_layout(axes = "collect")
  return(print(corp_together))
}

fig8<-Figure8()

#ggsave(here("figures", "Figure8.png"), fig8, dpi = 600, width = 6, height = 5, units  = "in")









