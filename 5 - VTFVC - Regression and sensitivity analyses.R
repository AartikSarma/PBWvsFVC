rm(list = ls())

dir.create("figures")
dir.create("results")

library(tidyverse)
library(ricu)
library(broom)
library(jtools)
library(mediation)
library(data.table)
library(stargazer)
library(patchwork)
library(lme4)

#Run this function to test for an association between VFR and mortality in 
# the MIMIC-IV and eICU datasets
run_analysis <- function(
  exclusion_criteria, #Additional exclusion criteria for sensitivity analysis
  analysis_description, #Folder name for output
  additional_covariates = ~., #Additional covariates for sensitivity analysis
  run_mediation = FALSE, #Set to TRUE to run mediation analyses
  regression_labs = c('VFR', 'Tidal volume\n(cc/kg)', 'SOFA', 'SaO2:FiO2 ratio*', 'Male', 'Age*', 'Other vs. White', 'Black vs. White')
){
  message("Loading MIMIC data...")
  miiv.included <- readRDS("miiv.included.Rds")
  
  message("Filtering MIMIC data...")
  miiv.included %>%
    as.data.frame() %>%
    filter(!!!rlang::parse_exprs(exclusion_criteria)) %>%
    mutate(age = age/10, safi = safi/10, race = relevel(as.factor(race), "WHITE")) %>%
    as.data.table -> 
    miiv.included
  
  miiv.avg.mortality <- mean(miiv.included$mortality)
  
  message("Loading eICU data...")
  eicu.included <- readRDS("eicu.included.Rds")
  
  eicu.included %>% 
    filter(!!!rlang::parse_exprs(exclusion_criteria)) %>%
    mutate(age = age/10, safi = safi/10, race = factor(race) %>% relevel(ref = "WHITE")) %>%
    as.data.table -> eicu.included
  
  formula.mortality <- mortality~ race + age + sex + safi + sofa + ccperkg + VFR
  formula.mortality <- update.formula(formula.mortality, new = additional_covariates)
  
  
  message("Running MIMIC regressions...")
  lm.miiv.mortality = glm(formula.mortality,
                          family=binomial(link='logit'),
                          miiv.included, control = list(maxit = 100, epsilon=1e-8))
  
  plot_coefs(lm.miiv.mortality, exp = T)+ 
    theme(aspect.ratio = 1, 
          panel.border = element_rect(size = 1, fill = NA)) + 
    labs(x = "Odds ratio") -> miiv.coef.plot
  
  effect_plot(lm.miiv.mortality, "VFR",interval = T, robust = F) + 
    theme(aspect.ratio = 1, 
          panel.border = element_rect(size = 1, fill = NA)) + 
    labs(y = "Adjusted mortality", x = "VT/FVC ratio (%)", title = "MIMIC-IV")+ 
    geom_hline(yintercept = miiv.avg.mortality, color = "red") -> miiv.vtfvc.effect.plot

  eicu.avg.mortality <- mean(eicu.included$mortality)
  
  message("Running eICU regressions...")
  
  formula.mortality <- update.formula(formula.mortality, new = ~. +(1|hospitalid))

  lm.eicu.mortality = glmer(formula.mortality,
                            family=binomial(link='logit'),
                            eicu.included)
  
  plot_coefs(lm.eicu.mortality, exp = T)+ 
    theme(aspect.ratio = 1, 
          panel.border = element_rect(size = 1, fill = NA)) + 
    labs(x = "Odds ratio") ->eicu.coef.plot
  
  effect_plot(lm.eicu.mortality, "VFR",interval = T, robust = T ,cluster = "hospitalid") + 
    theme(aspect.ratio = 1, 
          panel.border = element_rect(size = 1, fill = NA)) + 
    labs(y = "Adjusted mortality", x = "VT/FVC ratio (%)", title = "eICU") + 
    geom_hline(yintercept = eicu.avg.mortality, color = "red") -> eicu.vtfvc.effect.plot
  
  
  message("Saving figures and tables...")
  mimic.title <- paste("MIMIC-IV (N = ", nrow(miiv.included), ")", sep = "")
  eicu.title <- paste("eICU (N = ", nrow(eicu.included), ")", sep = "")
  
  file.prefix <- paste("VTFVC -", analysis_description)
  
  setwd("results")
  dir.create(file.prefix)
  setwd(file.prefix)

  save(list = ls(pattern = "^lm."),file = paste(file.prefix, "regression objects.Rda"))
  save(list = ls(pattern = ".plot"),file = paste(file.prefix, "plots.Rda"))
  
  pdf(file = paste(file.prefix, "regression figures.pdf"), height = 10, width = 16)
  miiv.coef.plot + labs(title = mimic.title) + 
    eicu.coef.plot + labs(title = eicu.title) & 
    scale_y_discrete(labels = regression_labs) & theme(axis.text.y = element_text(hjust = 1))& 
    theme(axis.text = element_text( size = 16),
          axis.text.y.left = element_text(size = 16),
          axis.title = element_text(size = 20), 
          plot.title = element_text(size = 24)) ->p 
  print(p)
  
  miiv.vtfvc.effect.plot + labs(title = mimic.title) + 
    eicu.vtfvc.effect.plot + labs(title = eicu.title) & 
    theme(axis.text = element_text( size = 16), 
          axis.title = element_text(size = 20), 
          plot.title = element_text(size = 24)) ->p 
  print(p)
  
  dev.off()
  
  stargazer(lm.miiv.mortality, lm.eicu.mortality,
            digits = 2,
            out = paste(file.prefix, "regression tables.html"),covariate.labels =
              c('Black vs. white', 'Other vs. white', 'Age*', 'Male', 'SaO2:FiO2 ratio*', 'SOFA', 'Tidal volume (cc/kg)', 'VFR'))

  setwd('..')
  setwd('..')
  
  if(run_mediation){
    formula.mediation <- VFR ~race + ccperkg+ age + sex + height
    
    lm.miiv.mediator =lm(formula.mediation,miiv.included)
    
    message("Running mediation analysis...")
    miiv.mediation.sex = mediate(model.m = lm.miiv.mediator,
                                 model.y =lm.miiv.mortality,
                                 treat = 'sex', 
                                 mediator='VFR',
                                 treat.value = "Female",
                                 control.value = "Male",
                                 boot=F)
    

    miiv.mediation.race.black = mediate(model.m = lm.miiv.mediator,
                                        model.y =lm.miiv.mortality,
                                        treat = 'race', 
                                        mediator='VFR',
                                        treat.value = "BLACK",
                                        control.value = "WHITE",
                                        boot=F)
    
    
    miiv.mediation.race.other = mediate(model.m = lm.miiv.mediator,
                                        model.y =lm.miiv.mortality,
                                        treat = 'race', 
                                        mediator='VFR',
                                        treat.value = "OTHER",
                                        control.value = "WHITE",
                                        boot=F)
    
    miiv.mediation.age = mediate(model.m = lm.miiv.mediator,
                                 model.y =lm.miiv.mortality,
                                 treat = 'age', 
                                 mediator='VFR',
                                 boot=F)

    
    
    message("Running eICU mediation analysis...")
    formula.mediation <- update.formula(formula.mediation, new = ~.+(1|hospitalid))
    lm.eicu.mediator =lmer(formula.mediation,eicu.included)
    
    eicu.mediation.sex = mediate(model.m = lm.eicu.mediator,
                                 model.y =lm.eicu.mortality,
                                 treat = 'sex', 
                                 mediator='VFR',
                                 treat.value = "Female",
                                 control.value = "Male",
                                 boot=F)
    
    eicu.mediation.race.black = mediate(model.m = lm.eicu.mediator,
                                        model.y =lm.eicu.mortality,
                                        treat = 'race', 
                                        mediator='VFR',
                                        treat.value = "BLACK",
                                        control.value = "WHITE",
                                        boot=F)
    
    eicu.mediation.race.other = mediate(model.m = lm.eicu.mediator,
                                        model.y =lm.eicu.mortality,
                                        treat = 'race', 
                                        mediator='VFR',
                                        treat.value = "OTHER",
                                        control.value = "WHITE",
                                        boot=F)
    
    eicu.mediation.age = mediate(model.m = lm.eicu.mediator,
                                 model.y =lm.eicu.mortality,
                                 treat = 'age', 
                                 mediator='VFR',
                                 boot=F)
    
    ls(pattern = "\\.mediation\\.")
    
    mediation.results <- NULL
    for(medres in ls(pattern = "\\.mediation\\.")){
      cohort <- str_extract(medres, "miiv|eicu")
      comparison <- str_remove(medres, ".+mediation\\.")
      med.obj <- get(medres)
      mediation.results <- rbind(
        mediation.results, 
        cbind(cohort = cohort, comparison = comparison, effect = "ACME", estimate  = med.obj$d.avg, ci.lo = med.obj$d.avg.ci[1], ci.hi = med.obj$d.avg.ci[2], pval = med.obj$d.avg.p),
        cbind(cohort = cohort, comparison = comparison, effect = "ADE", estimate = med.obj$z.avg, ci.lo = med.obj$z.avg.ci[1], ci.hi = med.obj$z.avg.ci[2], pval = med.obj$z.avg.p),
        cbind(cohort = cohort, comparison = comparison, effect = "Total effect", estimate  = med.obj$tau.coef, ci.lo = med.obj$tau.ci[1], ci.hi = med.obj$tau.ci[2], pval = med.obj$tau.p ),
        cbind(cohort = cohort, comparison = comparison, effect = "Proportion mediated", estimate  = med.obj$n.avg, ci.lo = med.obj$n.avg.ci[1], ci.hi = med.obj$n.avg.ci[2], pval =med.obj$n.avg.p)
      )
    }
    
    mediation.results %>%  as.data.frame %>%
      filter(!str_detect(effect, "Proportion")) %>%
      mutate(estimate = as.numeric(estimate), 
             ci.lo = as.numeric(ci.lo), 
             ci.hi = as.numeric(ci.hi)) %>%
      mutate(cohort = ifelse(cohort == "eicu", "E-ICU", "MIMIC-IV"), 
             comparison = case_when(
               comparison == "age" ~ "Age (10 years)", 
               comparison == "race.black" ~ "Race: Black vs. White", 
               comparison == "race.other" ~ "Race: Other vs. White", 
               comparison == "sex" ~ "Sex: Female vs. Male"
             ) %>%
               factor(levels = c("Age (10 years)","Sex: Female vs. Male","Race: Black vs. White","Race: Other vs. White" ) )) %>% 
      ggplot(aes(x = exp(estimate), y = effect, yend = effect)) + 
      facet_grid(cohort~comparison) + 
      geom_point(size = 4) + 
      geom_segment(aes(x = exp(ci.lo), xend = exp(ci.hi))) +
      labs(y = "", x = "Odds ratio for mortality") + 
      geom_vline(xintercept = 1, linetype = "dashed") + 
      theme_classic() + 
      theme(panel.border = element_rect(size = 1, fill = NA), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 14), 
            axis.line = element_blank(),
            axis.text = element_text(size = 12),
            aspect.ratio = 0.5) -> mediation.plot
    
    mediation.results %>% as.data.frame %>% 
      pivot_longer(cols = c(estimate, ci.lo, ci.hi), values_to = "value") %>%
      mutate(value = as.numeric(value)) %>%
      pivot_wider(id_cols = c(cohort, effect, name), names_from = comparison) %>%
      mutate(
        black.female = race.black + sex, 
        other.female = race.other + sex, 
      ) %>%
      dplyr::rename(
        white.female = sex, 
        black.male = race.black, 
        other.male = race.other, 
      ) %>% 
      pivot_longer(cols = !matches("cohort|effect|name"), names_to = "comparison") %>% 
      pivot_wider(id_cols = c(cohort, name, comparison), names_from = effect) %>%
      pivot_longer(cols = !matches("cohort|comparison|name"), names_to = "effect") %>%
      mutate(value = exp(value)) %>%
      filter(effect != "Proportion mediated") %>%
      pivot_wider(names_from = name, values_from = value) -> combined.mediation.effects
    
    combined.mediation.effects %>% 
      ggplot(aes(y = effect, x = estimate)) + 
      facet_grid(cohort~comparison) + 
      geom_point(size = 2.5) + 
      geom_vline(xintercept = 1) + 
      geom_segment(aes(x = ci.lo, xend = ci.hi, y = effect, yend = effect), linewidth = 1) + 
      theme(aspect.ratio = 0.5) + 
      theme(axis.text = element_text(size = 18), 
            strip.text =element_text(size = 18), 
            axis.title = element_text(size = 22)) + 
      theme_classic() + 
      theme(panel.border = element_rect(size = 1, fill = NA), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 14), 
            axis.line = element_blank(),
            axis.text = element_text(size = 12),
            aspect.ratio = 0.5) -> merged.mediation.plot
    
    setwd("results")
    setwd(file.prefix)
    
    save(list = ls(pattern = "mediation"),file = paste(file.prefix, "mediation objects.Rda"))
    
    pdf(file = paste(file.prefix, "mediation figures.pdf"), height = 10, width = 16)
    mediation.plot &
      theme(axis.text = element_text(size = 18), 
            strip.text =element_text(size = 18), 
            axis.title = element_text(size = 22)) ->p 
    print(p)
    print(merged.mediation.plot)
    dev.off()
    setwd("..")
    setwd("..")
  }
  
}

#Run the primary analysis and subgroup analyses
run_analysis(exclusion_criteria = "", analysis_description = "Primary analysis", run_mediation = T)
run_analysis(exclusion_criteria = "valid.dp", analysis_description = "Adjust for DP", additional_covariates = ~.+dp,
             regression_labs = c('Driving pressure', 'VFR', 'Tidal volume\n(cc/kg)', 'SOFA', 'SaO2:FiO2 ratio*', 'Male', 'Age*', 'Other vs. White', 'Black vs. White'))
run_analysis(exclusion_criteria = "valid.dp", analysis_description = "Adjust for Crs", additional_covariates = ~.+crs,
             regression_labs = c('Respiratory compliance', 'VFR', 'Tidal volume\n(cc/kg)', 'SOFA', 'SaO2:FiO2 ratio*', 'Male', 'Age*', 'Other vs. White', 'Black vs. White'))
run_analysis(exclusion_criteria = "safi <=235", analysis_description = "Low SF ratio",
             regression_labs = c('VFR', 'Tidal volume\n(cc/kg)', 'SOFA', 'SaO2:FiO2 ratio*', 'Male', 'Age*', 'Other vs. White', 'Black vs. White'))
run_analysis(exclusion_criteria = "ccperkg <= 6.5", analysis_description = "Low VT",
             regression_labs = c('VFR', 'Tidal volume\n(cc/kg)', 'SOFA', 'SaO2:FiO2 ratio*', 'Male', 'Age*', 'Other vs. White', 'Black vs. White'))

