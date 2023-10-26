library(tidyverse)
library(survival)
library(ricu)
library(lubridate)
library(data.table)

miiv.included <-readRDS("miiv.included.Rds")
censor_day <- 60


surv.time <- 
  miiv.included %>% 
  left_join(miiv$icustays  %>% as.data.table()) %>%
  left_join(miiv$patients %>%as.data.table) %>%
  mutate(survdays = dod - intime) %>% 
  dplyr::select(all_of(colnames(miiv.included)), survdays, dod, intime) %>%
  mutate(survdays = survdays/dminutes(1440)) %>% 
  mutate(age = age/10) %>%
  mutate(
    survdays = case_when(
      survdays > censor_day ~ censor_day,
      is.na(survdays) ~ censor_day, 
      TRUE ~ survdays
    ),
    surv01 = case_when(
    survdays < censor_day ~ 1, 
    TRUE ~ 0
  )) %>% 
  #Two subjects have data entry errors
  #Date of death many years prior to admission
  filter(survdays > -1) %>% #Because date of death does not include time of death, but intubation time is recorded by hour, R assumes time of death is 00:00 and calculates a small negative value for survival time 
  mutate(strata_var = Hmisc::cut2(VFR, g = 3)) 

surv.time %>%
  filter(survdays < 0) %>%
  dplyr::select(intime, dod, survdays)
# Fit stratified models
strata_names <- levels(surv.time$strata_var)

#Survival analysis for all subjects
coxph_fit0 <- coxph(Surv(survdays, surv01) ~VFR + ccperkg  + age + race + sex + safi + sofa, data = surv.time )

#Survival analysis for subjects stratified by VFR tercile
coxph_fit1 <- coxph(Surv(survdays, surv01) ~ccperkg + age + sex + race + safi + sofa, data = surv.time %>% filter(strata_var == strata_names[1])) 
coxph_fit2 <- coxph(Surv(survdays, surv01) ~ccperkg + age + sex + race + safi + sofa, data = surv.time %>% filter(strata_var == strata_names[2])) 
coxph_fit3 <- coxph(Surv(survdays, surv01) ~ccperkg + age + sex + race + safi + sofa, data = surv.time %>% filter(strata_var == strata_names[3])) 

#ANOVA comparing model fit with and without VFR
anova(
  coxph(Surv(survdays, surv01) ~ VFR +ccperkg + age + race + sex + safi + sofa, data = surv.time ), 
  coxph(Surv(survdays, surv01) ~ ccperkg + age + sex + race  + safi + sofa, data = surv.time )
)


library(broom)
coxph.results <- bind_rows(
  coxph_fit1 %>% survfit %>% tidy, 
  coxph_fit2%>% survfit%>%  tidy,
  coxph_fit3 %>% survfit%>% tidy,

  .id = "strata"
) %>% 
  mutate(strata = strata_names[as.numeric(strata)])


coxph.results$strata %>% unique

VFRLevels <- c("Low: 7.8-11.1%", "Middle: 11.1-12.9%", "High: 12.9-22.2%")
VFRLevels <- c("Low", "Middle", "High")

surv.plot <- ggplot(coxph.results, aes(x = time, y = estimate, fill = strata)) +
  geom_step(size = 1, aes( color = strata)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,censor_day)) +
  labs(
       x = "Days since intubation",  
       y = "Probability of survival", 
       color = "VFR") + 
  theme_minimal() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20), 
        aspect.ratio =  1/2, 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15),
        panel.border = element_rect(fill = NA, linewidth = 2)) + 
  scale_color_manual(values = c("green", "blue", "red"),labels =VFRLevels) + 
  scale_fill_manual(values = c("green", "blue", "red")) + 
  guides(fill = "none")


n.at.risk <- coxph.results %>% 
  mutate(time_strata = cut(time, breaks = c(-1, 20, 40, 59.9))) %>%
  arrange(time) %>%
  distinct(strata, time_strata, .keep_all = T) %>%
  mutate(strata = case_when(
    str_detect(strata, "7.82") ~ "Low - 7.8-11.1%",
    str_detect(strata, "22.2") ~ "High - 12.9-22.2%",
    TRUE ~ "Middle - 11.1-12.9%"
  ) %>%
    factor(levels = c("Low - 7.8-11.1%", "Middle - 11.1-12.9%","High - 12.9-22.2%") %>% rev)) %>%
 ggplot(aes(y = strata, x = time, label = n.risk)) + 
  geom_text(size = 5) +
  labs(x = "Number at risk") +
  theme_classic() + 
  theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size =14),
    axis.title.y = element_blank(),
    aspect.ratio = 1/6
    )

pdf("figures/survival.curve.pdf", width = 10, height =8)
surv.plot/n.at.risk
dev.off()
