
library(dplyr)
library(lubridate)
library(gombedatabase)
library(ggplot2)
library(lmerTest)
library(sjPlot) 
library(coxme)
library(survminer)
library(tidyr)
library(MuMIn)
library(performance)
library(DHARMa)
library(ggthemes)
library(lavaan)
library(lavaanPlot)


# ================ offspring survival to 1 ===============

infsurv_soc_1 = readRDS(file = "infsurv_1_anon.rds")

sapply(infsurv_soc_1, function(x) sum(is.na(x)))

bin_fit_soc = 
  glm(survto1 ~ 
        sex2 + 
        firstborn2 +
        mom_rank11_z +
        mom_age_z + 
        mom_imm_status2 + 
        mom_age_z2 +
        mom_siv_status2 + 
        nmales + 
        nfems 
      , data = infsurv_soc_1, family = binomial, na.action = na.fail)

check_collinearity(bin_fit_soc)

simout = simulateResiduals(bin_fit_soc)
plot(simout)

summary(bin_fit_soc)

allm_surv_soc_1 = dredge(bin_fit_soc, trace = 2, subset = !(mom_age_z2 && !(mom_age_z)))

allm_surv_soc_1 = subset(allm_surv_soc_1, !nested(.))

bestm_surv_soc_1 = get.models(allm_surv_soc_1, subset = delta < 6)

model.sel(bestm_surv_soc_1)

plot_model(bestm_surv_soc_1$`321`) + geom_hline(yintercept = 1, color = "gray50")

null_fit_soc = 
  glm(survto1 ~
        sex2 +
        nfems
      , data = infsurv_soc_1, family = binomial, na.action = na.fail)

fit_soc_csi_f_ratio = update(null_fit_soc, . ~ . + csi_f_ratio)
fit_soc_csi_m_ratio = update(null_fit_soc, . ~ . + csi_m_ratio)
fit_soc_csi_f_m_ratio = update(null_fit_soc, . ~ . + csi_m_ratio + csi_f_ratio)

AICc(null_fit_soc, fit_soc_csi_f_ratio, fit_soc_csi_m_ratio, fit_soc_csi_f_m_ratio) %>% 
  mutate(delta = AICc - AICc(null_fit_soc),
         weight = MuMIn::Weights(AICc)) %>% 
  arrange(AICc) %>%
  round(3)

# check model assumptions for best-fit model:
simout = DHARMa::simulateResiduals(fit_soc_csi_f_ratio) 
plot(simout)

# examine model results
plot_model(fit_soc_csi_f_ratio, type = "pred", terms = "csi_f_ratio [all]", show.data = T)


# results explained by differences in female observation time?
fit_soc_greg = update(null_fit_soc, . ~ . + time_obs_overall_ratio)

AICc(null_fit_soc, fit_soc_greg) %>% 
  mutate(delta = AICc - AICc(null_fit_soc)) %>% 
  arrange(AICc) %>%
  round(3)


# predictor effects in response scale:

newdf = data.frame(sex2 = 0,
                   nfems = mean(infsurv_soc_1$nfems),
                   csi_f_ratio = c(min(infsurv_soc_1$csi_f_ratio), 0.5, 1, 2, 3, max(infsurv_soc_1$csi_f_ratio)))

predres = data.frame(probs = predict(fit_soc_csi_f_ratio, newdata = newdf, type = "response"),
                     csi_f_ratio = newdf$csi_f_ratio)
predres

newdf = data.frame(sex2 = 0,
                   nfems = c(min(infsurv_soc_1$nfems), round(mean(infsurv_soc_1$nfems)), max(infsurv_soc_1$nfems)),
                   csi_f_ratio = mean(infsurv_soc_1$csi_f_ratio))

predres = data.frame(probs = predict(fit_soc_csi_f_ratio, newdata = newdf, type = "response"),
                     nfems = newdf$nfems)
predres

newdf = data.frame(sex2 = c(-1, 1),
                   nfems = mean(infsurv_soc_1$nfems),
                   csi_f_ratio = mean(infsurv_soc_1$csi_f_ratio))

predres = data.frame(probs = predict(fit_soc_csi_f_ratio, newdata = newdf, type = "response"),
                     nfems = newdf$sex2)
predres


# csi_f_ratio vs its component parts: 

fit_soc_time_f_ratio = update(null_fit_soc, . ~ . + time_obs_f_ratio)
fit_soc_grmrate_f_ratio = update(null_fit_soc, . ~ . + grmrate_f_ratio)
fit_soc_time_grmrate_f_ratio = update(null_fit_soc, . ~ . + time_obs_f_ratio + grmrate_f_ratio)


AICc(fit_soc_csi_f_ratio, fit_soc_time_f_ratio, fit_soc_grmrate_f_ratio, fit_soc_time_grmrate_f_ratio, null_fit_soc) %>% 
  arrange(AICc) %>% 
  mutate(delta = AICc - AICc(null_fit_soc))

set_theme(theme_bw())

# CSI-F vs kin presence
infsurv_soc_1 %>%
  mutate(totkin = daughters11 + grandmaingrp + replace_na(matsisters11, 0),
         kintype = case_when(daughters11 == 0 & grandmaingrp == 0 ~ "No daughters or mom",
                             daughters11 > 0 ~ "Adult daughter present",
                             grandmaingrp > 0 ~ "Mom's mom present")) %>% # 
  ggplot(aes(x = totkin, y = csi_f_ratio)) + 
  # geom_violin(draw_quantiles = 0.5) +
  geom_point(aes(color = kintype)) + 
  # geom_violin(aes(x = totyoungkin, y = csi_f), draw_quantiles = 0.5) +
  geom_smooth(method = "lm")

infsurv_soc_1 %>%
  mutate(kintype = case_when(daughters11 == 0 & grandmaingrp == 0 ~ "No daughters or mom",
                             daughters11 > 0 ~ "Adult daughter present",
                             grandmaingrp > 0 ~ "Mom's mom present"),
         died = survto1 == FALSE,
         survto1 = as.factor(survto1)) %>%
  ggplot(aes(x = survto1, y = csi_f_ratio)) + # csi_f
  geom_hline(yintercept = 1, color = "gray") + # 0
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(aes(shape = died, size = died, alpha = died, color = kintype), alpha = 0.7, width = 0.3, height = 0) + 
  scale_size_manual(values = c(2,4)) +
  scale_color_colorblind() +
  theme(legend.position = "none") + 
  xlab("Offspring survives to 1") + ylab("Maternal CSI with other females") +
  labs(tag = "A") + 
  labs(color = "Maternal kin presence", shape = "Offspring died", size = "Offspring died")

infsurv_soc_1 %>%
  mutate(kintype = case_when(daughters11 == 0 & grandmaingrp == 0 ~ "No daughters or mom",
                             daughters11 > 0 ~ "Adult daughter present",
                             grandmaingrp > 0 ~ "Mom's mom present"),
         died = survto1 == FALSE) %>%
  ggplot(aes(x = mom_age_z, y = csi_f_ratio)) + # basically same with csi_f_ratio & csi_f
  geom_hline(yintercept = 1, color = "gray") + # 0
  geom_point(aes(color = kintype, size = died, shape = died, alpha = died)) + 
  scale_shape_manual(values = c(16,17)) +
  scale_size_manual(values = c(2,4)) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_color_colorblind() +
  geom_smooth() + 
  xlab("Maternal age") + ylab("") +
  labs(tag = "B") + 
  labs(color = "Maternal kin present?", shape = "Offspring died", size = "Offspring died", alpha = "Offspring died")


# ------------------- sociality vs kin presence --------------------------

# basically, the model containing CSI-F fits better than any model I could 
#   come up with where I substituted kin presence for CSI-F

head(infsurv_soc_1)

summary(fit_soc_csi_f_ratio)

temp = 
  infsurv_soc_1 %>% 
  mutate(gmapres = grandmaingrp > 0,
         daughtpres = daughters11 > 0,
         matsispres = replace_na(matsisters11, 0) > 0,
         femkinpres = gmapres + daughtpres + matsispres > 0,
         totfemkin = gmapres + daughtpres + matsispres,
         gmapres_next = grandmaingrp_nextyr > 0,
         daughtpres_next = daughters11_nextyr > 0,
         matsispres_next = replace_na(matsisters11_nextyr, 0) > 0,
         femkinpres_next = gmapres_next + daughtpres_next + matsispres_next > 0,
         totfemkin_next = gmapres_next + daughtpres_next + matsispres_next)

temp %>% count(gmapres, grandmaingrp)
temp %>% count(daughtpres, daughters11)
temp %>% count(matsispres, matsisters11)

tempsocmod = glm(survto1 ~
                   sex2
                 + nfems
                 + csi_f_ratio   
                 , data = temp, family = binomial, na.action = na.fail)
AICc(tempsocmod, fit_soc_csi_f_ratio) # same as best model from main analysis

tempkintypemod = glm(survto1 ~
                       sex2
                     + nfems
                     + gmapres
                     + daughtpres
                     + matsispres
                     , data = temp, family = binomial, na.action = na.fail)

tempkintypetotmod = glm(survto1 ~
                          sex2
                        + nfems
                        + grandmaingrp
                        + daughters11
                        + replace_na(matsisters11, 0)
                        , data = temp, family = binomial, na.action = na.fail)

tempallkinmod = glm(survto1 ~
                      sex2
                    + nfems
                    + femkinpres
                    , data = temp, family = binomial, na.action = na.fail)

temptotkinmod = glm(survto1 ~
                      sex2
                    + nfems
                    + totfemkin
                    , data = temp, family = binomial, na.action = na.fail)


tempkintypemod_n = glm(survto1 ~
                         sex2
                       + nfems
                       + gmapres_next
                       + daughtpres_next
                       + matsispres_next
                       , data = temp, family = binomial, na.action = na.fail)

tempallkinmod_n = glm(survto1 ~
                        sex2
                      + nfems
                      + femkinpres_next
                      , data = temp, family = binomial, na.action = na.fail)

temptotkinmod_n = glm(survto1 ~
                        sex2
                      + nfems
                      + totfemkin_next
                      , data = temp, family = binomial, na.action = na.fail)

AICc(tempsocmod, 
     tempkintypemod, tempallkinmod, tempkintypetotmod, temptotkinmod, 
     tempkintypemod_n, tempallkinmod_n, temptotkinmod_n) %>% 
  arrange(AICc)

summary(tempallkinmod)
summary(tempallkinmod_n)

AICc(tempsocmod, 
     tempkintypemod, tempallkinmod, 
     tempkintypemod_n, tempallkinmod_n) %>% 
  arrange(AICc) %>% 
  mutate(delta = AICc - min(AICc))

# -------------------- path analysis survival to age 1 --------------------

head(infsurv_soc_1)

infsurv_soc_1 = 
  infsurv_soc_1 %>% 
  mutate(matkin = daughters11 + replace_na(matsisters11, 0) + grandmaingrp > 0) 
  # can change from binary to count variable, results don't change

mod_path <-'
csi_f_ratio ~ matkin + mom_age_z + mom_age_z2 + mom_rank11_z
mom_rank11_z ~ matkin + mom_age_z
survto1 ~ csi_f_ratio + matkin + mom_age_z + mom_age_z2 + mom_rank11_z
'

# run model: 
# (note that ordered = "sire" specifies that sire is a binary term)
fit_path <- sem(mod_path, data = infsurv_soc_1, ordered = "survto1")  

summary(fit_path, fit.measures = T, standardized = T, rsquare = T)

parameterestimates(fit_path, standardized = F)

labels = c(csi_f_ratio = "CSI-F", matkin = "Kin presence", mom_age_z = "Mat. age", 
           mom_age_z2 = "Mat. age^2", mom_rank11_z = "Mat. Elo", survto1 = "Survive to 1")

lavaanPlot(model = fit_path, coefs = T, sig = .05, labels = labels, 
                node_options = list( fontname = "Arial"), 
                edge_options = list(color = "grey", fontname = "Arial"))

# -------------------- CSI-F & CSI-M pre and post birth ----------------------------

# only offspring that survived, and for whom we have data in the year after birth
# (note we lack maternal association data for the entire year after birth for 
# a few offspring because they were born too recently)
temp = 
  infsurv_soc_1 %>% 
  filter(survto1 == T, !is.na(csi_f_ratio_nextyr)) 

# CSI-F pre and post, pearson correlation:
cor.test(temp$csi_f_ratio, temp$csi_f_ratio_nextyr)

# regression model:
mtemp = 
  temp %>% 
  mutate(femkin = ((daughters11_nextyr + 
                      replace_na(matsisters11_nextyr, 0) + 
                      grandmaingrp_nextyr) > 0) * 1) %>% 
  filter(survto1 == T) %>% 
  lmer(csi_f_ratio_nextyr ~ csi_f_ratio + femkin + (1|momid), data = .) 

summary(mtemp)
tab_model(mtemp, collapse.ci = T)

# plot: 
temp %>% 
  mutate(femkin = (daughters11 + replace_na(matsisters11, 0) + grandmaingrp > 0)*1, 
         femkin_next = (daughters11_nextyr + replace_na(matsisters11_nextyr, 0) + grandmaingrp_nextyr > 0)*1 #,
         # femkin_diff = femkin != femkin_next,
         # logcsif_next = log(csi_f_ratio_nextyr), 
         # logcsif = log(csi_f_ratio)
         ) %>% 
  ggplot(aes(x = csi_f_ratio, y = csi_f_ratio_nextyr)) + #   logcsif  logcsif_next
  geom_point(aes(color = femkin_next > 0)) + 
  geom_smooth(method = "lm") +
  xlab("CSI-F pre-birth") + ylab("CSI-F post-birth") + 
  ggtitle("CSI with females") +
  ggeasy::easy_add_legend_title("Fem. kin present") + 
  ggthemes::scale_color_colorblind()

# CSI-M pre and post, pearson correlation: 
cor.test(temp$csi_m_ratio, temp$csi_m_ratio_nextyr)

# regression model:
mtemp = 
  temp %>% 
  mutate(malekin = ((sons13 + replace_na(matbros13, 0)) > 0) * 1) %>% 
  lmer(csi_m_ratio_nextyr ~ csi_m_ratio + malekin  + (1|momid), data = .) # + mom_age_z + mom_age_z2 + mom_rank11_z

summary(mtemp)
tab_model(mtemp, collapse.ci = T)

# plot:
temp %>% 
  mutate(malekin = ((sons13 + replace_na(matbros13, 0)) > 0) * 1 #, 
         # logcsim_next = log(csi_m_ratio_nextyr), 
         # logcsim = log(csi_m_ratio)
         ) %>% 
  # filter(survto1 == T) %>% 
  ggplot(aes(x = csi_m_ratio, y = csi_m_ratio_nextyr)) + #   logcsif  logcsif_next
  geom_point(aes(color = malekin > 0)) + 
  geom_smooth(method = "lm") +
  xlab("CSI-M pre-birth") + ylab("CSI-M post-birth") + 
  ggtitle("CSI with males") +
  ggeasy::easy_add_legend_title("Male kin present") + 
  ggthemes::scale_color_colorblind()

# --------------------- dyad-level analysis -------------------------------

# dataset of each mother's top bond in the year leading up to giving birth:
readRDS("dyadic_data.rds")

# dyadic sociality index vs maternal characteristics and kinship status with the other female
dyad_mod = lmer(dsi_time_ratio ~ matkin + momelo + mom_age_z + mom_age_z2 + (1|id1), data = topmomdyads)
summary(dyad_mod)

tab_model(dyad_mod, collapse.ci = T)

# since there's more variance in DSI values (because dyadic bond strengths vary much more than do 
# individual social integration measures) we might want to Z transform the DSIs instead of standardizing
# by dividing by the mean.  Same results if you do it that way:
dyad_mod2 = lmer(dsi_time_z ~ matkin + momelo + mom_age_z + mom_age_z2 + (1|id1), data = topmomdyads)
summary(dyad_mod2)


# ===================== offspring survival to age 5 ========================

infsurv_soc_5 = readRDS(file = "infsurv_5_anon.rds")

table(infsurv_soc_5$survto5) # 38 died, 59 survived

bin_fit_soc = 
  glm(survto5 ~ 
        sex2 +
        firstborn2 +
        mom_rank11_z +
        mom_age_z + 
        mom_imm_status2 + 
        mom_age_z2 +
        mom_siv_status2 +
        nfems + 
        nmales
      , data = infsurv_soc_5, family = binomial, na.action = na.fail)

performance::check_collinearity(bin_fit_soc)

simout = DHARMa::simulateResiduals(bin_fit_soc)
plot(simout) # some issues with fit but we'll check the best fit model later on

summary(bin_fit_soc)

allm_surv_soc_5 = dredge(bin_fit_soc, trace = 2, subset = !(mom_age_z2 && !(mom_age_z))) 

allm_surv_soc_5 = subset(allm_surv_soc_5, !nested(.))

bestm_surv_soc_5 = get.models(allm_surv_soc_5, subset = delta < 6)

model.sel(bestm_surv_soc_5)

plot_model(bestm_surv_soc_5$`257`, vline.color = "darkgray")

null_fit_soc_5 = 
  glm(survto5 ~
        sex2 
      , data = infsurv_soc_5, family = binomial, na.action = na.fail)


fit_soc_5_csi_f_ratio = update(null_fit_soc_5, . ~ . + csi_f_ratio)
fit_soc_5_csi_m_ratio = update(null_fit_soc_5, . ~ . + csi_m_ratio)
fit_soc_5_csi_f_m_ratio = update(null_fit_soc_5, . ~ . + csi_f_ratio + csi_m_ratio)
fit_soc_5_greg = update(null_fit_soc_5, . ~ . + time_obs_overall_ratio)

AICc(null_fit_soc_5, fit_soc_5_greg, fit_soc_5_csi_f_ratio, fit_soc_5_csi_m_ratio, fit_soc_5_csi_f_m_ratio) %>%
  mutate(delta = AICc - AICc(null_fit_soc_5),
         weight = MuMIn::Weights(AICc)) %>% 
  arrange(AICc) %>%
  round(3)

simout = DHARMa::simulateResiduals(fit_soc_5_csi_f_ratio)
plot(simout) # OK no problem with fit with best-fitting model

plot_model(fit_soc_5_csi_f_ratio, type = "pred", terms = "csi_f_ratio [all]", show.data = T)


# CSI-F vs its component parts (survival to 5):
fit_soc_5_time_f_ratio = update(null_fit_soc_5, . ~ . + time_obs_f_ratio)
fit_soc_5_grmrate_f_ratio = update(null_fit_soc_5, . ~ . + grmrate_f_ratio)

AICc(fit_soc_5_csi_f_ratio, fit_soc_5_time_f_ratio, fit_soc_5_grmrate_f_ratio, null_fit_soc_5) %>%
  mutate(delta = AICc - AICc(null_fit_soc_5),
         weight = MuMIn::Weights(AICc)) %>% 
  arrange(AICc) %>%
  round(3) # once again model with CSI-F fits better than its component parts

# ----------------- path analysis survival to age 5 ----------------------------------

head(infsurv_soc_5)

infsurv_soc_5 = 
  infsurv_soc_5 %>% 
  mutate(matkin = daughters11 + replace_na(matsisters11, 0) + grandmaingrp > 0)   

mod_path2 <-'
csi_f_ratio ~ matkin + mom_age_z + mom_age_z2 + mom_rank11_z
mom_rank11_z ~ matkin + mom_age_z
survto5 ~ csi_f_ratio + matkin + mom_age_z + mom_age_z2 + mom_rank11_z
'

# run model: 
# (note that ordered = "sire" specifies that sire is a binary term)
fit_path2 <- sem(mod_path2, data = infsurv_soc_5, ordered = "survto5")  

summary(fit_path2, fit.measures = T, standardized = T, rsquare = T)

labels2 = c(csi_f_ratio = "CSI-F", matkin = "Kin presence", mom_age_z = "Mat. age", 
            mom_age_z2 = "Mat. age^2", mom_rank11_z = "Mat. Elo", survto5 = "Survive to 5")

lavaanPlot(model = fit_path2, coefs = T, sig = .05, labels = labels2, 
           node_options = list( fontname = "Arial"), 
           edge_options = list(color = "grey", fontname = "Arial"))
