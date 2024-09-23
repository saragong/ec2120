library(tidyverse)
library(R.matlab)
library(broom)

df <- readMat('/Users/saragong/ec2120/hw2/nls.mat') %>% 
  as_tibble() %>%
  mutate(
    log_uwe = log(uwe),
    # use the formula given in the homework
    exper = age80 - educ - 6,
    exper_sq = exper^2
  )
  
# Q1.1 (a) ------------
model_saturated <- lm(log_uwe ~ educ + exper + exper_sq + iq, data=df)
model_restricted <- lm(log_uwe ~ educ + exper + exper_sq, data=df)

model_saturated
model_restricted

# Q1.1 (b) ------------
# Let
# log_uwe ~ beta_0 + beta_educ*educ + beta_exper*exper + beta_exper_sq*exper_sq + beta_iq*iq
# log_uwe ~ alpha_0 + alpha_educ*educ + alpha_exper*exper + alpha_exper_sq*exper_sq 
# iq ~ gamma_0 + gamma_educ*educ + gamma_exper*exper + gamma_exper_sq*exper_sq


beta_educ_hat <- model_saturated %>% tidy() %>% filter(term=='educ') %>% pull(estimate)
beta_iq_hat <- model_saturated %>% tidy() %>% filter(term=='iq') %>% pull(estimate)
alpha_educ_hat <- model_restricted %>% tidy() %>% filter(term=='educ') %>% pull(estimate)

model_aux <- lm(iq ~ educ + exper + exper_sq, data=df)
model_aux
gamma_educ_hat <- model_aux %>% tidy() %>% filter(term=='educ') %>% pull(estimate)

# We can get the aux model's coef on education from the saturated and restricted models
cat("gamma_educ_hat from OVB formula:", round((alpha_educ_hat-beta_educ_hat)/beta_iq_hat, 10))
cat("gamma_educ_hat from aux regression:", round(gamma_educ_hat, 10))

# Q1.2 ------------
df <- df %>%
  mutate(
    iq_residualized = resid(model_aux)
  )

model_fwl <- lm(log_uwe ~ iq_residualized, data=df)

beta_iq_hat_from_fwl <- model_fwl %>% tidy() %>% filter(term=='iq_residualized') %>% pull(estimate)

cat("beta_iq_hat from full model:", beta_iq_hat)
cat("beta_iq_hat from FWL residualized regression:", beta_iq_hat_from_fwl)
