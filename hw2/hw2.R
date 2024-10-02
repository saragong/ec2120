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

# Q1.3.a ------------
lm(log_uwe ~ educ + exper + exper_sq + iq + fed + med, data=df)
#... in progress...

# Q1.4 --------------
# Question did not specify a model, so am using saturated model
X <- model.matrix(model_saturated) %>% t() # 815 x 5 design matrix
U_hat <- residuals(model_saturated) %>% as.vector()
Y <- df$log_uwe
n <- dim(X)[2]

# fyi - can get the OLS coefficients using the closed form expression
beta_hat <- solve((X %*% t(X))/n) %*% (X %*% Y)/n 

cat('Beta hat from the closed form expression:')
beta_hat

# compare to the beta hat we get from using lm off-the-shelf, it's the exact same
cat('Beta hat from using "lm" off the shelf:')
coef(model_saturated)


# Define this here for answering the problem
alpha_hat <- solve((X %*% t(X))/n)

# Following the formula on Note 3 Sec 9, this is going to create a list of 815 matrices,
# where each matrix is 5 x 5, then sum them elementwise and divide by N
Sigma_hat <- map(1:n, function(i) { 
  
  # slightly granular way to write this, but illustrative
  
  X_i <- as.matrix(X[,i])
  U_hat_i <- U_hat[i]
  
  return(U_hat_i^2 * (X_i %*% t(X_i)))
  
  }) %>% 
  Reduce(`+`, .)/n # sum the matrices elementwise

# This is our estimate of the asymptotic covariance matrix
Lambda_hat <- alpha_hat %*% Sigma_hat %*% t(alpha_hat)

# We can also get the asymptotic distribution and SE for the educ coefficient
l <- as.matrix(c(0, 1, 0, 0, 0))

t(l) %*% beta_hat # this is just the singular coefficient we want

# Standard error
se <- ((t(l) %*% Lambda_hat %*% l) / n)^(1/2)

cat('Estimated SE of beta_educ_hat:', se)

# Compare this to the SE we get off-the-shelf from lm - this is close ish but not exact
tidy(model_saturated) %>% filter(term=='educ') %>% pull(std.error) %>% cat('Estimated SE produced by using "lm" off the shelf:', .)

# (beta_educ_hat - beta_educ)/se converges to N(0,1), so
# beta_educ_hat converges to N(beta_educ, se^2)

cat('Estimated asymptotic 95% confidence interval of beta_educ_hat:')
ci_low <- beta_educ_hat+ se * qnorm(0.05/2)
ci_high <- beta_educ_hat+ se *qnorm(1-(0.05/2))
glue::glue('[{round(ci_low, 4)}-{round(ci_high, 4)}]')

# Compare this to the confidence interval we get off-the-shelf from LM - this is close ish but not exact
cat('Estimated 95% confidence interval of beta_educ_hat from using "lm" off the shelf')
ci_low <- tidy(model_saturated, conf.int = T, conf.level = 0.95) %>% filter(term=='educ') %>% pull(conf.low)
ci_high <- tidy(model_saturated, conf.int = T, conf.level = 0.95) %>% filter(term=='educ') %>% pull(conf.high)
glue::glue('[{round(ci_low, 4)}-{round(ci_high, 4)}]')

