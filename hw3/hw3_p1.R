library(tidyverse)
library(R.matlab)
library(broom)

df <- readMat('/Users/saragong/ec2120/hw3/twins.mat') %>% 
  as.data.frame()

# 1 -----------------
# Providing a consistent estimate for beta using the estimated covariance matrix.
# Note, the relevant reference is Note 5 pg 10.

n <- nrow(df)

# First I am estimating pi=[pi_1, pi_2], which are the coefficients in the unrestricted model.
# Then we can estimate the matrix Sigma=E[U_i U_i'].
pi_1_hat <- lm(lhrwage1 ~ educ1 + educ2, data=df) %>% coef()
pi_2_hat <- lm(lhrwage2 ~ educ1 + educ2, data=df) %>% coef()

# Slap it together into a 6 by 1 vector
pi_hat <- c(pi_1_hat, pi_2_hat) %>% matrix()

# Now make the Sigma_hat matrix which is a consistent estimator of Sigma by the WLLN
# Remember Sigma=E[U_i U_i'] so we are going to use 1/n sum_{i=1}^n U_i_hat U_i_hat'
# to make Sigma_hat.

# Similarly to estimating the asymptotic covariance matrix in pset 1, 
# I am going to create a list of 419 matrices, where each matrix is 2 x 2,
#  then sum the matrices elementwise and divide by n

# Note Sigma_hat is a 2 by 2 matrix.
Sigma_hat <- map(
  1:n, function(i) {
    
    # The outcome vector (2 by 1)
    Y_i <- df %>% select(lhrwage1, lhrwage2) %>% slice(i) %>% as.matrix() %>% t()
    
    # The predictor vector (intercept and two predictors, so 3 by 1)
    X_i <- df %>% 
      select(lhrwage1, lhrwage2) %>% 
      slice(i) %>% 
      as.numeric() %>%
      c(1, .) %>% # a 1 for the intercept term
      matrix()
    
    I <- diag(2) # 2 by 2 identity matrix
    
    # Use the formula to get the 2 x 1 vector of errors (one error for each twin)
    U_i_hat <- Y_i - kronecker(I, t(X_i)) %*% pi_hat
    
    return(U_i_hat %*% t(U_i_hat))
  }
)  %>%
  Reduce(`+`, .)/n # sum the matrices elementwise

# the weighting matrix Phi is the inverse of Sigma
# we construct Phi_hat as the inverse of Sigma_hat
# and it is consistent for Phi by Slutsky
Phi_hat <- solve(Sigma_hat)  

# Finally we construct beta_hat using the formula in Note 5 Section 5 pg 9,
# which by Claim 3 is a consistent estimator for beta

# This is the lefthand object in the formula (backticks to make the name compatible with R)
`R'PhiR_hat` <- map(
  1:n, function(i) {
    
    # R_i is a 2 by 2 matrix
    R_i <- df %>% 
      select(lhrwage1, lhrwage2) %>% 
      slice(i) %>% 
      as.matrix() %>% 
      t() %>%
      cbind(1, .) # add a column of 1 for intercept
    
    return(t(R_i) %*% Phi_hat %*% R_i)
  }
)  %>%
  Reduce(`+`, .)/n # sum the matrices elementwise

# This is the righthand object in the formula
`R'PhiY_hat` <- map(
  1:n, function(i) {
    
    # R_i is a 2 by 2 matrix
    R_i <- df %>% 
      select(lhrwage1, lhrwage2) %>% 
      slice(i) %>% 
      as.matrix() %>% 
      t() %>%
      cbind(1, .) # add a column of 1 for intercept
    
    # The outcome vector (2 by 1)
    Y_i <- df %>% select(lhrwage1, lhrwage2) %>% slice(i) %>% as.matrix() %>% t()
    
    return(t(R_i) %*% Phi_hat %*% Y_i)
  }
)  %>%
  Reduce(`+`, .)/n # sum the matrices elementwise

beta_hat <- solve(`R'PhiR_hat`) %*% `R'PhiY_hat`

cat('Consistent estimate for beta:')
beta_hat

# 2 -----------------

# The GLS predictor requires that:
# Y_1_tilde = beta1 + beta2*Z_1
# Y_2_tilde = beta1 + beta2*Z_2
# So, if X'=[1, Z_1, Z_2], then
# we need to find A_t such that Y_t = X' A_t beta is the model.
# In other words, using the notation on pg. 6, we need
# 1 = R_11 = X' a_11
# Z_1 = R_12 = X' a_12
# 1 = R_21 = X' a_21
# Z_2 = R_22 = X' a_22

# Following the notation of the notes, in this scenario q=3, M=2, and K=2.

a_11 <- matrix(c(1,0,0))
a_12 <- matrix(c(0,1,0))
a_21 <- matrix(c(1,0,0))
a_22 <- matrix(c(0,0,1))

A_1 <- cbind(a_11, a_12)
A_2 <- cbind(a_21, a_22)

# 6 x 2 matrix
A <- rbind(A_1, A_2)

# Now we plug this into that equation given in the homework question.

# Also following the notation of the notes, the choice variable is c, where
# c is a 2 by 1 vector


x <- df %>% select(educ1, educ2) %>% cbind(1, .) %>% as.matrix()

objective_function <- function(c, # c is the thing to optimize over
                               x, pi_hat, A, Phi_hat, n
) {
  
  # This formula is gross. but it's just copy-paste from the homework text.
  return(t((pi_hat - A %*% c)) %*% kronecker(Phi_hat, (t(x) %*% x)/n) %*% (pi_hat - A %*% c))
  
}

optim(
  par=c(0,0),
  fn=objective_function,
  x=x,
  pi_hat=pi_hat,
  A=A,
  Phi_hat=Phi_hat,
  n=n,
  method='BFGS'
)$par

# Just as a sanity check, try an off the shelf package
pivoted_df <- df %>%
  mutate(family=factor(row_number())) %>%
  pivot_longer(cols=-family) %>%
  mutate(
    twin = str_extract(name, "\\d+")  ,
    name = str_extract(name, "[a-zA-Z]+"), 
    ) %>%
  pivot_wider()

cat('Beta hat from off the shelf package lfe::felm')
lfe::felm(lhrwage ~ educ | family, data=pivoted_df) %>% tidy()

