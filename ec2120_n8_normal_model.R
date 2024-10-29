library(tidyverse)
options(scipen=999)

# Set seed for reproducibility
set.seed(123)

# Parameters
n <- 1000
k <- 4

beta_0 <- 1
beta_1 <- 2
beta_2 <- -0.5
beta_3 <- 1

beta <- c(beta_0, beta_1, beta_2, beta_3) %>% as.matrix()

sigma <- 1

# Generate the data
df <- tibble(
  x_1 = rnorm(n),
  x_2 = rnorm(n),
  x_3 = rnorm(n),
  v = rnorm(n, mean = 0, sd = sigma), # actually not observable, but include for illustration
  y = beta_0 + beta_1 * x_1 + beta_2 * x_2 + beta_3 * x_3 + v
)

model <- lm(y ~ x_1 + x_2 + x_3, data=df)

summary(model)

x <- model.matrix(model)
y <- df$y %>% as.matrix()
v <- df$v %>% as.matrix()

# Perform Singular Value Decomposition on the design matrix
svd_result <- svd(x, nu=n)

# Display the SVD components
q <- svd_result$u  # Left singular vectors
d <- diag(svd_result$d)  # Singular values
s <- svd_result$v  # Right singular vectors

dim(q)
dim(d)
dim(s)

# Confirm q is orthonormal
all.equal.numeric(q %*% t(q), diag(ncol(q)), tolerance = 1e-10) # should equal I_n up to rounding error
all.equal.numeric(s %*% t(s), diag(ncol(s)), tolerance = 1e-10) # should equal I_k up to rounding error

q1 <- q[1:n,1:k]

# Note that x = q1 %*% d %*% t(s) within rounding error

# these things are observable
y_star <- t(q) %*% y #
y_star_1 <- y_star[1:k] # = t(q1) %*% y
y_star_2 <- y_star[(k+1):n]

# these things are not observable
v_star <- t(q) %*% v
v_star_1 <- v_star[1:k]
v_star_2 <- v_star[(k+1):n]

# Solve for the beta and SSR
b <- solve(d %*% t(s)) %*% y_star_1

# The following three things should be equal exactly
as.numeric(coef(model))
as.numeric(b)
as.numeric(beta + sigma * (s %*% solve(d) %*% v_star_1))

# The following two things should be equal exactly
sum(residuals(model)^2)
sigma^2 * sum(v_star_2^2)




