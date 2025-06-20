# =======================
# LOAD DATA
# =======================

# SynNB: Synonymous rates data
SynNB <- read.table("Syn_rate_NB_v1", sep = " ", header = TRUE)

# Syn: Synonymous de novo counts (no header)
Syn <- read.table("syn_DD", sep = " ", header = FALSE)

# AM_rec: AlphaMissense recurrence data (skip first row)
AM_rec <- read.table("AM_rec", sep = " ", header = FALSE, skip = 1)

# =======================
# FIT POISSON MODEL
# =======================

# Fit Poisson GLM: V2 ~ log(V1) + offset(log(V3))
poisson_model <- glm(V2 ~ log(V1) + offset(log(V3)), 
                     family = poisson, 
                     data = Syn)

# Extract model coefficients
B0 <- coef(poisson_model)[1]  # Intercept
B1 <- coef(poisson_model)[2]  # Slope for log(V1)

# Compute expected rates (U)
Syn$U <- exp(B0 + B1 * log(Syn$V1))

# Combine expected rates, observed rate per V3, and original data
data_combined <- cbind(U = Syn$U, Roul = Syn$V2 / Syn$V3, Syn)

# =======================
# PREPARE FOR MERGING
# =======================

# Rename columns for merging
names(data_combined)[c(1, 3)] <- c('U', 'Roul')
names(AM_rec)[5] <- 'Roul'

# Merge with SynNB and AM_rec
SynNB <- merge(SynNB, data_combined[c("U", "Roul")])
AM_rec <- merge(AM_rec, data_combined[c("U", "Roul")])

# =======================
# LOG-LIKELIHOOD CALCULATIONS
# =======================

# Poisson log-likelihood
poisson_loglik <- sum(SynNB$mut * dpois(SynNB$rec, 
                                        lambda = SynNB$U, 
                                        log = TRUE))

# Negative binomial log-likelihood function with fixed k
neg_binom_loglik_fixed <- function(k, rec, mut, U) {
  r <- U / k
  mu <- U
  ll <- sum(mut * (lgamma(rec + r) - lgamma(r) - lgamma(rec + 1) +
                     r * log(r / (r + mu)) +
                     rec * log(mu / (r + mu))))
  return(-ll)  # Return negative for minimization
}

# =======================
# FIT NEGATIVE BINOMIAL MODEL
# =======================

# Initial value for k
k_init <- 1e-5

# Optimize to find best k
fit_fixed_k <- optim(par = k_init, 
                     fn = neg_binom_loglik_fixed, 
                     rec = SynNB$rec, 
                     mut = SynNB$mut, 
                     U = SynNB$U, 
                     method = "L-BFGS-B", 
                     lower = 1e-6)

# Extract k estimate and log-likelihood
k_est <- fit_fixed_k$par
nb_loglik_fixed <- -neg_binom_loglik_fixed(k_est, rec = SynNB$rec, mut = SynNB$mut, U = SynNB$U)

# =======================
# LIKELIHOOD RATIO TEST
# =======================

# Calculate LR statistic and p-value
lr_stat_fixed_poisson <- 2 * (nb_loglik_fixed - poisson_loglik)
p_value_fixed_poisson <- pchisq(lr_stat_fixed_poisson, df = 1, lower.tail = FALSE)

# Output results
cat("Fixed k LR Statistic (NB vs Poisson):", lr_stat_fixed_poisson, "\n")
cat("P-Value:", p_value_fixed_poisson, "\n")

# =======================
# P-VALUES FOR AM_rec
# =======================

# Compute p-values for AlphaMissense recurrence
AM_rec$p_value <- 1 - pnbinom(q = AM_rec$V6 - 1, 
                              size = 100 * AM_rec$U / k_est, 
                              mu = 100 * AM_rec$U)

# Adjust p-values using Benjamini-Hochberg correction
AM_rec$adjusted_pNB <- p.adjust(AM_rec$p_value * 49686008 / length(AM_rec$V1), 
                                method = "BH")