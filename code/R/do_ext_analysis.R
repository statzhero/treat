#' Extension for "Research on Corporate Transparency"
#' Author: Ulrich Atz
#' 
#' The code is adapted from 
#' Breuer, M., & Schütt, H. H. (2019). 
#' Accounting for Uncertainty: An Application of Bayesian Methods to Accruals Models 
#' (SSRN Scholarly Paper ID 3417406). 
#' Social Science Research Network. 
#' https://doi.org/10.2139/ssrn.3417406
#' https://github.com/hschuett/AccForUncertaintyCode/
#' ---------------------------------------------------------

library(tidyverse)
library(haven)
library(rstan)
options(mc.cores = parallel::detectCores())

# Load data
us_base_sample <- readRDS("data/generated/us_base_sample.rds")

# code that fits the model and saves output. Squeezed into one function:
run_model <- function(modelname,
                      filename,
                      parameters,
                      data) {
  # fit model
  fit <- stan(filename, data=input_data,
              iter=4000, warmup=1000, chains=2, seed=123,
              verbose=FALSE)
  
  # Save Results
  bayes_predicted <- as.matrix(fit, par="y_fit")
  saveRDS(bayes_predicted, paste0("output/modelfits/", modelname, "-ta-predicions.rds"))
  rm(bayes_predicted)
  
  bayes_predicted <- as.matrix(fit, par="y_fit_wo")
  saveRDS(bayes_predicted, paste0("output/modelfits/", modelname, "-ta-coef-only-preds.rds"))
  rm(bayes_predicted)
  
  ind_loklik <- as.matrix(fit, par="log_lik")
  saveRDS(ind_loklik, paste0("output/modelfits/", modelname, "-loklik.rds"))
  rm(ind_loklik)
  
  bayes_ind_coef <- as.matrix(fit, par="b")
  write.csv(bayes_ind_coef, gzfile(paste0("output/modelfits/", modelname, "-coef.csv.gz")))
  rm(bayes_ind_coef)
  
  bayes_results <- as.matrix(fit, par=parameters)
  write.csv(bayes_results, gzfile(paste0("output/modelfits/", modelname, "-results.csv.gz")))
  rm(bayes_results)
  
  sum_table <- summary(fit, pars=parameters)$summary
  write.csv(sum_table, paste0("output/modelfits/", modelname, "-summary.csv"))
  
  for (para in parameters){
    tp1 <- traceplot(fit, pars=c(para))
    ggsave(filename = paste0("output/modelfits/", modelname, "-trace-", para, ".pdf"), plot = tp1)
  }
  
  print(fit, par=parameters)
  return(fit)
}


# Industry-year-varying coefficients model -------------------------------------
# The Bayesian version of the most common approach that fits a regression based on industry
# year subsamples 
# Rather than doing that we use a hierarchical model 
# To pool industry-year
# coefficients and pool accross coefficients

min_obs <- 10

dta <- us_base_sample %>% 
  group_by(gvkey) %>% 
  arrange(gvkey, fyear) %>% 
  mutate(
    lagta = dplyr::lag(at),
    tacc = (ibc - oancf)/lagta,
    drev = (sale - dplyr::lag(sale) + recch)/lagta,
    inverse_a = 1/lagta,
    ppe = ppegt/lagta
  ) %>%
  drop_na(tacc, drev, ppe) %>% 
  select(gvkey, ff48_ind, fyear, tacc, drev, inverse_a, ppe) %>%
  group_by(ff48_ind, fyear) %>%
  filter(n() >= min_obs) %>%
  winsorize(drop = "fyear") %>% 
  fill(ff48_ind, .direction = "downup") %>% 
  tidylog::drop_na(ff48_ind) %>% 
  # ungroup %>% 
  #   tjmisc::sample_n_of(500, gvkey) %>% 
  group_by(fyear, ff48_ind) %>% 
  mutate(IndYearID = row_number())


modelname <- "mjones_hierarchical"
filename <- "code/stan/mjones_hierarchical.stan"
model_vars <- c("inverse_a", "drev", "ppe")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = dta$tacc,
                  x = dta[, model_vars],
                  N = nrow(dta),
                  J = max(dta$IndYearID),
                  K = length(model_vars),
                  IndYearID = dta$IndYearID)
fit <- run_model(modelname, filename, parameters, input_data)
