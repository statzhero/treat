#' Extension for "Research on Corporate Transparency"
#' Author: Ulrich Atz
#' 
#' The Bayesian/stan code is adapted from 
#' Breuer, M., & Schütt, H. H. (2019). 
#' Accounting for Uncertainty: An Application of Bayesian 
#' Methods to Accruals Models 
#' (SSRN Scholarly Paper ID 3417406). 
#' https://doi.org/10.2139/ssrn.3417406
#' https://github.com/hschuett/AccForUncertaintyCode/
#' ---------------------------------------------------------

library(tidyverse)
library(lme4)

# library(haven)
# library(rstan)
# options(mc.cores = parallel::detectCores())


# Functions
winsorize <- function(df, drop = NULL, ...) {
  if(is.null(drop)) ExPanDaR::treat_outliers(df, ...)
  else {
    vars <- !(names(df) %in% drop)
    ret <- df
    ret[vars] <- ExPanDaR::treat_outliers(ret[vars], ...)
    return(ret)
  }
}


# Load data
# Check ta !
source("code/R/theme_trr.R")
smp <- readRDS("data/generated/acc_sample.rds") 
smp_da <- smp %>%
  select(
    gvkey, fyear, ff12_ind, mj_da, dd_da, ln_ta, ln_mktcap, mtb,  
    ebit_avgta, sales_growth
  ) 
smp_da <- smp_da[is.finite(rowSums(smp_da %>% select(-gvkey, -ff12_ind))),]
# smp_da <- ExPanDaR::treat_outliers(smp_da, by = "fyear")

# Prep data
us_base_sample <- readRDS("data/generated/us_base_sample.rds")

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
  winsorize(drop = "fyear", truncate = FALSE) %>%
  fill(ff48_ind, .direction = "downup") %>%
  tidylog::drop_na(ff48_ind) %>%
  # ungroup %>%
  #   tjmisc::sample_n_of(100, gvkey) %>%
  group_by(fyear, ff48_ind) %>%
  mutate(IndYearID = interaction(ff48_ind, fyear),
         IndYearID = as.numeric(IndYearID))


# Hierarchical Jones model -----------------
#' ---------------------------------------------------------

# Model----------------------
mjh <- lmer(tacc ~ -1 + inverse_a + drev + ppe +
              (inverse_a + drev + ppe | IndYearID), 
            dta)

system.time(
mjh2 <- lmer(tacc ~ -1 + inverse_a + drev + ppe +
              (inverse_a + drev + ppe | gvkey), 
            dta)
)

# residuals(mjh2) %>% qqnorm
# residuals(mjh2) %>% qqline
# lattice::qqmath(mjh2)

# mjh2 <- lmer(tacc ~ -1 + inverse_a + drev + ppe +
#                (inverse_a + drev + ppe | ff48_ind) +
#                (inverse_a + drev + ppe | fyear), 
#              dta)

mjh_da <- dta %>% 
  ungroup %>% 
  mutate(tacc_hat = predict(mjh),
         mjh_da = tacc - tacc_hat) %>% 
  mutate(tacc_hat2 = predict(mjh2),
         mjh_da2 = tacc - tacc_hat2) %>% 
  select(gvkey, fyear, mjh_da, mjh_da2, 
         tacc, tacc_hat, tacc_hat2)


# Correlations----------------------
mjh_da %>% 
  left_join(mj, by = c("gvkey", "fyear")) %>%
  select(tacc, mj_da) %>% 
  mutate(ta_hat = tacc - mj_da) %>% 
  tidylog::drop_na() %>% 
  cor

cor(mjh_da$tacc, mjh_da$tacc_hat)
cor(mjh_da$tacc, mjh_da$mjh_da)

cor(mjh_da$tacc, mjh_da$tacc_hat2)
cor(mjh_da$tacc, mjh_da$mjh_da2)


# Figures----------------------
smp_ext <- smp_da %>% 
  left_join(mjh_da, by = c("gvkey", "fyear"))


fig_scatter_md_mjh <- ggplot(smp_ext, aes(x = mj_da, y = mjh_da)) +
  geom_bin2d(
    aes(fill = stat(log(count))),  bins = c(100, 100)
  ) +
  labs(x = "Modified Jones DA", y = "Hierarchical Jones DA (industry-year)") +
  scale_fill_trr266_c() + 
  theme_trr(axis_y_horizontal = FALSE)

fig_scatter_md_mjh2 <- ggplot(smp_ext, aes(x = mj_da, y = mjh_da2)) +
  geom_bin2d(
    aes(fill = stat(log(count))),  bins = c(100, 100)
  ) +
  labs(x = "Modified Jones DA", y = "Hierarchical Jones DA (firm-level)") +
  scale_fill_trr266_c() + 
  theme_trr(axis_y_horizontal = FALSE)

graph_boxplot <- function(df) {
  df <- df %>% select(fyear, mj_da, mjh_da) %>%
    pivot_longer(any_of(c("mj_da", "mjh_da")), names_to = "type", values_to = "da")
  
  ggplot(
    df, 
    aes(x = fyear, y = da, group = interaction(type, fyear), color = type)
  ) +
    geom_boxplot() +
    labs(x = "Fiscal year", y = NULL, color ="Type of discretionary accruals") +
    scale_color_trr266_d(labels = c("Modified Jones", "Hierarchical Jones (industry-year)")) + 
    theme_trr(legend = TRUE)
}

(fig_boxplot_hier <- graph_boxplot(smp_ext))


graph_boxplot2 <- function(df) {
  df <- df %>% select(fyear, mj_da, mjh_da2) %>%
    pivot_longer(any_of(c("mj_da", "mjh_da2")), names_to = "type", values_to = "da")
  
  ggplot(
    df, 
    aes(x = fyear, y = da, group = interaction(type, fyear), color = type)
  ) +
    geom_boxplot() +
    labs(x = "Fiscal year", y = NULL, color ="Type of discretionary accruals") +
    scale_color_trr266_d(labels = c("Modified Jones", "Hierarchical Jones (firm-level)")) + 
    theme_trr(legend = TRUE)
}

(fig_boxplot_hier2 <- graph_boxplot2(smp_ext))


graph_coef <- function(d, lab){
   ggplot() + 
    geom_boxplot(aes(x = d)) +
    labs(x = lab) + 
    scale_color_trr266_d() + 
    theme_trr(legend = TRUE)
}

# Industry-year
coef_mj <- list()
coef_mj[[1]] <- graph_coef(mj$mj_inverse_a, "Modified Jones INV_A coefficients")
coef_mj[[2]] <- graph_coef(mj$mj_drev, "Modified Jones DREV coefficients")
coef_mj[[3]] <- graph_coef(mj$mj_ppe, "Modified Jones PPE coefficients")

coef_mj[[4]] <- graph_coef(coef(mjh)$IndYearID$inverse_a, "Hier. mod. Jones INV_A coefficients")
coef_mj[[5]] <- graph_coef(coef(mjh)$IndYearID$drev, "Hier. mod. Jones DREV coefficients")
coef_mj[[6]] <- graph_coef(coef(mjh)$IndYearID$ppe, "Hier. mod. Jones PPE coefficients")


library(patchwork)
(fig_boxplot_coef <- coef_mj[[1]] + coef_mj[[4]] + coef_mj[[2]] +
                     coef_mj[[5]] + coef_mj[[3]] + coef_mj[[6]] +
  plot_layout(nrow = 3))

# Firm-level
coef_mj2 <- list()
coef_mj2[[1]] <- graph_coef(mj$mj_inverse_a, "Modified Jones INV_A coefficients")
coef_mj2[[2]] <- graph_coef(mj$mj_drev, "Modified Jones DREV coefficients")
coef_mj2[[3]] <- graph_coef(mj$mj_ppe, "Modified Jones PPE coefficients")

coef_mj2[[4]] <- graph_coef(coef(mjh2)$gvkey$inverse_a, "Hier. mod. Jones INV_A coefficients")
coef_mj2[[5]] <- graph_coef(coef(mjh2)$gvkey$drev, "Hier. mod. Jones DREV coefficients")
coef_mj2[[6]] <- graph_coef(coef(mjh2)$gvkey$ppe, "Hier. mod. Jones PPE coefficients")

(fig_boxplot_coef2 <- coef_mj2[[1]] + coef_mj2[[4]] + coef_mj2[[2]] +
    coef_mj2[[5]] + coef_mj2[[3]] + coef_mj2[[6]] +
    plot_layout(nrow = 3))



# Save outputs
save(
  list = c("smp_da", "var_names", ls(pattern = "fig_*"), ls(pattern = "tab_*")),
  file = "output/results.rda"
)


# NOT RUN

#' # Bayes Jones model -----------------
#' # code that fits the model and saves output. 
#' # Squeezed into one function:
#' run_model <- function(modelname,
#'                       filename,
#'                       parameters,
#'                       data) {
#'   # fit model
#'   fit <- stan(filename, data=input_data,
#'               iter=2000, warmup=1000, chains=2, seed=20210507,
#'               verbose=FALSE)
#'   
#'   # Save Results
#'   # bayes_predicted <- as.matrix(fit, par="y_fit_wo")
#'   # saveRDS(bayes_predicted, paste0("output/modelfits/", modelname, "-ta-coef-only-preds.rds"))
#' 
#'   bayes_ind_coef <- as.matrix(fit, par="b")[2000, ]
#'   write.csv(bayes_ind_coef, gzfile(paste0("output/modelfits/", modelname, "-coef.csv.gz")))
#' 
#'   bayes_results <- as.matrix(fit, par=parameters)
#'   write.csv(bayes_results, gzfile(paste0("output/modelfits/", modelname, "-results.csv.gz")))
#' 
#'   sum_table <- summary(fit, pars=parameters)$summary
#'   write.csv(sum_table, paste0("output/modelfits/", modelname, "-summary.csv"))
#'   
#'   for (para in parameters){
#'     tp1 <- traceplot(fit, pars=c(para))
#'     ggsave(filename = paste0("output/modelfits/", modelname, "-trace-", para, ".pdf"), plot = tp1)
#'   }
#'   
#'   print(fit, par=parameters)
#'   return(fit)
#' }
#' 
#' 
#' # Industry-year-varying coefficients model -----------------
#' #' ---------------------------------------------------------
#' # The Bayesian version allows the coefficients in the modified
#' # Jones model to vary across industry and year based on the 
#' # assumptions of the model. 

#' 
#' 
#' modelname <- "mjones_hierarchical"
#' filename <- "code/stan/mjones_hierarchical.stan"
#' model_vars <- c("inverse_a", "drev", "ppe")
#' parameters <- c("mu_b", "sigma", "L_Omega", "tau")
#' input_data = list(TA = dta$tacc,
#'                   x = dta[, model_vars],
#'                   N = nrow(dta),
#'                   J = max(dta$IndYearID),
#'                   K = length(model_vars),
#'                   IndYearID = dta$IndYearID)
#' fit <- run_model(modelname, filename, parameters, input_data)
