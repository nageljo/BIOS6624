######################
### Bayesian setup ###
######################
fit_brm <- function(formula, data, file = NULL) {
    # define priors and MCMC parameters
    priors <- c(
        set_prior("normal(0, 100)", class = "Intercept"),    # prior for intercept
        set_prior("normal(0, 100)", class = "b"),            # priors for fixed effects
        set_prior("normal(0, 100)", class = "sigma", lb = 0) # prior for sigma
    )
    chains <- 4
    iter <- 2000
    warmup <- 1000
    
    # create the required directories to store the model if they do not exist
    if (!is.null(file)) {
        out_dir <- dirname(file)
        if (!identical(out_dir, ".") && !dir.exists(out_dir)) {
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
    }
    
    # fit the model
    model <- brm(
        formula = formula,
        data = data,
        file = file,
        # prior = priors,
        chains = chains,
        iter = iter,
        warmup = warmup,
        file_refit = "on_change",
        family = gaussian(),
        seed = 4130,
        cores = 4,
        backend = "cmdstanr",
        refresh = 0,
        silent = 2
    )
    
    return(model)
}



####################################
### Helpers for Freq Bayes table ###
####################################

# Find the exact coefficient name for hard_drugs_baseline (handles e.g., 0/1 or factor coding)
find_term_name <- function(coef_names, var = "hard_drugs_baseline") {
    hits <- grep(paste0("^", var, "($|[^:])"), coef_names, value = TRUE)
    if (length(hits) == 0) {
        stop(sprintf("Could not find coefficient for '%s' in: %s",
                     var, paste(coef_names, collapse = ", ")))
    }
    if (length(hits) > 1) {
        message(sprintf("Multiple matches for '%s': using '%s'", var, hits[1]))
    }
    hits[1]
}

# Frequentist: estimate, 95% CI, p-value for hard_drugs_baseline
extract_lm_stats <- function(fit, var = "hard_drugs_baseline") {
    term <- find_term_name(names(coef(fit)), var)
    ci   <- confint(fit, level = 0.95)
    est  <- coef(fit)[term]
    lcl  <- ci[term, 1]
    ucl  <- ci[term, 2]
    pval <- summary(fit)$coefficients[term, "Pr(>|t|)"]
    list(estimate = est, lcl = lcl, ucl = ucl, p = pval)
}

# Bayesian: posterior mean & 95% HDI for hard_drugs_baseline
extract_brm_stats_hdi <- function(fit, var = "hard_drugs_baseline", cred = 0.95) {
    # Use posterior draws; column names are "b_<term>"
    draws <- as_draws_df(fit)
    # fixef names
    fx_names <- dimnames(brms::fixef(fit))[[1]]
    term     <- find_term_name(fx_names, var)
    col_name <- paste0("b_", term)
    
    if (!col_name %in% names(draws)) {
        stop(sprintf("Could not find draws column '%s' in posterior draws.", col_name))
    }
    
    est <- mean(draws[[col_name]])
    hdi_ci <- bayestestR::hdi(draws[[col_name]], ci = cred)
    lcl <- hdi_ci$CI_low
    ucl <- hdi_ci$CI_high
    list(estimate = est, lcl = lcl, ucl = ucl)
}

# Compute ΔLOOIC using prefit reduced model:
# ΔLOOIC = LOOIC(reduced) - LOOIC(full)
delta_loo_prefit <- function(full_fit, reduced_fit, reloo = TRUE) {
    # Add or recompute PSIS-LOO
    full_fit    <- brms::add_criterion(full_fit,    "loo", reloo = reloo)
    reduced_fit <- brms::add_criterion(reduced_fit, "loo", reloo = reloo)
    
    loo_full <- full_fit$criteria$loo
    loo_red  <- reduced_fit$criteria$loo
    
    # LOOIC = -2 * elpd_loo
    looic_full <- loo_full$estimates["looic","Estimate"]
    looic_red  <- loo_red$estimates["looic","Estimate"]
    
    delta    <- as.numeric(looic_red - looic_full)
    se_full  <- loo_full$estimates["looic","SE"]
    se_red   <- loo_red$estimates["looic","SE"]
    se_delta <- sqrt(se_full^2 + se_red^2)
    
    list(delta_looic = delta, se = se_delta)
}


# Helper to pretty format numbers
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)