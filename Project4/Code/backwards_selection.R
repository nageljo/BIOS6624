.create_formula <- function(lhs = "y", rhs = NULL) {
    if (is.null(rhs) || length(rhs) == 0) {
        formula_str <- paste(lhs, "~ 1")
    } else {
        formula_str <- paste(lhs, "~", paste(rhs, collapse = " + "))
    }
    return(as.formula(formula_str))
}

.get_nonzero_coef_names <- function(coef) {
    coef_names <- rownames(coef)[as.numeric(coef) != 0]
    return(coef_names[coef_names != "(Intercept)"])
}

.get_reduced_models_cv_glmnet <- function(simulated_data, alpha) {
    cv_fit <- glmnet::cv.glmnet(x = simulated_data$X, y = simulated_data$y, alpha = alpha)
    
    coef_lambda_min <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")
    coef_lambda_1se <- glmnet::coef.glmnet(cv_fit, s = "lambda.1se")
    
    names_lambda_min <- .get_nonzero_coef_names(coef_lambda_min)
    names_lambda_1se <- .get_nonzero_coef_names(coef_lambda_1se)
    
    df <- data.frame(y = simulated_data$y, simulated_data$X)
    reduced_model_min <- lm(.create_formula(rhs = names_lambda_min), data = df)
    reduced_model_1se <- lm(.create_formula(rhs = names_lambda_1se), data = df)
    
    return(list(reduced_model_min = reduced_model_min, reduced_model_1se = reduced_model_1se))
}

.get_reduced_model_p <- function(full_model, simulated_data, p_val) {
    reduced_model_p_ols <- olsrr::ols_step_backward_p(full_model, p_val = p_val)
    
    reduced_formula <- formula(reduced_model_p_ols$model)
    df <- data.frame(y = simulated_data$y, simulated_data$X)
    
    return(lm(formula = reduced_formula, data = df))
}

get_reduced_models <- function(simulated_data, alpha, p_val) {
    n_obs <- length(simulated_data$y)
    full_model <- lm(y ~ ., data = data.frame(y = simulated_data$y, simulated_data$X))
    
    lasso = .get_reduced_models_cv_glmnet(simulated_data, alpha = 1)
    elanet = .get_reduced_models_cv_glmnet(simulated_data, alpha = alpha)
    
    reduced_models = list(
        p = .get_reduced_model_p(full_model, simulated_data, p_val = p_val),
        aic = step(full_model, direction = "backward", trace = 0, k = 2),
        bic = step(full_model, direction = "backward", trace = 0, k = log(n_obs)),
        lasso_min = lasso$reduced_model_min,
        lasso_1se = lasso$reduced_model_1se,
        elanet_min = elanet$reduced_model_min,
        elanet_1se = elanet$reduced_model_1se
    )
    
    return(reduced_models)
}