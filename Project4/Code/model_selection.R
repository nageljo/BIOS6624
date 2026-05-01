#' Construct a linear model formula with no intercept
#'
#' @param lhs Character string giving the left-hand side variable name.
#' @param rhs Character vector of right-hand side variable names.
#'
#' @return An object of class \code{formula}.
#'
#' @keywords internal
.create_formula <- function(lhs = "y", rhs = NULL) {
    if (is.null(rhs) || length(rhs) == 0) {
        formula_str <- paste(lhs, "~ 0")
    } else {
        formula_str <- paste(lhs, "~ 0 + ", paste(rhs, collapse = " + "))
    }
    as.formula(formula_str)
}

#' Extract names of nonzero coefficients
#'
#' @param coef A coefficient matrix as returned by \code{coef.glmnet()}.
#'
#' @return A character vector of coefficient names with nonzero values.
#'
#' @keywords internal
.get_nonzero_coef_names <- function(coef) {
    rownames(coef)[as.numeric(coef) != 0]
}

#' Fit reduced linear models using cross-validated glmnet
#'
#' Uses cross-validation to select predictors via lasso or elastic net,
#' and refits ordinary least squares models using the selected variables.
#'
#' @param simulated_data A list containing:
#' \describe{
#'   \item{X}{A numeric model matrix of predictors.}
#'   \item{y}{A numeric response vector.}
#' }
#' @param alpha Elastic net mixing parameter.
#'   \code{alpha = 1} corresponds to lasso.
#'
#' @return A list with components:
#' \describe{
#'   \item{reduced_model_min}{OLS model using predictors at \code{lambda.min}.}
#'   \item{reduced_model_1se}{OLS model using predictors at \code{lambda.1se}.}
#' }
#'
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @keywords internal
.get_reduced_models_cv_glmnet <- function(simulated_data, alpha) {
    cv_fit <- glmnet::cv.glmnet(
        x = simulated_data$X,
        y = simulated_data$y,
        alpha = alpha,
        intercept = FALSE
    )
    
    coef_lambda_min <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")
    coef_lambda_1se <- glmnet::coef.glmnet(cv_fit, s = "lambda.1se")
    
    names_lambda_min <- .get_nonzero_coef_names(coef_lambda_min)
    names_lambda_1se <- .get_nonzero_coef_names(coef_lambda_1se)
    
    df <- data.frame(y = simulated_data$y, simulated_data$X)
    
    reduced_model_min <- lm(.create_formula(rhs = names_lambda_min), data = df)
    reduced_model_1se <- lm(.create_formula(rhs = names_lambda_1se), data = df)
    
    list(
        reduced_model_min = reduced_model_min,
        reduced_model_1se = reduced_model_1se
    )
}

#' Backward elimination based on p-values
#'
#' Iteratively removes the predictor with the largest p-value until all
#' remaining predictors have p-values below a specified threshold.
#'
#' @param full_model A fitted \code{lm} object containing all predictors.
#' @param p_val Numeric p-value cutoff for retention.
#'
#' @return A reduced \code{lm} object.
#'
#' @keywords internal
.get_reduced_model_p <- function(full_model, p_val) {
    current_model <- full_model
    
    repeat {
        current_terms <- attr(terms(current_model), "term.labels")
        if (length(current_terms) == 0) break
        
        test_results <- drop1(current_model, test = "F")
        if (!"Pr(>F)" %in% colnames(test_results)) break
        
        p_values <- test_results[["Pr(>F)"]][-1]
        if (is.null(p_values) || length(p_values) == 0) break
        names(p_values) <- current_terms
        
        max_p <- max(p_values, na.rm = TRUE)
        if (is.na(max_p) || max_p <= p_val) break
        
        var_to_drop <- names(p_values)[which.max(p_values)]
        current_terms <- attr(terms(current_model), "term.labels")
        new_terms <- setdiff(current_terms, var_to_drop)
        
        df <- model.frame(current_model)
        
        current_model <- lm(
            .create_formula(rhs = new_terms),
            data = df
        )
    }
    
    current_model
}

#' Fit multiple reduced linear models using different selection criteria
#'
#' Fits reduced models using p-value backward elimination, AIC, BIC,
#' lasso, and elastic net (via cross-validated glmnet).
#'
#' @param simulated_data A list containing:
#' \describe{
#'   \item{X}{A numeric model matrix of predictors.}
#'   \item{y}{A numeric response vector.}
#' }
#' @param alpha Elastic net mixing parameter passed to \code{glmnet}.
#' @param p_val Numeric p-value cutoff for backward elimination.
#'
#' @return A named list of \code{lm} objects with elements:
#' \itemize{
#'   \item \code{p} – p-value based backward elimination
#'   \item \code{aic} – AIC-based stepwise selection
#'   \item \code{bic} – BIC-based stepwise selection
#'   \item \code{lasso_min}, \code{lasso_1se}
#'   \item \code{elanet_min}, \code{elanet_1se}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- X[, 1] - X[, 3] + rnorm(100)
#' models <- get_reduced_models(
#'   list(X = X, y = y),
#'   alpha = 0.5,
#'   p_val = 0.05
#' )
#' }
#'
#' @export
get_reduced_models <- function(simulated_data, alpha, p_val) {
    n_obs <- length(simulated_data$y)
    df <- data.frame(y = simulated_data$y, simulated_data$X)
    
    full_model <- lm(y ~ 0 + ., data = df)
    
    lasso <- .get_reduced_models_cv_glmnet(simulated_data, alpha = 1)
    elanet <- .get_reduced_models_cv_glmnet(simulated_data, alpha = alpha)
    
    list(
        p = .get_reduced_model_p(full_model, p_val = p_val),
        aic = step(full_model, direction = "backward", trace = 0, k = 2),
        bic = step(full_model, direction = "backward", trace = 0, k = log(n_obs)),
        lasso_min = lasso$reduced_model_min,
        lasso_1se = lasso$reduced_model_1se,
        elanet_min = elanet$reduced_model_min,
        elanet_1se = elanet$reduced_model_1se
    )
}