options(progressr.enable = TRUE)

source("./Project4/Code/model_selection.R")
`%>%` <- magrittr::`%>%`

# define simulation parameters
params = list(
    replications = 10000,
    n = c(250, 500),
    p = 20,
    beta = c((1:5)/6, rep(0, 15)),
    rho = c(0, 0.35, 0.7),
    alpha = 0.5,
    p_val = 0.05,
    nvar = length(beta),
    family = "gaussian",
    correlation_structure = "exchangeable",
    seed = 4130,
    version = getRversion()
)
params$nvar <- length(params$beta)

# define simulation procedure
run_simulation <- function(n, rho, params) {
    simulated_data <- hdrm::gen_data(
        n = n,
        rho = rho,
        p = params$nvar,
        beta = params$beta,
        family = params$family,
        corr = params$correlation_structure
    )
    
    reduced_models <- get_reduced_models(
        simulated_data = simulated_data,
        alpha = params$alpha,
        p_val = params$p_val
    )
    
    (
        purrr::map(reduced_models, ~broom::tidy(., conf.int = TRUE))
        %>% dplyr::bind_rows(.id = "selection_method")
        %>% dplyr::select(-c(std.error, statistic))
    )
}

# run simulation
set.seed(params$seed)

old_plan <- future::plan(future::multisession, workers = parallel::detectCores() - 1)
on.exit(future::plan(old_plan), add = TRUE)

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

tictoc::tic()
progressr::with_progress({
    grid <- tidyr::expand_grid(
        n = params$n,
        rho = params$rho,
        replication = 1:params$replications
    )
    
    p <- progressr::progressor(along = seq_len(nrow(grid)))
    
    simulation_results <- (
        grid
        %>% dplyr::mutate(
            result = furrr::future_map2(
                n, rho,
                \(n, rho) {
                    out <- run_simulation(n, rho, params = params)
                    p()
                    out
                },
                .options = furrr::furrr_options(seed = TRUE)
            )
        )
        %>% tidyr::unnest(result)
        %>% janitor::clean_names()
    )
})
runtime <- tictoc::toc()
params$runtime <- runtime$callback_msg

# save simulation results
saveRDS(params, file = "./Project4/Data/simulation_parameters.rds")
readr::write_csv(simulation_results, "./Project4/Data/simulation_results.csv")
capture.output(sessionInfo(), file = "./Project4/Code/session_info_simulation.txt")