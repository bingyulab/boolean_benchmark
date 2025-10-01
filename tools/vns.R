#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(optparse)
  library(CellNOptR)
  library(MEIGOR)  # MEIGO R package for VNS optimization
  library(here)
  library(jsonlite)  # For saving results in JSON format
})

# Source helper functions
source(here::here("tools", "functions.R"))

# Define command line options
option_list <- list(
  make_option(c("-d", "--dataset"), type="character", default="toy",
              help="Dataset name (toy, apoptosis, dream, TCell) [default: %default]"),
  make_option(c("-k", "--kfold"), type="integer", default=5,
              help="Number of k-fold cross validation [default: %default]"),
  make_option(c("-c", "--change_percent"), type="double", default=0.0,
              help="Perturbation percentage for network [default: %default]"),
  make_option(c("-e", "--maxeval"), type="integer", default=1000,
              help="Maximum number of function evaluations for VNS [default: %default]"),
  make_option(c("-t", "--maxtime"), type="integer", default=3600,
              help="Maximum time in seconds for VNS optimization [default: %default]"),
  make_option(c("-l", "--use_local"), type="integer", default=1,
              help="Use local search in VNS [default: %default]"),
  make_option(c("-a", "--aggr"), type="integer", default=0,
              help="Aggregation strategy for VNS (0=none, 1=mean, 2=median) [default: %default]"),
  make_option(c("-s", "--local_search_type"), type="integer", default=1,
              help="Local search type for VNS (1=Nelder-Mead, 2=BFGS, 3=CG) [default: %default]"),
  make_option(c("-D", "--decomp"), type="integer", default=0,
              help="Decomposition strategy for VNS [default: %default]"),
  make_option(c("-m", "--maxdist"), type="integer", default=3,
              help="Maximum distance for neighborhood structures in VNS [default: %default]"),
  make_option(c("-i", "--iterprint"), type="integer", default=1,
              help="Print frequency for VNS iterations (0=no print, 1=print) [default: %default]"),
  make_option(c("-P", "--preprocessing"), type="logical", default=TRUE,
              help="Preprocess network before optimization [default: %default]"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE,
              help="Verbose output [default: %default]"),
  make_option(c("-o", "--output"), type="character", default="output/meigo",
              help="Output directory base path [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="MEIGO VNS Network Optimization Script")
opt <- parse_args(opt_parser)

# Dataset configuration mapping
get_dataset_config <- function(dataset_name) {
  dataset_configs <- list(
    "toy" = list(name="ToyModel", sif="ToyModel.sif", rdata="ToyModel.RData", 
                 midas="ToyModel.csv", bnet="ToyModel.bnet", default_time=10),
    "apoptosis" = list(name="Apoptosis", sif="Apoptosis.sif", rdata="Apoptosis.RData",
                       midas="Apoptosis.csv", bnet="Apoptosis.bnet", default_time=10),
    "dream" = list(name="DREAMmodel", sif="DreamModel.sif", rdata="DreamModel.RData",
                   midas="DreamModel.csv", bnet="DreamModel.bnet", default_time=30),
    "TCell" = list(name="TCell", sif="TCell.sif", rdata="TCell.RData",
                   midas="TCell.csv", bnet="TCell.bnet", default_time=10)
  )
  
  if (!dataset_name %in% names(dataset_configs)) {
    stop(paste("Unknown dataset:", dataset_name, 
               ". Available datasets:", paste(names(dataset_configs), collapse=", ")))
  }
  
  return(dataset_configs[[dataset_name]])
}

# Helper function to create objective function for VNS optimization
create_vns_objective_function <- function(cnolist, model) {
  cat("Creating VNS objective function...\n")
  
  # This function creates a closure that captures the model and data
  # The objective function will be called by MEIGO during optimization
  get_fobj <- function(cnolist, model) {
    f <- function(x, model1=model, cnolist1=cnolist) {
      # Prepare simulation list from the model
      simlist <- prep4sim(model1)
      
      # Compute the score for the current solution x
      # This represents how well the network with parameters x fits the data
      score <- computeScoreT1(cnolist1, model1, x)
      return(score)
    }
    return(f)
  }
  
  fobj <- get_fobj(cnolist, model)
  return(fobj)
}

# Helper function to save MEIGO results in multiple formats
save_meigo_results <- function(results_var_name, output_path, method_name) {
  cat("Saving MEIGO results in multiple formats...\n")
  
  # Save the R object
  rdata_file <- file.path(output_path, paste0(method_name, "_report.RData"))
  save(list = results_var_name, file = rdata_file, envir = .GlobalEnv)
  cat("R data saved to:", rdata_file, "\n")
  
  # Extract the results object for JSON conversion
  results_obj <- get(results_var_name, envir = .GlobalEnv)
  
  # Convert R object to a more JSON-friendly format
  results_list <- list()
  
  # Extract key components that are typically in MEIGO results
  if (!is.null(results_obj$xbest)) {
    results_list$best_solution <- as.vector(results_obj$xbest)
  }
  if (!is.null(results_obj$fbest)) {
    results_list$best_fitness <- as.numeric(results_obj$fbest)
  }
  if (!is.null(results_obj$numeval)) {
    results_list$num_evaluations <- as.numeric(results_obj$numeval)
  }
  if (!is.null(results_obj$time)) {
    results_list$optimization_time <- as.numeric(results_obj$time)
  }
  if (!is.null(results_obj$nfuneval)) {
    results_list$function_evaluations <- as.numeric(results_obj$nfuneval)
  }
  
  # Add any additional fields that might be present
  for (name in names(results_obj)) {
    if (!name %in% c("xbest", "fbest", "numeval", "time", "nfuneval")) {
      tryCatch({
        # Try to convert to a simple format
        if (is.vector(results_obj[[name]]) && length(results_obj[[name]]) < 1000) {
          results_list[[name]] <- as.vector(results_obj[[name]])
        } else if (is.matrix(results_obj[[name]]) && nrow(results_obj[[name]]) < 100) {
          results_list[[name]] <- as.matrix(results_obj[[name]])
        }
      }, error = function(e) {
        # Skip fields that can't be easily converted
        cat("Skipping field", name, "- conversion error\n")
      })
    }
  }
  
  # Save as JSON
  json_file <- file.path(output_path, paste0(method_name, "_report.json"))
  write_json(results_list, json_file, pretty = TRUE, auto_unbox = TRUE)
  cat("JSON results saved to:", json_file, "\n")
  
  return(results_list)
}

# Main VNS optimization function
run_vns_optimization <- function(dataset_name, vns_config, change_percent = 0.0, 
                                 output_base = "output/meigo", PreprocessingNetwork = TRUE) {
  
  cat("=== Starting VNS (MEIGO) Optimization ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Change Percent:", change_percent, "\n")
  cat("VNS Configuration:\n")
  for (param in names(vns_config)) {
    cat(sprintf("  %s: %s\n", param, vns_config[[param]]))
  }
  
  # Get dataset configuration
  dataset_config <- get_dataset_config(dataset_name)
  
  # Determine if perturbation is applied
  perturbation_applied <- change_percent > 0
  folder_name <- if (!perturbation_applied) "0_Modified" else paste0(change_percent * 100, "_Modified")
  
  # Set up file paths
  input_path <- file.path("data", dataset_config$name, folder_name)
  sif_file <- file.path(input_path, dataset_config$sif)
  midas_file <- file.path(input_path, dataset_config$midas)
  
  # Create output directory
  output_path <- file.path(output_base, dataset_config$name, folder_name, "VNS")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Validate input files
  if (!file.exists(sif_file)) {
    stop(paste("SIF file not found:", sif_file))
  }
  if (!file.exists(midas_file)) {
    stop(paste("MIDAS file not found:", midas_file))
  }
  
  cat("Input SIF file:", sif_file, "\n")
  cat("Input MIDAS file:", midas_file, "\n")
  cat("Output directory:", output_path, "\n")
  
  # Load network and data
  cat("\n=== Loading Network and Data ===\n")
  cat("Loading network from SIF file...\n")
  pknmodel <- readSIF(sif_file)
  
  cat("Loading experimental data from MIDAS file...\n")
  cnolist <- makeCNOlist(readMIDAS(midas_file, verbose = vns_config$verbose), subfield = FALSE)
  cnolist <- CNOlist(cnolist)  # Ensure proper CNOlist format for MEIGO
  
  # Save CNOlist plot
  plotCNOlistPDF(CNOlist = cnolist, 
                 filename = file.path(output_path, paste0(dataset_config$name, "_CNOlist.pdf")))
  
  # Preprocess network
  cat("\n=== Preprocessing Network ===\n")
  cat("Finding non-observable/non-controllable species and compressing model...\n")
  model <- preprocessing(data = cnolist, model = pknmodel)
  
  # Prepare VNS optimization problem
  cat("\n=== Setting Up VNS Optimization Problem ===\n")
  
  # Create objective function
  fobj <- create_vns_objective_function(cnolist, model)
  
  # Define problem dimensions based on the model
  nvar <- ncol(model$interMat)  # Number of variables (network edges/interactions)
  cat("Number of optimization variables:", nvar, "\n")
  
  # Define the optimization problem structure
  # This tells MEIGO the bounds and objective function
  problem <- list(
    f = fobj,                    # Objective function to minimize
    x_L = rep(0, nvar),         # Lower bounds (all interactions can be turned off)
    x_U = rep(1, nvar)          # Upper bounds (all interactions can be turned on)
  )
  
  # Configure VNS algorithm options
  opts <- list(
    maxeval = vns_config$maxeval,                    # Maximum function evaluations
    maxtime = vns_config$maxtime,                    # Maximum optimization time
    use_local = vns_config$use_local,                # Use local search refinement
    aggr = vns_config$aggr,                          # Aggregation strategy
    local_search_type = vns_config$local_search_type, # Type of local search
    decomp = vns_config$decomp,                      # Decomposition strategy
    maxdist = vns_config$maxdist,                    # Maximum neighborhood distance
    iterprint = vns_config$iterprint                 # Print iteration information
  )
  
  # Run VNS optimization
  cat("\n=== Running VNS Optimization ===\n")
  cat("Starting Variable Neighborhood Search algorithm...\n")
  cat("This may take several minutes depending on maxeval and maxtime settings...\n")
  
  # Record start time
  start_time <- Sys.time()
  
  # Run VNS optimization with timing
  optimization_time <- system.time({
    Results_VNS <- MEIGO(problem, opts, algorithm = "VNS")
  })
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("VNS optimization completed in %.2f seconds\n", optimization_time[["elapsed"]]))
  
  # Store results in global environment for consistency with original code
  assign("Results_VNS", Results_VNS, envir = .GlobalEnv)
  
  # Process optimization results
  cat("\n=== Processing VNS Results ===\n")
  cat("Best fitness achieved:", Results_VNS$fbest, "\n")
  cat("Number of function evaluations:", Results_VNS$nfuneval, "\n")
  
  # Cut the model using the best solution found by VNS
  optModel <- cutModel(model, Results_VNS$xbest)
  
  # Store optimized model in global environment
  assign("optModel", optModel, envir = .GlobalEnv)
  
  # Simulate the optimized network
  cat("Simulating optimized network...\n")
  simResults <- simulate_CNO(model = model, CNOlist = cnolist, bString = Results_VNS$xbest)
  
  # Save simulation results
  save(simResults, file = file.path(output_path, "evolSimRes.RData"))
  
  # Save MEIGO results in multiple formats
  results_summary <- save_meigo_results("Results_VNS", output_path, "VNSR")
  
  # Process and save optimized model
  cat("\n=== Saving Optimized Model ===\n")
  
  # Process the optimized model using helper functions
  optModel <- processOptimizedModel(optModel)
  
  # Save in different formats
  sif_output <- file.path(output_path, paste0("OPT_", dataset_config$name, ".sif"))
  rdata_output <- file.path(output_path, paste0("OPT_", dataset_config$name, ".RData"))
  bnet_output <- file.path(output_path, paste0("OPT_", dataset_config$name, ".bnet"))
  
  cat("Saving optimized model to:", sif_output, "\n")
  toSIF(optModel, file = sif_output, overwrite = TRUE)
  
  cat("Saving R data to:", rdata_output, "\n")
  save(optModel, file = rdata_output)
  
  cat("Saving BoolNet format to:", bnet_output, "\n")
  tryCatch({
    result <- writeBnetFromModel(optModel, bnet_output)
    verifyBoolNetConversion(optModel, bnet_output)
    conversion_success <- TRUE
  }, error = function(e) {
    cat("Warning: BoolNet conversion encountered an issue:", e$message, "\n")
    conversion_success <- FALSE
  })
  
  # Prepare comprehensive results summary
  results_data <- data.frame(
    dataset = dataset_name,
    method = "VNS",
    change_percent = round(change_percent, 4),
    mse = round(Results_VNS$fbest, 4),
    total_time = round(optimization_time[["elapsed"]], 4),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Save results summary to CSV
  results_file <- file.path(output_path, "optimization_results.csv")
  write.csv(results_data, results_file, row.names = FALSE)
  cat("Results summary saved to:", results_file, "\n")
  
  cat("\n=== VNS Optimization Complete ===\n")
  cat("Best Training Score:", round(Results_VNS$fbest, 4), "\n")
  cat("Total Time:", round(optimization_time[["elapsed"]], 4), "seconds\n")
  cat("Function Evaluations Used:", Results_VNS$nfuneval, "/", vns_config$maxeval, "\n")
  cat("Output directory:", output_path, "\n")
  
  return(results_data)
}

# Main execution function
main <- function() {
  tryCatch({
    # Prepare VNS configuration from command line arguments
    vns_config <- list(
      maxeval = opt$maxeval,
      maxtime = opt$maxtime,
      use_local = opt$use_local,
      aggr = opt$aggr,
      local_search_type = opt$local_search_type,
      decomp = opt$decomp,
      maxdist = opt$maxdist,
      iterprint = opt$iterprint,
      verbose = opt$verbose
    )
    
    cat("=== MEIGO VNS Network Optimization ===\n")
    cat("Starting optimization with the following configuration:\n")
    cat("Dataset:", opt$dataset, "\n")
    cat("Output base:", opt$output, "\n")
    cat("Change percent:", opt$change_percent, "\n")
    
    # Run VNS optimization
    results <- run_vns_optimization(
      dataset_name = opt$dataset,
      vns_config = vns_config,
      change_percent = opt$change_percent,
      output_base = opt$output,
      PreprocessingNetwork = opt$preprocessing
    )
    
    cat("\n=== Analysis Complete ===\n")
    cat("VNS optimization completed successfully!\n")
    cat("Check the output directory for detailed results and visualizations.\n")
    
  }, error = function(e) {
    cat("Error occurred during VNS optimization:\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Error traceback:\n")
    traceback()
    quit(status = 1)
  })
}

# Execute main function when script is run non-interactively
if (!interactive()) {
  main()
}