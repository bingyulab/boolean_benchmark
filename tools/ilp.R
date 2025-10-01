#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(optparse)
  library(CellNOptR)
  library(here)
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
  make_option(c("-x", "--cplex_path"), type="character", default="~/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex",
              help="Path to CPLEX executable [default: %default]"),
  make_option(c("-s", "--sizefac"), type="double", default=0.0001,
              help="Size factor for ILP [default: %default]"),
  make_option(c("-g", "--mipgap"), type="double", default=0.05,
              help="MIP gap for ILP [default: %default]"),
  make_option(c("-r", "--relgap"), type="double", default=0.0,
              help="Relative gap for ILP [default: %default]"),
  make_option(c("-t", "--timelimit"), type="integer", default=3600,
              help="Time limit in seconds for ILP [default: %default]"),
  make_option(c("-m", "--method"), type="character", default="quadratic",
              help="ILP method (linear or quadratic) [default: %default]"),
  make_option(c("-n", "--num_solutions"), type="integer", default=100,
              help="Number of solutions to find [default: %default]"),
  make_option(c("-l", "--limit_pop"), type="integer", default=500,
              help="Limit population in solution pool [default: %default]"),
  make_option(c("-i", "--pool_intensity"), type="integer", default=0,
              help="Solution pool search intensity [default: %default]"),
  make_option(c("-p", "--pool_replace"), type="integer", default=2,
              help="Solution pool replacement strategy [default: %default]"),
  make_option(c("-P", "--preprocessing"), type="logical", default=TRUE,
              help="Preprocess network before optimization [default: %default]"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE,
              help="Verbose output [default: %default]"),
  make_option(c("-o", "--output"), type="character", default="output/cellnopt",
              help="Output directory base path [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="CellNOpt ILP Optimization Script")
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

# Main ILP optimization function
run_ilp_optimization <- function(dataset_name, ilp_config, change_percent = 0.0, 
                                 output_base = "output/cellnopt", PreprocessingNetwork = TRUE) {
  
  cat("=== Starting ILP Optimization ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Change Percent:", change_percent, "\n")
  cat("ILP Configuration:\n")
  for (param in names(ilp_config)) {
    cat(sprintf("  %s: %s\n", param, ilp_config[[param]]))
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
  output_path <- file.path(output_base, dataset_config$name, folder_name, "ilp")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Validate input files
  if (!file.exists(sif_file)) {
    stop(paste("SIF file not found:", sif_file))
  }
  if (!file.exists(midas_file)) {
    stop(paste("MIDAS file not found:", midas_file))
  }
  
  # Check if CPLEX executable exists
  if (!file.exists(ilp_config$cplexPath)) {
    stop(paste("CPLEX executable not found at:", ilp_config$cplexPath))
  }
  
  cat("Input SIF file:", sif_file, "\n")
  cat("Input MIDAS file:", midas_file, "\n")
  cat("CPLEX path:", ilp_config$cplexPath, "\n")
  cat("Output directory:", output_path, "\n")
  
  # Load network and data
  cat("\n=== Loading Network and Data ===\n")
  cat("Loading network from SIF file...\n")
  pknmodel <- readSIF(sif_file)
  
  cat("Loading experimental data from MIDAS file...\n")
  md <- readMIDAS(midas_file)
  cnolist <- makeCNOlist(md, subfield = FALSE)

  # Save CNOlist plot
  plotCNOlistPDF(CNOlist = CNOlist(cnolist), 
                 filename = file.path(output_path, paste0(dataset_config$name, "_CNOlist.pdf")))
  
  # Preprocess network
  cat("\n=== Preprocessing Network ===\n")
  cat("Finding non-observable/non-controllable species and compressing model...\n")
  model <- preprocessing(data = cnolist, model = pknmodel)
  
  # Run ILP optimization
  cat("\n=== Running ILP Optimization ===\n")
  cat("Creating Linear Programming formulation and solving with CPLEX...\n")
  
  
  # Run ILP optimization with timing
  optimization_time <- system.time({
    opt_results <- ilpBinaryT1New(
      CNOlist(cnolist), 
      model, 
      md,
      cplexPath = ilp_config$cplexPath,
      sizeFac = ilp_config$sizeFac, 
      mipGap = ilp_config$mipGap,
      relGap = ilp_config$relGap, 
      timelimit = ilp_config$timelimit,
      method = ilp_config$method, 
      numSolutions = ilp_config$numSolutions,
      limitPop = ilp_config$limitPop, 
      poolIntensity = ilp_config$poolIntensity,
      poolReplace = ilp_config$poolReplace
    )
  })
  
  cat(sprintf("ILP optimization completed in %.2f seconds\n", optimization_time[["elapsed"]]))
  
  # Cut model and simulate results
  cat("\n=== Processing Results ===\n")
  optModel <- cutModel(model, opt_results$bitstringILP[[1]])
  simResults <- simulate_CNO(model = model, CNOlist = cnolist, bString = opt_results$bitstringILP[[1]])
  
  # Generate plots and save results
  cat("Generating plots and saving results...\n")
  
  # Plot cut model
  cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(opt_results$bitstringILP[[1]]), plotPDF = TRUE)
  
  # Plot model
  plotModel(model, cnolist, bString = opt_results$bitstringILP[[1]], output = "SVG", 
            filename = file.path(output_path, paste0(dataset_config$name, "_ilpFitT1.svg")))
  
  # Save simulation results
  save(simResults, file = file.path(output_path, paste0(dataset_config$name, "_ilpSimRes.RData")))
  
  # Process and save optimized model
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
  result <- writeBnetFromModel(optModel, bnet_output)
  verifyBoolNetConversion(optModel, bnet_output)
  
  # Prepare results summary
  results_summary <- data.frame(
    dataset = dataset_name,
    method = "ILP",
    change_percent = round(change_percent, 4),
    mse = round(opt_results$bScore, 4),
    total_time = round(optimization_time[["elapsed"]], 4),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Save results to CSV
  results_file <- file.path(output_path, "optimization_results.csv")
  write.csv(results_summary, results_file, row.names = FALSE)
  cat("Results summary saved to:", results_file, "\n")
  
  cat("\n=== ILP Optimization Complete ===\n")
  cat("Training Score:", round(opt_results$bScore, 4), "\n")
  cat("Total Time:", round(optimization_time[["elapsed"]], 4), "seconds\n")
  cat("Output directory:", output_path, "\n")
  
  return(results_summary)
}

# Main execution
main <- function() {
  tryCatch({
    # Prepare ILP configuration from command line arguments
    ilp_config <- list(
      cplexPath = opt$cplex_path,
      sizeFac = opt$sizefac,
      mipGap = opt$mipgap,
      relGap = opt$relgap,
      timelimit = opt$timelimit,
      method = opt$method,
      numSolutions = opt$num_solutions,
      limitPop = opt$limit_pop,
      poolIntensity = opt$pool_intensity,
      poolReplace = opt$pool_replace,
      verbose = opt$verbose
    )
    
    # Run ILP optimization
    results <- run_ilp_optimization(      
      dataset_name = opt$dataset,
      ilp_config = ilp_config,
      change_percent = opt$change_percent,
      output_base = opt$output,
      PreprocessingNetwork = opt$preprocessing
    )
    
    cat("\nILP optimization completed successfully!\n")
    
  }, error = function(e) {
    cat("Error occurred during ILP optimization:\n")
    cat(conditionMessage(e), "\n")
    quit(status = 1)
  })
}

# Execute main function
if (!interactive()) {
  main()
}