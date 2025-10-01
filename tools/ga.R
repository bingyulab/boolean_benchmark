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
  make_option(c("-g", "--maxgens"), type="integer", default=500,
              help="Maximum generations for GA [default: %default]"),
  make_option(c("-s", "--sizefac"), type="double", default=0.0001,
              help="Size factor for GA [default: %default]"),
  make_option(c("-p", "--popsize"), type="integer", default=50,
              help="Population size for GA [default: %default]"),
  make_option(c("-e", "--elitism"), type="integer", default=5,
              help="Elitism parameter for GA [default: %default]"),
  make_option(c("-m", "--stallgenmax"), type="integer", default=85,
              help="Maximum stall generations for GA [default: %default]"),
  make_option(c("-r", "--reltol"), type="double", default=0.1,
              help="Relative tolerance for GA [default: %default]"),
  make_option(c("-u", "--pmutation"), type="double", default=0.5,
              help="Mutation probability for GA [default: %default]"),
  make_option(c("-l", "--selpress"), type="integer", default=1.2,
              help="Selection pressure for GA [default: %default]"),
  make_option(c("-n", "--nafac"), type="double", default=1,
              help="NA factor for GA [default: %default]"),
  make_option(c("-t", "--maxtime"), type="integer", default=60,
              help="Maximum time in seconds for GA [default: %default]"),
  make_option(c("-P", "--preprocessing"), type="logical", default=TRUE,
              help="Preprocess network before optimization [default: %default]"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE,
              help="Verbose output [default: %default]"),
  make_option(c("-o", "--output"), type="character", default="output/cellnopt",
              help="Output directory base path [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="CellNOpt GA Optimization Script")
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

# Main GA optimization function
run_ga_optimization <- function(dataset_name, ga_config, change_percent = 0.0, 
                                output_base = "output/cellnopt", PreprocessingNetwork = TRUE) {
  
  cat("=== Starting GA Optimization ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Change Percent:", change_percent, "\n")
  cat("GA Configuration:\n")
  for (param in names(ga_config)) {
    cat(sprintf("  %s: %s\n", param, ga_config[[param]]))
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
  output_path <- file.path(output_base, dataset_config$name, folder_name, "ga")
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
  cnolist <- makeCNOlist(readMIDAS(midas_file, verbose = ga_config$verbose), subfield = FALSE)
  
  # Save CNOlist plot
  plotCNOlistPDF(CNOlist = CNOlist(cnolist), 
                 filename = file.path(output_path, paste0(dataset_config$name, "_CNOlist.pdf")))
  
  # Preprocess network
  if (PreprocessingNetwork) {
    cat("Preprocessing network...\n")
    cat("Finding non-observable/non-controllable species and compressing model...\n")
    model <- preprocessing(data = cnolist, model = pknmodel)
  } else {
    model <- pknmodel
  }

  cat("\n=== Preprocessing Network ===\n")
  cat("Finding non-observable/non-controllable species and compressing model...\n")
  model <- preprocessing(data = cnolist, model = pknmodel)
  
  # Run GA optimization
  cat("\n=== Running GA Optimization ===\n")
  cat("Starting genetic algorithm optimization...\n")
  
  # Calculate residual error and initialize bitstring
  resEcnolist <- residualError(cnolist)
  initBstring <- rep(1, length(model$reacID))
  
  # Record start time
  start_time <- Sys.time()
  
  # Run GA optimization with timing
  optimization_time <- system.time({
    opt_results <- gaBinaryT1(
      CNOlist = cnolist, 
      model = model, 
      initBstring = initBstring,
      maxGens = ga_config$maxGens, 
      sizeFac = ga_config$sizeFac,
      popSize = ga_config$popSize, 
      elitism = ga_config$elitism,
      stallGenMax = ga_config$stallGenMax, 
      relTol = ga_config$relTol,
      pMutation = ga_config$pMutation, 
      selPress = ga_config$selPress,
      NAFac = ga_config$NAFac,
      maxTime = ga_config$maxTime, 
      verbose = ga_config$verbose
    )
  })
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("GA optimization completed in %.2f seconds\n", optimization_time[["elapsed"]]))
  
  # Cut model and simulate results
  cat("\n=== Processing Results ===\n")
  optModel <- cutModel(model, opt_results$bString)
  simResults <- simulate_CNO(model = model, CNOlist = cnolist, bString = opt_results$bString)
  
  # Generate plots and save results
  cat("Generating plots and saving results...\n")
  
  # Plot cut model
  cutAndPlot(model = model, bStrings = list(opt_results$bString), CNOlist = cnolist, plotPDF = TRUE)
  
  # Plot fitness evolution
  pdf(file.path(output_path, paste0(dataset_config$name, "_evolFitT1.pdf")))
  plotFit(optRes = opt_results)
  dev.off()
  
  # Plot model
  plotModel(model, cnolist, bString = opt_results$bString, output = "SVG", 
            filename = file.path(output_path, paste0(dataset_config$name, "_mapback_evolFitT1_1.svg")))
  
  # Save simulation results
  save(simResults, file = file.path(output_path, paste0(dataset_config$name, "_evolSimRes.RData")))
  
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
    method = "GA",
    change_percent = round(change_percent, 4),
    mse = round(opt_results$bScore, 4),
    total_time = round(optimization_time[["elapsed"]], 4),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Save results to CSV
  results_file <- file.path(output_path, "optimization_results.csv")
  write.csv(results_summary, results_file, row.names = FALSE)
  cat("Results summary saved to:", results_file, "\n")
  
  cat("\n=== GA Optimization Complete ===\n")
  cat("Training Score:", round(opt_results$bScore, 4), "\n")
  cat("Total Time:", round(optimization_time[["elapsed"]], 4), "seconds\n")
  cat("Output directory:", output_path, "\n")
  
  return(results_summary)
}

# Main execution
main <- function() {
  tryCatch({
    # Prepare GA configuration from command line arguments
    ga_config <- list(
      maxGens = opt$maxgens,
      sizeFac = opt$sizefac,
      popSize = opt$popsize,
      elitism = opt$elitism,
      stallGenMax = opt$stallgenmax,
      relTol = opt$reltol,
      pMutation = opt$pmutation,
      selPress = opt$selpress,
      NAFac = opt$nafac,
      maxTime = opt$maxtime,
      verbose = opt$verbose
    )
    
    # Run GA optimization
    results <- run_ga_optimization(
      dataset_name = opt$dataset,
      ga_config = ga_config,
      change_percent = opt$change_percent,
      output_base = opt$output,
      PreprocessingNetwork = opt$preprocessing
    )
    
    cat("\nGA optimization completed successfully!\n")
    
  }, error = function(e) {
    cat("Error occurred during GA optimization:\n")
    cat(conditionMessage(e), "\n")
    quit(status = 1)
  })
}

# Execute main function
if (!interactive()) {
  main()
}