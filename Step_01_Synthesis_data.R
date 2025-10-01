library(BoolNet)
library(CellNOptR)
library(here)

source(here::here("tools", "data.R"))

# Goal: Based on prior knowledge network generate experiment data
# Input: SIF format prior knowledge network (loading via CellNOptR readSIF). At first my convert from sbml file.
# Output: Midas format data
# Step 1: Read data from sbml, convert into bnet and sif format
# Step 2: Choose Cues, Stimuli, Inhibitors and Signals
# Step 3: Generate CNOlist(following is example) based on attractors. Based on cues select network, 1 keep, 0 delete all its connections.
#         Then use function in BoolNet use  following function to simulate.  asynchronous,  set genesON and genesOFF based on experiment.
#             getAttractors(network, type = c("synchronous","asynchronous"), method = c("exhaustive", "sat.exhaustive", "sat.restricted", "random", "chosen"), startStates = list(), genesON = c(), genesOFF = c(), canonical = TRUE, randomChainLength = 10000, avoidSelfLoops = TRUE, geneProbabilities = NULL, maxAttractorLength = Inf, returnTable = TRUE).
#             Finally take average.
#             $namesCues
#             [1] "EGF"  "TNFa" "Raf"  "PI3K"

#             $namesStimuli
#             [1] "EGF"  "TNFa"

#             $namesInhibitors
#             [1] "Raf"  "PI3K"

#             $namesSignals
#             [1] "Akt"    "Hsp27"  "NFkB"   "Erk"    "p90RSK" "Jnk"    "cJun"

#             $timeSignals
#             [1]  0 10
#             $valueCues
#                 [,1] [,2] [,3] [,4]
#             [1,]    1    0    0    0
#             [2,]    0    1    0    0
#             [3,]    1    1    0    0
#             [4,]    1    0    1    0
#             [5,]    0    1    1    0
#             [6,]    1    1    1    0
#             [7,]    1    0    0    1
#             [8,]    0    1    0    1
#             [9,]    1    1    0    1

#             $valueInhibitors
#                 [,1] [,2]
#             [1,]    0    0
#             [2,]    0    0
#             [3,]    0    0
#             [4,]    1    0
#             [5,]    1    0
#             [6,]    1    0
#             [7,]    0    1
#             [8,]    0    1
#             [9,]    0    1

#             $valueStimuli
#                 [,1] [,2]
#             [1,]    1    0
#             [2,]    0    1
#             [3,]    1    1
#             [4,]    1    0
#             [5,]    0    1
#             [6,]    1    1
#             [7,]    1    0
#             [8,]    0    1
#             [9,]    1    1

#             $valueSignals
#             $valueSignals[[1]]
#                 [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#             [1,]    0    0    0    0    0    0    0
#             [2,]    0    0    0    0    0    0    0
#             [3,]    0    0    0    0    0    0    0
#             [4,]    0    0    0    0    0    0    0
#             [5,]    0    0    0    0    0    0    0
#             [6,]    0    0    0    0    0    0    0
#             [7,]    0    0    0    0    0    0    0
#             [8,]    0    0    0    0    0    0    0
#             [9,]    0    0    0    0    0    0    0

#             $valueSignals[[2]]
#                 [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#             [1,] 0.91  0.0 0.86  0.8 0.88 0.00  0.0
#             [2,] 0.82  0.7 0.90  0.0 0.00 0.25  0.4
#             [3,] 0.91  0.7 0.90  0.8 0.88 0.25  0.4
#             [4,] 0.91  0.0 0.86  0.0 0.00 0.00  0.0
#             [5,] 0.82  0.7 0.90  0.0 0.00 0.25  0.4
#             [6,] 0.91  0.7 0.90  0.0 0.00 0.25  0.4
#             [7,] 0.00  0.0 0.00  0.8 0.88 0.00  0.0
#             [8,] 0.00  0.7 0.90  0.0 0.00 0.25  0.4
#             [9,] 0.00  0.7 0.90  0.8 0.88 0.25  0.4


#             $valueVariances
#             $valueVariances[[1]]
#                 [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#             [1,]    0    0    0    0    0    0    0
#             [2,]    0    0    0    0    0    0    0
#             [3,]    0    0    0    0    0    0    0
#             [4,]    0    0    0    0    0    0    0
#             [5,]    0    0    0    0    0    0    0
#             [6,]    0    0    0    0    0    0    0
#             [7,]    0    0    0    0    0    0    0
#             [8,]    0    0    0    0    0    0    0
#             [9,]    0    0    0    0    0    0    0
# Step 4: Save data in midas format via writeMIDAS in CellNOptR.

# 1) Load SBML -> SIF -> model (readSIF output)
load_model_from_sbml <- function(sbml_file) {
    sif_file <- sub(".sbml$", ".sif", sbml_file)
    if (!file.exists(sif_file)) {
        message("Converting SBML to SIF: ", sbml_file)
        convertSBML(sbml_file) # CellNOptR
        if (!file.exists(sif_file)) {
            stop("Conversion failed; cannot find SIF: ", sif_file)
        }
    }
    model <- readSIF(sif_file)
    model
}

# 2) Detect node types (sources, sinks, internal, observable)
detect_node_types <- function(model) {
    # model$namesSpecies and model$interMat are expected
    if (is.null(model$namesSpecies) || is.null(model$interMat)) {
        stop("Model missing namesSpecies or interMat.")
    }
    names <- model$namesSpecies
    M <- model$interMat

    # degree computations (nonzero count)
    in_degree <- rowSums(M == 1)
    out_degree <- rowSums(M == -1)

    source_nodes <- names[in_degree == 0] # no incoming edges
    sink_nodes <- names[out_degree == 0] # no outgoing edges
    internal_nodes <- setdiff(names, union(source_nodes, sink_nodes))
    # observable nodes: readouts we can measure — pick sinks + a subset of internal nodes
    observable_nodes <- c(sink_nodes, internal_nodes) # conservative default

    list(
        names = names,
        indeg = in_degree,
        outdeg = out_degree,
        source_nodes = source_nodes,
        sink_nodes = sink_nodes,
        internal_nodes = internal_nodes,
        observable_nodes = unique(observable_nodes)
    )
}

# 3) Simple experimental design generator
# Choose cues/stimuli/inhibitors and create experiment matrix (binary treatment matrix)
design_experiments <- function(node_analysis, n_experiments = 30, n_timepoints = c(0, 30), n_readout = 8) {
    stimuli <- node_analysis$source_nodes
    readout <- node_analysis$sink_nodes

    inhibitor <- sample(node_analysis$internal_nodes, length((stimuli)))

    cues <- c(stimuli, inhibitor)

    readout <- c(readout, sample(setdiff(node_analysis$internal_nodes, inhibitor), n_readout - length(readout)))

    # Generate all subsets (powerset)
    all_subsets <- unlist(lapply(1:(length(cues) - 1), function(k) {
        combn(cues, k, simplify = FALSE)
    }), recursive = FALSE)

    # Convert each subset to binary representation
    binary_matrix <- t(sapply(all_subsets, function(subset) {
        as.integer(cues %in% subset)
    }))

    # Add row names to identify subsets
    rownames(binary_matrix) <- sapply(all_subsets, function(subset) {
        paste(subset, collapse = ",")
    })

    # Convert to data frame for readability
    binary_df <- as.data.frame(binary_matrix)
    set.seed(123)
    experiments <- sample(1:nrow(binary_df), n_experiments, replace = TRUE)
    experiments <- binary_df[experiments, , drop = FALSE]
    colnames(experiments) <- cues
    print(length(experiments))
    list(
        namesCues = cues,
        namesStimuli = stimuli,
        namesInhibitors = inhibitor,
        namesSignals = readout,
        timeSignals = n_timepoints,
        valueCues = experiments,
        valueStimuli = experiments[, stimuli, drop = FALSE],
        valueInhibitors = experiments[, inhibitor, drop = FALSE]
    )
}

# 4) Simulation engine (simple continuous propagation + logistic squashing)
# This generates valueSignals as a list of matrices, one per timepoint (experiments x signals)
simulate_network_response <- function(bnet_file, design, n_iterations = 20, simLength=50) {
    network <- loadNetwork(bnet_file)
    stimuli_transformed <- ifelse(design$valueStimuli == 0, -1, 1)
    inhibitor_transformed <- ifelse(design$valueInhibitors == 1, 0, -1)
    experiments <- cbind(stimuli_transformed, inhibitor_transformed)    
    # experiments <- design$valueCues
    readouts <- design$namesSignals
    timepoints <- design$timeSignals

    # Initialize valueSignals as list with matrices for each timepoint
    n_experiments <- nrow(experiments)
    n_readouts <- length(readouts)
    n_timepoints <- length(timepoints)

    valueSignals <- list()
    for (t in 1:n_timepoints) {
        valueSignals[[t]] <- matrix(NA_real_, nrow = n_experiments, ncol = n_readouts)
        colnames(valueSignals[[t]]) <- readouts
        rownames(valueSignals[[t]]) <- rownames(experiments)
    }

    for (exp_idx in 1:nrow(experiments)) {
        fixedOff <- names(experiments[exp_idx, ])[experiments[exp_idx, ] == 0]
        fixedOn <- names(experiments[exp_idx, ])[experiments[exp_idx, ] == 1]

        fixedmodel <- BoolNet::fixGenes(
            network,
            c(fixedOn, fixedOff),
            c(rep(1, length(fixedOn)), rep(0, length(fixedOff)))
        )

        timeSeries <- BoolNet::generateTimeSeries(
            network = fixedmodel,
            numSeries = 1000,
            numMeasurements = simLength,
            type = "asynchronous"
        )

        # Prepare results matrix
        numRuns <- length(timeSeries)
        numGenes <- length(fixedmodel$genes)
        avgActivationMat <- matrix(NA, nrow = numRuns, ncol = numGenes)
        colnames(avgActivationMat) <- fixedmodel$genes

        #Average gene activation per run
        for (i in seq_along(timeSeries)) {
            traj <- timeSeries[[i]]
            avgActivationMat[i, ] <- rowSums(traj) / ncol(traj)
        }
        
        # For each readout, assign values to timepoints
        for (j in seq_along(readouts)) {
            readout_name <- readouts[j]

            # Check if readout is present in the network (not NaN)
            readout_values <- avgActivationMat[, readout_name]
            is_present <- !all(is.nan(readout_values))

            # Timepoint 1 (t=0): 0 if node present, NaN if absent
            valueSignals[[1]][exp_idx, j] <- ifelse(is_present, 0, NaN)

            # Timepoint 2 (t>0): average value if present, NaN if absent
            if (is_present) {
                valueSignals[[2]][exp_idx, j] <- sample(readout_values, 1)
            } else {
                valueSignals[[2]][exp_idx, j] <- NaN
            }
        }
    }
    design$valueSignals <- valueSignals
    design
}

# Command-line interface
suppressPackageStartupMessages({
    library(optparse)
})

option_list <- list(
    make_option(c("-d", "--dataset"),
        type = "character", default = "tcell",
        help = "Dataset name or path to model file", metavar = "STRING"
    ),
    make_option(c("-t", "--timepoints"),
        type = "character", default = "0,30",
        help = "Timepoints (e.g. 0,10,30) for time-series MIDAS", metavar = "STRING"
    ),
    make_option(c("-p", "--proportion"),
        type = "integer", default = 60,
        help = "Proportion of simulated experiments", metavar = "INT"
    ),
    make_option(c("-s", "--seed"),
        type = "integer", default = 44,
        help = "Random seed [default %default]", metavar = "INT"
    )
)

parser <- OptionParser(
    option_list = option_list,
    description = "Generate MIDAS data from Boolean network models"
)

# Only parse arguments if script is run from command line
if (sys.nframe() == 0) {
    opt <- parse_args(parser)
    timepoints <- as.numeric(strsplit(opt$timepoints, ",")[[1]])

    message("=== MIDAS Data Generation Pipeline ===")
    message(sprintf("Dataset: %s", opt$dataset))
    message(sprintf("Timepoints: %s", timepoints))
    message(sprintf("Replicates: %d", opt$proportion))
    message(sprintf("Random seed: %d", opt$seed))

    set.seed(opt$seed)

    dataset_map <- list(
        toy       = c("ToyModel", "ToyModel.sif", "ToyModel.RData", "ToyModel.csv", "ToyModel.bnet"),
        apoptosis = c("Apoptosis", "Apoptosis.sif", "Apoptosis.RData", "Apoptosis.csv", "Apoptosis.bnet"),
        dream     = c("DREAMmodel", "DreamModel.sif", "DreamModel.RData", "DreamModel.csv", "DreamModel.bnet"),
        tcell     = c("TCell", "TCell.sif", "TCell.RData", "TCell.csv", "TCell.bnet")
    )
    # Extract file names for the specified dataset
    vals <- dataset_map[[opt$dataset]]
    base_name <- vals[1]
    sif_file <- file.path("data", base_name, vals[2])
    rdata_file <- file.path("data", base_name, vals[3])
    midas_file <- file.path("data", base_name, vals[4])
    bnet_file <- file.path("data", base_name, vals[5])

    # Load network based on dataset
    if (opt$dataset == "tcell") {
        message("Loading model from SBML/SIF ...")
        sbml_file <- sub(".bnet$", ".sbml", bnet_file)
        model <- load_model_from_sbml(sbml_file)
    }

    if (opt$dataset == "tcell") {
        node_analysis <- detect_node_types(model)
        cat("=== Network Analysis: ===\n")
        cat("- Source nodes (potential stimuli/inhibitors):", length(node_analysis$source_nodes), "\n")
        cat("- Sink nodes (potential readouts):", length(node_analysis$sink_nodes), "\n")
        cat("- Internal nodes:", length(node_analysis$internal_nodes), "\n")
        cat("- Total observable nodes:", length(node_analysis$observable_nodes), "\n")

        n_experiments = as.integer(opt$proportion * length(node_analysis$observable_nodes) / 100)
        experimental_design <- design_experiments(node_analysis,
            n_experiments = n_experiments,
            n_timepoints = timepoints,
            n_readout = 8
        )
    } else if (opt$dataset == "toy") {
        midas_file <- file.path("data", "ToyModel", "ToyModel_test.csv")
        experimental_design <- list()
        stimuli <- c("EGF", "TNFa")
        inhibitor <- c("Raf", "PI3K")
        cues <- c(stimuli, inhibitor)
        experimental_design$namesCues <- cues
        experimental_design$namesStimuli <- stimuli
        experimental_design$namesInhibitors <- inhibitor
        experimental_design$namesSignals <- c("Akt", "Hsp27", "NFkB", "Erk", "p90RSK", "Jnk", "cJun")
        experimental_design$timeSignals <- c(0, 10)
        experiments <- matrix(c(
            1, 0, 0, 0,
            0, 1, 0, 0,
            1, 1, 0, 0,
            1, 0, 1, 0,
            0, 1, 1, 0,
            1, 1, 1, 0,
            1, 0, 0, 1,
            0, 1, 0, 1,
            1, 1, 0, 1
        ), nrow = 9, ncol = 4, byrow = TRUE)
        colnames(experiments) <- cues
        experimental_design$valueCues <- experiments
        experimental_design$valueStimuli = experiments[, stimuli, drop = FALSE]
        experimental_design$valueInhibitors = experiments[, inhibitor, drop = FALSE]
    } else if (opt$dataset == "dream") {
        midas_file <- file.path("data", "DREAMmodel", "DreamModel_test.csv")
        experimental_design <- list()
        stimuli <- c("igf1", "il1a", "tgfa", "tnfa")
        inhibitor <- c("ikk", "mek12", "pi3k", "p38")
        cues <- c(stimuli, inhibitor)
        experimental_design$namesCues <- cues
        experimental_design$namesStimuli <- stimuli
        experimental_design$namesInhibitors <- inhibitor
        experimental_design$namesSignals <- c("akt", "erk12", "ikb", "jnk12", "p38", "hsp27", "mek12")
        experimental_design$timeSignals <- c(0, 30)
        experiments <- matrix(c(
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 1, 0, 1, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 0,
            1, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 1, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 1, 0, 0, 0, 1,
            0, 1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 1, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 1, 0, 0, 1, 0,
            0, 1, 0, 0, 0, 0, 1, 0,
            1, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 1, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 1, 0, 0, 0,
            0, 1, 0, 0, 1, 0, 0, 0,
            1, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 1, 0, 1, 0, 0, 0
        ), nrow = 25, ncol = 8, byrow = TRUE)
        colnames(experiments) <- cues
        experimental_design$valueCues <- experiments
        experimental_design$valueStimuli = experiments[, stimuli, drop = FALSE]
        experimental_design$valueInhibitors = experiments[, inhibitor, drop = FALSE]       
        matrix(c(
            0.91,0,0.86,0.8,0.88,0,0,
            0.82,0.7,0.9,0,0,0.25,0.4,
            0.91,0.7,0.9,0.8,0.88,0.25,0.4,
            0.91,0,0.86,0,0,0,0,
            0.82,0.7,0.9,0,0,0.25,0.4,
            0.91,0.7,0.9,0,0,0.25,0.4,
            0,0,0,0.8,0.88,0,0,
            0,0.7,0.9,0,0,0.25,0.4,
            0,0.7,0.9,0.8,0.88,0.25,0.4
        ), nrow = 6, ncol = 8, byrow = TRUE)
    }

    cat("\n=== Experimental Design: ===\n")
    cat("- Chosen stimuli:", paste(experimental_design$namesStimuli, collapse = ", "), "\n")
    cat("- Chosen inhibitors:", paste(experimental_design$namesInhibitors, collapse = ", "), "\n")
    cat("- Chosen readouts (signals):", paste(experimental_design$namesSignals, collapse = ", "), "\n")
    cat("- Number of experiments:", nrow(experimental_design$valueCues), "\n")
    cat("- Timepoints:", paste(experimental_design$timeSignals, collapse = ", "), "\n")

    message("Simulating network responses ...")
    sim <- simulate_network_response(bnet_file, experimental_design, n_iterations = 20, simLength = 240)

    sim$valueCues       <- as.matrix(sim$valueCues)
    sim$valueInhibitors <- as.matrix(sim$valueInhibitors)
    sim$valueStimuli    <- as.matrix(sim$valueStimuli)
    message("Building CNOlist and writing MIDAS file ...")
    writeMIDAS(sim, midas_file, overwrite=TRUE)
    
    cat("\n=== MIDAS Data Generation Complete! ===\n")
    cat("- Output file:", midas_file, "\n")
    cat("- Stimuli columns:", ifelse(length(sim$namesStimuli) > 0, paste(sim$namesStimuli, collapse = ", "), "(none)"), "\n")
    cat("- Inhibitors columns:", ifelse(length(sim$namesInhibitors) > 0, paste(sim$namesInhibitors, collapse = ", "), "(none)"), "\n")
    cat("- Signals:", paste(sim$namesSignals, collapse = ", "), "\n")
    cat("- Number of experiments:", nrow(sim$valueStimuli), "\n")
    cat("- Time points:", length(sim$timeSignals), "\n")

    cat("\nSample valueSignals (timepoint 1):\n")
}
