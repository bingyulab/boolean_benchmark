library(CellNOptR)

writeBnetFromModel <- function(model, filename) {
  # Extract species names - these will be our nodes
  species <- model$namesSpecies
  
  # Initialize rules list - each species will have a list of regulatory inputs
  rules <- setNames(vector("list", length(species)), species)
  
  # Process each reaction to build regulatory relationships
  for (j in seq_along(model$reacID)) {
    # Skip empty reactions
    if (model$reacID[j] == "") next
    
    # Find inputs (reactants) and outputs (products) for this reaction
    input_indices <- which(model$interMat[,j] == -1)
    output_indices <- which(model$interMat[,j] == 1)
    
    # Skip if no clear input-output relationship
    if (length(input_indices) == 0 || length(output_indices) == 0) next
    
    # Get the actual species names
    inputs <- rownames(model$interMat)[input_indices]
    outputs <- rownames(model$interMat)[output_indices]
    
    # Process each input-output pair
    for (input in inputs) {
      for (output in outputs) {
        # Check if this input acts as an inhibitor
        # Find the row index for this input species
        input_row_idx <- which(rownames(model$notMat) == input)
        
        # Determine if this is an inhibitory relationship
        is_inhibitor <- FALSE
        if (length(input_row_idx) > 0) {
          is_inhibitor <- model$notMat[input_row_idx, j] == 1
        }
        
        # Create the regulatory clause
        if (is_inhibitor) {
          clause <- paste0("!", input)
        } else {
          clause <- input
        }
        
        # Add this clause to the output's rule list
        rules[[output]] <- c(rules[[output]], clause)
      }
    }
  }
  
  # Identify input nodes (nodes with no incoming regulations)
  input_nodes <- species[sapply(rules, length) == 0]
  
  # Prepare BoolNet format lines
  bnet_lines <- c("targets, factors")
  
  # Add input nodes first (they regulate themselves)
  for (node in input_nodes) {
    bnet_lines <- c(bnet_lines, paste(node, ",", node))
  }
  
  # Add regulated nodes (combine multiple inputs with OR)
  regulated_nodes <- species[sapply(rules, length) > 0]
  for (node in regulated_nodes) {
    # Remove duplicates and combine with OR
    unique_rules <- unique(rules[[node]])
    rule_str <- paste(unique_rules, collapse = " | ")
    bnet_lines <- c(bnet_lines, paste(node, ",", rule_str))
  }
  # Write to file
  writeLines(bnet_lines, filename)
  
  # Return summary information for verification
  cat("Conversion Summary:\n")
  cat("Total species:", length(species), "\n")
  cat("Input nodes:", length(input_nodes), "-", paste(input_nodes, collapse = ", "), "\n")
  cat("Regulated nodes:", length(regulated_nodes), "\n")
  cat("Total reactions processed:", sum(model$reacID != ""), "\n")
  
  return(list(
    input_nodes = input_nodes,
    regulated_nodes = regulated_nodes,
    rules = rules
  ))
}

# Helper function to verify the conversion
verifyBoolNetConversion <- function(model, bnet_file) {
  # Read the generated BoolNet file
  bnet_lines <- readLines(bnet_file)
  
  cat("Generated BoolNet file contents:\n")
  cat(paste(bnet_lines, collapse = "\n"), "\n\n")
  
  # Parse and analyze the rules
  rule_lines <- bnet_lines[-1]  # Skip header
  
  cat("Analysis of generated rules:\n")
  for (line in rule_lines) {
    parts <- strsplit(line, ",")[[1]]
    target <- trimws(parts[1])
    factors <- trimws(parts[2])
    
    cat(sprintf("%-8s <- %s\n", target, factors))
  }
  
  # Check for potential issues
  cat("\nValidation checks:\n")
  
  # Check if all species are represented
  targets <- sapply(strsplit(rule_lines, ","), function(x) trimws(x[1]))
  missing_species <- setdiff(model$namesSpecies, targets)
  if (length(missing_species) > 0) {
    cat("WARNING: Missing species in output:", paste(missing_species, collapse = ", "), "\n")
  } else {
    cat("✓ All species are represented\n")
  }
  
  # Check for undefined regulators
  all_factors <- paste(sapply(strsplit(rule_lines, ","), function(x) trimws(x[2])), collapse = " ")
  # Remove logical operators and extract variable names
  factors_cleaned <- gsub("[!|&()]", " ", all_factors)
  unique_factors <- unique(trimws(unlist(strsplit(factors_cleaned, "\\s+"))))
  unique_factors <- unique_factors[unique_factors != ""]
  
  undefined_factors <- setdiff(unique_factors, model$namesSpecies)
  if (length(undefined_factors) > 0) {
    cat("WARNING: Undefined factors:", paste(undefined_factors, collapse = ", "), "\n")
  } else {
    cat("✓ All regulatory factors are defined species\n")
  }
}

write_boolnet_to_sif <- function(net, file) {
  sif <- data.frame()
  for (gene in names(net$interactions)) {
    expr <- net$interactions[[gene]]$expression
    # Remove parentheses and split by operators
    expr_clean <- gsub("[()]", "", expr)
    # Split by | or & (OR/AND), then trim whitespace
    regulators <- unlist(strsplit(expr_clean, "[|&]"))
    regulators <- trimws(regulators)
    regulators <- regulators[regulators != "" & regulators != gene]
    for (reg in regulators) {
      if (startsWith(reg, "!")) {
        interaction <- -1
        reg_clean <- substring(reg, 2)
      } else {
        interaction <- 1
        reg_clean <- reg
      }
      sif <- rbind(sif, data.frame(source=reg_clean, interaction=interaction, target=gene, stringsAsFactors=FALSE))
    }
  }
  write.table(sif, file=file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Usage:
# write_boolnet_to_sif(net, "network.sif")

SIFToBoolNet <- function(sifFile, boolnetFile, CNOlist, model=NULL, fixInputs=TRUE, preprocess=TRUE, ignoreAnds=TRUE)
{
  if (preprocess)
  {
    if (is.null(model))
      sif <- readSIF(sifFile)
    else
      sif = model
    sif <- preprocessing(data=CNOlist, model=sif)
    system("rm tmp.sif")
    writeSIF(sif,file="tmp.sif", overwrite=TRUE)
    sifFile <- "tmp.sif"
  }
    
  if ((class(CNOlist)=="CNOlist")==FALSE){
      CNOlist = CellNOptR::CNOlist(CNOlist)
  } 

  sif <- read.table(sifFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
  if (ncol(sif) != 3)
    sif <- read.table(sifFile,sep=" ",header=FALSE,stringsAsFactors=FALSE)
  
  genes <- unique(c(sif[,1],sif[,3],colnames(CNOlist@cues)))
  genes <- gsub("-","_",genes,fixed=TRUE)
  sif[,1] <- gsub("-","_",sif[,1],fixed=TRUE)
  sif[,3] <- gsub("-","_",sif[,3],fixed=TRUE)  
  
  andIdx <- grep("and",genes,fixed=TRUE)
  andIdx <- andIdx[order(genes[andIdx],decreasing=TRUE)]
  ands <- genes[andIdx]

  andStrings <- sapply(ands, function(gene)
  {
    facIdx <- which(sif[,3] == gene)
    if (ignoreAnds && length(facIdx) > 1) 
      return(NA)
      
    if (length(facIdx) > 0)
      paste("(",paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" & "),")",sep="")
  })
  
  if (ignoreAnds)
    ignoredAnds <- ands[is.na(andStrings)]
  else
    ignoredAnds <- c()
    
  geneStrings <- sapply(genes[-andIdx], function(gene)
  {
    facIdx <- which(sif[,3] == gene & !(sif[,1] %in% ignoredAnds))    
      
    if (length(facIdx) > 0)
      paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" | ")
    else
    if (fixInputs && (gene %in% colnames(CNOlist@cues)) && 
        !(gene %in% colnames(CNOlist@signals[[1]])))
      paste(gene, "[1.0]", sep="")
    else    
      gene
  })
  
  geneStrings <- sapply(geneStrings, function(s)
  {
    
    for (a in 1:length(andIdx))
    {
      s <- gsub(ands[a], andStrings[a], s, fixed=TRUE)
    }
    s
  })
  genes <- genes[-andIdx]
  print(genes)
  #print(outputStrings)
  sink(boolnetFile)
  cat("targets, factors\n")  
  for (i in 1:length(genes))
  {
    cat(genes[i],", ",geneStrings[i],"\n",sep="")
  }    
  sink()
}


# library(CellNOptR)
# CNOl <- makeCNOlist(readMIDAS(midasFile),subfield=FALSE)
# if ((class(CNOl)=="CNOlist")==FALSE){
#     CNOl = CellNOptR::CNOlist(CNOl)
# }
# CNOl <- makeCNOlist(readMIDAS(midasFile),subfield=FALSE)
# SIFToBoolNet(sifFile, "model_draft.txt", CNOl, fixInputs=fixInputs, preprocess=preprocess)
  

simulate_CNO <- function(model, bString, simList=NULL, CNOlist, indexList=NULL)
{

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
    if (is.null(simList)==TRUE){
        simList <- prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # keep simList and indxList for back compatibility ?
    modelCut <- cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim0[,indexList$inhibited] <- 1 - Sim0[,indexList$inhibited]
    simRes0 <- as.matrix(Sim0[,indexList$signals])

    # t1
    Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim[,indexList$inhibited] <- 1 - Sim[,indexList$inhibited]
    
    colnames(Sim) <- model$namesSpecies
    simRes <- as.matrix(Sim[,indexList$signals])

    sig <- #as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
              as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
    
    simResults <- list(input=CNOlist@cues,
                       t0=simRes0, t1=simRes, trueSig=CNOlist@signals[[2]])

    return(simResults)
}



fillNotMatFromReactions <- function(optModel) {
  # Initialize notMat if it doesn't exist or is empty
  if (is.null(optModel$notMat) || ncol(optModel$notMat) == 0) {
    n <- nrow(optModel$interMat)
    m <- ncol(optModel$interMat)
    optModel$notMat <- matrix(0, nrow = n, ncol = m)
    rownames(optModel$notMat) <- rownames(optModel$interMat)
    colnames(optModel$notMat) <- colnames(optModel$interMat)
    
    cat("Initialized empty notMat with dimensions:", n, "x", m, "\n")
  }
  
  # Process each reaction to identify negated components
  for (j in seq_along(optModel$reacID)) {
    reaction_id <- optModel$reacID[j]
    
    # Skip empty reactions
    if (is.na(reaction_id) || reaction_id == "") {
      next
    }
    
    cat("\nProcessing reaction", j, ":", reaction_id, "\n")
    
    # Parse the reaction to extract components
    reaction_info <- parseComplexReaction(reaction_id)
    
    if (is.null(reaction_info)) {
      cat("  WARNING: Could not parse reaction:", reaction_id, "\n")
      next
    }
    
    # Apply negation information to notMat
    for (component in reaction_info$source_components) {
      species_name <- component$species
      is_negated <- component$negated
      
      # Find the row index for this species
      species_row_idx <- which(rownames(optModel$notMat) == species_name)
      
      if (length(species_row_idx) == 0) {
        cat("  WARNING: Species", species_name, "not found in notMat rows\n")
        next
      }
      
      # Check if this species is actually an input in the interMat
      if (optModel$interMat[species_row_idx, j] == -1) {
        if (is_negated) {
          optModel$notMat[species_row_idx, j] <- 1
          cat("  Added negation:", species_name, "in reaction", reaction_id, "\n")
        } else {
          # Explicitly set to 0 for positive regulation (default, but good to be explicit)
          optModel$notMat[species_row_idx, j] <- 0
          cat("  Confirmed positive regulation:", species_name, "in reaction", reaction_id, "\n")
        }
      }
    }
  }
  
  return(optModel)
}


# Helper function to parse complex reaction strings
parseComplexReaction <- function(reaction_string) {
  # Split on "=" to separate source and target
  parts <- strsplit(reaction_string, "=")[[1]]
  
  if (length(parts) != 2) {
    return(NULL)  # Invalid reaction format
  }
  
  source_part <- parts[1]
  target_part <- parts[2]
  
  # Parse the source part to handle complex expressions
  source_components <- parseSourceComponents(source_part)
  
  return(list(
    source_components = source_components,
    target = target_part,
    original_string = reaction_string
  ))
}

# Helper function to parse source components (handles AND gates and negation)
parseSourceComponents <- function(source_string) {
  # Split on "+" to handle AND gates
  components <- strsplit(source_string, "\\+")[[1]]
  
  # Process each component to identify negation
  parsed_components <- list()
  
  for (i in seq_along(components)) {
    component <- trimws(components[i])  # Remove whitespace
    
    # Check for negation (starts with "!")
    if (startsWith(component, "!")) {
      species_name <- sub("^!", "", component)
      is_negated <- TRUE
    } else {
      species_name <- component
      is_negated <- FALSE
    }
    
    parsed_components[[i]] <- list(
      species = species_name,
      negated = is_negated
    )
  }
  
  return(parsed_components)
}

# Verification function to check the filled notMat
verifyNotMatFilling <- function(optModel) {
  cat("\n=== NotMat Verification ===\n")
  
  # Check each reaction
  for (j in seq_along(optModel$reacID)) {
    reaction_id <- optModel$reacID[j]
    
    if (is.na(reaction_id) || reaction_id == "") {
      next
    }
    
    cat("\nReaction", j, ":", reaction_id, "\n")
    
    # Find inputs and outputs in interMat
    inputs <- rownames(optModel$interMat)[optModel$interMat[, j] == -1]
    outputs <- rownames(optModel$interMat)[optModel$interMat[, j] == 1]
    
    cat("  Inputs:", paste(inputs, collapse = ", "), "\n")
    cat("  Outputs:", paste(outputs, collapse = ", "), "\n")
    
    # Check notMat entries for inputs
    for (input in inputs) {
      input_row_idx <- which(rownames(optModel$notMat) == input)
      if (length(input_row_idx) > 0) {
        not_value <- optModel$notMat[input_row_idx, j]
        regulation_type <- if (not_value == 1) "INHIBITORY" else "ACTIVATORY"
        cat("  ", input, "->", paste(outputs, collapse = ","), ":", regulation_type, "\n")
      }
    }
  }
  
  # Summary statistics
  total_reactions <- sum(optModel$reacID != "", na.rm = TRUE)
  total_negations <- sum(optModel$notMat == 1, na.rm = TRUE)
  
  cat("\n=== Summary ===\n")
  cat("Total reactions processed:", total_reactions, "\n")
  cat("Total negations identified:", total_negations, "\n")
  cat("Reactions with negation:", sum(colSums(optModel$notMat) > 0), "\n")
}

# Main function to process your optimized model
processOptimizedModel <- function(optModel) {
  cat("Processing optimized CellNOptR model...\n")
  
  # Fill the notMat
  optModel <- fillNotMatFromReactions(optModel)
  
  # Verify the results
  verifyNotMatFilling(optModel)
  
  return(optModel)
}


pkn2sif<-function(model,optimRes=NA,writeSif=FALSE, filename="Model"){
	
  if (is.na(optimRes[1])){
    BStimes<-rep(1,length(model$reacID))
  }else{
    BStimes<-optimRes$bString
  }	

	findOutput<-function(x){
		sp<-which(x == 1)
		sp<-model$namesSpecies[sp]
		}
		
	reacOutput<-apply(model$interMat,2,findOutput)
	
	findInput<-function(x){
		sp<-which(x == -1)
		sp<-model$namesSpecies[sp]
		}
		
	reacInput<-apply(model$interMat,2,findInput)
		
#if the class of reacInput is not a list, then there are no AND gates
	if(!is(reacInput,"list")){
	
		isNeg<-function(x){
			isNegI<-any(x == 1)
			return(isNegI)
			}
			
		inpSign<-apply(model$notMat,2,isNeg)
		inpSign<-!inpSign
		inpSign[inpSign]<-1
		inpSign[!inpSign]<--1
    
		sifFile<-cbind(reacInput,inpSign,reacOutput)
    # preserve matrix structure even if only one row matches
    sifFile <- sifFile[BStimes == 1, , drop = FALSE]
		colnames(sifFile)=NULL
    rownames(sifFile)=NULL
    
		}else{
		
#in this case there are AND gates and so we need to create dummy "and#" nodes
			sifFile<-matrix(ncol=3)
			nANDs<-1
			for(i in 1:length(reacOutput)){
			  if (BStimes[i]==1){

				  if(length(reacInput[[i]]) == 1){
            tmp<-matrix(0,nrow=1,ncol=3)
					  tmp[1,1]<-reacInput[[i]]
					  tmp[1,3]<-reacOutput[i]
					  tmp[1,2]<-ifelse(
						  any(model$notMat[,i] == 1),-1,1)
              sifFile<-rbind(sifFile,tmp)
          
					}else{
					
						for(inp in 1:length(reacInput[[i]])){
              tmp<-matrix(0,nrow=1,ncol=3)
							tmp[1,1]<-reacInput[[i]][inp]
							tmp[1,3]<-paste("and",nANDs,sep="")
              
              tmp[1,2]<-ifelse(
								model$notMat[which(reacInput[[i]][inp]==rownames(model$notMat)),i] == 1,
								-1,1)
              sifFile<-rbind(sifFile,tmp)
							}
						tmp<-matrix(0,nrow=1,ncol=3)	
						tmp[1,1]<-paste("and",nANDs,sep="")
						tmp[1,3]<-reacOutput[i]
						tmp[1,2]<-1
            sifFile<-rbind(sifFile,tmp)
						
            nANDs<-nANDs+1
            }
			    }
				}
				
      sifFile<-sifFile[2:dim(sifFile)[1],]
			
			}

#this is the edge attribute matrix that contains, for each edge, whether it is
#absent from the model (0), present at t1(1) or present at t2(2)
  if (writeSif==TRUE){
	filename<-paste(filename, ".sif", sep="")
    writeSIF(sifFile, filename=filename)
  }

  return(sifFile)
}


toSIF <- function(model, filename, overwrite=FALSE){

    # convert internal model structure to a SIF matrix
    sif = pkn2sif(model)

    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(sif[,1:3],file=filename,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
    else{
       if (overwrite==TRUE){
            write.table(sif[,1:3],file=filename,
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
        }
        else{
           stop(paste("File ", filename, "already exists.",  sep=""))
        }
    }


}


ilpBinaryT1New <- function(cnolist, 
                            model,
                            md,
                            cplexPath,
                            sizeFac = 0.0001, 
                            mipGap = 0, 
                            relGap = 0, 
                            timelimit = 3600, 
                            method = "quadratic",
                            numSolutions = 100, 
                            limitPop = 500, 
                            poolIntensity = 0, 
                            poolReplace = 2,
                            temp_dir = tempdir()){
  
  ## Initializing auxilliary objects for the analysis
  resILPAll <- list()
  exclusionList <- NULL
  cnolistReal <- cnolist
  
  if(!file.exists(cplexPath)){
    stop("User should provide a valid CPLEX path for using this function.")
  }
    
  ## Form of the objective function
  if(!(method%in%c("quadratic", "linear"))){
    print("'method' variable should be quadratic or linear. 
          Setting method = 'qadratic'")
    method <- "quadratic"
  }
  
  ## Penalty factors as decimals
  options(scipen = 10)

  startILP <- Sys.time()
  print("Creating LP file and running ILP...")
  print(paste("Cplex path:", cplexPath))
  print(paste("Method:", method))
  print(paste("Size factor:", sizeFac))
  print(paste("MIP gap:", mipGap))
  print(paste("Relative gap:", relGap))
  print(paste("Time limit:", timelimit, "seconds"))
  print(paste("Number of solutions:", numSolutions))
  print(paste("Limit population:", limitPop))
  print(paste("Pool intensity:", poolIntensity))
  print(paste("Pool replace:", poolReplace))

  resILP <- CellNOptR:::createAndRunILP(model, md, cnolistReal, accountForModelSize = TRUE, 
                           sizeFac = sizeFac, mipGap=mipGap, relGap=relGap, 
                           timelimit=timelimit, cplexPath = cplexPath, 
                           method = method, numSolutions = numSolutions, 
                           limitPop = limitPop, poolIntensity = poolIntensity, 
                           poolReplace = poolReplace)
  endILP <- Sys.time()
  
  CellNOptR:::cleanupILP()
  
  return(resILP)
  
}