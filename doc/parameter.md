# Parameter Documentation for Network Reconstruction Methods

## Executive Summary

This document provides comprehensive parameter documentation for CASPO, VNS, GA, and ILP optimization methods used in network reconstruction. Parameters are categorized by their utility for generalization control and network size management, enabling researchers to configure algorithms appropriately for their experimental objectives.

## CASPO (Constraint-Based Pathway Optimization)

### Core Parameters

**learn.Learner(zipped, dataset, length, discretization, factor)**

The CASPO learner accepts fundamental parameters that control the learning process and model constraints. 

- *zipped* parameter contains the network structure
- *dataset* provides experimental observations

**--threads T**
Controls parallel processing by specifying the number of threads for the clingo solver. Higher values accelerate computation but require adequate computational resources.
*Generalization Impact*: Moderate - affects computational efficiency rather than model generalization
*Network Size Impact*: None

**--conf C**
Configures thread behavior with default setting of "many" for optimal parallel performance. This parameter optimizes solver performance across different hardware configurations.
*Generalization Impact*: None
*Network Size Impact*: None

**--optimum O**
Accepts a logical network in CSV format to bypass learning and proceed directly to enumeration. When multiple networks are provided, the first network serves as the optimization target.
*Generalization Impact*: High - predetermined optimal networks may limit generalization to unseen data
*Network Size Impact*: High - directly constrains the resulting network structure

**--fit F**
Establishes tolerance over fitness with default value of 0. This parameter controls the acceptable deviation from perfect data fitting, allowing trade-offs between model accuracy and complexity.
*Generalization Impact*: High - higher tolerance may improve generalization by preventing overfitting
*Network Size Impact*: Moderate - relaxed fitness constraints may result in more parsimonious models

**--size S**
Sets tolerance over model size with default value of 0. This parameter directly controls the penalty for larger network structures during optimization.
*Generalization Impact*: High - size tolerance affects model complexity and generalization capacity
*Network Size Impact*: Very High - primary parameter for controlling network size

**--factor D**
Controls discretization over the range [0,D] with default value of 100. This parameter determines the granularity of continuous data discretization.
*Generalization Impact*: Moderate - affects how continuous data is represented in the discrete model
*Network Size Impact*: Low - primarily affects data representation rather than network structure

**--discretization T**
Specifies the discretization function: round, floor, or ceil, with default of round. This parameter determines how continuous values are converted to discrete states.
*Generalization Impact*: Moderate - different discretization methods may affect model performance on new data
*Network Size Impact*: None

**--length L**
Sets maximum conjunctions length (sources per hyperedge) with default value of 0 (unbounded). This parameter controls the complexity of logical relationships in the network.
*Generalization Impact*: High - longer conjunctions may lead to overfitting
*Network Size Impact*: High - directly controls the complexity of network interactions

## VNS (Variable Neighborhood Search)

### Core Configuration
**list(maxeval=2000, maxtime=30, use_local=1, aggr=0, local_search_type=1, decomp=1, maxdist=0.5)**

**maxeval**
Maximum number of function evaluations with default value of 1000. This parameter controls the computational budget for the optimization process.
*Generalization Impact*: Moderate - insufficient evaluations may result in suboptimal solutions with poor generalization
*Network Size Impact*: Low - affects optimization quality rather than network size directly

**maxtime**
Maximum CPU time in seconds with default value of 60. This parameter provides time-based termination criteria for the optimization process.
*Generalization Impact*: Moderate - insufficient time may prevent convergence to generalizable solutions
*Network Size Impact*: Low - affects optimization quality rather than network size directly

### Search Strategy Parameters

**maxdist**
Percentage of problem dimension perturbed in the furthest neighborhood, ranging from 0-1 with default value of 0.5. This parameter controls the exploration scope of the search algorithm.
*Generalization Impact*: High - appropriate exploration prevents premature convergence to local optima
*Network Size Impact*: Moderate - broader exploration may discover more compact network representations

**use_local**
Boolean parameter controlling local search activation with default value of 1. Local search improves solution quality through neighborhood exploitation.
*Generalization Impact*: High - local search refinement typically improves generalization
*Network Size Impact*: Moderate - local optimization may find more compact solutions

**aggr**
Aggressive search parameter with default value of 0. When set to 1, local search applies only when the best solution improves.
*Generalization Impact*: Moderate - aggressive search may improve efficiency but limit exploration
*Network Size Impact*: Low - affects search efficiency rather than network size

**local_search_type**
Specifies first improvement (1) or best improvement (2) scheme with default value of 1. This parameter controls the local search termination strategy.
*Generalization Impact*: Low - affects local search efficiency rather than generalization
*Network Size Impact*: None

**decomp**
Decomposition parameter with default value of 1. When activated, local search focuses on variables perturbed in the global phase.
*Generalization Impact*: Low - affects search efficiency rather than generalization
*Network Size Impact*: None

## GA (Genetic Algorithm)

### Population and Generation Parameters

**maxGens**
Maximum number of generations with default value of 1000. This parameter controls the evolutionary process duration.
*Generalization Impact*: Moderate - insufficient generations may prevent convergence to generalizable solutions
*Network Size Impact*: Low - affects optimization quality rather than network size directly

**stallGenMax**
Maximum stall generations with default value of Inf. This parameter terminates optimization when no improvement occurs for the specified number of generations.
*Generalization Impact*: Moderate - premature termination may prevent finding better generalizable solutions
*Network Size Impact*: Low - affects optimization quality rather than network size directly

**popSize**
Population size with default value of 200. This parameter controls the diversity of the evolutionary search.
*Generalization Impact*: High - larger populations maintain diversity and improve generalization
*Network Size Impact*: Moderate - diverse populations may explore different network sizes

**elitism**
Number of best individuals preserved across generations with default value of popSize/10. This parameter maintains high-quality solutions during evolution.
*Generalization Impact*: High - elitism preserves good solutions that may generalize well
*Network Size Impact*: Moderate - elite solutions may influence network size distribution

### Objective Function Parameters

**sizeFac**
Size factor with default value of 1e-04. This parameter penalizes increased model size in the objective function.
*Generalization Impact*: High - appropriate size penalties prevent overfitting and improve generalization
*Network Size Impact*: Very High - primary parameter for controlling network size through regularization

**NAFac**
NA factor with default value of 1. This parameter penalizes unresolved states in the model output.
*Generalization Impact*: High - penalizing unresolved states improves model reliability
*Network Size Impact*: Moderate - may influence network structure to reduce unresolved states

**numStarts**
Number of optimization starts with default value of 5. This parameter controls the number of independent optimization runs.
*Generalization Impact*: High - multiple starts improve the likelihood of finding generalizable solutions
*Network Size Impact*: Moderate - multiple runs may identify different network sizes

## ILP (Integer Linear Programming)

### Solution Quality Parameters

**accountForModelSize**
Boolean parameter with default value of TRUE. When activated, the objective function includes model size considerations.
*Generalization Impact*: High - size accounting typically improves generalization by preventing overfitting
*Network Size Impact*: Very High - directly controls whether network size affects optimization

**sizeFac**
Size factor with default value of 0.0001. This parameter controls the penalty weight for model size in the objective function.
*Generalization Impact*: High - appropriate size penalties improve generalization
*Network Size Impact*: Very High - primary parameter for network size control

**mipGap**
Mixed integer programming gap parameter with default value of 0. This parameter controls the optimality tolerance for the solver.
*Generalization Impact*: Moderate - tighter gaps may improve solution quality
*Network Size Impact*: Low - affects solution quality rather than network size

**relGap**
Relative gap parameter controlling solution quality compared to optimal solutions. Looser gaps enable retrieval of more solutions.
*Generalization Impact*: High - gap tolerance affects the diversity of solutions explored
*Network Size Impact*: Moderate - solution diversity may include different network sizes

### Computational Parameters

**timelimit**
Time limit in seconds with default value of 3600. This parameter controls the maximum computation time for optimization.
*Generalization Impact*: Moderate - sufficient time enables finding better generalizable solutions
*Network Size Impact*: Low - affects optimization quality rather than network size directly

**numSolutions**
Number of solutions to retrieve. This parameter controls the diversity of solutions returned by the solver.
*Generalization Impact*: High - multiple solutions enable selection of models with better generalization
*Network Size Impact*: High - diverse solutions may exhibit different network sizes

**limitPop**
Population limit with default value of 500. This parameter controls the maximum number of solutions maintained during optimization.
*Generalization Impact*: Moderate - larger populations may improve solution quality
*Network Size Impact*: Moderate - population diversity may include different network sizes

**poolIntensity**
Pool intensity with default value of 0. This parameter controls the intensity of solution pool management.
*Generalization Impact*: Low - affects solver efficiency rather than generalization
*Network Size Impact*: None

**poolReplace**
Pool replacement strategy with default value of 2. This parameter controls how solutions are replaced in the solution pool.
*Generalization Impact*: Low - affects solver efficiency rather than generalization
*Network Size Impact*: None

## Recommendations for Parameter Configuration

### For Improved Generalization
Configure parameters that control model complexity and prevent overfitting. Key parameters include sizeFac (GA, ILP), fit tolerance (CASPO), length constraints (CASPO), and population diversity (GA). Moderate size penalties combined with sufficient optimization time typically yield models with better generalization capacity.

### For Network Size Control
Focus on parameters that directly influence model complexity. Primary parameters include sizeFac (GA, ILP), size tolerance (CASPO), length constraints (CASPO), and accountForModelSize (ILP). Increasing size penalties generally results in more parsimonious networks, while reducing these penalties allows larger, more complex models.

### For Experimental Robustness
Configure parameters to ensure consistent performance across different perturbation levels. Use multiple optimization starts (numStarts in GA), diverse solution pools (numSolutions in ILP), and appropriate termination criteria (maxeval, maxtime) to maintain algorithm stability under varying network conditions.