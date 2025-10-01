#!/usr/bin/env python3
"""
CellNOpt Example using RPy2
This script demonstrates how to use CellNOpt with RPy2 for Boolean network optimization.
"""

import os
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from tools.config import dataset_map
from tools.comparison import limit_float, AttractorAnalysis
from tools.topology_analysis import NetworkTopologyAnalyzer, BooleanNetworkGraph

        
class CellNOptAnalyzer:
    """
    A Python class to interface with CellNOpt R package using RPy2
    """

    def __init__(self, dataset="toy", manager=None, kflod=5):
        """Initialize the CellNOpt analyzer"""
        self.cellnoptr = importr('CellNOptR')
        self.PERTUB = True if manager.change_percent > 0 else False
        self.ChangePct = manager.change_percent
        self.ilp_config = manager.get_ilp_config()
        self.ga_config = manager.get_ga_config()
        self.nfold = kflod
        self.parse(dataset)        
        
        # Initialize results storage
        self.cv_results = []
        self.fold_performances = {}
        print(f"CellNOpt Cross-Validator initialized for {dataset} with {self.nfold}-fold CV")
    
    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name, _ = dataset_map[dataset]
        
        # Core paths
        self.dataset = base_name
        self.filename = "0_Modified" if not self.PERTUB else f"{self.ChangePct * 100:.0f}_Modified"
        
        # Main data directory (where the perturbed network and CV splits are stored)
        self.input_path = os.path.join("data", base_name, self.filename)
        
        # Output directory for cross-validation results
        self.cv_output_path = os.path.join("output/cellnopt", base_name, self.filename, "cross_validation")
        
        # Original model files
        self.sif_file = os.path.join(self.input_path, sif_name)
        self.data_file = os.path.join(self.input_path, rdata_name)
        self.midas_file = os.path.join(self.input_path, midas_name)
        self.output_file = os.path.join("output/cellnopt", self.dataset, self.filename)       
        
        self.GD_MODEL = os.path.join("data", base_name, bnet_name)
        self.GD_MODEL_SIF = os.path.join("data", base_name, sif_name)
        
        if not os.path.exists(self.output_file):
            os.makedirs(self.output_file, exist_ok=True)
    
    def load_network(self):
        """
        Load a network from file
        Args:
            file (str): Path to file
        Returns:
            R object containing the network
        """
        print(f"Loading network from {self.sif_file}")

        r(f'pknmodel <- readSIF("{self.sif_file}")')
        return r('pknmodel')
    
    def load_data(self):
        """
        Load experimental data from MIDAS file
        Args:
            midas_file (str): Path to MIDAS file
        Returns:
            R object containing the experimental data
        """        
        output_file = os.path.join(self.output_file, f"{self.dataset}_CNOlist.pdf")
        print(f"Loading data from {self.midas_file}, and saving to {output_file}")
        r(f'cnolist <- makeCNOlist(readMIDAS("{self.midas_file}", verbose=TRUE), subfield=FALSE)')
        r(f'plotCNOlistPDF(CNOlist=CNOlist(cnolist), filename="{output_file}")') 
        return r('cnolist')
        
    def preprocess_network(self):
        """
        Preprocess the network and data
        Finding and cutting the non observable and non controllable species
        Compressing the model
        Expanding the gates
        Returns:
            R object containing the preprocessed model
        """
        print("Preprocessing network...")
        r('model <- preprocessing(data = cnolist, model = pknmodel)')
        return r('model')

    def optimize_network(self, output_file, PLOT=True, method='ga', fold_idx=None):
        """
        Run Boolean network optimization
        Args:
            output_file (str): Directory to save output files within method directory, different for self.output_file
            PLOT (bool): Whether to plot the results
            method (str): Optimization method to use ('ga' or 'ilp')
        Returns:
            R object containing optimization results
        """
        print(f"Running optimization (method = {method}...")
        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')   

        if method == 'ga':  
            r(f'''
                resEcnolist <- residualError(cnolist);
                initBstring <- rep(1,length(model$reacID));
                t <- system.time(opt_results<-gaBinaryT1(
                    CNOlist=cnolist, model=model, initBstring=initBstring,
                    maxGens={self.ga_config["maxGens"]}, sizeFac={self.ga_config["sizeFac"]}, 
                    popSize={self.ga_config["popSize"]},  elitism={self.ga_config["elitism"]}, 
                    stallGenMax={self.ga_config["stallGenMax"]}, relTol={self.ga_config["relTol"]},
                    pMutation={self.ga_config["pMutation"]}, selPress={self.ga_config["selPress"]},
                    NAFac={self.ga_config["NAFac"]},
                    maxTime={self.ga_config["maxTime"]}, verbose={self.ga_config["verbose"]})
                )
                cat("Time taken for optimization:", t[3], "seconds\\n");
                optModel <- cutModel(model, opt_results$bString);
                simResults <- simulate_CNO(model=model,#_orig,
                            CNOlist=cnolist,
                            bString=opt_results$bString)                
            ''')
            if PLOT:
                r(f'''
                    cutAndPlot(model=model, bStrings=list(opt_results$bString),
                        CNOlist=cnolist,plotPDF=TRUE);
                    pdf("{output_file}/{self.dataset}_evolFitT1.pdf");
                    plotFit(optRes=opt_results);
                    dev.off();
                    plotModel(model, cnolist, bString=opt_results$bString, output="SVG", filename="{output_file}/{self.dataset}_mapback_evolFitT1_1.svg");
                    save(simResults,file=paste("{output_file}/{self.dataset}_evolSimRes.RData",sep=""))                    
                ''')
            # Extract results from R
            training_score = limit_float(r('opt_results$bScore'), nbit=4)
        elif method == 'ilp':  
            print("Creating LP file and running ILP...")
            r(f'''     
                t <- system.time(opt_results <- ilpBinaryT1New(
                    CNOlist(cnolist), model, cplexPath = "{self.ilp_config['cplexPath']}",
                    sizeFac = {self.ilp_config['sizeFac']}, mipGap={self.ilp_config['mipGap']}, 
                    relGap={self.ilp_config['relGap']}, timelimit={self.ilp_config['timelimit']}, 
                    method = "{self.ilp_config['method']}", numSolutions = {self.ilp_config['numSolutions']}, 
                    limitPop = {self.ilp_config['limitPop']}, poolIntensity = {self.ilp_config['poolIntensity']}, 
                    poolReplace = {self.ilp_config['poolReplace']})
                    );              
            ''')
            r(f'''                             
                optModel <- cutModel(model, opt_results$bitstringILP[[1]]);
                simResults <- simulate_CNO(model=model,#_orig,
                            CNOlist=cnolist,
                            bString=opt_results$bitstringILP[[1]])
            ''')
            if PLOT:
                r(f'''
                    cutAndPlot(CNOlist = cnolist, model = model, bStrings = list(opt_results$bitstringILP[[1]]), plotPDF=TRUE)
                    plotModel(model, cnolist, bString=opt_results$bitstringILP[[1]], output="SVG", filename="{output_file}/{self.dataset}_ilpFitT1.svg");
                    save(simResults,file=paste("{output_file}/{self.dataset}_ilpSimRes.RData",sep=""))  
                ''')            
            # Extract results from R
            training_score = limit_float(r('opt_results$bScore'), nbit=4)

        # Save the optimized model for this fold
        if fold_idx is not None:
            self.save_fold_model(output_file, fold_idx)
        
        return {
            'training_score': training_score,
            'method': method
        }       

    def save_results(self, output_dir, fold_idx=None):
        """
        Save optimization results
        Args:
            output_dir (str): Output directory
        """
        print(f"Saving results to {output_dir}/")        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')
        r(f'''
            optModel <- processOptimizedModel(optModel)
        ''')

        if fold_idx is None:
            sif_fname = os.path.join(output_dir, f"OPT_{self.dataset}.sif")
            rdata_fname = os.path.join(output_dir, f"OPT_{self.dataset}.RData")
            boolnet_fname = os.path.join(output_dir, f"OPT_{self.dataset}.bnet")
        else:
            sif_fname = os.path.join(output_dir, f"fold_{fold_idx}_optimized.sif")
            rdata_fname = os.path.join(output_dir, f"fold_{fold_idx}_optimized.RData")
            boolnet_fname = os.path.join(output_dir, f"fold_{fold_idx}_optimized.bnet")
            
        
        print(f"Saving SIF to {sif_fname}...")
        
        r(f'''
            toSIF(optModel, file="{sif_fname}", overwrite=TRUE)
        ''')
        print(f"Saving RData to {rdata_fname}...")
        r(f'''
            save(optModel, file="{rdata_fname}")
        ''')
        print(f"Saving BoolNet to {boolnet_fname}...")
        r(f'''            
            # SIFToBoolNet(sifFile     = "{sif_fname}",
            #             boolnetFile = "{boolnet_fname}",
            #             CNOlist     = cnolist,
            #             model       = optModel,
            #             fixInputs   = FALSE,
            #             preprocess  = TRUE,
            #             ignoreAnds  = TRUE)
            
            result <- writeBnetFromModel(optModel, "{boolnet_fname}")

            verifyBoolNetConversion(optModel, "{boolnet_fname}")
        ''')
        
    def evaluate_model(self, output_dir):
        """
        Evaluate the model using the CellNOptR package
        """
        print("Evaluating model...")
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(output_dir, fname)
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")
        
        AA = AttractorAnalysis(self.GD_MODEL, opt_fname)
        results = AA.comparison()
        
        if 'ilp' in output_dir:
            results['method']     = "ILP"
        else:
            results['method']     = "GA"
        # print("Total time taken for evaluation:", r('t[["elapsed"]]'))
        results['total_time']         = limit_float(r('t[["elapsed"]]'))
        results['change_percent']     = limit_float(self.ChangePct)
        # results.to_csv(os.path.join(output_dir, "results.csv"), index=False)

        opt_sif = os.path.join(output_dir, f"OPT_{self.dataset}.sif")
        print(f"Comparing original topology {self.GD_MODEL} with optimized topology {opt_sif}")
        sif_net1 = BooleanNetworkGraph.read_sif(self.GD_MODEL_SIF)
        sif_net2 = BooleanNetworkGraph.read_sif(opt_sif)

        print(f"SIF Network 1: {sif_net1.number_of_nodes()} nodes, {sif_net1.number_of_edges()} edges")
        print(f"SIF Network 2: {sif_net2.number_of_nodes()} nodes, {sif_net2.number_of_edges()} edges")
        
        # Quick similarity check
        sif_analyzer = NetworkTopologyAnalyzer(sif_net1, sif_net2)
        jaccard = sif_analyzer.jaccard_similarity()

        print(f"Jaccard similarity between SIF networks: {jaccard}")
        results['jaccard_topology'] = jaccard

        return results

    def run_full_analysis(self, method="ga"):
        # Load network and data
        model = self.load_network()
        cnolist = self.load_data()
        
        # Preprocess
        self.preprocess_network()

        # Optimize
        output_file = os.path.join(self.output_file, method)                
        os.makedirs(output_file, exist_ok=True)
        print(f"Optimizing with method: {method}")
        opt_results = self.optimize_network(PLOT=True, method=method, output_file=output_file)

        # Save results
        self.save_results(output_file, fold_idx=None)
        print(f"Results saved to {output_file}/")
        
        # Evaluate and cross-validate
        print("Evaluating and cross-validating the model...")
        results = self.evaluate_model(output_file)        
        print("Analysis completed successfully!")
        return model, cnolist, results
                
        
def main(dataset="toy", method="ga", manager=None):
    """Example usage of CellNOptAnalyzer"""
    
    # file = "data/apoptosis.xml"
    # For MIDAS file, you'll need to create or provide one
    # This is just an example - you'll need actual experimental data
    # midas_file = "data/experimental_data.csv"  # You need to create this

    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer(
        dataset=dataset,
        manager=manager
    )

    model, cnolist, results = analyzer.run_full_analysis(
        method=method
    )
    return model, cnolist, results

if __name__ == "__main__":
    # model, cnolist, results = main(Test=True, method="ilp")
    from tools.config import NetworkPerturbationConfig, AdaptiveParameterManager
    config = NetworkPerturbationConfig(
        change_percent=0.,
        size_adaptation_strength=2.0,
        generalization_focus=True
    )
    manager = AdaptiveParameterManager(config)
    model, cnolist, results = main(dataset="toy", method="ga", manager=manager)