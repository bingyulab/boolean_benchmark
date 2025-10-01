import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os
import shutil
import json
from tools.config import dataset_map
from tools.comparison import limit_float, AttractorAnalysis
from tools.topology_analysis import NetworkTopologyAnalyzer, BooleanNetworkGraph


class MEIGOOptimizer:
    def __init__(self, dataset="toy", manager=None):
        """Initialize the MEIGO analyzer"""
        self.cellnoptr = importr('MEIGOR')
        self.PERTUB = True if manager.change_percent > 0 else False
        self.ChangePct = manager.change_percent
        self.vns_config = manager.get_vns_config()
        self.parse(dataset)
        print("MEIGO loaded successfully")
        
    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name, _ = dataset_map[dataset]
        self.dataset = base_name
        self.filename = "0_Modified" if not self.PERTUB else f"{self.ChangePct * 100:.0f}_Modified"
        self.input_path = os.path.join("data", base_name, self.filename)
        self.sif_file = os.path.join(self.input_path, sif_name)
        self.data_file = os.path.join(self.input_path, rdata_name)
        self.midas_file = os.path.join(self.input_path, midas_name)
        self.output_file = os.path.join("output/meigo", self.dataset, self.filename)       
        
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
        r('cnolist <- CNOlist(cnolist)')
        r(f'plotCNOlistPDF(CNOlist=cnolist, filename="{output_file}")') 
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
    
    @staticmethod
    def convert(rdf):
        """
        Convert R DataFrame to pandas DataFrame
        """
        with localconverter(ro.default_converter + pandas2ri.converter):
            output = {}
            for name, obj in rdf.items():
                if hasattr(obj, 'to_dict'):
                    # pandas.DataFrame → list of record‐dicts
                    output[name] = obj.to_dict(orient='records')
                elif hasattr(obj, 'tolist'):
                    # pandas.Series or numpy array → list
                    output[name] = obj.tolist()
                else:
                    # fallback for scalars or other R objects
                    try:
                        output[name] = obj.tolist()
                    except Exception:
                        output[name] = obj
        return output
    
    def _save_results(self, method):
        # save RData file into json.
        # eSSR VNSR 
        print(f'Saving results to {method}_report.RData')
        filepath = f'{method}_report.RData'
        rdf = ro.r['load'](filepath)
        
        output = MEIGOOptimizer.convert(rdf)
        # 3. Write the combined JSON
        new_file = f'{self.output_file}/{method}_report.json'
        with open(new_file, 'w') as f:
            json.dump(output, f, indent=2, default=str)

        # Move file to new directory (keeping original name)
        shutil.move(filepath, os.path.join(self.output_file, f'{method}_report.RData'))

    def run_vns(self):        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')   

        r(f'''
            get_fobj <- function(cnolist, model){{
                f <- function(x, model1=model, cnolist1=cnolist){{
                    simlist = prep4sim(model1)
                    score = computeScoreT1(cnolist1, model1, x)
                    return(score)
                }}
                return(f)
            }}
            fobj <- get_fobj(cnolist, model)
            nvar <- ncol(model$interMat)
            problem <- list(f=fobj, x_L=rep(0, nvar), x_U=rep(1, nvar))
            opts <- list(
                maxeval={self.vns_config['maxeval']}, maxtime={self.vns_config['maxtime']}, 
                use_local={self.vns_config['use_local']}, aggr={self.vns_config['aggr']}, 
                local_search_type={self.vns_config['local_search_type']}, 
                decomp={self.vns_config['decomp']}, maxdist={self.vns_config['maxdist']},
                iterprint={self.vns_config['iterprint']})
            t <- system.time(Results_VNS <- MEIGO(problem, opts, "VNS"))
            optModel <- cutModel(model, Results_VNS$xbest)
            # plotModel(optModel,cnolist)
            assign("optModel_VNS", Results_VNS, envir = .GlobalEnv)
            simResults <- simulate_CNO(model=model,#_orig,
                                CNOlist=cnolist,
                                bString=optModel$xbest)
            save(simResults,file="{self.output_file}/evolSimRes.RData")                    
        ''')
        print("VNS optimization completed.")
        

    def run_ess(self):
        r('''
        initial_pars=createLBodeContPars(model, LB_n = 1, LB_k = 0.09,
	        LB_tau = 0.1, UB_n = 5, UB_k = 0.95, UB_tau = 10, random = TRUE)
        f_hepato <- getLBodeContObjFunction(cnolist, model, initial_pars, indices=NULL,
            time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 1e-03,
            maxStepSize = Inf, maxNumSteps = 1e4, maxErrTestsFails = 50, nan_fac = 5)
        n_pars <- length(initial_pars$LB)
        problem <- list(f=f_hepato, x_L=initial_pars$LB[initial_pars$index_opt_pars],
            x_U=initial_pars$UB[initial_pars$index_opt_pars], x_0=initial_pars$LB[initial_pars$index_opt_pars])
        opts <- list(maxeval=100, local_solver=0, ndiverse=10, dim_refset=6, save_results=1)
        Results_ESS <- MEIGO(problem, opts, algorithm="ESS")
        opt_pars <- initial_pars;
        opt_pars$parValues <- Results_ESS$xbest;
        simData <- plotLBodeFitness(cnolist, model,opt_pars,
        reltol = 1e-05, atol = 1e-03, maxStepSize = 0.01)
        assign("optModel_ESS", Results_ESS, envir = .GlobalEnv)
        ''')
        print("ESS optimization completed.")

    def evaluate_model(self, output_dir):
        print("Evaluating model...")
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(output_dir, fname)
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")
        
        AA = AttractorAnalysis(self.GD_MODEL, opt_fname)
        results = AA.comparison()

        results['total_time']         = limit_float(r('t[["elapsed"]]'))
        results['method']             = "VNS"
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

    def save_results(self, output_dir):
        print(f"Saving results to {output_dir}/")        
        r('library(here)')
        r('source(here::here("tools", "functions.R"))')
        r(f'''
            optModel <- processOptimizedModel(optModel)
        ''')
                        
        sif_fname = os.path.join(output_dir, f"OPT_{self.dataset}.sif")
        rdata_fname = os.path.join(output_dir, f"OPT_{self.dataset}.RData")
        boolnet_fname = os.path.join(output_dir, f"OPT_{self.dataset}.bnet")
        
        print(f"Saving SIF to {sif_fname}...")
        
        r(f'''
            toSIF(optModel, file="{sif_fname}", overwrite=TRUE)
        ''')
        print(f"Saving RData to {rdata_fname}...")
        r(f'''
            save(optModel, file="{rdata_fname}")
        ''')
        print(f"Saving BoolNet to {boolnet_fname}...")
        
        try:
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
            
            # Check if the output file was actually created and is not empty
            if os.path.exists(boolnet_fname) and os.path.getsize(boolnet_fname) > 0:
                conversion_success = True
            else:
                print(f"SIFToBoolNet completed but output file {boolnet_fname} was not created or is empty")
                conversion_success = False
                
        except Exception as e:
            # Handle any R errors or Python exceptions
            print(f"Error during SIFToBoolNet conversion: {e}")
            conversion_success = False
        return conversion_success
    
    def get_vns_results(self):
        # Retrieve VNS results from R
        return r('Results_VNS')

    def get_ess_results(self):
        # Retrieve ESS results from R
        return r('Results_ESS')
    
    def run_full_analysis(self, method="VNS"):
        """
        Run complete CellNOpt analysis pipeline
        Args:
            file (str): Path to network file
            midas_file (str): Path to MIDAS data file
            output_dir (str): Output directory for results
            numSol (int): Number of solutions for ILP optimization
            relGap (float): Relative gap for ILP optimization
        """
        # Load network and data
        model = self.load_network()
        cnolist = self.load_data()
        
        # Preprocess
        self.preprocess_network()

        # Optimize
        output_file = os.path.join(self.output_file, method)                
        os.makedirs(output_file, exist_ok=True)
        print(f"Optimizing with method: {method}")
        if method == "VNS":
            self.run_vns()
            # Save results
            self._save_results("VNSR")            
            
            conversion_success = self.save_results(output_file)
            print(f"Results saved to {output_file}/")
            
            if not conversion_success:
                print(f"Error: SIFToBoolNet conversion failed. Check the output directory {output_file} for details.")
                return model, cnolist, None
            # Evaluate and cross-validate
            print("Evaluating and cross-validating the model...")
            results = self.evaluate_model(output_file)        
            
            print("Analysis completed successfully!")
            return model, cnolist, results
        
        elif method == "ESS":
            self.run_ess()
            # Save results
            self._save_results("eSSR")
            return model, cnolist, self.get_ess_results()
        else:
            raise ValueError(f"Unknown method: {method}")


# Example usage:
if __name__ == "__main__":
    from tools.config import NetworkPerturbationConfig, AdaptiveParameterManager
    config = NetworkPerturbationConfig(
        change_percent=0.0,
        size_adaptation_strength=2.0,
        generalization_focus=True
    )
    manager = AdaptiveParameterManager(config)
    
    optimizer = MEIGOOptimizer(dataset="toy", manager=manager)
    model, cnolist, vns_results = optimizer.run_full_analysis("VNS")
    model, cnolist, ess_results = optimizer.run_full_analysis("ESS")
    print("VNS and ESS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    print("ESS Results:", ess_results)
    
    # Test for general data
    optimizer = MEIGOOptimizer(dataset="toy", manager=manager)
    model, cnolist, vns_results = optimizer.run_full_analysis("VNS")
    print("VNS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    # print("ESS Results:", ess_results)
    