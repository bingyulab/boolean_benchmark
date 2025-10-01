from dataclasses import dataclass
from typing import Dict, Any, Sequence
import json, os, sys


dataset_map = {
    "toy": ("ToyModel", "ToyModel.sif", "ToyModel.RData", "ToyModel.csv", "ToyModel.bnet", 10),
    "apoptosis": ("Apoptosis", "Apoptosis.sif", "Apoptosis.RData", "Apoptosis.csv", "Apoptosis.bnet", 10),
    "dream": ("DREAMmodel", "DreamModel.sif", "DreamModel.RData", "DreamModel.csv", "DreamModel.bnet", 30),
    "TCell": ("TCell", "TCell.sif", "TCell.RData", "TCell.csv", "TCell.bnet", 30),
}


@dataclass
class NetworkPerturbationConfig:
    """
    Configuration class for adaptive parameter adjustment based on network perturbation levels.
    
    This class implements a strategy where parameters are adjusted based on the perturbation
    percentage to maintain network size and improve generalization across different modification levels.
    """
    change_percent: float  # Perturbation level from 0.0 to 1.0
    base_size_factor: float = 1e-4  # Base size penalty factor
    size_adaptation_strength: float = 2.0  # How strongly to adapt size penalties
    generalization_focus: bool = True  # Whether to prioritize generalization
    
    def __post_init__(self):
        """Validate configuration parameters after initialization."""
        if not 0.0 <= self.change_percent <= 1.0:
            raise ValueError("change_percent must be between 0.0 and 1.0")
        if self.base_size_factor <= 0:
            raise ValueError("base_size_factor must be positive")

class AdaptiveParameterManager:
    """
    Manages parameter adaptation across different optimization methods based on network perturbation levels.
    
    The core principle is that as perturbation increases, we need to adjust parameters to:
    1. Maintain reasonable network sizes (prevent over-shrinking)
    2. Improve generalization by controlling overfitting
    3. Ensure consistent performance across perturbation levels
    """
    
    def __init__(self, config: NetworkPerturbationConfig):
        self.config = config
        self.change_percent = config.change_percent
        
        # Calculate adaptation factors based on perturbation level
        self._calculate_adaptation_factors()
    
    def _calculate_adaptation_factors(self):
        """
        Calculate how much to adjust parameters based on perturbation level.
        
        The logic here is:
        - As perturbation increases, we reduce size penalties to maintain network size
        - We increase population diversity to improve robustness
        - We adjust termination criteria to allow more exploration
        """
        # Size factor adjustment: decrease penalties as perturbation increases
        # This prevents networks from becoming too small at high perturbation levels
        self.size_factor_multiplier = 1.0 - (self.change_percent * 0.8)
        
        # Exploration factor: increase exploration as perturbation increases
        # This helps algorithms find better solutions in perturbed landscapes
        self.exploration_multiplier = 1.0 + (self.change_percent * 0.5)
        
        # Tolerance factor: increase tolerance as perturbation increases
        # This allows more flexibility in fitting to potentially noisy data
        self.tolerance_multiplier = 1.0 + (self.change_percent * 0.3)
        
        # Population factor: increase population diversity as perturbation increases
        self.population_multiplier = 1.0 + (self.change_percent * 0.4)
    
    def get_caspo_config(self) -> Dict[str, Any]:
        """
        Generate CASPO configuration parameters adapted for the current perturbation level.
        
        Key adaptations:
        - Adjust size tolerance to maintain network size
        - Modify length constraints based on perturbation level
        - Adapt fit tolerance to handle potentially noisier data
        """
        # Base configuration
        
        base_threads = 16
        base_factor = 100
        
        config = {
            "change_percent": self.change_percent,
            'threads': 1,  # Use multiple threads for efficiency
            'conf': 'many',  # Standard configuration
            'fit': 0.04,  # tolerance over fitness (Default to 0)
            'size': 0,  # tolerance over size (Default to 0)
            'factor': 100,  # discretization range [0, factor]
            'discretization': 'round',  # discretization function: round, floor, ceil (Default to round)
            'length': 0,  # max conjunction length (0 = unbounded)
        }

        return config
    
    def _calculate_caspo_length(self) -> int:
        """
        Calculate optimal conjunction length for CASPO based on perturbation level.
        
        The logic here balances model expressiveness with overfitting risk:
        - At low perturbation: Allow moderate complexity (length 3-4)
        - At high perturbation: Reduce complexity to prevent overfitting (length 2-3)
        - Use exploration_multiplier to fine-tune based on search capacity
        """
        if self.change_percent == 0:
            return 0  # Unbounded for original network
        elif self.change_percent < 0.3:
            # Low perturbation: allow moderate complexity
            return max(2, int(4 * self.exploration_multiplier))
        elif self.change_percent < 0.7:
            # Medium perturbation: balanced complexity
            return max(2, int(3 * self.exploration_multiplier))
        else:
            # High perturbation: constrain complexity to prevent overfitting
            return max(2, int(2 * self.exploration_multiplier))
        
    def get_vns_config(self) -> Dict[str, Any]:
        """
        Generate VNS configuration parameters adapted for the current perturbation level.
        
        Key adaptations:
        - Increase exploration distance for higher perturbation levels
        - Adjust termination criteria to allow more thorough search
        - Configure local search parameters for robustness
        """
        base_maxeval = 1000
        base_maxtime = 60
        
        config = {
            "change_percent": self.change_percent,
            'maxeval': int(base_maxeval * self.exploration_multiplier),  # total function evaluations
            'maxtime': int(base_maxtime * self.exploration_multiplier),
            'use_local': 1,  # Always use local search
            'aggr': 1,  # Non-aggressive search for better exploration
            'local_search_type': 2,  # 1 = first improvement, 2 = best
            'decomp': 1,  # focus local search on perturbed vars
            'maxdist': min(1.0, self.exploration_multiplier),  # furthest neighborhood (0–1)
            'iterprint': 0,  # Reduce output for cleaner logs
        }
        
        return config
    
    def get_ga_config(self) -> Dict[str, Any]:
        """
        Generate GA configuration parameters adapted for the current perturbation level.
        
        Key adaptations:
        - Adjust size factor to maintain network size
        - Increase population size for better exploration
        - Modify termination criteria based on perturbation level
        """
        base_pop_size = 200
        base_max_gens = 1000
        
        config = {
            "change_percent": self.change_percent,
            # Population and generation parameters
            'popSize': int(base_pop_size * self.population_multiplier), # population size
            'maxGens': int(base_max_gens * self.exploration_multiplier), # maximum generations
            'stallGenMax': max(50, int(100 * self.exploration_multiplier)), # stop if no improvement (None = ∞)
            'elitism': None,  # number of elites (None = popSize//10)
            
            # Objective function parameters - critical for network size control
            'sizeFac': self.config.base_size_factor, # penalty weight for network size
            'NAFac': 1.0,  # penalty for unresolved states
            
            # Termination and search parameters
            'maxTime': max(60, int(60 * self.exploration_multiplier)),
            
            # Genetic algorithm specific parameters
            'pMutation': min(0.7, 0.5 + self.change_percent * 0.2),  # Increase mutation for higher perturbation
            'selPress': 1.2,  # Standard selection pressure
            'relTol': 0.1 * self.tolerance_multiplier,  # Increase tolerance for higher perturbation
            'verbose': "FALSE",  # Reduce output, we need pass to R script
        }
        
        # Calculate elitism based on population size
        config['elitism'] = max(5, int(config['popSize'] / 10))
        
        return config
    
    def get_ilp_config(self) -> Dict[str, Any]:
        """
        Generate ILP configuration parameters adapted for the current perturbation level.
        
        Key adaptations:
        - Adjust size factor to maintain network size
        - Modify gap parameters for solution quality vs. speed trade-off
        - Configure solution diversity parameters
        """
        num_solutions = 3
        based_time_limit = 3600  # Base time limit in seconds
        based_pop_size = 500
        
        if sys.platform.startswith("darwin"):
            cplex_path = "~/CPLEX_Studio2211/cplex/bin/x86-64_osx/cplex"
        elif sys.platform.startswith("linux"):
            cplex_path = "~/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex"
        else:
            cplex_path = "~/CPLEX_Studio2211/cplex/bin/x86-64_windows/cplex"

        config = {
            # Model size control - most important parameter
            'accountForModelSize': "TRUE",  # include size in objective
            'sizeFac': self.config.base_size_factor, # penalty weight for size

            # Solution quality parameters
            'mipGap': 0,  # absolute MIP gap tolerance, Exact solutions preferred
            'relGap': 0.05,  # relative gap tolerance, Increase tolerance for higher perturbation
            
            # Computational parameters
            'timelimit': max(based_time_limit, int(based_time_limit * self.exploration_multiplier)),

            'cplexPath': os.path.expanduser(cplex_path),
            'method': 'quadratic',  # Standard method
            
            # Solution diversity parameters
            'numSolutions': max(num_solutions, int(num_solutions * self.exploration_multiplier)), # how many solutions to retrieve
            'limitPop': max(based_pop_size, int(based_pop_size * self.population_multiplier)), # max solutions in the pool
            'poolIntensity': 0,  # pool management intensity
            'poolReplace': 2,  # pool replacement strategy
        }
        
        return config
    
    def get_all_configs(self) -> Dict[str, Dict[str, Any]]:
        """
        Generate all configuration dictionaries for all optimization methods.
        
        Returns a dictionary with keys 'caspo', 'vns', 'ga', 'ilp' containing
        their respective configuration parameters.
        """
        return {
            'caspo': self.get_caspo_config(),
            'vns': self.get_vns_config(),
            'ga': self.get_ga_config(),
            'ilp': self.get_ilp_config()
        }
    
    def save_config(self, filepath: str):
        """Save all configurations to a JSON file for reproducibility."""
        config_data = {
            'perturbation_level': self.change_percent,
            'adaptation_factors': {
                'size_factor_multiplier': self.size_factor_multiplier,
                'exploration_multiplier': self.exploration_multiplier,
                'tolerance_multiplier': self.tolerance_multiplier,
                'population_multiplier': self.population_multiplier,
            },
            'configurations': self.get_all_configs()
        }
        
        with open(filepath, 'w') as f:
            json.dump(config_data, f, indent=2)
    
    def print_adaptation_summary(self):
        """Print a summary of how parameters are being adapted."""
        print(f"Parameter Adaptation Summary for {self.change_percent:.1%} perturbation:")
        print(f"  Size Factor Multiplier: {self.size_factor_multiplier:.3f} (reduces size penalties)")
        print(f"  Exploration Multiplier: {self.exploration_multiplier:.3f} (increases search effort)")
        print(f"  Tolerance Multiplier: {self.tolerance_multiplier:.3f} (increases fitting tolerance)")
        print(f"  Population Multiplier: {self.population_multiplier:.3f} (increases population diversity)")
        print()


# Example usage and testing functions
def create_experiment_configs(perturbation_levels: Sequence[float], base_size_factor: float = 1e-4) -> Dict[float, AdaptiveParameterManager]:
    """
    Create configuration managers for multiple perturbation levels.
    
    This function is designed for your experimental setup where you test
    perturbation levels from 0.1 to 0.9.
    """
    configs = {}
    
    for level in perturbation_levels:
        config = NetworkPerturbationConfig(
            change_percent=level,
            base_size_factor=base_size_factor,
            size_adaptation_strength=2.0,
            generalization_focus=True
        )
        configs[level] = AdaptiveParameterManager(config)
    
    return configs


def demonstrate_usage():
    """
    Demonstrate how to use the adaptive parameter system in your experimental workflow.
    """
    print("=== Adaptive Parameter Configuration System Demo ===\n")
    
    # Define your experimental perturbation levels
    perturbation_levels = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]
    
    # Create configuration managers for each perturbation level
    experiment_configs = create_experiment_configs(perturbation_levels)
    
    # Show how parameters adapt across perturbation levels
    print("Parameter Adaptation Across Perturbation Levels:")
    print("=" * 50)
    
    for level in perturbation_levels:
        manager = experiment_configs[level]
        manager.print_adaptation_summary()
        
        # Show key parameter changes
        ga_config = manager.get_ga_config()
        caspo_config = manager.get_caspo_config()
        
        print(f"  GA sizeFac: {ga_config['sizeFac']:.2e}")
        print(f"  GA popSize: {ga_config['popSize']}")
        print(f"  CASPO size_tolerance: {caspo_config['size']}")
        print(f"  CASPO fit_tolerance: {caspo_config['fit']:.3f}")
        print("-" * 30)
    
    # Example: Using configurations in your experimental loop
    print("\nExample: Experimental Loop Integration")
    print("=" * 40)
    
    # Simulate your experimental workflow
    for perturbation_level in [0.1, 0.5, 0.9]:
        print(f"\nProcessing perturbation level: {perturbation_level:.1%}")
        
        # Get the appropriate configuration manager
        manager = experiment_configs[perturbation_level]
        
        # Get configurations for all methods
        all_configs = manager.get_all_configs()
        
        # Show how you would use these in your code
        print(f"  VNS maxeval: {all_configs['vns']['maxeval']}")
        print(f"  GA sizeFac: {all_configs['ga']['sizeFac']:.2e}")
        print(f"  ILP sizeFac: {all_configs['ilp']['sizeFac']:.2e}")
        
        # Save configuration for reproducibility
        config_filename = f"config_perturbation_{perturbation_level:.1f}.json"
        manager.save_config(config_filename)
        print(f"  Configuration saved to: {config_filename}")


if __name__ == "__main__":
    demonstrate_usage()
