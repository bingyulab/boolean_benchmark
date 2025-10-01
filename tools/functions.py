#!/usr/bin/env python3
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from collections import defaultdict
import re
import logging

def setup_logger(name="network_analysis", log_file="network_analysis.log", level=logging.INFO):
    """Setup logger with file and console handlers, preventing duplicates."""
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # If logger already has handlers, don't add more (prevents duplicates on re-runs)
    if logger.handlers:
        return logger
    
    class ColoredFormatter(logging.Formatter):
        RESET = "\033[0m"
        COLORS = {
            logging.DEBUG: "\033[36m",     # Cyan
            logging.INFO: "\033[32m",      # Green
            logging.WARNING: "\033[33m",   # Yellow
            logging.ERROR: "\033[31m",     # Red
            logging.CRITICAL: "\033[41m",  # Red background
        }
        
        def format(self, record):
            msg = super().format(record)
            color = self.COLORS.get(record.levelno, "")
            return f"{color}{msg}{self.RESET}" if color else msg
    
    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_formatter = ColoredFormatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(level)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Prevent propagation to root logger (avoids duplicate messages)
    logger.propagate = False
    
    return logger

def parse_namedlist_to_dataframes(named_list):
    """
    Parse rpy2 NamedList object containing attractor comparison results
    """
    results = {}
    
    # Get the names (keys) from the NamedList
    names = list(named_list.names)
    print(f"Found {len(names)} files: {names}")
    
    for name in names:
        print(f"\nProcessing: {name}")
        
        # Access the R DataFrame for this name
        r_dataframe = named_list[name]
        print(f"R DataFrame type: {type(r_dataframe)}")
        
        # Convert R DataFrame to pandas DataFrame
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.conversion.rpy2py(r_dataframe)
                results[name] = df
                print(f"Successfully converted {name} to pandas DataFrame")
                print(f"Shape: {df.shape}")
                print(f"Columns: {list(df.columns)}")
                print(f"Index: {list(df.index)}")
        except Exception as e:
            print(f"Conversion failed for {name}: {e}")
    
    return results


def parse_rpy2_results(r_list_vector):
    """
    Parse RPy2 ListVector containing attractor comparison results using pandas2ri
    """
    results = {}
    
    # Enable pandas conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        # Convert the R list to Python dict
        named_list = ro.conversion.rpy2py(r_list_vector)
        names = list(named_list.names())
        print(f"Found {len(names)} files: {names}")
        
        for i, name in enumerate(names):
            print(f"\nProcessing: {name}")
            
            # Access the R DataFrame for this name
            r_dataframe = named_list[i]
            print(f"R DataFrame type: {type(r_dataframe)}")
            
            df = ro.conversion.rpy2py(r_dataframe)
            results[name] = df
            print(f"Successfully converted {name} to pandas DataFrame")
    return results

# Convert R list to Python dict
def rlist_to_pydict(rlist):
    py_dict = {}
    for name in rlist.names:
        value = rlist.rx2(name)
        # Convert R vectors to Python lists
        if hasattr(value, 'tolist'):
            py_dict[name] = value.tolist()
        else:
            py_dict[name] = value
    return py_dict
    
    
def analyze_attractor_performance(df):
    """
    Analyze attractor performance metrics
    """
    print("\nPerformance Analysis:")
    print(f"Best accuracy: {df['accuracy'].max():.4f} (Attractor: {df['accuracy'].idxmax()})")
    print(f"Best F1 score: {df['f1'].max():.4f} (Attractor: {df['f1'].idxmax()})")
    print(f"Average hamming distance: {df['hamming'].mean():.2f}")
    print(f"Average normalized hamming: {df['norm_hamming'].mean():.4f}")
    
    # Count valid precision/recall values
    valid_precision = df['precision'].notna().sum()
    valid_f1 = df['f1'].notna().sum()
    print(f"Attractors with valid precision: {valid_precision}/{len(df)}")
    print(f"Attractors with valid F1: {valid_f1}/{len(df)}")
    

def caspo_to_boolnet(caspo_results, output_file=None, frequency_threshold=0.5):
    """
    Convert CASPO training results to BoolNet format.
    
    Parameters:
    - caspo_results: Can be either:
        1. A string with comma-separated rules (first format)
        2. A pandas DataFrame with columns: mapping, frequency, inclusive, exclusive
        3. A CSV file path containing the DataFrame format
    - output_file: Path to save the BoolNet file (optional)
    - frequency_threshold: Minimum frequency to include a rule (default: 0.5)
    
    Returns:
    - Boolean indicating success/failure
    - String containing the BoolNet content (if successful)
    """
    
    try:
        # Parse input based on format
        if isinstance(caspo_results, str):
            if caspo_results.endswith('.csv'):
                # File path
                df = pd.read_csv(caspo_results)
                rules = df[df['frequency'] >= frequency_threshold]['mapping'].tolist()
            elif ',' in caspo_results and '<-' in caspo_results:
                # Comma-separated rules format
                rules = [rule.strip() for rule in caspo_results.split(',')]
            else:
                raise ValueError("Invalid string format")
        elif isinstance(caspo_results, pd.DataFrame):
            # DataFrame format
            rules = caspo_results[caspo_results['frequency'] >= frequency_threshold]['mapping'].tolist()
        else:
            raise ValueError("Unsupported input format")
        
        # Group rules by target gene
        gene_rules = defaultdict(list)
        
        for rule in rules:
            if '<-' not in rule:
                continue
                
            target, regulators = rule.split('<-')
            target = target.strip()
            regulators = regulators.strip()
            
            # Convert regulators to Boolean format
            # Handle negation (!) and conjunction (+)
            boolean_expr = regulators.replace('+', ' & ').replace('!', '!')
            
            gene_rules[target].append(boolean_expr)
        
        # Generate BoolNet content
        boolnet_lines = ["targets, factors"]
        
        for target, expressions in gene_rules.items():
            if len(expressions) == 1:
                boolean_rule = expressions[0]
            else:
                # Multiple rules for same target are combined with OR
                boolean_rule = ' | '.join(f"{expr}" for expr in expressions)
            
            boolnet_lines.append(f"{target}, {boolean_rule}")
        
        # Add any genes that appear as regulators but not as targets
        all_genes = set(gene_rules.keys())
        for rule in rules:
            if '<-' in rule:
                regulators = rule.split('<-')[1].strip()
                # Extract gene names (remove operators)
                gene_names = re.findall(r'[A-Za-z0-9_]+', regulators)
                all_genes.update(gene_names)
        
        # Add missing genes as constants (they maintain their state)
        for gene in all_genes:
            if gene not in gene_rules:
                boolnet_lines.append(f"{gene}, {gene}")
        
        boolnet_content = '\n'.join(boolnet_lines)
        
        # Save to file if specified
        if output_file:
            with open(output_file, 'w') as f:
                f.write(boolnet_content)
            print(f"BoolNet file saved to: {output_file}")
        
        return True, boolnet_content
        
    except Exception as e:
        print(f"Error converting CASPO to BoolNet: {e}")
        return False, None

# Example usage functions for both formats
def convert_caspo_string_format(caspo_string, output_file=None):
    """Convert CASPO string format to BoolNet."""
    return caspo_to_boolnet(caspo_string, output_file)

def convert_caspo_csv_format(csv_file_or_df, output_file=None, frequency_threshold=0.5):
    """Convert CASPO CSV format to BoolNet."""
    return caspo_to_boolnet(csv_file_or_df, output_file, frequency_threshold)


def caspo_to_sif(caspo_results, output_file=None, frequency_threshold=0.5):
    """
    Convert CASPO training results to SIF format.

    Parameters:
    - caspo_results: Can be either:
        1. A string with comma-separated rules (first format)
        2. A pandas DataFrame with columns: mapping, frequency, inclusive, exclusive
        3. A CSV file path containing the DataFrame format
    - output_file: Path to save the SIF file (optional)
    - frequency_threshold: Minimum frequency to include a rule (default: 0.5)

    Returns:
    - Boolean indicating success/failure
    - String containing the SIF content (if successful)
    """
    try:
        # Load or parse input
        if isinstance(caspo_results, str):
            if caspo_results.endswith('.csv'):
                df = pd.read_csv(caspo_results)
                rules = df[df['frequency'] >= frequency_threshold]['mapping'].tolist()
            elif ',' in caspo_results and '<-' in caspo_results:
                rules = [rule.strip() for rule in caspo_results.split(',')]
            else:
                raise ValueError("Invalid string format for CASPO results")
        elif isinstance(caspo_results, pd.DataFrame):
            df = caspo_results
            rules = df[df['frequency'] >= frequency_threshold]['mapping'].tolist()
        else:
            raise ValueError("Unsupported input format for CASPO results")

        sif_lines = []

        for rule in rules:
            if '<-' not in rule:
                continue

            target, regulators = rule.split('<-')
            target = target.strip()
            # split on '+' but keep '!' attached
            for part in regulators.split('+'):
                part = part.strip()
                sign = +1
                if part.startswith('!'):
                    sign = -1
                    regulator = part[1:]
                else:
                    regulator = part
                sif_lines.append(f"{regulator}\t{sign}\t{target}")

        sif_content = "\n".join(sif_lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(sif_content)
            print(f"SIF file saved to: {output_file}")

        return True, sif_content

    except Exception as e:
        print(f"Error converting CASPO to SIF: {e}")
        return False, None
