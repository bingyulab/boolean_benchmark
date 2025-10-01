import pyboolnet.file_exchange as FileExchange
import pyboolnet.basins_of_attraction as Basins
from pyboolnet.trap_spaces import compute_steady_states
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from rpy2.robjects.vectors import FloatVector
import logging


def limit_float(x, nbit=4):
    if isinstance(x, FloatVector):
        x = float(x[0])
    else:
        x = float(x)
    s = str(x)
    if '.' in s and len(s.split('.')[1]) > nbit:
        return round(x, nbit)
    return x


def get_regulators(node_rules):
    """Extract set of regulators from prime rules of a node."""
    regs = set()
    for conjuncts in node_rules:
        for term in conjuncts:
            regs.update(term.keys())
    return regs

def jaccard_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Calculate Jaccard similarity between two vectors.
    NaN values are treated as non-matching positions.
    
    Jaccard = |A ∩ B| / |A ∪ B|
    For binary vectors: intersection = both are 1, union = at least one is 1
    """
    
    # For binary vectors, calculate intersection and union
    intersection = np.sum((vec1 == 1) & (vec2 == 1))
    union = len(vec1) - np.sum((vec1 == 0) | (vec2 == 0))

    if union == 0:
        return 1.0  # Both vectors are all zeros in valid positions
    
    return intersection / union

def levenshtein_distance(s1, s2):
    """
    Compute Levenshtein (edit) distance between two sequences.
    Handles NaN values by treating them as a special character.
    """
    # Convert to strings, treating NaN as special symbol
    str1 = ''.join([ 'N' if np.isnan(x) else str(int(x)) for x in s1])
    str2 = ''.join([str(int(x)) for x in s2 if not np.isnan(x)])

    m, n = len(str1), len(str2)
    
    # Create DP matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i-1] == str2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],    # deletion
                    dp[i][j-1],    # insertion
                    dp[i-1][j-1]   # substitution
                )
    return dp[m][n]


def hamming_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """
    Calculate Hamming similarity between two vectors.
    NaN values are treated as non-matching positions.
    Only considers positions that are valid in both vectors.
    
    Similarity = (total_positions - hamming_distance) / total_positions
    """
    if len(vec1) != len(vec2):
        raise ValueError("Vectors must have the same length for Hamming similarity")
    # Calculate matches
    matches = np.sum(vec1 == vec2)
    return matches / len(vec1)


def longest_common_subsequence(s1, s2):
    """
    Compute Longest Common Subsequence length between two sequences.
    Handles NaN values by treating them as a special character.
    """
    # Convert to strings, treating NaN as special symbol
    str1 = ''.join([ 'N' if np.isnan(x) else str(int(x)) for x in s1])
    str2 = ''.join([str(int(x)) for x in s2 if not np.isnan(x)])
    
    m, n = len(str1), len(str2)
    
    # Create DP matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i-1] == str2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    return dp[m][n]

class AttractorAnalysis:
    def __init__(self, ori_bnet: str, compared_bnet):
        """
        Initialize AttractorAnalysis with original and compared Boolean network files.

        Parameters:
        -----------
        ori_bnet : str
            Path to the original Boolean network file.
        compared_bnet : str or list of str
            Path(s) to the compared Boolean network file(s).
        """
        self.ori_bnet = ori_bnet
        self.compared_bnet = compared_bnet

    @staticmethod
    def load_bnet(bnet_file):
        primes = FileExchange.bnet2primes(bnet_file)
        # primes = AttractorAnalysis.refine_bnet(primes)
        return primes

    @staticmethod
    def refine_bnet(prime):
        # Filter out nodes with only self as regulator
        cleaned_prime = {}
        for node, rules in prime.items():
            regulators = get_regulators(rules)
            if regulators == {node}:  # only self-loop
                # check if node appears in other nodes' regulators
                appears_elsewhere = any(
                    node in get_regulators(other_rules)
                    for other_node, other_rules in prime.items()
                    if other_node != node
                )
                if appears_elsewhere:
                    cleaned_prime[node] = rules  # keep because it's used elsewhere
                else:
                    print(f"Removing isolated self-loop node: {node}")
            else:
                cleaned_prime[node] = rules
        return cleaned_prime
    
    @staticmethod
    def compute_attractors(primes):  
        logging.disable(logging.INFO)
        # result = Attractors.compute_attractors(primes, "synchronous") 
        result = compute_steady_states(primes)
        logging.disable(logging.NOTSET)
        return result

    def get_attractors(self, bnet_file):
        primes = AttractorAnalysis.load_bnet(bnet_file)
        attrs = AttractorAnalysis.compute_attractors(primes)
        # return primes, [x['state'] for x in attrs['attractors']]
        return primes, attrs
    
    @staticmethod
    def get_basin_sizes(primes, state):   
        logging.disable(logging.INFO)
        weak = Basins.weak_basin(primes, "asynchronous", state)
        logging.disable(logging.NOTSET)
        return weak.get('perc', 0.0)
    
    def compare_attractors(self, ori_primes, ori_attrs, com_primes, comp_attrs):
        comparator = AttractorComparison(ori_primes, ori_attrs, com_primes, comp_attrs)
        results = comparator.comprehensive_comparison()
        return results
    
    def compare_multiple_attractors(self, ori_primes, ori_attrs, com_primes, comp_attrs):
        results = []
        for i, reconstructed_attractors in enumerate(comp_attrs):
            comparator = AttractorComparison(ori_primes, ori_attrs, com_primes[i], reconstructed_attractors)
            comparison_result = comparator.comprehensive_comparison()
            comparison_result['reconstruction_id'] = i
            results.append(comparison_result)                    
        return results
    
    def comparison(self):
        ori_primes, ori_attrs = self.get_attractors(self.ori_bnet)
        if isinstance(self.compared_bnet, str):
            comp_primes, comp_attrs = self.get_attractors(self.compared_bnet)
            results = self.compare_attractors(ori_primes, ori_attrs, comp_primes, comp_attrs)
        else:
            tmp_res = [self.get_attractors(bnet) for bnet in self.compared_bnet]
            com_primes = [x[0] for x in tmp_res]
            comp_attrs = [x[1] for x in tmp_res]
            results = self.compare_multiple_attractors(ori_primes, ori_attrs, com_primes, comp_attrs)
        return results
            

class AttractorComparison:
    """
    A comprehensive toolkit for comparing attractors from different network models,
    handling variable network sizes and attractor counts.
    """

    def __init__(self, ori_primes, original_attractors, recon_primes, 
                 reconstructed_attractors, basin_weight=0.3, strategy="random"):
        self.ori_primes = ori_primes
        self.recon_primes = recon_primes
        self.basin_weight = basin_weight  # Weight for basin size consideration        
        
        self.recon_total = len(reconstructed_attractors)
        self.original_attractors = original_attractors
        if strategy is None:
            self.reconstructed_attractors = reconstructed_attractors
        else:
            self.reconstructed_attractors = self._select_attractors_for_comparison(
                reconstructed_attractors, len(original_attractors), strategy
            )
            
        # Get node sets
        self.orig_nodes = self._get_nodes_from_attractors(self.original_attractors)
        self.recon_nodes = self._get_nodes_from_attractors(self.reconstructed_attractors)
        self.common_nodes = self.orig_nodes.intersection(self.recon_nodes)
        self.all_nodes = sorted(list(self.orig_nodes.union(self.recon_nodes)))
        # Build enhanced representations
        self._build_enhanced_representations() 
        
    def _compute_basin_sizes(self, primes, attractors):
        weight = [AttractorAnalysis.get_basin_sizes(primes, att['attr']) for att in attractors]        
        weights = np.array(weight)
        # Normalize to prevent scale issues
        return weights / np.sum(weights) if np.sum(weights) > 0 else weights
    
    def _select_attractors_for_comparison(self, attractors, max_count, strategy):
        if len(attractors) <= max_count:
            return attractors
        top_indices = []
        if strategy == "top_k":
            # Select top K attractors by basin size
            basin_sizes = self._compute_basin_sizes(self.recon_primes, attractors)
            #  Normalize basin sizes to [0,1] to prevent overflow
            basin_sizes = basin_sizes / np.max(basin_sizes) if np.max(basin_sizes) > 0 else basin_sizes
        # Get indices of top K attractors
            top_indices = np.argsort(basin_sizes)[::-1][:max_count]
        elif strategy == "random":
            # Select random attractors
            top_indices = np.random.choice(len(attractors), size=max_count, replace=False)
        return [attractors[i] for i in top_indices]  
    
    def _get_nodes_from_attractors(self, attractors):
        """Get all unique nodes from attractors."""
        if not attractors:
            return set()
        nodes = set()
        for att in attractors:
            nodes.update(att.keys())
        return nodes
    
    def _build_enhanced_representations(self):
        """
        Build enhanced representations that handle missing nodes properly.
        Uses three-state representation: 0, 1, NaN (missing)
        """
        self.orig_enhanced = []
        self.recon_enhanced = []
        
        # Build for original attractors
        for att in self.original_attractors:
            vector = []
            for node in self.all_nodes:
                if node in att:
                    vector.append(float(att[node]))
                else:
                    vector.append(np.nan)  # Missing node
            self.orig_enhanced.append(vector)
        
        # Build for reconstructed attractors
        for att in self.reconstructed_attractors:
            vector = []
            for node in self.all_nodes:
                if node in att:
                    vector.append(float(att[node]))
                else:
                    vector.append(np.nan)  # Missing node
            self.recon_enhanced.append(vector)
        
        self.orig_enhanced = np.array(self.orig_enhanced)
        self.recon_enhanced = np.array(self.recon_enhanced)

    def _compute_pairwise_similarity(self, vec1, vec2, similarity_type='jaccard'):
        """
        Compute similarity between two vectors handling NaN values properly.
        """        
        
        if similarity_type == 'jaccard':
            return jaccard_similarity(vec1, vec2)

        elif similarity_type == 'hamming':
            return hamming_similarity(vec1, vec2)
        
        elif similarity_type == 'levenshtein':
            # Levenshtein similarity (1 - edit_distance/max_length)
            edit_dist = levenshtein_distance(vec1, vec2)
            max_len = len(vec1)
            return 1 - (edit_dist / max_len) if max_len > 0 else 0.0
        
        elif similarity_type == 'lcs':
            # LCS similarity (lcs_length / max_length)
            lcs_len = longest_common_subsequence(vec1, vec2)
            max_len = max(len(vec1), len(vec2))
            return lcs_len / max_len if max_len > 0 else 0.0

    def _compute_similarity_matrix(self, similarity_type='jaccard'):
        """Compute similarity matrix between all pairs of attractors."""
        n_orig = len(self.orig_enhanced)
        n_recon = len(self.recon_enhanced)
        
        if n_orig == 0 or n_recon == 0:
            return np.array([]).reshape(n_orig, n_recon)
        
        sim_matrix = np.zeros((n_orig, n_recon))
        
        for i in range(n_orig):
            for j in range(n_recon):
                # print(f"Computing {similarity_type} similarity between original attractor {i} and reconstructed attractor {j}")
                sim_matrix[i, j] = self._compute_pairwise_similarity(
                    self.orig_enhanced[i], self.recon_enhanced[j], similarity_type)
        
        return sim_matrix

    def compute_optimal_matching_metrics(self):
        """
        Compute metrics based on optimal bipartite matching.
        This addresses the third problem by finding best matches first.
        """        
        # Compute similarity matrices for different metrics
        jaccard_matrix = self._compute_similarity_matrix('jaccard')
        hamming_matrix = self._compute_similarity_matrix('hamming')
        # lcs_matrix = self._compute_similarity_matrix('lcs')
        # levenshtein_matrix = self._compute_similarity_matrix('levenshtein')

        n_orig = len(self.original_attractors)
        n_recon = len(self.reconstructed_attractors)
        
        if n_orig == 0 or n_recon == 0:
            return self._empty_results()
        
        # Use Jaccard for matching (most robust for categorical data)
        cost_matrix = 1 - jaccard_matrix
        
        # Handle rectangular matrices for Hungarian algorithm
        if n_orig != n_recon:
            # Pad with high costs to handle different sizes
            max_dim = max(n_orig, n_recon)
            padded_cost = np.full((max_dim, max_dim), 2.0)  # Cost > 1
            padded_cost[:n_orig, :n_recon] = cost_matrix
            row_ind, col_ind = linear_sum_assignment(padded_cost)
            
            # Filter out padding matches
            valid_matches = [(i, j) for i, j in zip(row_ind, col_ind) 
                           if i < n_orig and j < n_recon]
        else:
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            valid_matches = list(zip(row_ind, col_ind))
        
        # Compute metrics for matched pairs
        matched_jaccard = []
        matched_hamming = []
        # matched_lcs = []
        # matched_levenshtein = []

        for i, j in valid_matches:
            matched_jaccard.append(jaccard_matrix[i, j])
            matched_hamming.append(hamming_matrix[i, j])
            # matched_lcs.append(lcs_matrix[i, j])
            # matched_levenshtein.append(levenshtein_matrix[i, j])

        true_positives = 0
        false_positives = 0
        false_negatives = 0

        for i, j in valid_matches:
            vec_true = self.orig_enhanced[i]
            vec_pred = self.recon_enhanced[j]
            
            # Iterate over all nodes in the union space
            for node_idx in range(len(self.all_nodes)):
                true_val = vec_true[node_idx]
                pred_val = vec_pred[node_idx]
                
                # Valid comparison is possible only if the original value is not NaN
                if not np.isnan(true_val):
                    # True Positive: Both are 1 and match
                    if true_val == 1 and pred_val == 1:
                        true_positives += 1
                    # False Positive: Original is 0, but prediction is 1
                    elif true_val == 0 and pred_val == 1:
                        false_positives += 1
                    # False Negative: Original is 1, but prediction is 0 or NaN (missing)
                    elif true_val == 1 and (pred_val == 0 or np.isnan(pred_val)):
                        false_negatives += 1

        # Calculate precision, recall, and F1 score from the counts
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0.0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0.0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        return {
            'jaccard': np.mean(matched_jaccard) if matched_jaccard else 0.0,
            'hamming': np.mean(matched_hamming) if matched_hamming else 0.0,
            # 'lcs': np.mean(matched_lcs) if matched_lcs else 0.0,
            # 'levenshtein': np.mean(matched_levenshtein) if matched_levenshtein else 0.0,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'total_matches': len(valid_matches),
            'orig_count': n_orig,
            'recon_count': n_recon,
            'recon_total': self.recon_total
        }
    
    def _empty_results(self):
        """Return empty results when comparison is not possible."""
        return {
            'jaccard': 0.0, 'hamming': 0.0, 
            # 'levenshtein': 0.0, 'lcs': 0.0,
            'precision': 0.0, 'recall': 0.0, 'f1_score': 0.0,
            'total_matches': 0,
            'orig_count': len(self.original_attractors),
            'recon_count': len(self.reconstructed_attractors),
            'recon_total': self.recon_total
        }
    
    def comprehensive_comparison(self, return_df=True):
        """
        Perform comprehensive comparison with improved metrics.
        """
        # Optimal matching metrics
        matching_results = self.compute_optimal_matching_metrics()
        
        
        # Node coverage analysis
        node_results = {
            'common_nodes': len(self.common_nodes),
            'orig_nodes': len(self.orig_nodes),
            'recon_nodes': len(self.recon_nodes),
            'node_coverage': len(self.common_nodes) / len(self.orig_nodes) if self.orig_nodes else 0.0
        }
        
        # Combine all results
        results = {**matching_results, **node_results}

        # Compute composite score with balanced weighting
        composite_score = (
            0.35 * results['jaccard'] +
            0.35 * results['hamming'] +
            0.3 * results['f1_score']
        )
        
        results['composite_score'] = composite_score
        
        if return_df:
            return pd.DataFrame([results])
        return results     
    
    
if __name__ == "__main__":
    # ori_primes = "data/ToyModel/ToyModel.bnet"
    ori_primes = "output/cellnopt/ToyModel/0_Modified/ga/OPT_ToyModel.bnet"
    # ori_primes = "output/meigo/ToyModel/0_Modified/VNS/OPT_ToyModel.bnet"
    # ori_primes = "output/caspo/ToyModel/0_Modified/OPT_ToyModel.bnet"
    # compared_bnet = "data/ToyModel/ToyModel.bnet"
    compared_bnet = "output/cellnopt/ToyModel/80_Modified/ga/OPT_ToyModel.bnet"
    # compared_bnet = "output/caspo/ToyModel/50_Modified/OPT_ToyModel.bnet"
    # compared_bnet = "output/caspo/ToyModel/10_Modified/OPT_ToyModel.bnet"
    # The code you provided is a Python code snippet that seems to be commented out. It appears to be
    # assigning a file path to a variable `ori_primes`. The file path seems to be related to a
    # Bayesian network file in the context of a cell signaling network model. However, since the code
    # is commented out, it is not actively doing anything in the program.
    # compared_bnet = "output/meigo/ToyModel/60_Modified/VNS/OPT_ToyModel.bnet"
    analysis = AttractorAnalysis(ori_primes, compared_bnet)
    results = analysis.comparison()
    print(results)
    
    results_dict = results.iloc[0].to_dict()

    print(f"\nSimilarity Metrics Comparison:")
    print(f"Jaccard: {results_dict['jaccard']:.3f}")
    print(f"Hamming: {results_dict['hamming']:.3f}")
    # print(f"LCS: {results_dict['lcs']:.3f}")
    # print(f"Levenshtein: {results_dict['levenshtein']:.3f}")