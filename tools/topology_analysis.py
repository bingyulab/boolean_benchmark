import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from tools.comparison import limit_float


class BooleanNetworkGraph(nx.MultiDiGraph):
    """
    Boolean network graph extending networkx.MultiDiGraph for topology analysis
    Built upon the caspo Graph class structure
    """

    @classmethod
    def from_tuples(cls, tuples):
        """
        Creates a graph from tuples (source, target, sign)        
        Parameters
        ----------
        tuples : iterable[(str,str,int)]
            Tuples describing signed directed edges            
        Returns
        -------
        BooleanNetworkGraph
            Created graph instance
        """
        return cls((source, target, {'sign': sign}) for source, target, sign in tuples)

    @classmethod
    def read_sif(cls, path):
        """
        Creates a graph from SIF (Simple Interaction Format) file
        
        The SIF format expects: source interaction_type target
        We map interaction types to signs: activation=1, inhibition=-1
        
        Parameters
        ----------
        path : str
            Path to SIF file
            
        Returns
        -------
        BooleanNetworkGraph
            Created graph instance
        """
        df = pd.read_csv(path, sep=r'\s+', names=['source', 'interaction', 'target']).drop_duplicates()
        
        # Map interaction types to signs
        sign_mapping = {
            'activates': 1, 'activation': 1, '+': 1, '1': 1,
            'inhibits': -1, 'inhibition': -1, '-': -1, '-1': -1,
            'unknown': 0, '0': 0
        }
        
        # Convert interaction types to numerical signs
        df['sign'] = df['interaction'].map(lambda x: sign_mapping.get(str(x).lower(), 0))
        
        edges = [(source, target, {'sign': sign, 'interaction': interaction}) 
                for _, source, interaction, target, sign in df.itertuples()]
        
        return cls(edges)

    def to_simple_graph(self):
        """
        Convert to simple directed graph (no multiple edges, no signs)
        Useful for topology-only analysis
        
        Returns
        -------
        nx.DiGraph
            Simple directed graph
        """
        simple_graph = nx.DiGraph()
        simple_graph.add_nodes_from(self.nodes())
        
        # Add edges without considering signs or multiple connections
        for source, target in self.edges():
            if not simple_graph.has_edge(source, target):
                simple_graph.add_edge(source, target)
                
        return simple_graph

    def get_adjacency_matrix(self):
        """
        Get adjacency matrix representation
        
        Returns
        -------
        numpy.ndarray
            Adjacency matrix
        """
        simple_graph = self.to_simple_graph()
        return nx.adjacency_matrix(simple_graph).toarray()

    def compute_centralities(self):
        """
        Compute various centrality measures for network analysis
        
        Returns
        -------
        dict
            Dictionary containing different centrality measures
        """
        simple_graph = self.to_simple_graph()
        
        centralities = {
            'degree': dict(simple_graph.degree()),
            'in_degree': dict(simple_graph.in_degree()),
            'out_degree': dict(simple_graph.out_degree()),
            'betweenness': nx.betweenness_centrality(simple_graph),
            'closeness': nx.closeness_centrality(simple_graph),
            # 'eigenvector': nx.eigenvector_centrality(simple_graph, max_iter=1000),
            'pagerank': nx.pagerank(simple_graph)
        }
        
        return centralities

    def get_topology_features(self):
        """
        Extract comprehensive topological features
        
        Returns
        -------
        dict
            Dictionary of topological features
        """
        simple_graph = self.to_simple_graph()
        
        features = {
            'num_nodes': simple_graph.number_of_nodes(),
            'num_edges': simple_graph.number_of_edges(),
            'density': nx.density(simple_graph),
            # 'is_connected': nx.is_weakly_connected(simple_graph),
            'num_strongly_connected_components': nx.number_strongly_connected_components(simple_graph),
            'num_weakly_connected_components': nx.number_weakly_connected_components(simple_graph),
            # 'average_clustering': nx.average_clustering(simple_graph.to_undirected()),
            'transitivity': nx.transitivity(simple_graph.to_undirected())
        }
        
        # Degree distribution statistics
        degrees = [d for n, d in simple_graph.degree()]
        if degrees:
            features.update({
                'avg_degree': np.mean(degrees),
                'std_degree': np.std(degrees),
                'max_degree': max(degrees),
                'min_degree': min(degrees)
            })        
        return features
    

class NetworkTopologyAnalyzer:
    """
    Class for comparing topological similarity between boolean networks
    """
    
    def __init__(self, network1: BooleanNetworkGraph, network2: BooleanNetworkGraph):
        """
        Initialize analyzer with two networks
        
        Parameters
        ----------
        network1, network2 : BooleanNetworkGraph
            Networks to compare
        """
        self.network1 = network1
        self.network2 = network2
        self.simple1 = network1.to_simple_graph()
        self.simple2 = network2.to_simple_graph()
        
    def jaccard_similarity(self):
        """
        Calculate Jaccard similarity of edges
        
        This measures what fraction of edges are shared between networks.
        Higher values (closer to 1) indicate more similar topologies.
        
        Returns
        -------
        float
            Jaccard similarity coefficient (0-1)
        """
        edges1 = set(self.simple1.edges())
        edges2 = set(self.simple2.edges())
        
        intersection = len(edges1.intersection(edges2))
        union = len(edges1.union(edges2))

        return limit_float(intersection / union) if union > 0 else 0.0

    def node_overlap_similarity(self):
        """
        Calculate node overlap similarity
        
        Returns
        -------
        float
            Node Jaccard similarity (0-1)
        """
        nodes1 = set(self.network1.nodes())
        nodes2 = set(self.network2.nodes())
        
        intersection = len(nodes1.intersection(nodes2))
        union = len(nodes1.union(nodes2))
        
        return intersection / union if union > 0 else 0.0

    def degree_distribution_similarity(self):
        """
        Compare degree distributions using Kolmogorov-Smirnov test
        
        This statistical test tells us if the degree distributions
        come from the same underlying distribution.
        
        Returns
        -------
        dict
            KS test results with statistic and p-value
        """
        degrees1 = [d for n, d in self.simple1.degree()]
        degrees2 = [d for n, d in self.simple2.degree()]
        
        if not degrees1 or not degrees2:
            return {'statistic': 1.0, 'pvalue': 0.0, 'similar': False}
        
        ks_stat, p_value = ks_2samp(degrees1, degrees2)
        
        return {
            'statistic': ks_stat,  # Lower values indicate more similar distributions
            'pvalue': p_value,     # Higher p-values indicate similar distributions
            'similar': p_value > 0.05  # Common significance threshold
        }

    def centrality_correlation(self):
        """
        Calculate correlations between centrality measures of common nodes
        
        Returns
        -------
        dict
            Correlation coefficients for different centrality measures
        """
        # Get common nodes
        common_nodes = set(self.network1.nodes()).intersection(set(self.network2.nodes()))
        
        if len(common_nodes) < 2:
            return {}
        
        # Compute centralities
        cent1 = self.network1.compute_centralities()
        cent2 = self.network2.compute_centralities()
        
        correlations = {}
        
        for centrality_type in cent1.keys():
            values1 = [cent1[centrality_type][node] for node in common_nodes 
                      if node in cent1[centrality_type] and node in cent2[centrality_type]]
            values2 = [cent2[centrality_type][node] for node in common_nodes 
                      if node in cent1[centrality_type] and node in cent2[centrality_type]]
            
            if len(values1) > 1 and len(values2) > 1:
                correlation = np.corrcoef(values1, values2)[0, 1]
                correlations[centrality_type] = correlation if not np.isnan(correlation) else 0.0
        
        return correlations

    def graph_edit_distance(self, normalized=True):
        """
        Calculate graph edit distance (approximate)
        
        This measures the minimum number of operations needed to transform
        one graph into another. Lower values indicate more similar graphs.
        
        Parameters
        ----------
        normalized : bool
            Whether to normalize by graph size
            
        Returns
        -------
        float
            Graph edit distance
        """
        try:
            # Use NetworkX's graph edit distance (this can be slow for large graphs)
            distance = nx.graph_edit_distance(self.simple1, self.simple2, timeout=30)
            
            if normalized and distance is not None:
                max_size = max(self.simple1.number_of_nodes(), self.simple2.number_of_nodes())
                distance = distance / max_size if max_size > 0 else 0
                
            return distance if distance is not None else float('inf')
        except:
            # Fallback: approximate based on edge and node differences
            nodes1, nodes2 = set(self.simple1.nodes()), set(self.simple2.nodes())
            edges1, edges2 = set(self.simple1.edges()), set(self.simple2.edges())
            
            node_ops = len(nodes1.symmetric_difference(nodes2))
            edge_ops = len(edges1.symmetric_difference(edges2))
            
            distance = node_ops + edge_ops
            
            if normalized:
                max_size = max(len(nodes1), len(nodes2), len(edges1), len(edges2))
                distance = distance / max_size if max_size > 0 else 0
                
            return distance

    def comprehensive_similarity_report(self):
        """
        Generate comprehensive similarity analysis
        
        Returns
        -------
        dict
            Complete similarity analysis results
        """
        report = {
            'basic_metrics': {
                'jaccard_edge_similarity': self.jaccard_similarity(),
                'node_overlap_similarity': self.node_overlap_similarity(),
                # 'graph_edit_distance': self.graph_edit_distance(normalized=True)
            },
            'degree_distribution': self.degree_distribution_similarity(),
            'centrality_correlations': self.centrality_correlation(),
            'network1_features': self.network1.get_topology_features(),
            'network2_features': self.network2.get_topology_features()
        }

        return report

    def visualize_comparison(self, save_path=None):
        """
        Create visualization comparing the two networks
        
        Parameters
        ----------
        save_path : str, optional
            Path to save the visualization
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Network 1 structure
        pos1 = nx.spring_layout(self.simple1, seed=42)
        nx.draw(self.simple1, pos=pos1, ax=axes[0,0], with_labels=True, 
                node_color='lightblue', node_size=500, font_size=8)
        axes[0,0].set_title("Network 1 Topology")
        
        # Network 2 structure  
        pos2 = nx.spring_layout(self.simple2, seed=42)
        nx.draw(self.simple2, pos=pos2, ax=axes[0,1], with_labels=True,
                node_color='lightcoral', node_size=500, font_size=8)
        axes[0,1].set_title("Network 2 Topology")
        
        # Degree distributions
        degrees1 = [d for n, d in self.simple1.degree()]
        degrees2 = [d for n, d in self.simple2.degree()]
        
        axes[1,0].hist([degrees1, degrees2], bins=max(5, max(max(degrees1 + [0]), max(degrees2 + [0]))), 
                      alpha=0.7, label=['Network 1', 'Network 2'])
        axes[1,0].set_xlabel('Degree')
        axes[1,0].set_ylabel('Frequency')
        axes[1,0].set_title('Degree Distribution Comparison')
        axes[1,0].legend()
        
        # Centrality comparison (if common nodes exist)
        common_nodes = set(self.network1.nodes()).intersection(set(self.network2.nodes()))
        if common_nodes:
            cent1 = self.network1.compute_centralities()
            cent2 = self.network2.compute_centralities()
            
            # Plot degree centrality correlation
            common_list = list(common_nodes)
            deg1 = [cent1['degree'][node] for node in common_list if node in cent1['degree']]
            deg2 = [cent2['degree'][node] for node in common_list if node in cent2['degree']]
            
            if deg1 and deg2:
                axes[1,1].scatter(deg1, deg2, alpha=0.7)
                axes[1,1].plot([min(deg1 + deg2), max(deg1 + deg2)], 
                              [min(deg1 + deg2), max(deg1 + deg2)], 'r--', alpha=0.5)
                axes[1,1].set_xlabel('Network 1 Degree')
                axes[1,1].set_ylabel('Network 2 Degree')
                axes[1,1].set_title('Degree Centrality Correlation')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


if __name__ == "__main__":
    
    print("\n=== Reading from SIF Files ===")
    # Test reading from SIF files
    # sif_net1 = BooleanNetworkGraph.read_sif('data/ToyModel/ToyModel.sif')
    # sif_net1 = BooleanNetworkGraph.read_sif('output/caspo/ToyModel/0_Modified/OPT_ToyModel.sif')
    sif_net1 = BooleanNetworkGraph.read_sif('output/meigo/ToyModel/0_Modified/VNS/OPT_ToyModel.sif')
    # sif_net1 = BooleanNetworkGraph.read_sif('data/DREAMmodel/DreamModel.sif')
    # sif_net2 = BooleanNetworkGraph.read_sif('data/DREAMmodel/30_Modified/DreamModel.sif')
    # sif_net2 = BooleanNetworkGraph.read_sif('output/caspo/ToyModel/60_Modified/OPT_ToyModel.sif')
    sif_net2 = BooleanNetworkGraph.read_sif('output/meigo/ToyModel/80_Modified/VNS/OPT_ToyModel.sif')
    print(f"SIF Network 1: {sif_net1.number_of_nodes()} nodes, {sif_net1.number_of_edges()} edges")
    print(f"SIF Network 2: {sif_net2.number_of_nodes()} nodes, {sif_net2.number_of_edges()} edges")
    
    # Quick similarity check
    sif_analyzer = NetworkTopologyAnalyzer(sif_net1, sif_net2)
    jaccard = sif_analyzer.jaccard_similarity()
    print(f"SIF networks Jaccard similarity: {jaccard:.3f}")
    
    report = sif_analyzer.comprehensive_similarity_report()
    
    print("\n=== Topology Similarity Analysis ===")
    print(f"Jaccard Edge Similarity: {report['basic_metrics']['jaccard_edge_similarity']:.3f}")
    print(f"Node Overlap Similarity: {report['basic_metrics']['node_overlap_similarity']:.3f}")
    
    print(f"\nDegree Distribution Similar: {report['degree_distribution']['similar']}")
    print(f"Degree Distribution p-value: {report['degree_distribution']['pvalue']:.3f}")
    
    print(f"\nCentrality Correlations:")
    for centrality, correlation in report['centrality_correlations'].items():
        print(f"  {centrality}: {correlation:.3f}")
    
    # Show individual network features
    print(f"\n=== Network 1 Features ===")
    for feature, value in report['network1_features'].items():
        print(f"{feature}: {value}")
    
    print(f"\n=== Network 2 Features ===")
    for feature, value in report['network2_features'].items():
        print(f"{feature}: {value}")
        
    # Demonstrate centrality analysis
    print(f"\n=== Centrality Analysis Example ===")
    centralities1 = sif_net1.compute_centralities()
    print("Network 1 degree centralities:")
    for node, degree in centralities1['degree'].items():
        print(f"  {node}: {degree}")
    
    # Create visualization
    print(f"\n=== Creating Visualization ===")
    sif_analyzer.visualize_comparison("output/visualization_output.png")