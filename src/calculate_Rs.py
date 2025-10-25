# src/calculate_Rs.py

import networkx as nx
import numpy as np

def calculate_laplacian_pseudoinverse(G):
    """
    Calculates the Moore-Penrose pseudoinverse of the graph's Laplacian matrix.
    This is a key step for calculating effective resistance.

    Args:
        G (nx.Graph): A connected NetworkX graph.

    Returns:
        np.ndarray: The pseudoinverse of the Laplacian matrix.
    """
    if not nx.is_connected(G):
        raise ValueError("Graph must be connected to calculate its Laplacian pseudoinverse.")
    
    # Get the graph Laplacian matrix as a dense NumPy array
    L = nx.laplacian_matrix(G).toarray()
    
    # Calculate the Moore-Penrose pseudoinverse
    L_plus = np.linalg.pinv(L)
    
    return L_plus

def calculate_all_pairs_effective_resistance(G, L_plus=None):
    """
    Calculates the effective resistance (Omega_uv) for all edges in the graph.
    The formula used is: Omega_uv = L+_uu + L+_vv - 2*L+_uv

    Args:
        G (nx.Graph): A connected NetworkX graph.
        L_plus (np.ndarray, optional): Pre-computed pseudoinverse of the Laplacian.
                                       If None, it will be computed automatically.

    Returns:
        dict: A dictionary where keys are edge tuples (u, v) and values are
              the effective resistance Omega_uv.
    """
    if L_plus is None:
        L_plus = calculate_laplacian_pseudoinverse(G)
        
    # Create a mapping from node labels to matrix indices
    nodes = list(G.nodes())
    node_map = {node: i for i, node in enumerate(nodes)}
    
    resistances = {}
    for u, v in G.edges():
        i, j = node_map[u], node_map[v]
        omega_uv = L_plus[i, i] + L_plus[j, j] - 2 * L_plus[i, j]
        # Resistance cannot be negative; floating point errors can make it tiny negative
        resistances[(u, v)] = max(0, omega_uv)
        
    return resistances

def calculate_structural_penalty_weights(G):
    """
    Calculates the structural heterogeneity penalty factor (w_uv) for each edge.
    Formula: w_uv = 1 - (|d_u - d_bar| + |d_v - d_bar|) / (2 * N * d_bar)

    Args:
        G (nx.Graph): A NetworkX graph.

    Returns:
        dict: A dictionary where keys are edge tuples (u, v) and values are
              the structural penalty weights w_uv.
    """
    N = G.number_of_nodes()
    if N == 0 or G.number_of_edges() == 0:
        return {} # Return empty dict for empty graph
        
    degrees = dict(G.degree())
    d_bar = sum(degrees.values()) / N # Average degree
    
    if d_bar == 0: # Handle graph with no edges
        return {edge: 1.0 for edge in G.edges()}

    weights = {}
    denominator = 2 * N * d_bar
    if denominator == 0: # Should not happen if there are edges, but as a safeguard
        return {edge: 1.0 for edge in G.edges()}
        
    for u, v in G.edges():
        du, dv = degrees[u], degrees[v]
        numerator = abs(du - d_bar) + abs(dv - d_bar)
        w_uv = 1 - (numerator / denominator)
        weights[(u, v)] = w_uv
        
    return weights

def calculate_Rs(G):
    """
    Calculates the Structurally-Weighted Resistance Entropy (Rs) for a given graph.
    This is the main function that orchestrates the calculation.
    Rs(G) = sum_{u,v in E} [ p_uv * log(1/p_uv) * w_uv ]

    Args:
        G (nx.Graph): A connected NetworkX graph.

    Returns:
        float: The calculated Rs value for the graph.
               Returns 0 if the graph has no edges.
    """
    if not nx.is_connected(G):
        # The Rs metric as defined relies on effective resistance, which is for connected graphs.
        raise ValueError("Rs calculation requires a connected graph.")

    if G.number_of_edges() == 0:
        return 0.0

    # Step 1: Calculate effective resistances (Omega_uv)
    edge_resistances = calculate_all_pairs_effective_resistance(G)
    
    # Step 2: Calculate the edge probability distribution (p_uv)
    total_resistance = sum(edge_resistances.values())
    if total_resistance == 0:
        return 0.0 # Avoid division by zero
        
    edge_probabilities = {edge: res / total_resistance for edge, res in edge_resistances.items()}

    # Step 3: Calculate structural penalty weights (w_uv)
    structural_weights = calculate_structural_penalty_weights(G)

    # Step 4: Calculate Rs by summing the contributions of each edge
    Rs_value = 0.0
    for edge in G.edges():
        p_uv = edge_probabilities.get(edge)
        w_uv = structural_weights.get(edge)

        # The term is p_uv * log(1/p_uv). If p_uv is 0, this term is 0.
        if p_uv is not None and w_uv is not None and p_uv > 0:
            spectral_flow_term = p_uv * np.log(1 / p_uv)
            Rs_value += spectral_flow_term * w_uv
            
    return Rs_value
