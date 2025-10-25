# src/clinical_analysis.py

import networkx as nx
import numpy as np
import itertools
from .calculate_Rs import calculate_Rs

def reconstruct_vascular_network(adj_matrix):
    """
    Reconstructs a vascular network from an adjacency matrix.

    In a real-world scenario, this adjacency matrix would be derived from
    DCE-MRI data, where nodes are vessel bifurcations and edges represent
    vessel segments.

    Args:
        adj_matrix (np.ndarray): A square NumPy array representing the
                                 adjacency matrix of the vascular network.

    Returns:
        nx.Graph: A NetworkX graph object representing the network.
    """
    G = nx.from_numpy_array(adj_matrix)
    
    # Ensure the graph is connected. In a real scenario, we would only
    # analyze the largest connected component of the tumor vasculature.
    if not nx.is_connected(G):
        # Get the largest connected component subgraph
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
        print("Warning: Input matrix resulted in a disconnected graph. "
              "Using the largest connected component.")
              
    return G

def find_optimal_targets_greedy(G, num_targets=1, budget=None):
    """
    Finds the optimal set of edges (vessel segments) to remove in order to
    maximize the increase in Rs of the remaining network.

    This implements the greedy search algorithm for Rs-guided therapy. It
    iteratively finds the edge whose removal causes the largest increase (or smallest decrease)
    in Rs, removes it, and repeats.

    Args:
        G (nx.Graph): The initial vascular network graph.
        num_targets (int): The number of edges to target for removal.
        budget (float, optional): A constraint on the total "cost" of removed edges.
                                  For this example, we assume cost is uniform (1 per edge).
                                  If specified, the algorithm stops when the budget is
                                  exhausted or `num_targets` is reached.

    Returns:
        tuple: A tuple containing:
            - best_targets (list): A list of edge tuples to be targeted.
            - delta_rs_history (list): A list of the change in Rs (Delta_Rs) at each step.
            - final_G (nx.Graph): The graph after removing the targeted edges.
    """
    if not nx.is_connected(G):
        raise ValueError("Input graph must be connected.")

    # Create a copy to work with
    current_G = G.copy()
    
    # Calculate initial Rs
    try:
        initial_rs = calculate_Rs(current_G)
    except ValueError as e:
        print(f"Could not calculate initial Rs: {e}")
        return [], [], G

    best_targets = []
    delta_rs_history = []
    
    # Determine the number of iterations
    iterations = num_targets
    if budget is not None:
        iterations = min(num_targets, int(budget)) # Assuming cost=1 per edge

    print(f"Starting greedy search for {iterations} targets. Initial Rs: {initial_rs:.4f}")

    for i in range(iterations):
        candidate_edges = list(current_G.edges())
        best_edge_to_remove = None
        max_delta_rs = -np.inf
        next_rs = initial_rs

        if not candidate_edges:
            print("No more edges to remove.")
            break

        # Iterate through all possible edges to find the best one to remove
        for edge in candidate_edges:
            temp_G = current_G.copy()
            temp_G.remove_edge(*edge)

            # Important: The Rs metric is defined for a CONNECTED graph.
            # The goal of the therapy is to fragment the network.
            # We calculate Rs on the largest remaining connected component.
            if not nx.is_connected(temp_G):
                largest_cc_nodes = max(nx.connected_components(temp_G), key=len)
                largest_cc_graph = temp_G.subgraph(largest_cc_nodes)
            else:
                largest_cc_graph = temp_G

            # If the removal shatters the graph into tiny pieces, Rs might be 0
            if largest_cc_graph.number_of_edges() == 0:
                rs_after_removal = 0.0
            else:
                try:
                    rs_after_removal = calculate_Rs(largest_cc_graph)
                except ValueError:
                    # This can happen if a component becomes invalid for Rs calculation
                    continue
            
            delta_rs = rs_after_removal - initial_rs
            
            if delta_rs > max_delta_rs:
                max_delta_rs = delta_rs
                best_edge_to_remove = edge
                next_rs = rs_after_removal

        if best_edge_to_remove is None:
            print("No beneficial edge removal found. Stopping.")
            break

        # Permanently remove the best edge found in this iteration
        current_G.remove_edge(*best_edge_to_remove)
        
        # Update state
        initial_rs = next_rs # The new baseline Rs is the one from the new graph state
        best_targets.append(best_edge_to_remove)
        delta_rs_history.append(max_delta_rs)
        
        print(f"Step {i+1}/{iterations}: Targeted edge {best_edge_to_remove}. "
              f"New Rs: {initial_rs:.4f} (Delta Rs: {max_delta_rs:+.4f})")

    return best_targets, delta_rs_history, current_G
