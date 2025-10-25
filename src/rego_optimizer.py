# src/rego_optimizer.py

import networkx as nx
import random
import copy
from .calculate_Rs import calculate_Rs

def rego_optimizer(G_initial, n_iter=1000, verbose=True):
    """
    Optimizes a network's topology to maximize its Rs value using the
    REGO (Resistance Entropy Gradient Optimization) algorithm.

    The algorithm performs degree-preserving edge swaps. At each step, it
    attempts a random swap and accepts it only if it increases the Rs value.

    Args:
        G_initial (nx.Graph): The initial connected NetworkX graph.
        n_iter (int): The number of iterations (potential swaps) to perform.
        verbose (bool): If True, prints progress updates.

    Returns:
        tuple: A tuple containing:
            - G_optimized (nx.Graph): The final, optimized graph.
            - rs_history (list): A list of Rs values at each accepted swap.
    """
    if not nx.is_connected(G_initial):
        raise ValueError("Initial graph must be connected for REGO optimization.")

    # Create a deep copy to avoid modifying the original graph
    G = copy.deepcopy(G_initial)
    
    # Initial Rs calculation
    try:
        current_rs = calculate_Rs(G)
    except ValueError as e:
        print(f"Error calculating initial Rs: {e}")
        return G_initial, []

    rs_history = [current_rs]
    
    if verbose:
        print(f"Starting REGO optimization. Initial Rs: {current_rs:.4f}")

    # Get a list of edges to sample from
    edges = list(G.edges())
    n_edges = len(edges)
    if n_edges < 2:
        print("Optimization requires at least 2 edges. Returning initial graph.")
        return G, rs_history

    accepted_swaps = 0
    for i in range(n_iter):
        # 1. Select two random edges (u, v) and (x, y) for a potential swap.
        # Ensure the edges do not share any nodes.
        if n_edges < 2: break
        
        # Try to find a valid pair of edges for swapping
        attempts = 0
        while attempts < 100: # Add a limit to prevent infinite loops
            edge1_idx, edge2_idx = random.sample(range(n_edges), 2)
            u, v = edges[edge1_idx]
            x, y = edges[edge2_idx]
            
            # Check if they share any nodes
            if len({u, v, x, y}) == 4:
                break # Found valid edges to swap
            attempts += 1
        else:
            # Could not find swappable edges after 100 attempts
            if verbose: print("Could not find swappable edges. Continuing...")
            continue
            
        # 2. Propose a swap: (u, v), (x, y) -> (u, y), (v, x)
        # Check if the new edges would create multi-edges.
        if G.has_edge(u, y) or G.has_edge(v, x):
            continue

        # 3. Temporarily perform the swap
        G.remove_edge(u, v)
        G.remove_edge(x, y)
        G.add_edge(u, y)
        G.add_edge(v, x)

        # 4. Check connectivity. If disconnected, revert.
        if not nx.is_connected(G):
            G.remove_edge(u, y)
            G.remove_edge(v, x)
            G.add_edge(u, v)
            G.add_edge(x, y)
            continue

        # 5. Calculate new Rs and decide
        try:
            new_rs = calculate_Rs(G)
        except ValueError:
            G.remove_edge(u, y)
            G.remove_edge(v, x)
            G.add_edge(u, v)
            G.add_edge(x, y)
            continue
            
        if new_rs > current_rs:
            # Accept the swap
            current_rs = new_rs
            rs_history.append(current_rs)
            edges = list(G.edges()) # Update edge list
            n_edges = len(edges)
            accepted_swaps += 1
            if verbose and (accepted_swaps % 10 == 0):
                print(f"Iter {i+1}/{n_iter} | Accepted: {accepted_swaps} | Rs: {current_rs:.4f}")
        else:
            # Reject and revert
            G.remove_edge(u, y)
            G.remove_edge(v, x)
            G.add_edge(u, v)
            G.add_edge(x, y)

    if verbose:
        print(f"\nOptimization finished.")
        print(f"Total iterations: {n_iter}, Total accepted swaps: {accepted_swaps}.")
        print(f"Final Rs: {current_rs:.4f}")

    return G, rs_history
