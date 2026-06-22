import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pycirclize import Circos
from collections import Counter
import pandas as pd
import sys 

"""Code to load the results of NBS and visualize the significant edges in a circular plot.
    Morevoer it creates particpation counts for each network and network pair, and identifies top nodes by strength.
    Results loaded here come from Matlab code, so they are in .mat format."""

def save_summary_tables(summary, out_prefix):
    print(out_prefix)
    pd.DataFrame(
        [(k, v) for k, v in summary["network_edge_counts"].items()],
        columns=["network", "edge_participation"]
    ).to_csv(f"{out_prefix}_network_edge_counts.csv", index=False)

    pd.DataFrame(
        [(f"{a}-{b}", c) for (a, b), c in summary["network_pair_counts"].items()],
        columns=["network_pair", "count"]
    ).to_csv(f"{out_prefix}_network_pair_counts.csv", index=False)

    pd.DataFrame(
        summary["top_nodes_by_strength"],
        columns=["label", "network", "nodal_strength"]
    ).to_csv(f"{out_prefix}_top_nodes_by_strength.csv", index=False)


def summarize_nbs_component(binary_mat, tval_mat, labels, network_labels, top_n=10):
    """
    Summarize a significant NBS component.

    Parameters
    ----------
    binary_mat : np.ndarray
        Binary adjacency matrix of significant edges.
    tval_mat : np.ndarray
        Matrix of edgewise test statistics.
    labels : list of str
        Node labels.
    network_labels : list of str
        Network label for each node.
    top_n : int
        Number of top nodes/pairs to return.

    Returns
    -------
    summary : dict
        Dictionary with counts and rankings.
    """

    binary_mat = np.asarray(binary_mat)
    tval_mat = np.asarray(tval_mat)
    weighted_mat = binary_mat * tval_mat

    n = binary_mat.shape[0]

    # collect significant edges 
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if binary_mat[i, j] != 0:
                edges.append((i, j, weighted_mat[i, j]))

    n_edges = len(edges)

    #  within vs between
    within_count = 0
    between_count = 0
    pair_counter = Counter()
    network_edge_counter = Counter()
    network_node_counter = Counter()

    involved_nodes = set()

    for i, j, w in edges:
        ni = network_labels[i]
        nj = network_labels[j]

        involved_nodes.update([i, j])

        pair = tuple(sorted((ni, nj)))
        pair_counter[pair] += 1

        if ni == nj:
            network_edge_counter[ni] += 1   # one within-network edge → count once
            within_count += 1
        else:
            network_edge_counter[ni] += 1   # between-network edge → count once per side
            network_edge_counter[nj] += 1
            between_count += 1

    # unique nodes per network (not edge-endpoint occurrences)
    network_node_counter = Counter(network_labels[i] for i in involved_nodes)

    nodal_strength = weighted_mat.sum(axis=1)

    top_nodes = sorted(
        [(labels[i], network_labels[i], nodal_strength[i]) for i in range(n) if i in involved_nodes],
        key=lambda x: abs(x[2]),
        reverse=True
    )[:top_n]

    summary = {
        "n_edges": n_edges,
        "within_edges": within_count,
        "between_edges": between_count,
        "within_pct": 100 * within_count / n_edges if n_edges > 0 else 0,
        "between_pct": 100 * between_count / n_edges if n_edges > 0 else 0,
        "network_edge_counts": dict(network_edge_counter.most_common()),
        "network_node_counts": dict(network_node_counter.most_common()),
        "network_pair_counts": dict(pair_counter.most_common()),
        "top_nodes_by_strength": top_nodes,
    }

    return summary

def threshold_to_edges(matrix, abs_threshold=0, binarize=False):
    """
    Convert a matrix into a list of edges.

    Parameters
    ----------
    matrix : np.ndarray
        Square matrix.
    abs_threshold : float
        Keep only edges with abs(weight) > abs_threshold.
    binarize : bool
        If True, all kept edges get weight 1.

    Returns
    -------
    edges : list of tuple
        Each tuple is (i, j, weight).
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input matrix must be square.")

    n = matrix.shape[0]
    edges = []

    for i in range(n):
        for j in range(i + 1, n):
            w = matrix[i, j]
            if np.abs(w) > abs_threshold:
                if binarize:
                    w = 1.0
                edges.append((i, j, float(w)))

    return edges


def plot_connectivity_circos_weighted(
    adj_matrix,
    labels,
    network_labels,
    network_colors,
    title="Circular connectivity plot",
    edge_threshold=0.0,
    use_abs_threshold=True,
    only_upper_triangle=True,
    figsize=(10, 10),
    label_size=12,
    min_lw=0.3,
    max_lw=4.0,
    fig_name=None
):
    """
    Plot edges from a weighted 2D adjacency matrix on a circular plot.
    Node labels are color-coded by network.
    Edge thickness increases with edge weight.

    Parameters
    ----------
    adj_matrix : np.ndarray
        Square 2D adjacency matrix.
    labels : list of str
        Node labels, one per row/column of adj_matrix.
    network_labels : list of str
        Network membership of each node, same length as labels.
    network_colors : dict
        Dictionary mapping network name to color.
    title : str
        Figure title.
    edge_threshold : float
        Minimum threshold to keep an edge.
    use_abs_threshold : bool
        If True, threshold uses abs(weight).
        If False, threshold uses raw weight.
    only_upper_triangle : bool
        If True, draw only upper triangle to avoid duplicates.
    figsize : tuple
        Figure size.
    label_size : int
        Font size for node labels.
    min_lw : float
        Minimum line width for weakest retained edge.
    max_lw : float
        Maximum line width for strongest retained edge.
    """

    adj_matrix = np.asarray(adj_matrix)

    if adj_matrix.ndim != 2 or adj_matrix.shape[0] != adj_matrix.shape[1]:
        raise ValueError("adj_matrix must be a square 2D matrix.")

    n = adj_matrix.shape[0]

    if len(labels) != n:
        raise ValueError("labels must have the same length as matrix dimension.")

    if len(network_labels) != n:
        raise ValueError("network_labels must have the same length as matrix dimension.")

    # Build list of edges to draw
    edges = []
    for i in range(n):
        j_range = range(i + 1, n) if only_upper_triangle else range(n)

        for j in j_range:
            if i == j:
                continue

            w = adj_matrix[i, j]

            keep = abs(w) > edge_threshold if use_abs_threshold else w > edge_threshold
            if keep:
                edges.append((i, j, w))

    if len(edges) == 0:
        raise ValueError("No edges survived the threshold.")

    # For weighted linewidths
    abs_weights = np.array([abs(w) for _, _, w in edges], dtype=float)
    wmin = abs_weights.min()
    wmax = abs_weights.max()

    def scale_linewidth(w):
        aw = abs(w)
        if wmax == wmin:
            return (min_lw + max_lw) / 2
        return min_lw + (aw - wmin) / (wmax - wmin) * (max_lw - min_lw)

    # One sector per node
    sectors = {label: 1 for label in labels}
    circos = Circos(sectors, space=1)

    # Add outer track and colored labels
    for i, sector in enumerate(circos.sectors):
        net = network_labels[i]
        node_color = network_colors.get(net, "lightgray")

        track = sector.add_track((92, 100))
        track.axis(fc=node_color, ec="white", lw=0.5)

        sector.text(
            sector.name,
            r=106,
            size=label_size,
            color=node_color,
            orientation="vertical",
            fontweight="bold",
        )

    # Draw weighted edges
    for i, j, w in edges:
        lw = scale_linewidth(w)

        # Optional: positive and negative edges in different colors
        edge_color = "black" #"red" if w > 0 else "blue"

        circos.link(
            (labels[i], 0.85, 0.85),
            (labels[j], 0.85, 0.85),
            color=edge_color,
            lw=lw,
            alpha=0.7,
            
        )

    fig = circos.plotfig(figsize=figsize)
    fig.suptitle(title, y=0.98)
    fig.savefig(fig_name, dpi=300)
    plt.show()

 
if __name__ == "__main__":
    do_between = True  
    do_within = False
    if do_between == True and do_within == False: 
        file_path = "/home/riccardo/codici_progetti/Salerno/structural_similarity/matlab_stuff/NBS_results/EULER_TH_120/protocol_yes/TH_3.1_PERM10000_ALPHA0.05/between_groups"
        file_list = ["MSA_gt_HC_NBS_results.mat", "HC_gt_MSA_NBS_results.mat"]
    elif do_between == False and do_within == True:  
        file_path = "/home/riccardo/codici_progetti/Salerno/structural_similarity/matlab_stuff/NBS_results/EULER_TH_120/protocol_yes/TH_3.1_PERM10000_ALPHA0.05/within_group"
        file_list = ["MSA_Disease_Duration_NBS_results.mat", "MSA_CerebGM_NBS_results.mat", "MSA_Pons_NBS_results.mat"]
    for f in file_list:
        mat_file = f"{file_path}/{f}"
        f_name = f.split(".")[0]
        mat = loadmat(mat_file, struct_as_record=False, squeeze_me=True)
        res = mat["res"]
        
        # binary matrix in sparse form 
        binary_mat = res.NBS.con_mat.toarray()  # Convert sparse matrix to dense array 
        # t values 
        tval_mat = np.array(res.NBS.test_stat)
        # labels 
        labels_old = res.NBS.node_label.tolist()  # Convert MATLAB cell array to Python list  
        labels_new = [s.replace("rh_", "").replace("lh_", "").replace("7Networks","").replace("_LH","LH").replace("_RH","RH").replace("Default","DMN") \
                    .replace("Vis", "VIS").replace("SomMot", "SOMMOT").replace("DorsAttn", "DAN").replace("SalVentAttn", "SAL").replace("Limbic", "LIMB").replace("Cont", "CONT") for s in labels_old]
        
        count_vis = sum("VIS" in s for s in labels_new)
        count_sm = sum("SOMMOT" in s for s in labels_new)
        count_dan = sum("DAN" in s for s in labels_new)
        count_sal = sum("SAL" in s for s in labels_new)
        count_limb = sum("LIMB" in s for s in labels_new)
        count_fp = sum("CONT" in s for s in labels_new)
        count_dmn = sum("DMN" in s for s in labels_new)

        print(count_vis, count_sm, count_dan, count_sal, count_limb, count_fp, count_dmn)

        network_labels = (
                        ["VIS"] * 17 +
                        ["SOMMOT"] * 14 +
                        ["DAN"] * 15 +
                        ["SAL"] * 12 +
                        ["LIMB"] * 5 +
                        ["CONT"] * 13 +
                        ["DMN"] * 24
                        )
        network_colors = {
                "VIS": "#1f77b4",
                "SOMMOT":  "#d62728",
                "DAN": "#2ca02c",
                "SAL": "#ff7f0e", 
                "LIMB": "#9467bd",
                "CONT": "#17becf", 
                "DMN":  "#8c564b",  # keep only if you actually use SAL separately
            }
        
        #print(tval_mat*binary_mat)

        #plot_connectivity_circos_weighted(
        #    adj_matrix=tval_mat*binary_mat,
        #    labels=labels_new,
        #    network_labels=network_labels,
        #    network_colors=network_colors,
        #    title=None, #"Weighted connectivity",
        #    edge_threshold=0,
        #    label_size=12,
        #    min_lw=0.5,
        #    max_lw=5,
        #    fig_name = f"{file_path}/{f_name}_weighted_circos.png") 
        
        #plot_connectivity_circos_weighted(
        #    adj_matrix=binary_mat,
        #    labels=labels_new,
        #    network_labels=network_labels,
        #    network_colors=network_colors,
        #    title=None, #"Binary connectivity",
        #    edge_threshold=0,
        #    label_size=12,
        #    min_lw=0.5,
        #    max_lw=5,
        #    fig_name = f"{file_path}/{f_name}_binary_circos.png")

        summary = summarize_nbs_component(
            binary_mat=binary_mat,
            tval_mat=tval_mat,
            labels=labels_new,
            network_labels=network_labels,
            top_n=20
        )

        #print(f"\n===== {f_name} =====")
        #print(f"Number of significant edges: {summary['n_edges']}")
        print(
            f"Within-network edges: {summary['within_edges']} "
            f"({summary['within_pct']:.1f}%)"
        )
        print(
            f"Between-network edges: {summary['between_edges']} "
            f"({summary['between_pct']:.1f}%)"
        )

        #print("\nMain networks involved (edge participation):")
        for net, count in summary["network_edge_counts"].items():
            print(f"  {net}: {count}")

        #print("\nMost represented network pairs:")
        for pair, count in summary["network_pair_counts"].items():
            print(f"  {pair[0]}-{pair[1]}: {count}")

        #print("\nTop nodes by nodal strength:")
        for lab, net, val in summary["top_nodes_by_strength"]:
            print(f"  {lab} ({net}): {val:.3f}")

        save_summary_tables(summary, f"{file_path}/{f_name}")
        #print(f"{file_path}/{f_name}")
        #print(f"finished {f_name}")


