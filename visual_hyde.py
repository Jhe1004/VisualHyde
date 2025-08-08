# coding: utf-8
import os
import sys
import argparse
import csv
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.colorbar import ColorbarBase
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces, random_color, TextFace
from PIL import Image, ImageEnhance
import warnings
import colorsys

# --- Dependencies and Description ---
matplotlib.use('Agg')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

Description = (
    '''
    -----------------------------------------------------------------------------
    |                               visual_hyde.py                              |
    -----------------------------------------------------------------------------

    Created by Jian He (j.he930724@gamail.com)
    Modified for enhanced logging and advanced algorithms by the Assistant.
    The node calculation now uses a 'difference consensus' algorithm.
    The script now filters results based on P-value instead of Z-score.
    Added an optional false positive filter based on null-hypothesis simulations,
    with distinct logic for leaf and node modes.
    Added robust error checking for all input file name consistencies.

    Dependencies:
    python3
    ete3 (conda install -c etetoolkit ete3 ete_toolchain)
    matplotlib (conda install matplotlib)
    numpy (conda install numpy)
    pandas (conda install pandas)
    Pillow (conda install Pillow or pip install Pillow)
    colorsys (standard library)
    csv (standard library)
    '''
)

# --- Utility Functions ---
def get_safe_leaf_name(leaf_name):
    """Converts a leaf name into a filesystem-safe string."""
    return "".join(c if c.isalnum() else "_" for c in str(leaf_name))

def apply_relative_filter(real_heatmap, null_heatmap, threshold):
    """
    Applies the relative strength filter to a heatmap DataFrame.
    Returns the filtered DataFrame.
    """
    if null_heatmap is None:
        return real_heatmap
    
    filtered_heatmap = real_heatmap.copy()
    filtered_count = 0
    for p1 in filtered_heatmap.index:
        for p2 in filtered_heatmap.columns:
            gamma_real = filtered_heatmap.loc[p1, p2]
            if pd.notna(gamma_real):
                try:
                    gamma_null = null_heatmap.loc[p1, p2]
                    if pd.notna(gamma_null):
                        s_gamma_real = min(gamma_real, 1 - gamma_real)
                        s_gamma_null = min(gamma_null, 1 - gamma_null)

                        if s_gamma_real > 1e-9:
                            ratio = s_gamma_null / s_gamma_real
                            if ratio >= threshold:
                                filtered_heatmap.loc[p1, p2] = np.nan
                                filtered_count += 1
                        elif s_gamma_null > 1e-9:
                            filtered_heatmap.loc[p1, p2] = np.nan
                            filtered_count += 1
                except KeyError:
                    continue
    if filtered_count > 0:
        print(f"    ...filtered out {filtered_count} potential false positive signals.")
    return filtered_heatmap


# --- Data Parsing and Preparation Functions ---
def make_predefined_clade_file(tree_file, max_leaves_per_clade=5):
    """
    Automatically generates a simple predefined clade file based on the input tree.
    """
    try:
        t = Tree(tree_file)
        if len(t.children) != 2:
            print("Warning: Tree does not appear to be bifurcating or root position is unclear, automatic clade definition might be inaccurate.")
            if len(t.children) > 0:
                ingroup_node = max(t.children, key=len)
                print(f"Warning: Attempting to use subtree with {len(ingroup_node)} leaf nodes as ingroup for automatic clades.")
            else:
                ingroup_node = t
        elif len(t.children[0]) > len(t.children[1]):
            ingroup_node = t.children[0]
        else:
            ingroup_node = t.children[1]

        output_filename = "Predefined_clade.txt"
        with open(output_filename, "w") as write_file:
            if len(ingroup_node) <= max_leaves_per_clade and not ingroup_node.is_leaf():
                leaf_names = ingroup_node.get_leaf_names()
                write_file.write(",".join(leaf_names) + ",\n")
            else:
                def return_species_name_clade(node):
                    child_node = node.children
                    for each_node in child_node:
                        if not each_node.is_leaf() and len(each_node) == 0: continue
                        node_leaves = each_node.get_leaf_names()
                        if len(node_leaves) <= max_leaves_per_clade:
                            for each_species_name in node_leaves:
                                write_file.write(each_species_name + ",")
                            write_file.write("\n")
                        else:
                            return_species_name_clade(each_node)
                return_species_name_clade(ingroup_node)

        print(f"Automatically generated predefined clade file: {output_filename}")
        return output_filename
    except Exception as e:
        print(f"Error: Could not automatically generate predefined clade file: {e}")
        return None

def parse_tree(tree_file, predefined_clade_file):
    """
    Parses the species tree file and the predefined clade file.
    """
    try:
        t = Tree(tree_file)
        leaf_labels = []
        if len(t.children) == 2:
            child1_len = len(t.children[0])
            child2_len = len(t.children[1])
            if child1_len > child2_len: ingroup_node = t.children[0]
            elif child2_len > child1_len: ingroup_node = t.children[1]
            else:
                print("Warning: Tree root children have equal size. Using all leaf nodes as labels.")
                ingroup_node = t
            leaf_labels = ingroup_node.get_leaf_names()
        else:
            print("Warning: Tree root does not appear standard, using all leaf nodes as labels.")
            leaf_labels = t.get_leaf_names()

        if not leaf_labels:
            print("Error: Could not extract leaf node labels from the tree.")
            return None

        subtrees = []
        try:
            with open(predefined_clade_file, "r") as read_file:
                for each_line in read_file:
                    line_list_orig = [name.replace("\n","").replace(" ","") for name in each_line.split(",")]
                    line_list = [name for name in line_list_orig if name and name in leaf_labels]
                    if line_list:
                        try:
                            first_leaf_index = leaf_labels.index(line_list[0])
                            subtrees.append(line_list + [first_leaf_index])
                        except ValueError:
                            print(f"Warning: Name '{line_list[0]}' in predefined clade not found in leaf labels, skipping.")
        except FileNotFoundError:
            print(f"Error: Predefined clade file not found: {predefined_clade_file}")
            return None

        subtrees.sort(key=lambda elem: elem[-1])
        highlight_subtrees = [elem[:-1] for elem in subtrees]
        name_len_list = [len(str(name)) for name in leaf_labels]
        name_len = (min(name_len_list) if name_len_list else 0, max(name_len_list) if name_len_list else 0)
        return t, leaf_labels, highlight_subtrees, name_len
    except FileNotFoundError:
        print(f"Error: Tree file not found - {tree_file}")
        return None
    except Exception as e:
        print(f"Error: Error parsing tree or clade file: {e}")
        return None

def make_hotmap_table_gamma(leaf_labels, hypothesis_hybrid_species, hyde_data_df, pvalue_threshold, 
                            null_heatmap_data=None, filter_ratio_threshold=0.05):
    """
    Generates the heatmap data table (DataFrame) for the specified hybrid.
    Optionally filters the results based on a null-hypothesis dataset.
    """
    hypothesis_hybrid_species_str = str(hypothesis_hybrid_species)
    leaf_labels_str = [str(lbl) for lbl in leaf_labels]
    hyde_data_df_copy = hyde_data_df.copy()
    for col in ['Hybrid', 'P1', 'P2']:
        if col in hyde_data_df_copy.columns:
            hyde_data_df_copy[col] = hyde_data_df_copy[col].astype(str)

    sub_table = hyde_data_df_copy[
        (hyde_data_df_copy["Hybrid"] == hypothesis_hybrid_species_str) &
        (hyde_data_df_copy["Pvalue"].astype(float) < pvalue_threshold)
    ]

    def get_gamma(hyde_out_table, each_index, each_column):
        idx_str = str(each_index)
        col_str = str(each_column)
        sub_table_p1p2 = hyde_out_table[(hyde_out_table["P1"] == idx_str) & (hyde_out_table["P2"] == col_str)]
        if not sub_table_p1p2.empty:
            return pd.to_numeric(sub_table_p1p2["Gamma"].iloc[0], errors='coerce')
        else:
            sub_table_p2p1 = hyde_out_table[(hyde_out_table["P1"] == col_str) & (hyde_out_table["P2"] == idx_str)]
            if not sub_table_p2p1.empty:
                gamma_p2p1 = pd.to_numeric(sub_table_p2p1["Gamma"].iloc[0], errors='coerce')
                if pd.notna(gamma_p2p1): return 1.0 - gamma_p2p1
                else: return np.nan
            else: return np.nan

    df = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
    for each_index in leaf_labels_str:
        for each_column in leaf_labels_str:
            if each_index == each_column: continue
            gamma = get_gamma(sub_table, each_index, each_column)
            if pd.notna(gamma):
                if 0 <= gamma <= 1: df.at[each_index, each_column] = gamma
                else: pass

    if null_heatmap_data:
        safe_leaf_name_key = get_safe_leaf_name(hypothesis_hybrid_species_str)
        null_heatmap = null_heatmap_data.get(safe_leaf_name_key)
        if null_heatmap is not None:
            print(f"  Applying null hypothesis filter for leaf: {hypothesis_hybrid_species_str}...")
            df = apply_relative_filter(df, null_heatmap, filter_ratio_threshold)

    return df

def calculate_all_node_heatmaps(tree, leaf_labels, hyde_data_df, pvalue_threshold, 
                                gamma_diff_threshold=0.2, 
                                diff_consensus_threshold=0.05,
                                diff_consensus_ratio=0.5):
    """
    Calculates node heatmaps using the two-pass 'difference consensus' algorithm.
    NOTE: This function intentionally does NOT apply the null-hypothesis filter during its
    internal calculations to ensure the consensus algorithm works on unfiltered data.
    Filtering should be applied as a post-processing step in the main loop.
    """
    leaf_labels_str = [str(lbl) for lbl in leaf_labels]
    all_leaf_labels_set_str = set(leaf_labels_str)
    
    # --- STAGE 1: Postorder traversal to generate DRAFT and DIFFERENCE heatmaps ---
    print(f"Stage 1 (Postorder): Generating draft values and difference maps with loose threshold ({gamma_diff_threshold})...")
    draft_heatmaps = {}
    diff_heatmaps = {}
    
    all_nodes_postorder = list(tree.traverse("postorder"))
    total_nodes = len(all_nodes_postorder)
    processed_count = 0

    for node in all_nodes_postorder:
        processed_count += 1
        if processed_count % 20 == 0 or processed_count == total_nodes:
            percentage = (processed_count / total_nodes) * 100
            print(f"  ...processed {processed_count}/{total_nodes} nodes ({percentage:.1f}%)")

        if node.is_leaf():
            node_name_str = str(node.name)
            if node_name_str in all_leaf_labels_set_str:
                # *** MODIFIED: Node mode calculation IGNORES leaf-level filters ***
                draft_heatmaps[node] = make_hotmap_table_gamma(leaf_labels_str, node_name_str, hyde_data_df, pvalue_threshold)
            else:
                draft_heatmaps[node] = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
        else: # Internal node
            children = node.children
            if len(children) == 2:
                child1, child2 = children
                child1_heatmap = draft_heatmaps.get(child1)
                child2_heatmap = draft_heatmaps.get(child2)

                current_draft_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                current_diff_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)

                if child1_heatmap is not None and child2_heatmap is not None:
                    node_descendants_str = set(str(n) for n in node.get_leaf_names())
                    for p1_name in leaf_labels_str:
                        if p1_name in node_descendants_str: continue
                        for p2_name in leaf_labels_str:
                            if p2_name in node_descendants_str or p1_name == p2_name: continue
                            try:
                                gamma1 = child1_heatmap.loc[p1_name, p2_name]
                                gamma2 = child2_heatmap.loc[p1_name, p2_name]
                                if pd.notna(gamma1) and pd.notna(gamma2):
                                    diff = abs(gamma1 - gamma2)
                                    current_diff_heatmap.loc[p1_name, p2_name] = diff
                                    if diff < gamma_diff_threshold:
                                        avg_gamma = (gamma1 + gamma2) / 2.0
                                        current_draft_heatmap.loc[p1_name, p2_name] = avg_gamma
                            except (KeyError, TypeError):
                                continue
                
                draft_heatmaps[node] = current_draft_heatmap
                diff_heatmaps[node] = current_diff_heatmap
            else:
                draft_heatmaps[node] = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                diff_heatmaps[node] = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
    
    print("Stage 1 completed. Draft and difference maps generated for all nodes.")

    # --- STAGE 2: Preorder traversal for FINAL heatmaps using 'difference consensus' ---
    final_heatmaps = {}
    final_logs = {}

    print(f"\nStage 2 (Preorder, 'difference consensus'): Refining with diff threshold={diff_consensus_threshold}, ratio={diff_consensus_ratio}...")

    for node in tree.traverse("preorder"):
        if node.is_leaf():
            final_heatmaps[node] = draft_heatmaps.get(node)
            continue

        final_heatmap_for_node = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
        current_node_log_data = [] 
        
        subtree_nodes = [n for n in node.traverse() if not n.is_leaf()]
        node_descendants_str = set(str(n) for n in node.get_leaf_names())

        for p1_name in leaf_labels_str:
            if p1_name in node_descendants_str: continue
            for p2_name in leaf_labels_str:
                if p2_name in node_descendants_str or p1_name == p2_name: continue
                
                diff_values_in_subtree = [
                    diff_heatmaps.get(sub_node).loc[p1_name, p2_name]
                    for sub_node in subtree_nodes
                    if sub_node in diff_heatmaps and pd.notna(diff_heatmaps.get(sub_node).loc[p1_name, p2_name])
                ]
                
                total_internal_nodes = len(diff_values_in_subtree)
                if total_internal_nodes == 0: continue

                low_diff_count = sum(1 for diff in diff_values_in_subtree if diff < diff_consensus_threshold)
                ratio_found = low_diff_count / total_internal_nodes if total_internal_nodes > 0 else 0
                
                result = 'FAIL'
                final_gamma_val = 'N/A'

                if ratio_found > diff_consensus_ratio:
                    result = 'PASS'
                    final_gamma = draft_heatmaps.get(node).loc[p1_name, p2_name]
                    if pd.notna(final_gamma):
                        final_heatmap_for_node.loc[p1_name, p2_name] = final_gamma
                        final_gamma_val = f"{final_gamma:.4f}"

                diff_values_str = ";".join(f"{d:.4f}" for d in diff_values_in_subtree)

                log_row = [
                    p1_name, p2_name,
                    f"{diff_consensus_ratio:.2f}", f"{ratio_found:.2f}",
                    low_diff_count, total_internal_nodes,
                    f"< {diff_consensus_threshold}",
                    final_gamma_val, result, diff_values_str
                ]
                current_node_log_data.append(log_row)

        final_heatmaps[node] = final_heatmap_for_node
        final_logs[node] = current_node_log_data

    print("Stage 2 completed. All node heatmaps have been finalized.")
    return final_heatmaps, final_logs

def check_name_consistency(tree, hyde_df):
    """
    Strictly checks if all names in the HyDe file are present in the tree.
    Exits with an error if any mismatch is found.
    """
    leaves_name_in_species_tree = set()
    try:
        if len(tree.children) == 2 and len(tree.children[0]) != len(tree.children[1]):
            ingroup_node = max(tree.children, key=len)
        else: ingroup_node = tree
        leaves_name_in_species_tree = set(str(leaf.name) for leaf in ingroup_node.get_leaves() if leaf.name is not None)
    except Exception as e:
        print(f"FATAL ERROR: Could not extract leaf names from tree: {e}"); sys.exit(1)
    
    if not leaves_name_in_species_tree:
        print("FATAL ERROR: Could not extract any valid leaf node names from the provided tree file."); sys.exit(1)
    
    hyde_names = set()
    try:
        for col in ['P1', 'Hybrid', 'P2']:
            if col in hyde_df.columns:
                valid_names = {str(name) for name in hyde_df[col].dropna().unique() if name}
                hyde_names.update(valid_names)
    except Exception as e: 
        print(f"FATAL ERROR: Could not extract names from HyDe data file: {e}"); sys.exit(1)

    if not hyde_names: 
        print("FATAL ERROR: No valid names (P1/Hybrid/P2) could be extracted from the HyDe file."); sys.exit(1)

    missing_in_tree = hyde_names - leaves_name_in_species_tree
    if not missing_in_tree:
        print("Name consistency check passed: All HyDe participants found in tree leaves.")
        missing_in_hyde = leaves_name_in_species_tree - hyde_names
        if missing_in_hyde: 
            print(f"  (Note: These tree leaves are not participants in any HyDe triplet: {', '.join(missing_in_hyde)})")
        return True
    else:
        print("\nFATAL ERROR: Name mismatch detected. The script cannot continue.")
        print(f"The following names were found in your HyDe input file but are MISSING from your species tree:")
        for name in sorted(list(missing_in_tree)):
            print(f"  - {name}")
        print("Please correct the names in your input files to ensure they match perfectly.")
        sys.exit(1)

# --- NEW CONSOLIDATED PLOTTING FUNCTION ---
def generate_hyde_visualization(fig_prefix, hyde_output_array, tree_file, node_num, highlight_target, clade_definitions, name_len, colorbar_brightness_factor):
    """
    Consolidates the logic of drawing the heatmap, tree, colorbar, and combining them into a single figure.
    This function creates and deletes temporary files as part of its execution.
    Returns True on success, False on failure.
    """
    temp_hotmap_file = None
    temp_colorbar_file = None
    temp_tree_file = None

    # --- Part 1: Draw Hotmap (from draw_hotmap) ---
    try:
        hyde_output_array.to_csv(fig_prefix + ".csv", index=True, sep=',')
    except Exception as e:
        print(f"Warning: Could not save heatmap data to CSV '{fig_prefix}.csv': {e}")
    
    num_leaves = hyde_output_array.shape[0]
    if num_leaves == 0:
        print("Warning: Heatmap data is empty, cannot draw heatmap.")
        return False

    inches_per_leaf = 0.4
    base_inches = 2.0
    max_inches = 30.0
    fig_inches = min(max_inches, base_inches + num_leaves * inches_per_leaf)
    
    fig_hotmap = plt.figure(figsize=(fig_inches, fig_inches))
    border_width = 0.00001
    ax_size = [0 + border_width, 0 + border_width, 1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig_hotmap.add_axes(ax_size)
    
    cmap_list = []
    color_val = 1.0; lucency = 0.0; steps = 5000
    for _ in range(steps):
        cmap_list.append([0, 0, color_val, lucency])
        color_val -= (1.0 / steps); lucency += (1.0 / steps)
    color_val = 0.0; lucency = 1.0
    for _ in range(steps):
        cmap_list.append([color_val, 0, 0, lucency])
        color_val += (1.0 / steps); lucency -= (1.0 / steps)
    
    cmap_array = np.clip(np.array(cmap_list), 0, 1)
    original_cmap = ListedColormap(cmap_array)
    original_cmap.set_bad(color='white', alpha=0)
    norm = Normalize(vmin=0, vmax=1)
    
    ax.imshow(hyde_output_array.astype(float), norm=norm, cmap=original_cmap, interpolation='nearest', aspect='equal')
    grid_linewidth = max(0.5, min(2, fig_inches / 10.0))
    ax.set_xticks(np.arange(hyde_output_array.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(hyde_output_array.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=grid_linewidth)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xticklabels([]); ax.set_yticklabels([])
    
    temp_hotmap_file = fig_prefix + "_hotmap_temp.png"
    try:
        plt.savefig(temp_hotmap_file, dpi=200)
    except Exception as e:
        print(f"Error: Error saving heatmap: {e}")
        plt.cla(); plt.clf(); plt.close(fig_hotmap); plt.close("all")
        return False
    finally:
        plt.cla(); plt.clf(); plt.close(fig_hotmap); plt.close("all")

    # --- Part 2: Draw Colorbar ---
    darker_colors_rgba = []
    cmap_for_colorbar = None
    if original_cmap:
        heatmap_colors_rgba = original_cmap.colors
        for r, g, b, a in heatmap_colors_rgba:
            try:
                h, s, v = colorsys.rgb_to_hsv(r, g, b)
                v_new = max(0.0, min(1.0, v * colorbar_brightness_factor))
                new_r, new_g, new_b = colorsys.hsv_to_rgb(h, s, v_new)
                darker_colors_rgba.append([new_r, new_g, new_b, a])
            except Exception as e:
                darker_colors_rgba.append([r, g, b, a])
        cmap_for_colorbar = ListedColormap(darker_colors_rgba)
        cmap_for_colorbar.set_bad(color='white', alpha=0)

    if cmap_for_colorbar and norm:
        fig_cbar = plt.figure(figsize=(2, 10))
        ax_cbar = fig_cbar.add_axes([0.2, 0.05, 0.2, 0.9])
        temp_colorbar_file = f"{fig_prefix}_colorbar_temp.png"
        try:
            cb = ColorbarBase(ax_cbar, cmap=cmap_for_colorbar, norm=norm, orientation='vertical')
            cb.ax.tick_params(labelsize=40)
            plt.savefig(temp_colorbar_file, dpi=200, bbox_inches='tight', pad_inches=0.1)
        except Exception as e:
            print(f"Error: Error saving colorbar image '{temp_colorbar_file}': {e}")
            plt.cla(); plt.clf(); plt.close(fig_cbar); plt.close("all")
            if os.path.exists(temp_hotmap_file): os.remove(temp_hotmap_file)
            return False
        finally:
            plt.cla(); plt.clf(); plt.close(fig_cbar); plt.close("all")
    else:
        if os.path.exists(temp_hotmap_file): os.remove(temp_hotmap_file)
        return False

    # --- Part 3: Draw Tree ---
    try:
        t = Tree(tree_file)
    except Exception as e:
        print(f"Error: Could not reload tree file {tree_file} in draw_tree: {e}")
        if os.path.exists(temp_hotmap_file): os.remove(temp_hotmap_file)
        if os.path.exists(temp_colorbar_file): os.remove(temp_colorbar_file)
        return False
    
    def node_layout_background(node):
        ns = NodeStyle(); ns["hz_line_width"]=4.5; ns["vt_line_width"]=4.5; ns["size"]=0; node.set_style(ns)
    def node_layout(node, color):
        if node.is_leaf():
            node_name_str = str(node.name) if node.name is not None else ""
            padding_dots = "·" * ((name_len[1] - len(node_name_str)) * 2 + 10)
            leaf_name_text = node_name_str + padding_dots
            descFace = faces.TextFace(leaf_name_text, fsize=30, fgcolor=color if color else 'black')
            descFace.margin_top = 3; descFace.margin_bottom = 3; descFace.border.margin = 1
            node.add_face(descFace, column=1, position='aligned')
        ns = NodeStyle(); ns["hz_line_width"]=4.5; ns["vt_line_width"]=4.5
        ns["vt_line_color"]=color; ns["hz_line_color"]=color; ns["size"]=0; node.set_style(ns)
    
    all_clade_leaves = set(str(name) for clade in clade_definitions for name in clade if name)
    for each_node in t.traverse():
        node_name_str = str(each_node.name) if each_node.name is not None else ""
        if each_node.is_leaf():
            highlight_target_str = str(highlight_target) if highlight_target is not None else None
            if isinstance(highlight_target, str) and node_name_str == highlight_target_str:
                node_face = TextFace("o", fsize=40); node_face.background.color = "red"
                each_node.add_face(node_face, column=0, position="branch-right")
            if node_name_str not in all_clade_leaves: pass
        else: node_layout_background(each_node)
    
    h = 0
    for each_subtree in clade_definitions:
        color = random_color(h, s=0.9, l=0.4); h = h + 0.58
        try:
            subtree_str = {str(name) for name in each_subtree if name}
            if len(subtree_str) == 1:
                leaf_name = list(subtree_str)[0]
                matching_leaves = [leaf for leaf in t.get_leaves() if str(leaf.name) == leaf_name]
                if matching_leaves: node_layout(matching_leaves[0], color)
            elif len(subtree_str) > 1:
                nodes_in_subtree = [leaf for leaf in t.get_leaves() if str(leaf.name) in subtree_str]
                if len(nodes_in_subtree) >= 1:
                    Highlight_node = t.get_common_ancestor(nodes_in_subtree)
                    for node in Highlight_node.traverse():
                        if node.is_leaf():
                            if str(node.name) in subtree_str: node_layout(node, color)
                        else: node_layout(node, color)
        except Exception as e: print(f"Warning: Error processing predefined clade {each_subtree}: {e}")
    
    if node_num and isinstance(highlight_target, list) and highlight_target:
        try:
            target_leaves_str = {str(name) for name in highlight_target}
            target_nodes = [leaf for leaf in t.get_leaves() if str(leaf.name) in target_leaves_str]
            if len(target_nodes) >= 1:
                target_ancestor = t.get_common_ancestor(target_nodes)
                if target_ancestor:
                    node_face = TextFace(str(node_num), fsize=40); node_face.background.color = "LightGreen"
                    target_ancestor.add_face(node_face, column=0, position="branch-right")
        except Exception as e: print(f"Warning: Could not label node {node_num}: {e}")
    
    ts = TreeStyle(); ts.scale = 40; ts.draw_guiding_lines = True
    ts.show_leaf_name = False; ts.force_topology = True; ts.show_scale = False
    temp_tree_file = "tree_temp.png"
    try:
        num_leaves_tree = len(t.get_leaves())
        render_height = max(100, num_leaves_tree * 20)
        t.render(temp_tree_file, h=render_height, tree_style=ts, dpi=200)
    except Exception as e:
        print(f"Error: Error rendering tree: {e}")
        if os.path.exists(temp_hotmap_file): os.remove(temp_hotmap_file)
        if os.path.exists(temp_colorbar_file): os.remove(temp_colorbar_file)
        return False

    # --- Part 4: Combine Figures ---
    try:
        treepic_orig = Image.open(temp_tree_file)
        hotpic_orig = Image.open(temp_hotmap_file)
        colorbarpic_orig = Image.open(temp_colorbar_file)
    except (FileNotFoundError, Exception) as e:
        print(f"Error: Cannot find/open temporary image file for combining: {e}.")
        return False
    
    try:
        treepic_rotate_orig = treepic_orig.rotate(90, expand=True)
        tree_w_orig, tree_h_orig = treepic_orig.size
        hot_w_orig, hot_h_orig = hotpic_orig.size
        
        target_h, target_w = hot_h_orig, hot_w_orig
        
        new_tree_w = int(target_h * (tree_w_orig / tree_h_orig))
        treepic = treepic_orig.resize((new_tree_w, target_h), Image.LANCZOS)
        tree_w, tree_h = treepic.size
        
        new_rot_h = int(target_w * (treepic_rotate_orig.height / treepic_rotate_orig.width))
        treepic_rotate = treepic_rotate_orig.resize((target_w, new_rot_h), Image.LANCZOS)
        
        new_cbar_w = int(target_h * (colorbarpic_orig.width / colorbarpic_orig.height))
        colorbarpic = colorbarpic_orig.resize((new_cbar_w, target_h), Image.LANCZOS)
        cbar_w, cbar_h = colorbarpic.size
        
        padding, padding_between = 20, 10
        total_width = tree_w + target_w + cbar_w + 4 * padding
        total_height = max(tree_h, target_h + padding_between + new_rot_h) + 2 * padding
        
        combine = Image.new("RGB", (total_width, total_height), "#FFFFFF")
        
        combine.paste(treepic, (padding, padding))
        combine.paste(hotpic_orig, (padding + tree_w + padding, padding))
        combine.paste(treepic_rotate, (padding + tree_w + padding, padding + target_h + padding_between))
        combine.paste(colorbarpic, (padding + tree_w + padding + target_w + padding, padding))
        
        output_filename = fig_prefix + "_combined.png"
        combine.save(output_filename)
        print(f"Final image saved: {output_filename}")
        return True
    except Exception as e:
        print(f"Error during image combination or saving: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        for img in ['treepic_orig', 'hotpic_orig', 'colorbarpic_orig', 'treepic_rotate_orig', 'treepic', 'colorbarpic', 'treepic_rotate']:
            try: locals()[img].close()
            except (KeyError, NameError): pass
        for f_path in [temp_hotmap_file, temp_tree_file, temp_colorbar_file]:
            try:
                if f_path and os.path.exists(f_path): os.remove(f_path)
            except OSError as e: print(f"Warning: Could not delete temp file {f_path}: {e}")

# --- Main Program ---
def main():
    parser = argparse.ArgumentParser(
        description="Options for visual_hyde.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=Description
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-i', '--infile', action="store", metavar='FILE', type=str, required=True, help="HyDe output file")
    required.add_argument('-t', '--treefile', action="store", metavar='FILE', type=str, required=True, help="Species tree file (Newick)")
    
    node_args = parser.add_argument_group("Arguments for Node Mode")
    node_args.add_argument('-n', '--node', action="store_true", default=False, help='Activate Node model (uses difference consensus algorithm)')
    node_args.add_argument('-gdt', '--gamma_diff_threshold', action="store", type=float, default=0.2, metavar='F', help='LOOSE Gamma diff threshold for Stage 1 draft calculation (def: 0.2)')
    node_args.add_argument('--diff_consensus_threshold', action="store", type=float, default=0.1, metavar='F', help="Threshold for a gamma difference to be considered 'low' (def: 0.1)")
    node_args.add_argument('--diff_consensus_ratio', action="store", type=float, default=0.8, metavar='F', help="Minimum ratio of nodes with 'low' difference to form consensus (def: 0.8)")

    filter_args = parser.add_argument_group("Optional False Positive Filtering")
    filter_args.add_argument('--filter_dir', action="store", metavar='DIR', type=str, help='Directory with null-hypothesis CSVs for filtering')
    filter_args.add_argument('--filter_ratio_threshold', action="store", type=float, default=0.05, metavar='F', help='Ratio threshold for null-hypothesis filter (def: 0.05)')

    other = parser.add_argument_group("Other arguments")
    other.add_argument('-c', '--preclade', action="store", metavar='FILE', type=str, help='Predefined clades file')
    other.add_argument('-l', '--leaves', action="store", metavar='LEAF', type=str, help='(Leaf mode only) Process single leaf')
    other.add_argument('-p', '--pvalue', action="store", type=float, default=0.05, metavar='F', help='P-value threshold for significant results (def: 0.05)')
    args = parser.parse_args()

    # --- File checks and data loading ---
    if not os.path.exists(args.treefile): print(f"FATAL ERROR: Tree file not found: {args.treefile}"); sys.exit(1)
    if not os.path.exists(args.infile): print(f"FATAL ERROR: HyDe output file not found: {args.infile}"); sys.exit(1)
    
    try:
        input_tree = Tree(args.treefile)
        if not input_tree.get_leaves(): print("FATAL ERROR: Input tree file contains no leaves."); sys.exit(1)
    except Exception as e:
        print(f"FATAL ERROR: Error parsing tree file '{args.treefile}': {e}"); sys.exit(1)
    
    print(f"Reading HyDe output file: {args.infile} ...")
    try:
        hyde_data_df = pd.read_csv(args.infile, sep="\t")
        required_cols = ["P1", "Hybrid", "P2", "Pvalue", "Gamma"]
        if not all(col in hyde_data_df.columns for col in required_cols):
            missing = [col for col in required_cols if col not in hyde_data_df.columns]
            print(f"FATAL ERROR: HyDe file missing columns: {missing}"); sys.exit(1)
        hyde_data_df['Gamma'] = pd.to_numeric(hyde_data_df['Gamma'], errors='coerce')
        hyde_data_df['Pvalue'] = pd.to_numeric(hyde_data_df['Pvalue'], errors='coerce')
        initial_rows = len(hyde_data_df)
        hyde_data_df.dropna(subset=['Gamma', 'Pvalue', 'P1', 'Hybrid', 'P2'], inplace=True)
        for col in ['P1', 'Hybrid', 'P2']:
            hyde_data_df = hyde_data_df[hyde_data_df[col].astype(str).str.strip() != '']
        if len(hyde_data_df) < initial_rows:
            print(f"Warning: Dropped {initial_rows - len(hyde_data_df)} rows from HyDe data due to invalid/missing values.")
        if hyde_data_df.empty:
            print("FATAL ERROR: No valid data rows remaining in HyDe file after cleaning."); sys.exit(1)
        print("File reading and basic validation completed.")
    except Exception as e:
        print(f"FATAL ERROR: Error reading or processing HyDe file '{args.infile}': {e}"); sys.exit(1)
    
    if not check_name_consistency(input_tree, hyde_data_df):
        sys.exit(1)
    
    actual_predefined_clade_file = None
    if args.preclade:
        if os.path.exists(args.preclade):
            actual_predefined_clade_file = args.preclade
        else:
            print(f"Warning: Clade file '{args.preclade}' not found, generating automatically.")
            actual_predefined_clade_file = make_predefined_clade_file(args.treefile)
    else:
        print("Clade file not specified, generating automatically.")
        actual_predefined_clade_file = make_predefined_clade_file(args.treefile)
    
    if not actual_predefined_clade_file or not os.path.exists(actual_predefined_clade_file):
        print("FATAL ERROR: Could not find or generate clade file."); sys.exit(1)
    
    print(f"Using clade file: {actual_predefined_clade_file}")
    parse_result = parse_tree(args.treefile, actual_predefined_clade_file)
    if parse_result is None:
        print("FATAL ERROR: Failed to parse tree or clade file."); sys.exit(1)
    
    _, leaf_labels_ingroup, clade_defs, name_len_range = parse_result
    parsed_tree = input_tree
    leaf_labels_ingroup_str = [str(lbl) for lbl in leaf_labels_ingroup]
    if not leaf_labels_ingroup_str :
        print("FATAL ERROR: No ingroup labels found after parsing."); sys.exit(1)

    # --- Load and VALIDATE Null Hypothesis Data (if provided) ---
    null_heatmap_data = {}
    if args.filter_dir:
        if not os.path.isdir(args.filter_dir):
            print(f"FATAL ERROR: Filter directory not found: {args.filter_dir}")
            sys.exit(1)
        print(f"Loading and validating null hypothesis filter data from: {args.filter_dir}...")
        
        for filename in os.listdir(args.filter_dir):
            if filename.endswith(".csv"):
                safe_name = os.path.splitext(filename)[0]
                try:
                    null_df = pd.read_csv(os.path.join(args.filter_dir, filename), index_col=0)
                    if set(null_df.index) != set(leaf_labels_ingroup_str):
                        print(f"\nFATAL ERROR: Name mismatch inside filter file '{filename}'.")
                        print("The row/column names in this file do not exactly match the leaf names in the main species tree.")
                        sys.exit(1)
                    null_heatmap_data[safe_name] = null_df
                except Exception as e:
                    print(f"Warning: Could not load or parse null data file {filename}: {e}")
        print(f"Loaded and validated {len(null_heatmap_data)} filter files.")

    # --- Main Processing ---
    if args.node:
        node_heatmaps_cache, node_calculation_logs = calculate_all_node_heatmaps(
            parsed_tree, leaf_labels_ingroup_str, hyde_data_df, args.pvalue, 
            gamma_diff_threshold=args.gamma_diff_threshold, 
            diff_consensus_threshold=args.diff_consensus_threshold,
            diff_consensus_ratio=args.diff_consensus_ratio
        )
        print(f"\nRunning Node Mode using 'difference consensus' algorithm...")

        all_nodes_postorder = list(parsed_tree.traverse("postorder"))
        plottable_node_indices = []
        for i, n in enumerate(all_nodes_postorder):
            if not n.is_leaf() and not n.is_root():
                node_leaves_in_ingroup = [leaf for leaf in n.get_leaves() if str(leaf.name) in leaf_labels_ingroup_str]
                if len(node_leaves_in_ingroup) >= 2 and len(n.children) == 2:
                    plottable_node_indices.append(i)

        total_plottable_nodes = len(plottable_node_indices)
        print(f"Found {total_plottable_nodes} potentially plottable internal nodes (explicitly excluding the root node).")

        processed_nodes = 0
        current_plot_index = 0
        for node_idx, each_node in enumerate(all_nodes_postorder):
            if node_idx not in plottable_node_indices: continue

            current_plot_index += 1
            node_name_str = str(each_node.name) if hasattr(each_node, 'name') and each_node.name else f"Internal_{current_plot_index}"
            fig_prefix = f"Node_{current_plot_index}"
            print(f"\nProcessing Node {current_plot_index}/{total_plottable_nodes} ({node_name_str})...")

            calculation_log = node_calculation_logs.get(each_node)
            if calculation_log:
                # *** MODIFIED: Change extension to .tsv and use tab delimiter ***
                log_filename = f"{fig_prefix}_calculation_details.tsv"
                print(f"  Writing calculation details to {log_filename}...")
                try:
                    with open(log_filename, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f, delimiter='\t')
                        header = [
                            "Node_ID", "Parent_1", "Parent_2", "Consensus_Threshold_Ratio", "Low_Diff_Ratio_Found",
                            "Low_Diff_Points", "Total_Internal_Nodes_in_Subtree", "Low_Diff_Threshold",
                            "Final_Gamma_Value", "Result", "Subtree_Gamma_Diffs"
                        ]
                        writer.writerow(header)
                        for row in calculation_log:
                            writer.writerow([f"Node_{current_plot_index}"] + row)
                except Exception as e:
                    print(f"  Warning: Could not write calculation TSV log file: {e}")

            hyde_node_array = node_heatmaps_cache.get(each_node)
            if hyde_node_array is None:
                print(f"  Skipping plot: Heatmap data calculation failed or node not found in cache.")
                continue
            
            if args.filter_dir:
                null_node_heatmap = null_heatmap_data.get(fig_prefix)
                if null_node_heatmap is not None:
                    print(f"  Applying null hypothesis filter for node: {fig_prefix}...")
                    hyde_node_array = apply_relative_filter(hyde_node_array, null_node_heatmap, args.filter_ratio_threshold)
                else:
                    print(f"  Warning: No null hypothesis filter file found for {fig_prefix}.csv. Proceeding without filtering this node.")


            if hyde_node_array.isnull().all().all():
                print(f"  Note: No significant data found for this node. Plotting with empty heatmap.")
            
            node_ingroup_leaves = [str(leaf.name) for leaf in each_node.get_leaves() if str(leaf.name) in leaf_labels_ingroup_str]
            
            success = generate_hyde_visualization(
                fig_prefix, hyde_node_array, args.treefile, current_plot_index, 
                node_ingroup_leaves, clade_defs, name_len_range, 0.6
            )

            if success:
                processed_nodes += 1
            else:
                print(f"  Warning: Image generation failed for node {node_name_str}. Skipping.")

        print(f"\nNode mode processing completed. Generated images for {processed_nodes} internal nodes.")

    else: # Leaf Mode
        leaves_to_process = []
        hyde_hybrid_names = {str(name) for name in hyde_data_df['Hybrid'].unique()}
        if args.leaves:
            single_leaf_str = str(args.leaves)
            if single_leaf_str in leaf_labels_ingroup_str:
                if single_leaf_str in hyde_hybrid_names:
                    leaves_to_process = [single_leaf_str]
                else:
                    print(f"Warning: Specified leaf '{single_leaf_str}' ... is not found as a 'Hybrid' ... No plot generated.")
                    sys.exit(0)
            else:
                print(f"FATAL ERROR: Specified leaf '{single_leaf_str}' not found in the ingroup leaves..."); sys.exit(1)
        else:
            leaves_to_process = [leaf for leaf in leaf_labels_ingroup_str if leaf in hyde_hybrid_names]
            skipped_leaves = set(leaf_labels_ingroup_str) - set(leaves_to_process)
            if not leaves_to_process:
                print(f"FATAL ERROR: None of the ingroup leaves were found as 'Hybrid' in the HyDe file..."); sys.exit(1)
            print(f"Running Leaf Mode for {len(leaves_to_process)} leaves...")
            if skipped_leaves:
                print(f"  (Skipping {len(skipped_leaves)} ingroup leaves not listed as 'Hybrid')")

        processed_leaves_count = 0
        for i, current_leaf in enumerate(leaves_to_process):
            safe_leaf_name = get_safe_leaf_name(current_leaf)
            fig_prefix = f"{safe_leaf_name}"
            print(f"Plotting leaf node: {current_leaf} ({i+1}/{len(leaves_to_process)})...")
            
            hyde_leaf_array = make_hotmap_table_gamma(leaf_labels_ingroup_str, current_leaf, hyde_data_df, args.pvalue,
                                                      null_heatmap_data, args.filter_ratio_threshold)
            if hyde_leaf_array is None:
                print(f"  Warning: Heatmap data calculation failed for {current_leaf}. Skipping.")
                continue
            if hyde_leaf_array.isnull().all().all():
                print(f"  Note: No significant data (P < {args.pvalue}) found for {current_leaf}. Plotting with empty heatmap.")
            
            success = generate_hyde_visualization(
                fig_prefix, hyde_leaf_array, args.treefile, "", 
                current_leaf, clade_defs, name_len_range, 0.6
            )

            if success:
                processed_leaves_count += 1
            else:
                print(f"  Warning: Image generation failed for leaf {current_leaf}. Skipping.")

        print(f"Leaf mode processing completed. Generated images for {processed_leaves_count} leaves.")

if __name__ == "__main__":
    main()
