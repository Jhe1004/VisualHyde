# coding: utf-8
import os
import sys
import argparse
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, Normalize # Added Normalize
from matplotlib.colorbar import ColorbarBase # Added ColorbarBase
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces, random_color, TextFace
from PIL import Image, ImageEnhance # Added ImageEnhance if needed later, Image is main one
import warnings # Keep warnings import
import colorsys # <-- Needed for brightness adjustment

# Set Matplotlib backend and Qt platform variable (if needed)
# Using 'Agg' is generally good for scripts without a GUI
matplotlib.use('Agg')
# This QT variable might not be necessary with 'Agg' but leave it for now
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

Description = (
    '''
    -----------------------------------------------------------------------------
    |                           visual_hyde.py                                  |
    -----------------------------------------------------------------------------

    Created by Jian He (j.he930724@gamail.com)

    Dependencies:
    python3
    ete3 (conda install -c etetoolkit ete3 ete_toolchain)
    matplotlib (conda install matplotlib)
    numpy (conda install numpy)
    pandas (conda install pandas)
    Pillow (conda install Pillow or pip install Pillow) # PIL is part of Pillow
    colorsys (standard library)

    '''
)

# --- Helper Functions ---

def make_predefined_clade_file(tree_file, max_leaves_per_clade=5):
    """
    Automatically generates a simple predefined clade file based on the input tree.
    Each automatically defined clade contains no more than max_leaves_per_clade leaf nodes.

    Args:
        tree_file (str): Path to the species tree file.
        max_leaves_per_clade (int): Maximum number of leaf nodes per automatically defined clade.
    """
    try:
        t = Tree(tree_file)
        # Find ingroup node (assuming outgroup is a single leaf or smaller subtree)
        if len(t.children) != 2:
            print("Warning: Tree does not appear to be bifurcating or root position is unclear, automatic clade definition might be inaccurate.")
            if len(t.children) > 0:
                ingroup_node = max(t.children, key=len)
                print(f"Warning: Attempting to use subtree with {len(ingroup_node)} leaf nodes as ingroup for automatic clades.")
            else:
                ingroup_node = t # If no children, use the whole tree
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

    Args:
        tree_file (str): Path to the species tree file (Newick format).
        predefined_clade_file (str): Path to the predefined clade file.

    Returns:
        tuple: Containing the ETE3 Tree object, list of leaf node names (in tree order),
               list of highlight clades, and the range of name lengths. Returns None on error.
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
        name_len = (min(name_len_list) if name_len_list else 0,
                    max(name_len_list) if name_len_list else 0)

        return t, leaf_labels, highlight_subtrees, name_len
    except FileNotFoundError:
        print(f"Error: Tree file not found - {tree_file}")
        return None
    except Exception as e:
        print(f"Error: Error parsing tree or clade file: {e}")
        return None

def make_hotmap_table_gamma(leaf_labels, hypothesis_hybrid_species, hyde_data_df, zscore):
    """Generates the heatmap data table (DataFrame) for the specified hybrid."""
    hypothesis_hybrid_species_str = str(hypothesis_hybrid_species)
    leaf_labels_str = [str(lbl) for lbl in leaf_labels]
    hyde_data_df_copy = hyde_data_df.copy()
    for col in ['Hybrid', 'P1', 'P2']:
         if col in hyde_data_df_copy.columns:
             hyde_data_df_copy[col] = hyde_data_df_copy[col].astype(str)

    sub_table = hyde_data_df_copy[
        (hyde_data_df_copy["Hybrid"] == hypothesis_hybrid_species_str) &
        (hyde_data_df_copy["Zscore"].astype(float) > zscore)
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
    return df

def calculate_all_node_heatmaps(tree, leaf_labels, hyde_data_df, zscore_threshold, gamma_diff_threshold=0.2):
    """Calculates and caches the heatmap for each node using postorder traversal (v1.8 logic)."""
    node_heatmap_cache = {}
    leaf_labels_str = [str(lbl) for lbl in leaf_labels]
    all_leaf_labels_set_str = set(leaf_labels_str)

    print(f"Starting calculation of all node heatmaps (postorder traversal, gamma diff threshold={gamma_diff_threshold})...")
    processed_node_count = 0
    total_nodes = len(list(tree.traverse()))

    for node in tree.traverse("postorder"):
        processed_node_count += 1
        if processed_node_count % 10 == 0 or processed_node_count == total_nodes:
            node_id_str = str(node.name) if hasattr(node, 'name') and node.name else f"InternalNode_{processed_node_count}"
            if node.is_leaf(): node_id_str = str(node.name)
            print(f"  Processing node {processed_node_count}/{total_nodes} ({node_id_str})...")

        node_descendants_str = set(str(n) for n in node.get_leaf_names())

        if node.is_leaf():
            node_name_str = str(node.name)
            if node_name_str in all_leaf_labels_set_str:
                leaf_heatmap = make_hotmap_table_gamma(leaf_labels_str, node_name_str, hyde_data_df, zscore_threshold)
                node_heatmap_cache[node] = leaf_heatmap
            else:
                nan_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                node_heatmap_cache[node] = nan_heatmap
        else:
            children = node.children
            if len(children) == 2:
                child1, child2 = children
                child1_heatmap = node_heatmap_cache.get(child1)
                child2_heatmap = node_heatmap_cache.get(child2)

                if child1_heatmap is None or child2_heatmap is None:
                    node_name_str = str(node.name) if hasattr(node,'name') and node.name else 'Unnamed'
                    print(f"Warning: Child heatmap not found for internal node {node_name_str}, skipping.")
                    nan_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                    node_heatmap_cache[node] = nan_heatmap
                    continue

                current_node_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                for p1_name in leaf_labels_str:
                    if p1_name in node_descendants_str: continue
                    for p2_name in leaf_labels_str:
                        if p2_name in node_descendants_str or p1_name == p2_name: continue
                        try:
                            gamma1 = child1_heatmap.loc[p1_name, p2_name]
                            gamma2 = child2_heatmap.loc[p1_name, p2_name]
                        except KeyError: gamma1, gamma2 = np.nan, np.nan

                        avg_gamma = np.nan
                        try:
                            if pd.notna(gamma1) and pd.notna(gamma2):
                                if abs(gamma1 - gamma2) < gamma_diff_threshold:
                                    avg_gamma = (gamma1 + gamma2) / 2.0
                        except TypeError: avg_gamma = np.nan

                        if pd.notna(avg_gamma):
                            current_node_heatmap.loc[p1_name, p2_name] = avg_gamma

                node_heatmap_cache[node] = current_node_heatmap
            else:
                node_name_str = str(node.name) if hasattr(node,'name') and node.name else 'Unnamed'
                print(f"Warning: Internal node {node_name_str} is not binary ({len(children)} children), skipping.")
                nan_heatmap = pd.DataFrame(np.nan, index=leaf_labels_str, columns=leaf_labels_str, dtype=float)
                node_heatmap_cache[node] = nan_heatmap

    print("All node heatmap calculations completed.")
    return node_heatmap_cache

# --- Plotting Functions ---
# --- Draw heatmap WITHOUT colorbar (Dynamic Size) ---
def draw_hotmap(pic_name, hyde_output_array):
    """
    Draws the heatmap WITHOUT the colorbar, using a dynamic size based on
    the number of leaves and the original saturation.
    Returns the path to the saved heatmap image, the ORIGINAL colormap,
    and the normalization object.
    """
    # Save data to CSV
    try:
        hyde_output_array.to_csv(pic_name + ".csv", index=True, sep=',')
    except Exception as e:
        print(f"Warning: Could not save heatmap data to CSV '{pic_name}.csv': {e}")

    # --- Calculate dynamic figure size ---
    num_leaves = hyde_output_array.shape[0] # Get number of rows (leaves)
    if num_leaves == 0:
        print("Warning: Heatmap data is empty, cannot draw heatmap.")
        return None, None, None

    # --- Parameters for dynamic sizing (ADJUST THESE AS NEEDED) ---
    inches_per_leaf = 0.4  # How many inches per leaf in the heatmap
    base_inches = 2.0      # A minimum base size in inches
    max_inches = 30.0      # Maximum size in inches (original value)
    # --- End Parameters ---

    # Calculate the figure size, ensuring it doesn't exceed max_inches
    fig_inches = min(max_inches, base_inches + num_leaves * inches_per_leaf)
    print(f"Calculated heatmap figsize: {fig_inches:.2f} x {fig_inches:.2f} inches for {num_leaves} leaves.")

    # Create figure and axes with dynamic size
    fig = plt.figure(figsize=(fig_inches, fig_inches)) # <-- USE DYNAMIC SIZE
    border_width = 0.00001 # Keep border calculation relative
    ax_size = [0 + border_width, 0 + border_width,
               1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig.add_axes(ax_size)

    # --- Create ORIGINAL colormap (original logic) ---
    # (Colormap creation logic remains the same)
    cmap_list = []
    color_val = 1.0
    lucency = 0.0
    steps = 5000 # Keep high resolution for smoothness
    for _ in range(steps):
        cmap_list.append([0, 0, color_val, lucency]) # Blue to transparent
        color_val -= (1.0 / steps)
        lucency += (1.0 / steps)

    color_val = 0.0
    lucency = 1.0
    for _ in range(steps):
        cmap_list.append([color_val, 0, 0, lucency]) # Transparent to red
        color_val += (1.0 / steps)
        lucency -= (1.0 / steps)

    cmap_array = np.clip(np.array(cmap_list), 0, 1)
    original_cmap = ListedColormap(cmap_array)
    original_cmap.set_bad(color='white', alpha=0) # Set NaN color

    # Create normalization object
    norm = Normalize(vmin=0, vmax=1)

    # Draw the image using the ORIGINAL colormap
    im = ax.imshow(hyde_output_array.astype(float),
                   norm=norm,
                   cmap=original_cmap,
                   interpolation='nearest', # 'nearest' is good for discrete blocks
                   aspect='equal')

    # Add grid lines - Adjust linewidth based on figure size? Optional.
    grid_linewidth = max(0.5, min(2, fig_inches / 10.0)) # Optional scaling
    ax.set_xticks(np.arange(hyde_output_array.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(hyde_output_array.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=grid_linewidth) # Use scaled linewidth

    # Hide axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Save heatmap to temporary file
    temp_hotmap_file = pic_name + "_hotmap_temp.png"
    try:
        # Keep DPI consistent, pixel size will now vary with fig_inches
        plt.savefig(temp_hotmap_file, dpi=200)
        print(f"Saved temporary heatmap: {temp_hotmap_file} ({(fig_inches*200):.0f}x{(fig_inches*200):.0f} px)")
    except Exception as e:
        print(f"Error: Error saving heatmap: {e}")
        temp_hotmap_file = None
    finally:
        # Ensure matplotlib resources are released
        plt.cla()
        plt.clf()
        plt.close(fig)
        plt.close("all") # Close any other potential figures

    # Return temporary filename, the ORIGINAL colormap, and normalization object
    return temp_hotmap_file, original_cmap, norm

# --- Draw only the colorbar ---
def draw_colorbar_only(colorbar_pic_name, cmap, norm):
    """Draws only the colorbar to a separate image file."""
    fig = plt.figure(figsize=(2, 10)) # Tall and narrow
    ax = fig.add_axes([0.2, 0.05, 0.2, 0.9]) # Position for the bar
    try:
        cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
        cb.ax.tick_params(labelsize=40) # Label size
        plt.savefig(colorbar_pic_name, dpi=200, bbox_inches='tight', pad_inches=0.1)
        return colorbar_pic_name
    except Exception as e:
        print(f"Error: Error saving colorbar image '{colorbar_pic_name}': {e}")
        return None
    finally:
        plt.cla()
        plt.clf()
        plt.close(fig)
        plt.close("all")

# --- Draw tree (Increased Height Factor) ---
def draw_tree(tree_file, node_num, highlight_target, clade_definitions, name_len, picture_size): # picture_size unused currently
    """Draws the tree (original logic + increased height factor)."""
    try:
        t = Tree(tree_file)
    except Exception as e:
        print(f"Error: Could not reload tree file {tree_file} in draw_tree: {e}")
        return None

    def node_layout_background(node):
        ns = NodeStyle(); ns["hz_line_width"]=4.5; ns["vt_line_width"]=4.5; ns["size"]=0; node.set_style(ns)

    def node_layout(node, color):
        if node.is_leaf():
            node_name_str = str(node.name) if node.name is not None else ""
            padding_dots = "·" * ((name_len[1] - len(node_name_str)) * 2 + 10)
            leaf_name_text = node_name_str + padding_dots
            # Consider reducing fsize if height is still an issue
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
            else: pass
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
        num_leaves = len(t.get_leaves())
        # *** INCREASED HEIGHT MULTIPLIER (adjust 60 if needed) ***
        render_height = max(100, num_leaves * 20)
        print(f"Rendering tree with estimated height: {render_height}")
        t.render(temp_tree_file, h=render_height, tree_style=ts, dpi=200)
        return temp_tree_file
    except Exception as e:
        print(f"Error: Error rendering tree: {e}")
        return None

# --- Combine Figure (Rotated Tree Below Heatmap, Backward Compatible Resize) ---
def combine_fig(fig_name, tree_img_path, hotmap_img_path, colorbar_img_path):
    """
    Combines the species tree image (left), the heatmap image (top-middle),
    a rotated tree (bottom-middle), and the separate colorbar image (right).
    Uses older Pillow resize syntax for backward compatibility (Pillow < 9.0).
    """
    try:
        treepic_orig = Image.open(tree_img_path)
        hotpic_orig = Image.open(hotmap_img_path)
        colorbarpic_orig = Image.open(colorbar_img_path)
    except FileNotFoundError as e:
        print(f"Error: Cannot find temporary image file for combining: {e}. Cannot combine images.")
        return
    except Exception as e:
        print(f"Error: Error opening temporary images: {e}")
        return

    try:
        # --- Rotate tree ---
        treepic_rotate_orig = treepic_orig.rotate(90, expand=True)

        # --- Get original dimensions ---
        tree_w_orig, tree_h_orig = treepic_orig.size
        hot_w_orig, hot_h_orig = hotpic_orig.size
        rot_w_orig, rot_h_orig = treepic_rotate_orig.size
        cbar_w_orig, cbar_h_orig = colorbarpic_orig.size

        # --- Resizing Strategy: Heatmap drives alignment ---
        target_h = hot_h_orig
        target_w = hot_w_orig

        # Resize Left Tree - match heatmap height
        aspect_ratio_tree = tree_w_orig / tree_h_orig
        new_tree_w = int(target_h * aspect_ratio_tree)
        treepic = treepic_orig.resize((new_tree_w, target_h), Image.LANCZOS) # Use older LANCZOS
        tree_w, tree_h = treepic.size
        print(f"Resized left tree to {tree_w}x{tree_h}")

        # Resize Rotated Tree - match heatmap width
        aspect_ratio_rot = rot_h_orig / rot_w_orig
        new_rot_h = int(target_w * aspect_ratio_rot)
        treepic_rotate = treepic_rotate_orig.resize((target_w, new_rot_h), Image.LANCZOS) # Use older LANCZOS
        rot_w, rot_h = treepic_rotate.size
        print(f"Resized rotated tree to {rot_w}x{rot_h}")

        # Resize Colorbar - match heatmap height
        aspect_ratio_cbar = cbar_w_orig / cbar_h_orig
        new_cbar_w = int(target_h * aspect_ratio_cbar)
        colorbarpic = colorbarpic_orig.resize((new_cbar_w, target_h), Image.LANCZOS) # Use older LANCZOS
        cbar_w, cbar_h = colorbarpic.size
        print(f"Resized colorbar to {cbar_w}x{cbar_h}")

        # Keep original heatmap
        hotpic = hotpic_orig
        hot_w, hot_h = hotpic.size

        # --- Define padding ---
        padding = 20 # General padding
        padding_between = 10 # Between heatmap and rotated tree

        # --- Calculate total canvas size ---
        total_width = tree_w + hot_w + cbar_w + 4 * padding
        middle_col_h = hot_h + padding_between + rot_h
        total_height = max(tree_h, middle_col_h) + 2 * padding

        # Create canvas
        combine = Image.new("RGB", (total_width, total_height), "#FFFFFF")
        print(f"Creating canvas: {total_width}x{total_height}")

        # --- Calculate Paste Coordinates ---
        paste_x_hot = padding + tree_w + padding
        paste_y_hot = padding
        print(f"Pasting heatmap at ({paste_x_hot}, {paste_y_hot}) size {hotpic.size}")
        combine.paste(hotpic, (paste_x_hot, paste_y_hot))

        paste_x_tree = padding
        paste_y_tree = padding
        print(f"Pasting tree at ({paste_x_tree}, {paste_y_tree}) size {treepic.size}")
        combine.paste(treepic, (paste_x_tree, paste_y_tree))

        paste_x_rot = paste_x_hot
        paste_y_rot = paste_y_hot + hot_h + padding_between
        print(f"Pasting rotated tree at ({paste_x_rot}, {paste_y_rot}) size {treepic_rotate.size}")
        combine.paste(treepic_rotate, (paste_x_rot, paste_y_rot))

        paste_x_cbar = paste_x_hot + hot_w + padding
        paste_y_cbar = padding
        print(f"Pasting colorbar at ({paste_x_cbar}, {paste_y_cbar}) size {colorbarpic.size}")
        combine.paste(colorbarpic, (paste_x_cbar, paste_y_cbar))

        # Save image
        output_filename = fig_name + ".png"
        combine.save(output_filename)
        print(f"Final image saved: {output_filename}")

    except Exception as e:
        print(f"Error during image combination or saving: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Clean up temporary files
        try: treepic_orig.close()
        except: pass
        try: hotpic_orig.close()
        except: pass
        try: colorbarpic_orig.close()
        except: pass
        try: treepic_rotate_orig.close()
        except: pass
        try: treepic.close()
        except: pass
        try: hotpic.close() # hotpic is just hotpic_orig now
        except: pass
        try: colorbarpic.close()
        except: pass
        try: treepic_rotate.close()
        except: pass

        for f_path in [hotmap_img_path, tree_img_path, colorbar_img_path]:
            try:
                if f_path and os.path.exists(f_path): os.remove(f_path)
            except OSError as e: print(f"Warning: Could not delete temp file {f_path}: {e}")

# --- Name Consistency Check Function ---
def check_name_consistency(tree, hyde_df):
    """Checks name consistency between tree leaves and HyDe participants."""
    leaves_name_in_species_tree = set()
    try:
        # Simplified ingroup logic for checking (adjust if needed)
        if len(tree.children) == 2 and len(tree.children[0]) != len(tree.children[1]):
             ingroup_node = max(tree.children, key=len)
        else: ingroup_node = tree # Check all if ambiguous/not bifurcating
        leaves_name_in_species_tree = set(str(leaf.name) for leaf in ingroup_node.get_leaves() if leaf.name is not None)
    except Exception as e:
        print(f"Error extracting leaf names from tree: {e}"); return False
    if not leaves_name_in_species_tree:
        print("Error: Could not extract any valid leaf node names from tree."); return False

    hyde_names = set()
    try:
        for col in ['P1', 'Hybrid', 'P2']:
            if col in hyde_df.columns:
                valid_names = {str(name) for name in hyde_df[col].dropna().unique() if name}
                hyde_names.update(valid_names)
    except Exception as e: print(f"Error extracting names from HyDe data: {e}"); return False
    if not hyde_names: print("Error: No valid names (P1/Hybrid/P2) from HyDe file."); return False

    if leaves_name_in_species_tree.issuperset(hyde_names): # Allow tree to have extra leaves (outgroup)
        print("Name consistency check passed: All HyDe participants found in tree leaves.")
        missing_in_hyde = leaves_name_in_species_tree - hyde_names
        if missing_in_hyde: print(f"  (Tree leaves not in HyDe: {missing_in_hyde})")
        return True
    else:
        print("Warning: Name mismatch!")
        missing_in_tree = hyde_names - leaves_name_in_species_tree
        if missing_in_tree: print(f"  Names in HyDe but not in tree: {missing_in_tree}")
        if not leaves_name_in_species_tree.intersection(hyde_names):
             print("Error: No common names found. Cannot proceed."); return False
        print("Proceeding with common names, but be aware of mismatches.")
        return True # Allow proceeding with partial match






# --- Main Program ---
def main():
    parser = argparse.ArgumentParser(
        description="Options for visual_hyde.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # epilog=Description # Make sure 'Description' variable is defined if used
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-i', '--infile', action="store", metavar='FILE', type=str, required=True, help="HyDe output file")
    required.add_argument('-t', '--treefile', action="store", metavar='FILE', type=str, required=True, help="Species tree file (Newick)")

    additional = parser.add_argument_group("Additional arguments")
    additional.add_argument('-n', '--node', action="store_true", default=False, help='Node model')
    additional.add_argument('-gdt', '--gamma_diff_threshold', action="store", type=float, default=0.2, metavar='F', help='Gamma diff threshold for node avg (def: 0.2)')
    additional.add_argument('-c', '--preclade', action="store", metavar='FILE', type=str, help='Predefined clades file')
    additional.add_argument('-l', '--leaves', action="store", metavar='LEAF', type=str, help='(Leaf mode only) Process single leaf')
    additional.add_argument('-z', '--zscore', action="store", type=float, default=3.0, metavar='F', help='Z-score threshold (def: 3.0)')
    # additional.add_argument('-cb', '--cbright', action="store", type=float, default=0.6, metavar='F', help='Colorbar brightness factor (def: 0.6)') # Example if brightness becomes an argument

    args = parser.parse_args()

    csv_file_name = args.infile
    tree_file = args.treefile
    node_model = args.node
    gamma_diff_threshold = args.gamma_diff_threshold
    predefined_clade_file_arg = args.preclade
    single_leaf_to_process = args.leaves
    zscore_threshold = args.zscore
    # colorbar_brightness_factor = args.cbright # Use if it becomes an arg
    colorbar_brightness_factor = 0.6 # Define the brightness factor

    # --- File Checks ---
    if not os.path.exists(tree_file):
        print(f"Error: Tree file not found: {tree_file}"); sys.exit(1)
    if not os.path.exists(csv_file_name):
        print(f"Error: HyDe output file not found: {csv_file_name}"); sys.exit(1)

    # --- Parse Tree (initial for checks) ---
    try:
        input_tree = Tree(tree_file)
        if not input_tree.get_leaves():
            print("Error: Input tree file contains no leaves."); sys.exit(1)
    except Exception as e:
        print(f"Error: Error parsing tree file '{tree_file}': {e}"); sys.exit(1)

    # --- Read HyDe Data ---
    print(f"Reading HyDe output file: {csv_file_name} ...")
    try:
        hyde_data_df = pd.read_csv(csv_file_name, sep="\t")
        required_cols = ["P1", "Hybrid", "P2", "Zscore", "Gamma"]
        if not all(col in hyde_data_df.columns for col in required_cols):
            missing = [col for col in required_cols if col not in hyde_data_df.columns]
            print(f"Error: HyDe file missing columns: {missing}"); sys.exit(1)

        # Convert necessary columns and handle NaNs/empty strings
        hyde_data_df['Gamma'] = pd.to_numeric(hyde_data_df['Gamma'], errors='coerce')
        hyde_data_df['Zscore'] = pd.to_numeric(hyde_data_df['Zscore'], errors='coerce')
        initial_rows = len(hyde_data_df)
        hyde_data_df.dropna(subset=['Gamma', 'Zscore', 'P1', 'Hybrid', 'P2'], inplace=True)
        for col in ['P1', 'Hybrid', 'P2']:
            # Ensure names are treated as strings and remove rows with empty strings
            hyde_data_df = hyde_data_df[hyde_data_df[col].astype(str).str.strip() != '']
        if len(hyde_data_df) < initial_rows:
            print(f"Warning: Dropped {initial_rows - len(hyde_data_df)} rows from HyDe data due to invalid/missing values.")
        if hyde_data_df.empty:
            print("Error: No valid data rows remaining in HyDe file after cleaning."); sys.exit(1)
        print("File reading and basic validation completed.")
    except Exception as e:
        print(f"Error reading or processing HyDe file '{csv_file_name}': {e}"); sys.exit(1)

    # --- Check Names ---
    if not check_name_consistency(input_tree, hyde_data_df):
        print("Exiting due to name inconsistency."); sys.exit(1)

    # --- Handle Clade File ---
    actual_predefined_clade_file = None
    if predefined_clade_file_arg:
        if os.path.exists(predefined_clade_file_arg):
            actual_predefined_clade_file = predefined_clade_file_arg
        else:
            print(f"Warning: Clade file '{predefined_clade_file_arg}' not found, generating automatically.")
            actual_predefined_clade_file = make_predefined_clade_file(tree_file)
    else:
        print("Clade file not specified, generating automatically.")
        actual_predefined_clade_file = make_predefined_clade_file(tree_file)

    if not actual_predefined_clade_file or not os.path.exists(actual_predefined_clade_file):
        print("Error: Could not find or generate clade file."); sys.exit(1)
    print(f"Using clade file: {actual_predefined_clade_file}")


    # --- Parse Tree and Clades (main) ---
    parse_result = parse_tree(tree_file, actual_predefined_clade_file)
    if parse_result is None:
        print("Error: Failed to parse tree or clade file."); sys.exit(1)
    _, leaf_labels_ingroup, clade_defs, name_len_range = parse_result
    parsed_tree = input_tree # Use the tree loaded earlier
    leaf_labels_ingroup_str = [str(lbl) for lbl in leaf_labels_ingroup]
    if not leaf_labels_ingroup_str :
        print("Error: No ingroup labels found after parsing."); sys.exit(1)

    # --- Ensure colorsys is available ---
    try:
        colorsys
    except NameError:
        print("Error: colorsys module not imported!"); sys.exit(1)

    # --- Main Processing ---
    if node_model:
        # --- Node Mode Loop ---
        node_heatmaps_cache = calculate_all_node_heatmaps(parsed_tree, leaf_labels_ingroup_str, hyde_data_df, zscore_threshold, gamma_diff_threshold)
        print(f"Running Node Mode ...")

        all_nodes_postorder = list(parsed_tree.traverse("postorder"))
        # Identify structurally suitable nodes first, EXCLUDING the root
        plottable_node_indices = []
        for i, n in enumerate(all_nodes_postorder):
            # *** MODIFIED CONDITION TO EXCLUDE ROOT ***
            # Check if it's not a leaf AND explicitly not the root node
            if not n.is_leaf() and not n.is_root():
                node_leaves_in_ingroup = [leaf for leaf in n.get_leaves() if str(leaf.name) in leaf_labels_ingroup_str]
                # Criteria: Internal node (not root), >=2 ingroup descendants, binary structure assumed by calculation logic
                if len(node_leaves_in_ingroup) >= 2 and len(n.children) == 2:
                    plottable_node_indices.append(i)

        total_plottable_nodes = len(plottable_node_indices)
        # Updated print statement for clarity
        print(f"Found {total_plottable_nodes} potentially plottable internal nodes (explicitly excluding the root node).")

        processed_nodes = 0
        current_plot_index = 0 # Tracks the index among plottable nodes
        for node_idx, each_node in enumerate(all_nodes_postorder):
            # Skip nodes not in our pre-filtered list (which now excludes the root)
            if node_idx not in plottable_node_indices:
                continue

            # --- No extra root check needed here anymore ---

            current_plot_index += 1 # Increment only for non-root, plottable nodes
            node_name_str = str(each_node.name) if hasattr(each_node, 'name') and each_node.name else f"Internal_{current_plot_index}"
            fig_prefix = f"Node_{current_plot_index}"
            # Slightly adjusted print statement to reflect the filtered count
            print(f"Processing Node {current_plot_index}/{total_plottable_nodes} ({node_name_str})...")

            # Retrieve calculated heatmap data for this node
            hyde_node_array = node_heatmaps_cache.get(each_node)

            # --- MODIFIED CHECK FOR NODE MODE ---
            # Skip ONLY if data calculation failed (is None), otherwise proceed
            if hyde_node_array is None:
                print(f"  Skipping plot: Heatmap data calculation failed or node not found in cache.")
                continue

            # Check if data is empty, print message but DO NOT SKIP
            is_empty_heatmap = hyde_node_array.isnull().all().all()
            if is_empty_heatmap:
                print(f"  Note: No significant data found for this node. Plotting with empty heatmap.")
                # No 'continue' here - proceed to plot
            # --- END MODIFIED CHECK ---


            # --- Plotting Steps (executed even if heatmap is empty) ---

            # 1. Draw heatmap (using dynamic sizing version, handles empty array)
            temp_hotmap, cmap_for_heatmap, norm = draw_hotmap(fig_prefix, hyde_node_array)
            if not temp_hotmap:
                print(f"  Warning: Heatmap drawing failed for {node_name_str}. Skipping combine.")
                continue # Skip if drawing itself failed

            # 2. Create Darker Colormap
            print(f"  Creating darker colorbar (brightness factor: {colorbar_brightness_factor})...")
            darker_colors_rgba = []
            cmap_for_colorbar = None # Initialize
            if cmap_for_heatmap: # Check if base cmap was created
                heatmap_colors_rgba = cmap_for_heatmap.colors
                for r, g, b, a in heatmap_colors_rgba:
                    try:
                        h, s, v = colorsys.rgb_to_hsv(r, g, b)
                        v_new = max(0.0, min(1.0, v * colorbar_brightness_factor))
                        new_r, new_g, new_b = colorsys.hsv_to_rgb(h, s, v_new)
                        darker_colors_rgba.append([new_r, new_g, new_b, a])
                    except Exception as e:
                        print(f"  Warning: Color darkening failed for one color ({r},{g},{b}): {e}. Using original.")
                        darker_colors_rgba.append([r, g, b, a])
                cmap_for_colorbar = ListedColormap(darker_colors_rgba)
                cmap_for_colorbar.set_bad(color='white', alpha=0)
            else:
                print("  Warning: Cannot create darker colormap because base cmap failed. Skipping colorbar.")

            # 3. Draw colorbar (only if cmap is valid)
            temp_colorbar = None
            if cmap_for_colorbar and norm:
                temp_colorbar = draw_colorbar_only(f"{fig_prefix}_colorbar_temp.png", cmap_for_colorbar, norm)

            if not temp_colorbar:
                print(f"  Warning: Colorbar drawing failed or skipped for {node_name_str}. Skipping combine.")
                # Clean up heatmap if it exists, then skip combining
                if temp_hotmap and os.path.exists(temp_hotmap): os.remove(temp_hotmap)
                continue

            # 4. Draw tree
            # Get leaves belonging to this node to highlight the clade
            node_ingroup_leaves = [str(leaf.name) for leaf in each_node.get_leaves() if str(leaf.name) in leaf_labels_ingroup_str]
            temp_tree = draw_tree(tree_file, current_plot_index, node_ingroup_leaves, clade_defs, name_len_range, 0)
            if not temp_tree:
                print(f"  Warning: Tree drawing failed for {node_name_str}. Skipping combine.")
                if temp_hotmap and os.path.exists(temp_hotmap): os.remove(temp_hotmap)
                if temp_colorbar and os.path.exists(temp_colorbar): os.remove(temp_colorbar)
                continue

            # 5. Combine
            print(f"  Combining images for {fig_prefix}...")
            combine_fig(fig_prefix + "_combined", temp_tree, temp_hotmap, temp_colorbar)
            processed_nodes += 1 # Count successful combinations

        print(f"Node mode processing completed. Generated images for {processed_nodes} internal nodes (root excluded).")

    else: # Leaf Mode (remains unchanged from previous version - plots all specified leaves)
        # --- Leaf Mode Loop ---
        leaves_to_process = []
        hyde_hybrid_names = {str(name) for name in hyde_data_df['Hybrid'].unique()}

        if single_leaf_to_process:
            single_leaf_str = str(single_leaf_to_process)
            if single_leaf_str in leaf_labels_ingroup_str:
                if single_leaf_str in hyde_hybrid_names:
                    leaves_to_process = [single_leaf_str]
                else:
                    print(f"Warning: Specified leaf '{single_leaf_str}' exists in tree ingroup but is not found as a 'Hybrid' in the HyDe file '{csv_file_name}'. No plot generated.")
                    sys.exit(0)
            else:
                print(f"Error: Specified leaf '{single_leaf_str}' not found in the determined ingroup leaves of the tree '{tree_file}'."); sys.exit(1)
        else:
            leaves_to_process = [leaf for leaf in leaf_labels_ingroup_str if leaf in hyde_hybrid_names]
            skipped_leaves = set(leaf_labels_ingroup_str) - set(leaves_to_process)
            if not leaves_to_process:
                print(f"Error: None of the {len(leaf_labels_ingroup_str)} ingroup leaves were found as 'Hybrid' in the HyDe file '{csv_file_name}'. Cannot generate plots in Leaf Mode."); sys.exit(1)
            print(f"Running Leaf Mode for {len(leaves_to_process)} leaves (found as 'Hybrid' in HyDe data)...")
            if skipped_leaves:
                print(f"  (Skipping {len(skipped_leaves)} ingroup leaves not listed as 'Hybrid' in HyDe file)")


        total_leaves_to_process = len(leaves_to_process)
        processed_leaves_count = 0
        for i, current_leaf in enumerate(leaves_to_process):
            fig_num_prefix = f"{i+1}_" if not single_leaf_to_process else ""
            safe_leaf_name = "".join(c if c.isalnum() else "_" for c in current_leaf)
            fig_prefix = f"{fig_num_prefix}{safe_leaf_name}"
            print(f"Plotting leaf node: {current_leaf} ({i+1}/{total_leaves_to_process})...")

            # 0. Calculate heatmap data
            hyde_leaf_array = make_hotmap_table_gamma(leaf_labels_ingroup_str, current_leaf, hyde_data_df, zscore_threshold)

            if hyde_leaf_array is None:
                print(f"  Warning: Heatmap data calculation failed unexpectedly for {current_leaf}. Skipping plot.")
                continue

            # Check if data is empty, print message but DO NOT SKIP
            if hyde_leaf_array.isnull().all().all():
                print(f"  Note: No significant data (Z > {zscore_threshold}) found for {current_leaf}. Plotting with empty heatmap.")
                # No 'continue' here

            # 1. Draw heatmap
            temp_hotmap, cmap_for_heatmap, norm = draw_hotmap(fig_prefix, hyde_leaf_array)
            if not temp_hotmap:
                print(f"  Warning: Heatmap drawing failed for {current_leaf}. Skipping combine.")
                continue

            # 2. Create Darker Colormap
            print(f"  Creating darker colorbar (brightness factor: {colorbar_brightness_factor})...")
            darker_colors_rgba = []
            cmap_for_colorbar = None
            if cmap_for_heatmap:
                heatmap_colors_rgba = cmap_for_heatmap.colors
                for r, g, b, a in heatmap_colors_rgba:
                    try:
                        h, s, v = colorsys.rgb_to_hsv(r, g, b)
                        v_new = max(0.0, min(1.0, v * colorbar_brightness_factor))
                        new_r, new_g, new_b = colorsys.hsv_to_rgb(h, s, v_new)
                        darker_colors_rgba.append([new_r, new_g, new_b, a])
                    except Exception as e:
                        print(f"  Warning: Color darkening failed for one color ({r},{g},{b}): {e}. Using original.")
                        darker_colors_rgba.append([r, g, b, a])
                cmap_for_colorbar = ListedColormap(darker_colors_rgba)
                cmap_for_colorbar.set_bad(color='white', alpha=0)
            else:
                print("  Warning: Cannot create darker colormap because base cmap failed. Skipping colorbar.")

            # 3. Draw colorbar
            temp_colorbar = None
            if cmap_for_colorbar and norm:
                temp_colorbar = draw_colorbar_only(f"{fig_prefix}_colorbar_temp.png", cmap_for_colorbar, norm)

            if not temp_colorbar:
                print(f"  Warning: Colorbar drawing failed or skipped for {current_leaf}. Skipping combine.")
                if temp_hotmap and os.path.exists(temp_hotmap): os.remove(temp_hotmap)
                continue

            # 4. Draw tree
            temp_tree = draw_tree(tree_file, "", current_leaf, clade_defs, name_len_range, 0)
            if not temp_tree:
                print(f"  Warning: Tree drawing failed for {current_leaf}. Skipping combine.")
                if temp_hotmap and os.path.exists(temp_hotmap): os.remove(temp_hotmap)
                if temp_colorbar and os.path.exists(temp_colorbar): os.remove(temp_colorbar)
                continue

            # 5. Combine
            print(f"  Combining images for {fig_prefix}...")
            combine_fig(fig_prefix + "_combined", temp_tree, temp_hotmap, temp_colorbar)
            processed_leaves_count += 1

        print(f"Leaf mode processing completed. Generated images for {processed_leaves_count} leaves.")


if __name__ == "__main__":
    main()
