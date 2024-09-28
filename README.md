
# visual_hyde.py

Created by Jian He (j.he930724@gmail.com)                      
This is version 1.0 created on 4th May 2021                                                                                                                            
                                                                                        
## Dependencies
- Python 3.x
- ete3 (conda install -c etetoolkit ete3 ete_toolchain)
- matplotlib (conda install matplotlib)
- numpy (conda install numpy)
- pandas (conda install pandas)
  
 
## Description

**HyDe** (Blischak et al., 2018) is a Python package used to detect hybridization events using phylogenetic invariants based on the coalescent model. To visualize the HyDe results, this script generates heatmaps (Figure a) for each potential hybrid sample ("H"). In these heatmaps, the parental taxa ("P1" and "P2") are placed on the x and y axes, respectively, with all taxa arranged phylogenetically. Each square in the heatmap represents a HyDe test conducted in a 4-taxon tree, “((P1, H, P2), O)”. The color of each square corresponds to different γ values, representing the genetic contribution from the parent on the x or y axis (Figure a).

However, a positive HyDe result for "H" being a hybrid of "P1" and "P2" may not solely indicate a recent hybridization event but could also reflect an ancient hybridization event at deeper phylogenetic nodes (Blischak et al., 2018). Specifically, if a deeper node before “P1” was involved in a hybridization event, all descendants of that node (including "P1") would detect the same hybridization signal, and the corresponding γ values would be similar. Therefore, in the heatmap, if squares with similar colors (indicating similar γ values) form a larger block, and the corresponding taxa on the axes form a monophyletic clade, the most recent common ancestor (MRCA) of this clade is likely to be one of the true parents of "H" (Figure a). In practice, when multiple squares in a heatmap show similar γ values and form a larger block, the average γ value of these squares is calculated. This average γ value reflects the contribution from the MRCA of the clade associated with the larger block.

Furthermore, if the MRCA of a clade was itself the result of a hybridization event, all taxa within this clade should also exhibit the same hybridization signal. To illustrate this, stacked heatmaps were generated for each clade by combining individual heatmaps of all samples within that clade (Figure b). In these stacked heatmaps, only squares that were consistently light across all individual heatmaps remain light. Thus, light squares in the stacked heatmap indicate that the hybridization event occurred at the MRCA of the clade (Figures c and d). 

In practice, the final γ value representing the genetic contribution from the MRCA of the clade is calculated by summing the average γ values across all small squares in all samples within the clade, using the formula:

![Gamma Formula](https://github.com/Jhe1004/VisualHyde/blob/main/gamma_formula.png)

where:
where M is the number of samples within the clade, Nj is the number of small squares for the j-th sample, and γij is the γ value of the i-th square for the j-th sample. By averaging these values, the final γ value accurately represents the genetic contribution from the MRCA of the clade.
  

 ![Image text](https://github.com/Jhe1004/VisualHyde/blob/main/figure_1.jpg)
> Figure Captions: Schematic diagram of hybridization detection using HyDe by drawing heatmaps. (a) A heatmap for a potential hybrid sample is shown. The offspring sample (Hybrid Clade Species 1) is marked with a red star on the cladogram. The small colored squares in the heatmap indicate the corresponding hybridization signal detected by HyDe. The parents of this sample (Hybrid Clade Species 1) are represented by the corresponding samples on the horizontal (Parental Clade 1 Species 1-2) and vertical (Parental Clade 2 Species 1-5) axes. The shade of the squares indicates the inheritance probabilities of the corresponding taxa on the left axis. If small squares with similar shades form a larger square in the heatmap, and their corresponding taxa on the coordinate axis form a monophyletic clade, then the most recent common ancestor (MRCA) of this clade may be indicated as one of the parents. (b) Stacking all the heatmaps for each sample (Hybrid Clade Species 1-3) of the clade. Only if all the stacked squares are colorful will the corresponding square after stacking be colorful. (c) The stacked heatmap. The colorful squares indicate that the hybridization event occurred at the MRCA of these samples (Hybrid Clade Species 1-3) stacked in (b). (d) The schematic diagram of a phylogenetic network drawn based on (c). The numbers to the right of the red branches in the plot indicate the inheritance probabilities (γ).

## Features

- **Leaf Mode**: Generates a heatmap for each leaf species that shows potential hybridization events.
- **Node Mode**: Stacks up heatmaps for monophyletic clades to highlight shared hybridization events.
- **Customizable clades**: Option to specify predefined clades, or allow automatic grouping of up to five leaves per clade.
- **Output**: Combined heatmap and phylogenetic tree in a single image.

## Requirements

- Python 3.x
- **ete3**: Toolkit for phylogenetic tree analysis. Install with `conda install -c etetoolkit ete3 ete_toolchain`
- **matplotlib**: Plotting library for heatmap generation. Install with `conda install matplotlib`
- **numpy**: For numerical operations. Install with `conda install numpy`
- **pandas**: For reading and handling tabular data. Install with `conda install pandas`

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/visual_hyde.git
    ```

2. Install the required dependencies:
    ```bash
    conda install -c etetoolkit ete3 ete_toolchain
    conda install matplotlib numpy pandas
    ```

## Usage

```bash
python visual_hyde.py -i <hyde_output.tsv> -t <species_tree.newick> -z <Z-score threshold>
```

### Command Line Arguments

- `-i` / `--infile`: **(Required)** Name of the HyDe output file (TSV format).
- `-t` / `--treefile`: **(Required)** Name of the phylogenetic tree file (Newick format).
- `-n` / `--node`: **(Optional)** Run in node mode. This stacks up heatmaps for each monophyletic clade.
- `-c` / `--preclade`: **(Optional)** Specify a predefined clade file. If not provided, the script will auto-define clades with up to five leaves.
- `-l` / `--leaves`: **(Optional)** Specify the name of a leaf for which to generate a heatmap. If not specified, the script generates heatmaps for all leaves.
- `-s` / `--picturesize`: **(Optional)** Specify the size of the output figure. Default is 4000.
- `-z` / `--zscore`: **(Optional)** Set a Z-score threshold for filtering HyDe output. Default is 3.

### Example

#### Leaf Mode (Default)

Generate heatmaps for all leaves based on the HyDe output:

```bash
python visual_hyde.py -i hyde_output.tsv -t species_tree.newick -z 3
```

Generate a heatmap for a specific leaf:

```bash
python visual_hyde.py -i hyde_output.tsv -t species_tree.newick -l SpeciesA -z 3
```

#### Node Mode

Run in node mode to generate heatmaps showing shared hybridization events for clades:

```bash
python visual_hyde.py -i hyde_output.tsv -t species_tree.newick -n -z 3
```

## Output

The script generates:
1. **Heatmaps**: Shows hybridization signals for each leaf or node. The heatmaps are saved as PNG files.


## Visualization Explanation

- **Leaves Model**: For each leaf, the x and y axes represent potential parental species, and each cell in the heatmap reflects the strength of hybridization (based on Z-score or Gamma value).
  
- **Node Model**: The heatmap highlights hybridization events shared across all taxa in a node, providing a summary of common hybridization signals in a clade.

## License

This project is licensed under the MIT License.

## Contact

For any issues or questions, please contact Jian He at j.he930724@gmail.com.
