# VisualHyde
## Dependancies: 
  python3 (Recommended Anaconda)
  
  ete3 (http://etetoolkit.org/;  install: conda install -c etetoolkit ete3 ete_toolchain)
  
  matplotlib (Install: conda install matplotlib)
  
  numpy (Install: conda install numpy)
  
  pandas (Install: conda install pandas)
  

## Description:

  HyDe (Blischak et al., 2018) is a python package to detect hybridization using phylogenetic invariants based on the coalescent model. Similar to Patterson’s D-statistic (Patterson et al. 2012), HyDe considers a rooted, 4-taxon tree “((P1, H, P2), O)” that consisting of an outgroup “O” and a triplet of ingroup. The hybrid taxa “H” is modeled as a hybridization between two parental taxa, “P1” and “P2”. According to Meng and Davey’s hybrid model,the hybrid taxa is either sister to “P1” with probability (γ) or sister to “P2” with probability (1–γ), The null hypothesis was that when hybridization was absent, the expected value of γ should be 0. HyDe perform a formal statistical test of γ = 0 versus γ > 0 using Z-test. The higher the Z-score, the more reliable of a hybrid event.

  To visualize the HyDe results, we have drawn a heatmap for each potential hybrid sample “H” (Figure a). For the heatmaps, the parents (“P1” and “P2”) were placed in the heatmap on “x” or “y” axes respectively. all the taxa on the axes were arranged phylogenetically. Each square in the heatmap was represents a HyDe test which preform in a 4-taxon tree “((P1, H, P2), O)”. The different color of the squares corresponds different "γ" values that represent genetic contribution of sample on the left axes. 
  
  However, when a HyDe test shows that the "H" was detected as a hybrid of "P1"   and "P2", it does not only represent a recent hybridization event that related to these three samples, but also were possible represent the ancient hybridization event at the deeper nodes. Specifically, if the deeper nodes before of the "P1" or "P2" was involved in a hybridization event, then all the descendants after the nodes should be able to detect the same hybridization event too. Therefore, in a heatmap, if some squares with similar colors formed a larger one in the heatmap, and their corresponding taxa on the axis form a monophyletic clade, then the most recent common ancestor (MRCA) of the monophyletic clade may be indicated as the true parent of "H" (Figure a). Further, if the MRCA of a clade was origin from a hybridization event, then all the taxa on this clade should also detect the same hybridization event signal. Therefore, for each clade, we stacked all the heatmaps which was drawn for each sample in the clade (Figure b). For the heatmap after stacking, only if all the stacked squares were light, the corresponding square would be light. In other words, the light squares in the stacked heatmap means that the hybridization event occurred on the MRCA of these samples which involved in stacking (Figure c; d)

 ![Image text](https://github.com/Jhe1004/VisualHyde/blob/main/figure_1.jpg)

## Useage:

python3 visual_hyde.py -h
