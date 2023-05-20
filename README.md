# VisualHyde
## Dependancies: 
  python3
  
  ete3 (http://etetoolkit.org/;  install: conda install -c etetoolkit ete3 ete_toolchain)
  
  matplotlib (install: conda install matplotlib)
  
  numpy (install: conda install numpy)
  
  pandas (install: conda install pandas)
  

## Description:

  HyDe (Blischak et al., 2018) is a python package to detect hybridization using phylogenetic invariants based on the coalescent model. Similar to Patterson’s D-statistic (Patterson et al. 2012), HyDe considers a rooted, 4-taxon tree “((P1, H, P2), O)” that consisting of an outgroup “O” and a triplet of ingroup. The hybrid taxa “H” is modeled as a hybridization between two parental taxa, “P1” and “P2”. According to Meng and Davey’s hybrid model,the hybrid taxa is either sister to “P1” with probability (γ) or sister to “P2” with probability (1–γ), The null hypothesis was that when hybridization was absent, the expected value of γ should be 0. HyDe perform a formal statistical test of γ = 0 versus γ > 0 using Z-test. The higher the Z-score, the more reliable of a hybrid event.

  To visualize the HyDe results, we made a heatmap for each potential hybridization sample (“PH-sample”). In the heatmap, the other samples were assumed as parents (“P1” and “P2”) and placed in the heat map at “x” and “y” axes respectively. the samples in the axes were arranged phylogenetically (species tree were generated by concatenated super-alignment dataset). Thus, each square in the heatmap was represents a HyDe test which preform in a 4-taxon tree “((P1, PH-sample, P2), O)” and lighter yellows indicated higher γ for this test while darker indicate lower. However, what needs to be clearly understood here is that the light squares in heatmap does not necessarily indicate the “PH-sample” was a resent mixture between “P1” and “P2”, it may also indicate older hybridization events that involved the ancestors of the PH-sample, “P1” and “P2” in the phylogeny (Appelhans et al., 2020). In the other words, if the common ancestor of a clade was one of the parents of “PH-sample”, then all taxa in the clade will be successfully detected as one of parents of “PH-sample” and the γ of those tests should be relatively similar. Thus, if a few lighted up squares (with similar lighting) in the heatmap formed a larger square, and this larger square correspond to two monophyletic clades in the “x” and “y” axes respectively, the parents “P1” and “P2” were most likely the last common ancestors of two clades respectively.

  In addition (-n Node model), if the common ancestor of a clade was arisen through a hybridization event, this event should be detected also in all taxa in this clade. Thus, we stack up all the heatmaps for each monophyletic clade respectively. For a certain square in the superimposed heatmap, only the squares in all heatmaps were light, the square after superimposition will be light. Therefore, in the superimposed heatmap, a light square means this hybridization event was shared by all “PH-sample” of this clade.   
