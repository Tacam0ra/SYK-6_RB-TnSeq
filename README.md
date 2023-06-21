# SYK-6_RB-TnSeq
Code repository for "Multiplexed fitness profiling by RB-TnSeq elucidates pathways for lignin-related aromatic catabolism in Sphingobium sp. SYK-6"

DOI: 10.5281/zenodo.8066380

BarSeqProc.py: performs strain and gene fitness calculations as described in Wetmore et al., mBio (2015) (doi: 10.1128/mBio.00306-15). This is a custom Python script translated from the original BarSeqR.pl code.

align_allAnalyzedGenes_OneSample_triplicate.py: reads in the outputs of BarSeqProc.py (list of analyzed genes) from each of three biological triplicate fitness experiments. Outputs genes that are common to all replicates, along with their normalized gene fitness scores. Performs a one-sample, two-sided t-test with the null hypothesis that the mean gene fitness of the three replicates is equal to zero. Adjusts p-values with the false discovery rate correction (FDR, Benjamini-Hochberg), outputting a q-value for each gene.

align_allAnalyzedGenes_TwoSammple_triplicate.py: reads in the normalized gene fitness scores from two sets of align_AllAnalyzedGenes_OneSample_triplicate.py. Outputs genes that are common to each set of three replicates, along with their normalized gene fitness scores. Performs a two-sample, two-sided t-test with the null hypothesis that the mean gene fitness of the two sets of triplicate enrichment conditions (both compared to the same baseline in BarSeqProc.py) are equal. Adjusts p-values with the false discovery rate correction (FDR, Benjamini-Hochberg), outputting a q-value for each gene.

meanNGF_comparisonPlotter_TwoSample_triplicate.py: reads in all output files from align_allAnalyzedGenes_TwoSample_triplicate.py (normalized gene fitness scores, t-like test statistics, and p/q values). Plots of mean normalized gene fitness scores from one set (enrichment1-vs-baseline) on one axis and normalized gene fitness scores from the other set (enrichment2-vs-baseline). Uses q-value to mark significant genes.
