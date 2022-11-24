# SYK-6_RB-TnSeq
Code repository for "Multiplexed fitness profiling by RB-TnSeq elucidates pathways for lignin-related aromatic catabolism in Sphingobium sp. SYK-6"

BarSeqProc.py: performs strain and gene fitness calculations as described in Wetmore et al., mBio (2015) (doi: 10.1128/mBio.00306-15). This is a custom Python script translated from the original BarSeqR.pl code.

align_allAnalyzedGenes_OneSample_triplicate.py: reads in the outputs of BarSeqProc.py (list of analyzed genes) from each of three biological triplicate fitness experiments. Outputs genes that are common to all replicates, along with their normalized gene fitness scores. Performs a one-sample, two-sided t-test with the null hypothesis that the mean gene fitness of the three replicates is equal to zero. Adjusts p-values with the false discovery rate correction (FDR, Benjamini-Hochberg), outputting a q-value for each gene.

