# READS IN THE NGF SCORES FROM OUTPUT OF 20220202_align_AllAnalyzedGenes_OneSample_triplicate.py 
# OUTPUTS GENES THAT ARE COMMON TO THE THREE SETS, ALONG WITH SCORES
# T-TEST USED IS A TWO SAMPLE, TWO SIDED TEST 
# NULL HYPOTHESIS: MEANS OF THREE REPLICATE SETS (TEST + BaselineB, BOTH COMPARED TO BaselineA IN BARSEQ) ARE EQUAL
# UPDATE, 2/10/22: ADJUSTS P-VALUES USING THE FDR CORRECTION (B-H PROCEDURE)

import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import matplotlib.pylab as plb

# load data from BarSeqProc
    
array1 = np.loadtxt('./fitness/stat_analysis/10-vs-1_1sampleTest_STATS.csv', delimiter=',', usecols=(0,4,5,6,7,8,9,10), skiprows=1, dtype=str) # this is your BLB-vs-BLA sample set
array2 = np.loadtxt('./fitness/stat_analysis/8-vs-1_1sampleTest_STATS.csv', delimiter=',', usecols=(0,4,5,6,7,8,9,10), skiprows=1, dtype=str) # this is your test-vs-BLA sample set


# identify genes for which both sets (test-vs-BLA and BLB-vs-BLA) contain data (i.e., gene not eliminated on count basis) 

same = np.intersect1d((array1[:,0].astype(str)), (array2[:,0].astype(str)))
same.sort()
comparison = same.astype(str)


# build an array that contains gene # and NGF values for both sets (test-vs-BLA and BLB-vs-BLA) of three replicates

sv1 = []
sv2 = []

for i in comparison:
    if (i in array1[:,0]):
        pos1 = np.where(array1[:,0] == i)
        sv1.append(array1[pos1])
    if (i in array2[:,0]):
        pos2 = np.where(array2[:,0] == i)
        sv2.append(array2[pos2])

sv1 = np.concatenate(sv1)
sv2 = np.concatenate(sv2)

sharedValues = np.concatenate((sv1[:,0:7],sv2[:,1:8]), axis=1)


# conduct the two-sample, two-sided t-test for significance, using null hypothesis of mu_blb-vs-bla = mu_test-vs-bla

alignedNGFs = []
for i in range(0,len(sharedValues)):
    alignedNGFs.append(np.asarray((sharedValues[i,1], sharedValues[i,3], sharedValues[i,5], sharedValues[i,7], sharedValues[i,9], sharedValues[i,11]), dtype=float))

tscores = []
pvalues = []
for r in range(0,len(alignedNGFs)):
    foo = alignedNGFs[r]
    ts, p = ttest_ind(foo[0:3], foo[3:6], equal_var=True)
    tscores.append(ts)
    pvalues.append(p)

joined = np.asarray(list(zip(tscores,pvalues)), dtype=str)

sharedValues_wStats = np.concatenate((joined, sharedValues), axis=1)


# apply the False Discovery Rate (FDR), aka Benjamini-Hochberg procedure
# this allows you to determine the SUBSET of p-values for which you are EVEN MORE CONFIDENT that there is no Type 1 (false positive) error
# formula: fdr_adj_p-val = p_val*len_p_val_list/rank

p_val_list = sharedValues_wStats[:,(1,2)]         # index the list by position in sharedValues_wStats
p_val_list = p_val_list.tolist()

p_val_list.sort(key=lambda x:x[0])                       # sort the list of p-values from smallest to largest, maintaining position index
len_p_val_list = len(p_val_list)

rank = 1
p_adj_list = []
for p in range(0,len_p_val_list,1):
    fdr_adj_p_val = ((float(p_val_list[p][0]))*len_p_val_list)/rank
    rank += 1
    p_adj_list.append([fdr_adj_p_val, p_val_list[p][0], p_val_list[p][1]])
# put this in the "append" term if you want to make a dictionary:
#            {
#                    "p_val": p_val_list[p][0]
#                   "fdr_adj_p_val": fdr_adj_p_val
#                   "index": p_val_list[p][1]
#            }
#    )

p_adj_list.sort(key=lambda x:x[2])                  # re-sort according to index from sharedValues_wStats

qvalues = np.asarray(p_adj_list)[:,0]               # pull out the index-sorted q-values into their own array

sharedValues_ALLstats = np.concatenate((qvalues[:,None],sharedValues_wStats), axis=1)


# save output 
np.savetxt(('./fitness/stat_analysis/8-vs-10_2sampleTest_STATS.csv'), sharedValues_ALLstats[:,(3,0,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16)], delimiter=',', fmt='%s', header='gene, 2samp_qValue, 2samp_tStat, 2samp_pValue, NGF_BLB-vs-BLA_a, tLike_BLB-vs-BLA_a, NGF_BLB-vs-BLA_b, tLike_BLB-vs-BLA_b, NGF_BLB-vs-BLA_c, tLike_BLB-vs-BLA_c, NGF_test-vs-BLA_a, tLike_test-vs-BLA_a, NGF_test-vs-BLA_b, tLike_test-vs-BLA_b, NGF_test-vs-BLA_c, tLike_test-vs-BLA_c, description', comments='')
