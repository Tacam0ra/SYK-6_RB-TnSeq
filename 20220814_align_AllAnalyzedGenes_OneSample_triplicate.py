# READS IN THE LISTS OF ANALYZED GENES FROM BIOLOGICAL TRIPLICATE FITNESS EXPERIMENTS
# OUTPUTS GENES THAT ARE COMMON TO ALL REPLICATES, ALONG WITH SCORES
# T-TEST USED IS A ONE SAMPLE, TWO SIDED TEST WITH NULL HYPOTHESIS THAT MEAN OF THREE REPLICATES IS EQUAL TO ZERO
# UPDATE, 2/10/22: ADJUSTS P-VALUES USING THE FDR CORRECTION (B-H PROCEDURE)

import numpy as np
from scipy.stats import ttest_1samp
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from mpl_toolkits.mplot3d import axes3d, Axes3D

# load data from BarSeqProc

array1 = np.loadtxt('./fitness/7A_vs_13A_allAnalyzedGenes.csv', delimiter=',', skiprows=1, usecols=(0,1,2), dtype=str)
array2 = np.loadtxt('./fitness/7B_vs_13B_allAnalyzedGenes.csv', delimiter=',', skiprows=1, usecols=(0,1,2), dtype=str)
array3 = np.loadtxt('./fitness/7C_vs_13C_allAnalyzedGenes.csv', delimiter=',', skiprows=1, usecols=(0,1,2,3), dtype=str)


# identify genes for which all replicates contain data (i.e., gene not eliminated on count basis), then resort in order of gene rank

same1 = np.intersect1d((array1[:,0].astype(str)), (array2[:,0].astype(str)))
same2 = np.intersect1d(same1, (array3[:,0].astype(str)))
same2.sort()
comparison = same2.astype(str)


# build an array that contains gene #, NGF value, t-like test stat, gene identifier, and description for each of the three replicates

sv1 = []
sv2 = []
sv3 = []

for i in comparison:
    if (i in array1[:,0]):
        pos1 = np.where(array1[:,0] == i)
        sv1.append(array1[pos1])
    if (i in array2[:,0]):
        pos2 = np.where(array2[:,0] == i)
        sv2.append(array2[pos2])
    if (i in array3[:,0]):
        pos3 = np.where(array3[:,0] == i)
        sv3.append(array3[pos3])

sv1 = np.concatenate(sv1)
sv2 = np.concatenate(sv2)
sv3 = np.concatenate(sv3)

sharedValues = np.concatenate((sv1,sv2,sv3), axis=1)


# conduct the one-sample, two-sided t-test for significance, using null hypothesis of mu=0

alignedNGFs = []
for i in range(0,len(sharedValues)):
    alignedNGFs.append(np.asarray((sharedValues[i,1], sharedValues[i,4], sharedValues[i,7]), dtype=float))

tscores = []
pvalues = []
for x in alignedNGFs:
    ts, p = ttest_1samp(x, popmean=0)
    tscores.append(ts)
    pvalues.append(p)

joined = np.asarray(list(zip(tscores,pvalues)), dtype=str)

sharedValues_wStats = np.concatenate((joined,sharedValues[:,(0,1,2,4,5,7,8,9)]), axis=1)


# apply the False Discovery Rate (FDR), aka Benjamini-Hochberg procedure
# this allows you to determine the SUBSET of p-values for which you are EVEN MORE CONFIDENT that there is no Type 1 (false positive) error
# formula: fdr_adj_p-val = p_val*len_p_val_list/rank

p_val_list = sharedValues_wStats[:,1:3]        # index the list by position in sharedValues_wStats
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

# plot NGFs as a 3D plot

x = np.asarray(sharedValues_ALLstats[:,4], dtype=float)
y = np.asarray(sharedValues_ALLstats[:,6], dtype=float)
z = np.asarray(sharedValues_ALLstats[:,8], dtype=float)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

scatter = ax.scatter(x, y, z, c='b', alpha=0.5)
#plb.plot(x, p(x), 'm-')
ax.set_xlabel('NGF1')
ax.set_ylabel('NGF2')
ax.set_zlabel('NGF3')
plt.title('NGFs_7-vs-13_1sampleTest')


# save outputs

plt.savefig(('./fitness/stat_analysis/NGFs_7-vs-13_1sampleTest.png'),format='png', dpi=400, bbox_inches='tight', edgecolor='none')
np.savetxt(('./fitness/stat_analysis/7-vs-13_1sampleTest_STATS.csv'), sharedValues_ALLstats[:,(3,0,1,2,4,5,6,7,8,9,10)], delimiter=',', fmt='%s', header='gene, 1samp_qValue, 1samp_tStat, 1samp_pValue, NGF_A, tLike_A, NGF_B, tLike_B, NGF_C, tLike_C, description', comments='')
