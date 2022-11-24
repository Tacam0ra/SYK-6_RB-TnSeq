# INPUTS: ALL OUTPUT FILES (NGF SCORES, T-LIKE STATS, AND STATS) FROM TWO-SAMPLE ANALYSIS: 
# C:\Users\ableem\Documents\3-PROJECTS\tn-seq\sequencing\SCRIPTS\20220202_align_AllAnalyzedGenes_TwoSample_triplicate.py 
# OUTPUTS: PLOTS OF MEAN NGF SCORE (test-vs-BLA) VS. MEAN NGF SCORE (BLB-vs-BLA)
# UPDATE 2/10/22: USE Q-VALUE INSTEAD OF P-VALUE AS CUTOFF FOR "STARRED" MARKER IN PLOTS


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc, font_manager
import glob


# set directory for file inputs/outputs
workDir = './fitness/stat_analysis'

# load two-sample analysis files
filenames = glob.glob(workDir+'/*_2sampleTest_STATS.csv')


# take the mean of the two baseline NGF scores and the two test NGF scores; output as a single array along with the name of each condition
baseline_meanNGF = []
test_meanNGF = []
baseline_meanTS = []
test_meanTS = []
conditions = []
qValues = []
pValues = []
gene = []
desc = []
for f in filenames:
    print('Now reading:'+f)
    data = np.loadtxt(fname=f, skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), delimiter=',', dtype=float)
    names = np.loadtxt(fname=f, skiprows=1, usecols=(0,16), delimiter=',', dtype=str)
    qValues.append(data[:,0])
    pValues.append(data[:,2])
    gene.append(names[:,0])
    desc.append(names[:,1])
    av_BLB_NGF = data[:,(3,5,7)].mean(axis=1)
    av_test_NGF = data[:,(9,11,13)].mean(axis=1)
    av_BLB_TS = data[:,(4,6,8)].mean(axis=1)
    av_test_TS = data[:,(10,12,14)].mean(axis=1)
    cond_name = f.replace('./fitness/stat_analysis\\','')
    cond_name = cond_name.replace('_2sampleTest_STATS.csv','')
    baseline_meanNGF.append(av_BLB_NGF)
    test_meanNGF.append(av_test_NGF)
    baseline_meanTS.append(av_BLB_TS)
    test_meanTS.append(av_test_TS)
    conditions.append(cond_name)


# construct a list of arrays, which each array contains the baseline NGF mean, test NGF mean, 2-sample p-value, and FDR adjusted p-value (q-value), for a given exp
variablesList = []
for i in range(0,len(baseline_meanNGF),1):
    variablesList.append(np.asarray(list(zip(baseline_meanNGF[i], test_meanNGF[i], baseline_meanTS[i], test_meanTS[i], pValues[i], qValues[i], gene[i], desc[i]))))
    
for i in range(0, len(variablesList), 1):
    np.savetxt((workDir+'/'+conditions[i]+'_means-pVals-qVals.csv'), variablesList[i], delimiter=',', fmt='%s', header='baseline_meanNGF, test_meanNGF, baseline_meanTS, test_meanTS, pValue, qValue, gene, description', comments='')


## split each set of variables into two lists: one for those genes that have p < 0.05, and one for those that hae p > 0.05
#signif = []
#nonsignif = [] 
#for i in range(0, len(variablesList), 1):
#    list1 = []
#    list2 = []
#    for j in range(0, len(variablesList[i]), 1):
#        if (float(variablesList[i][j,4]) < 0.05):
#            list1.append(variablesList[i][j,:])
#        else:
#            list2.append(variablesList[i][j,:])
#    signif.append(np.asarray(list1))
#    nonsignif.append(np.asarray(list2))
#
#for i in range(0, len(signif), 1):      # rearrange so that each list element is a single array
#    signif[i] = np.asarray(signif[i])
#    nonsignif[i] = np.asarray(nonsignif[i])
    
# split each set of variables into two lists: one for those genes that have q < 0.1, and one for those that have q >= 0.1
signif = []
nonsignif = [] 
for i in range(0, len(variablesList), 1):
    list1 = []
    list2 = []
    for j in range(0, len(variablesList[i]), 1):
        if (float(variablesList[i][j,5]) < 0.1):
            list1.append(variablesList[i][j,:])
        else:
            list2.append(variablesList[i][j,:])
    signif.append(np.asarray(list1))
    nonsignif.append(np.asarray(list2))

for i in range(0, len(signif), 1):      # rearrange so that each list element is a single array
    signif[i] = np.asarray(signif[i])
    nonsignif[i] = np.asarray(nonsignif[i])


# plot properties
fontProperties = {'family':'sans-serif','sans-serif':['Arial'],
            'weight' : 'normal', 'size' : 14}

ticks_font = font_manager.FontProperties(family='Arial', style='normal',
            size=14, weight='bold', stretch='normal')

rc('font',**fontProperties)


# plot 
x = np.linspace(-50,50,num=100)
for u in range(0, len(variablesList), 1):
#    bl_signif = signif[u][:,0] # plot significant variables separately from nonsignificant variables
#    test_signif = signif[u][:,1]
#    bl_nonsignif = nonsignif[u][:,0]
#    test_nonsignif = nonsignif[u][:,1]
    title = conditions[u]
    fig = plt.figure(figsize=(4,8))
    
    ax1 = fig.add_subplot(211)
    plt.style.use('seaborn-whitegrid')
    ax1.plot(x, x + 0, '-b', alpha=0.5, label='_nolegend_')  # solid blue line to indicate x=y
    ax1.plot(x, x - 1.5, '--r', alpha=0.5) # dashed red lines to indicate |x-y| > 1.5
    ax1.plot(x, x + 1.5, '--r', alpha=0.5, label='_nolegend_') # dashed red lines to indicate |x-y| > 1.5
    ax1.scatter((signif[u][:,0]).astype(float), (signif[u][:,1]).astype(float), c='#00cc00', s=50, marker='*', alpha=1)
    ax1.scatter((nonsignif[u][:,0]).astype(float), (nonsignif[u][:,1]).astype(float), c='#505050', s=10, alpha=0.5, label='_nolegend_')
    ax1.set_xlim([-6,6])
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(3))
    ax1.set_ylim([-6,6])
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(3))
    ax1.set_xlabel('Average Baseline NGF')
    ax1.set_ylabel('Average Test NGF')
    ax1.set_title(title)
#    ax1.legend(('|x-y| > 1.5','p-value < 0.05'), bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
    ax1.legend(('|x-y| > 1.5','q-value < 0.1'), bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
    
    ax2 = fig.add_subplot(212)
    plt.style.use('seaborn-whitegrid')
    ax2.plot(x, x + 0, '-b', alpha=0.5, label='_nolegend_')  # solid black line to indicate x=y
    ax2.plot(x, x - 5, '--r', alpha=0.5) # dashed red lines to indicate |x-y| > 1.5
    ax2.plot(x, x + 5, '--r', alpha=0.5, label='_nolegend_') # dashed red lines to indicate |x-y| > 1.5
    ax2.scatter((signif[u][:,2]).astype(float), (signif[u][:,3]).astype(float), c='#00cc00', s=50, marker='*', alpha=1)
    ax2.scatter((nonsignif[u][:,2]).astype(float), (nonsignif[u][:,3]).astype(float), c='#505050', s=10, alpha=0.5, label='_nolegend_')
    ax2.set_xlim([0,28])
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(4))
    ax2.set_ylim([0,28])
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(4))
    ax2.set_xlabel('Average Baseline T-like Stat')
    ax2.set_ylabel('Average Test T-like Stat')
    ax2.set_title(title)
#    ax2.legend(('|x-y| > 3','p-value < 0.05'), bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)  
    ax2.legend(('|x-y| > 5','q-value < 0.1'), bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)  
    
    fig.subplots_adjust(hspace=0.4)
    
    plt.savefig(workDir+'/'+title+'_2sample_means-and-qVals_PLOT.png',format='png', dpi=400, bbox_inches='tight', edgecolor='none')

    
    
#    plt.style.use('seaborn-whitegrid')
#    plt.plot(x, x + 0, '-b', alpha=0.5, label='_nolegend_')  # solid blue line to indicate x=y
#    plt.plot(x, x - 1.5, '--r', alpha=0.5) # dashed red lines to indicate |x-y| > 1.5
#    plt.plot(x, x + 1.5, '--r', alpha=0.5, label='_nolegend_') # dashed red lines to indicate |x-y| > 1.5
#    plt.scatter(signif[u][:,0], signif[u][:,1], c='#00cc00', s=50, marker='*', alpha=1)
#    plt.scatter(nonsignif[u][:,0], nonsignif[u][:,1], c='#505050', s=10, alpha=0.5, label='_nolegend_')
#    plt.xlim([-15,6])
#    plt.xticks(np.arange(-15,6,step=3))
#    plt.ylim([-15,6])
#    plt.yticks(np.arange(-15,6,step=3))
#    plt.xlabel('Average Baseline NGF')
#    plt.ylabel('Average Test NGF')
#    plt.title(title)
#    plt.legend(('|x-y| > 1.5','p-value < 0.05'), bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
#    plt.savefig('C:/Users/ableem/Documents/3-PROJECTS/tn-seq/sequencing/NGS_july2021/BarSeq/barseqR/fitness/two_mean_plots/'+title+'_2sample_NGFmeansPlot.png',format='png', dpi=400, bbox_inches='tight', edgecolor='none')
#
#    
#    
