#!/users/bleema/anaconda3/bin/python
import numpy as np
import pandas as pd
from matplotlib import rc, font_manager
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import copy
import os


# sysnames
baseline = str(input('What is the baseline condition (e.g. 1A)?'))
condition = str(input('What is the test condition (e.g. 2A)?'))

# set directories for file inputs/outputs, and for genes table
workDir = '.'
genesDir = '../../../../genesets/syk-6/'


#################################################################################
##### sort all.poolcount according to gene and elminate non-qualified genes #####
################################################################################# 

# import all.poolcount for given exp (BarSeqR output)
allCounts = np.loadtxt(workDir+'/all.poolcount', skiprows=0, delimiter='\t', dtype=str) 


# use pandas to make into dataframe; rename columns to shorter sample name
# note that you will need to build a "poolcount_column_names" text file separately, with titles for each experiment column 
allCounts = allCounts[1:len(allCounts),:]
new_cols = pd.read_table(workDir+'/poolcount_column_names.txt', sep='\t') #rename columns to shorter abbrev's
col_mapping = [f"{c[1]}" for c in enumerate(new_cols.columns)] 
cols = np.asarray(col_mapping)
aCobj = pd.DataFrame(allCounts, columns=cols)


# convert numbered strings to numeric values in dataframe
numcols = cols[4:len(cols)]
del1 = np.where(numcols==('locusId'))
numcols = np.delete(numcols, del1)
del2 = np.where(numcols==('f'))
numcols = np.delete(numcols, del2)
for f in numcols:
    aCobj[f] = pd.to_numeric(aCobj[f])


# use only columns that are defined in the set (baseline/condition)
aCobj_trim1 = aCobj[['barcode','rcbarcode','scaffold','strand','pos','locusId','f',baseline]]
aCobj_trim2 = aCobj[condition]
frames = [aCobj_trim1, aCobj_trim2,]
dataFrame = pd.concat(frames, axis=1, join='inner')


#select only rows that are within genes (ingenes); write output to csv    
ingenes = dataFrame[dataFrame['locusId'] != '']
#ingenes.to_csv(r'counts-in-genes.csv', index=False, header=True)
    
    
# select only rows that have >= 3 reads/strain in T0 (in3genes); write output to csv
in3genes = ingenes[ingenes[baseline] >= 3]
in3genes.to_csv(r''+workDir+'/fitness/'+baseline+'_in3genes.csv', index=False, header=True)


# reload, group strains as list of arrays where each array corresponds to one gene (locusID)
raw = np.loadtxt(workDir+'/fitness/'+baseline+'_in3genes.csv', skiprows=0, delimiter=',', dtype=str)
rawLabels = raw[0,:] 
size=raw.shape[1]
raw =raw[1:len(raw),5:size]
groupedByGene = [raw[raw[:,0]==k] for k in raw[:,0][np.sort(np.unique(raw[:,0], return_index=True)[1])]]
    

# ID genes that have <30 reads per gene in Time0
lowSums = []
for i in range(0,len(groupedByGene)):
#    if float(sum(groupedByGene[i][:,2]) < float(30)):
    if float(sum((groupedByGene[i][:,2]).astype('float64')) < float(30)):
        lowSums.append(groupedByGene[i][0,0])
        
print(("the following genes contain insufficient reads (<30) in baseline and will be removed:",lowSums))


# remove genes that have <30 reads per gene in Time0; write genes used to GENES_USED_FOR_ANALYSIS.csv
remove = []
for k in lowSums:
    for j in range(0,len(groupedByGene)):
        if groupedByGene[j][0,0] == k:
            remove.append(j)

GenesUsed = np.delete(groupedByGene, remove)

usedForAnalysis = []
for p in range(0,len(GenesUsed)):
    usedForAnalysis.append(GenesUsed[p][0,0])

np.savetxt(workDir+'/fitness/'+baseline+'_GenesUsedforAnalysis.csv', usedForAnalysis, delimiter=',', fmt='%s')
 


######################################################################
##### now perform fitness calcs for the desired exp versus Time0 #####
######################################################################   
   
# write list of experiments to be evaluated
expList = np.ndarray.tolist(rawLabels[8:len(rawLabels)])
print("The experiment to be evaluated is:", expList)
expIndex = int(3)
BLIndex = int(2)


# obtain counts for T0 and Exp
BLCounts = []
strainCounts = []
strainLoc = []
BLSums = [] 
strainSums = []
for u in range(0,len(GenesUsed)):
    BLCounts.append((GenesUsed[u][:,BLIndex]).astype('float64'))
    strainCounts.append((GenesUsed[u][:,expIndex]).astype('float64'))
    strainLoc.append((GenesUsed[u][:,1]).astype('float64'))
    BLSums.append(sum((GenesUsed[u][:,BLIndex]).astype('float64')))
    strainSums.append(sum((GenesUsed[u][:,expIndex]).astype('float64')))
    
BLSums = np.asarray(BLSums)
strainSums = np.asarray(strainSums)


# eliminate genes if strainSums = 0
remove2 = []
for i in range(0,len(strainSums)):
        if strainSums[i] == 0: remove2.append([i])

GenesUsed = np.delete(GenesUsed, remove2)
BLCounts = np.delete(BLCounts, remove2)
strainCounts = np.delete(strainCounts, remove2)
strainLoc = np.delete(strainLoc, remove2)
BLSums = np.delete(BLSums, remove2)
strainSums = np.delete(strainSums, remove2)

usedForAnalysis = []
for p in range(0,len(GenesUsed)):
    usedForAnalysis.append(GenesUsed[p][0,0])

#np.savetxt('GenesUsedforAnalysis_noZeroSums.csv', usedForAnalysis, delimiter=',', fmt='%s')

print(len(remove2), "additional genes contain zero reads in the pool and will be removed.")


# perform preliminary strain/gene fitness calculations and set median prelim gene fitness to 0
readRatio = strainSums/BLSums
prelimStrainF = []
prelimGeneF = []
for u in range(0,len(GenesUsed)):
    psf = np.log2(strainCounts[u] + np.sqrt(readRatio[u])) - np.log2(BLCounts[u] + 1/(np.sqrt(readRatio[u])))
    prelimStrainF.append(psf)
    pgf = np.median(psf)
    prelimGeneF.append(pgf)

prelimGeneF = np.asarray(prelimGeneF)
prelimGeneFnorm = prelimGeneF - np.median(prelimGeneF) 
#NOTE: prelimGeneF, NOT prelimGeneFnorm, is  used in the calculation of psi below. This is in accordance with the psi calculation in FEBA.R code from Wetmore et al. 


# determine new strain fitness values with pseudocount (psi in Wetmore et al. calcs) 
# psi = (2^prelimGeneF)*readRatio for genes with >=3 strains. otherwise, psi = readRatio
# avoids noise in low-count estimates
psi = []
for u in range(0,len(GenesUsed)):
    if len(GenesUsed[u]) >= 3:
        psi1 = (2**prelimGeneF[u])*readRatio[u]
    else:
        psi1 = readRatio[u]
    psi.append(psi1)

pseudoCount = np.asarray(psi)

strainPseudoCount = np.sqrt(pseudoCount)
BL_PseudoCount = 1/np.sqrt(pseudoCount)

strainFitness = []
for u in range(0,len(GenesUsed)):
    sf = np.log2(strainPseudoCount[u] + strainCounts[u]) - np.log2(BL_PseudoCount[u] + BLCounts[u])
    strainFitness.append(sf)
    
    
# calculate strain variance, use it to weight gene fitness calc's (weight inversely prop. to strain variance)
strainVar = []                          
for u in range(0,len(GenesUsed)):
    sv = ((1/(1+strainCounts[u])) + (1/(1+BLCounts[u])))/((np.log(2))**2)               # naive strain variance
    strainVar.append(sv)

maxWt = ((2/21)/((np.log(2))**2))**(-1)     # strains with low variance are favored, but not too much

strainWt = []
sw2 = np.array([])
for u in range(0,len(GenesUsed)):
    for i in range(0,len(GenesUsed[u])):
        sw1 = min(maxWt, (1/(strainVar[u][i])))
        sw2 = np.append(sw2, [sw1])
    strainWt.append(sw2)
    sw2 = np.array([])      # clears after each i loop since want a new array for each gene


# calculate unnormalized gene fitness, uGeneFitness. save as .csv
uGeneFitness = []
num = np.array([])
denom = np.array([])
for u in range(0,len(GenesUsed)):
    for i in range(0,len(GenesUsed[u])):
        m1 = strainWt[u][i]*strainFitness[u][i]
        m2 = strainWt[u][i]
        num = np.append(num, [m1])
        denom = np.append(denom, [m2])
    uGeneFitness.append(sum(num)/sum(denom))
    num = np.array([])
    denom = np.array([])

uGeneFitness = np.asarray(uGeneFitness)
#np.savetxt('unnormGeneFit.csv', uGeneFitness, delimiter=',', fmt='%s')



##################################
##### normalize gene fitness #####
##################################

# normalize gene fitness using median of 251-gene window (125 genes on either side of GOI)
normGeneFitness = []

# break this up into three chunks to account for the fact that the genome is circular
normGeneFitness_middle = []
for l in range(125,(len(GenesUsed)-124)):
    m = l-125
    n = l+126
    low = min(m,n)
    high = max(m,n)
    nGF_mid = uGeneFitness[l] - np.median(uGeneFitness[low:high])
    normGeneFitness_middle.append(nGF_mid)
    
normGeneFitness_first125 = []   
ceiling = len(GenesUsed)
for l in range(0,125):
    m2 = l-125
    if m2 < 0: m2_new = len(GenesUsed) - 125 - (-l)
    else: m2_new = m
    n2 = l+126
    range1 = uGeneFitness[m2_new:ceiling]
    range2 = uGeneFitness[0:n2]
    full = np.concatenate([range1,range2])
    nGF_first = uGeneFitness[l] - np.median(full)
    normGeneFitness_first125.append(nGF_first)
    
normGeneFitness_last125 = []   
ceiling = len(GenesUsed)
for l in range((len(GenesUsed)-124), len(GenesUsed)):
    m3 = l-125
    n3 = l+126
    if n3 > len(GenesUsed): n3_new = n3 - len(GenesUsed)
    else: n3_new = n3
    range3 = uGeneFitness[m3:ceiling]
    range4 = uGeneFitness[0:n3_new]
    fuller = np.concatenate([range3,range4])
    nGF_last = uGeneFitness[l] - np.median(fuller)
    normGeneFitness_last125.append(nGF_last)
    
# compile the median-smoothed nGFs from each of the three chunks into a single array of normalized gene fitness values    
normGeneFitness = np.concatenate([normGeneFitness_first125, normGeneFitness_middle, normGeneFitness_last125])




###########################################
##### calculate t-like test statistic #####
###########################################
    
# calculate sumSq, weighted sum of squared differences of strain fitness for the gene 
sumSq = []             
numss = np.array([])
denomss = np.array([])
for u in range(0,len(GenesUsed)):
    for i in range(0,len(GenesUsed[u])):
        mm1 = strainWt[u][i]*((strainFitness[u][i]-uGeneFitness[u])**2)
        mm2 = strainWt[u][i]
        numss = np.append(numss, [mm1])
        denomss = np.append(denomss, [mm2])
    sumSq.append(sum(numss)/sum(denomss))
    numss = np.array([])
    denomss = np.array([])  

# for the calculation of typical gene variance, Vt, exclude genes without at least 15 time-zero reads in each half of the gene.
# first, define a list, T0count_Vg, of the time-zero strain counts for each gene
T0count_Vg = []                          
for u in range(0,len(GenesUsed)):
    tzc = GenesUsed[u][:,2].astype(float)              
    T0count_Vg.append(tzc)

# then, sort T0 counts, strain fitness, and strain weights into two lists/each, based on position in gene (1st or 2nd half)
tzc1 = copy.deepcopy(T0count_Vg)
tzc2 = copy.deepcopy(T0count_Vg)
sF1 = copy.deepcopy(strainFitness)
sF2 = copy.deepcopy(strainFitness)
sW1 = copy.deepcopy(strainWt)
sW2 = copy.deepcopy(strainWt)
for u in range(0,len(GenesUsed)):                # this makes each half have the same list sizing as the full gene lists and fills them with zeros
    tzc1[u].fill(0)
    tzc2[u].fill(0)
    sF1[u].fill(0)
    sF2[u].fill(0) 
    sW1[u].fill(0)
    sW2[u].fill(0)                              

for u in range(0,len(GenesUsed)):                # splits T0 counts, strain fitness values, and strain weights into two lists, cutoff at 50% of gene length
    for i in range(0,len(GenesUsed[u])):
        if strainLoc[u][i] <= 0.5: 
            tzc1[u][i] = T0count_Vg[u][i]
            sF1[u][i] = strainFitness[u][i]
            sW1[u][i] = strainWt[u][i]
        else: 
            tzc2[u][i] = T0count_Vg[u][i]
            sF2[u][i] = strainFitness[u][i] 
            sW2[u][i] = strainWt[u][i]

# sum the T0 counts for each half of each gene. 
# If the sum of T0 counts in either half is less than 15, eliminate that gene from the list.
delete_tzc = []
for u in range(0,len(GenesUsed)):
    if np.sum(tzc1[u]) < 15 or np.sum(tzc2[u]) < 15:
        delete_tzc.append(u)

tzc1_trim = np.delete(tzc1, delete_tzc)                # redefine each list, with genes removed if they don't meet the 15-per-half T0 count threshold
tzc2_trim = np.delete(tzc2, delete_tzc)
sF1_trim = np.delete(sF1, delete_tzc)
sF2_trim = np.delete(sF2, delete_tzc)
sW1_trim = np.delete(sW1, delete_tzc)
sW2_trim = np.delete(sW2, delete_tzc)

# calculate uGF1 and uGF2, unnormalized gene fitness values for each gene half, using the same formlua as for uGeneFitness, but including only strain fitness/weights from genes that meet the 15-per-half T0 count threshold
uGF1 = []
uGF2 = []
num1 = np.array([])
denom1 = np.array([])
num2 = np.array([])
denom2 = np.array([])
for u in range(0,len(tzc1_trim)):
    for i in range(0,len(tzc1_trim[u])):
        num1 = np.append(num1, [sW1_trim[u][i]*sF1_trim[u][i]])
        denom1 = np.append(denom1, [sW1_trim[u][i]])
        num2 = np.append(num2, [sW2_trim[u][i]*sF2_trim[u][i]])
        denom2 = np.append(denom2, [sW2_trim[u][i]])
    uGF1.append(sum(num1)/sum(denom1))
    uGF2.append(sum(num2)/sum(denom2))
    num1 = np.array([])
    denom1 = np.array([])
    num2 = np.array([])
    denom2 = np.array([])
    
# in rare instances, it is possible that uGF1 or uGF2 is equal to zero. if true, remove that gene
deleteGF = []
for j in range(0,len(tzc1_trim)):
    if np.isnan(uGF1[j]) or np.isnan(uGF2[j]):
        deleteGF.append(j)

uGF1 = np.delete(uGF1, deleteGF)
uGF2 = np.delete(uGF2, deleteGF)


# calculate typVar (Vt), variance in the typical gene, based on median absolute difference between the two halves
mad12 = np.array([])                     

for b in range(0,len(uGF1)):
    mabsdiff = np.absolute(uGF1[b] - uGF2[b])
    mad12 = np.append(mad12, mabsdiff)                  # median absolute difference between half1, half2

qnorm = 0.674                                           # qnorm(0.75) = 0.674, 75th percentile of normal distribution
typVar = (np.median(mad12)**2)/((2*qnorm)**2)           # variance in typical gene!


# calculate pseudovar (Vg), prior estimate of variance in gene fitness
# from Wetmore et al., Vg = Vt * [Vn/median (Vn)]^2, where the median is over genes used to estimate Vt (uGF1,2)
naiveGeneVar = []
for u in range(0,len(GenesUsed)):
    Vn = ((1/(1+strainSums[u])) + (1/(1+BLSums[u])))/((np.log(2))**2)
    naiveGeneVar.append(Vn)

strainSums_halves = np.delete(strainSums, delete_tzc)         # get count sums for just the genes used to estimate Vt
BL_sums_halves = np.delete(BLSums, delete_tzc)
    
naiveGeneVar_halves = []                                    # calculate Vn for just the genes used to estimate Vt
for b in range(0,len(strainSums_halves)):
    VnH = ((1/(1+strainSums_halves[b])) + (1/(1+BL_sums_halves[b])))/((np.log(2))**2)
    naiveGeneVar_halves.append(VnH)

medVn_halves = np.median(naiveGeneVar_halves)               # take median over just the genes used to estimate Vt

pseudovar = []                                              # calculate pseudovariance, Vg, for each gene, where the denominator is the median taken only over genes used to estimate Vt
for u in range(0,len(GenesUsed)):
    Vg = typVar*((naiveGeneVar[u]/medVn_halves)**2)
    pseudovar.append(Vg)


# calculate estimated variance for each gene. From Wetmore et al., Ve = (sumSq + Vg) / totStrains
totStrains = []
for u in range(0,len(GenesUsed)):
    tS = len(GenesUsed[u])
    totStrains.append(tS)

estGeneVar = []
for u in range(0,len(GenesUsed)):
    Ve = (sumSq[u] + pseudovar[u])/totStrains[u]
    estGeneVar.append(Ve)
    
    
# calculate t-like test stat! 
# From Wetmore et al., t = normGeneFitness/sqrt(sigma^2 + max(Ve,Vn))
# sigma is a small constant to represent uncertainty in normalization for small fitness values. Set to 0.1
sigma = 0.1
tStat = []
for u in range(0,len(GenesUsed)):
    tStat.append(normGeneFitness[u]/(np.sqrt((sigma**2) + max(estGeneVar[u], naiveGeneVar[u]))))

tStat = np.asarray(tStat)



#######################
#### volcano plots ####
#######################

# find lowest normGeneFitness value, assoc gene, and tStat
lowestnGF = min(normGeneFitness)
#lowestIndex = normGeneFitness.index(lowestnGF)
lowestIndex = np.where(normGeneFitness == lowestnGF)[0][0]
print("The gene with the lowest normalized fitness is #", usedForAnalysis[lowestIndex])
print("It has a value of", lowestnGF)
print("and a t-like test stat of", np.absolute(tStat[lowestIndex]))

# plot normd gene fitness vs absolute t score
tStat_abs = np.absolute(tStat)
minT = np.zeros((len(GenesUsed),))
minT.fill(3)

fig2 = plt.figure(figsize=(5,4.7))
#fig2 = plt.figure(figsize=(10,10))
ax2 = fig2.add_subplot(111)
ax2.scatter(normGeneFitness, tStat_abs, c='#57B7C8', marker='o', alpha=0.5, edgecolors='#256874')
ax2.axhline(y=3, linewidth=1, linestyle='--', color='#A9A9A9')
ax2.axvline(x=0, linewidth=1, linestyle='--', color='#A9A9A9')
ax2.set_ylabel('absolute t-score', fontweight='bold')
ax2.set_xlabel('normGeneFitness', fontweight='bold')
ax2.set_ylim([-0.5,(max(tStat_abs)+1)])
#ax2.set_ylim([3,14])
ax2.set_xlim([(min(normGeneFitness)-1),(max(normGeneFitness)+1)])
#ax2.set_xlim([-4.5,2])
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(4))
ax2.set_title((condition+'_vs_'+baseline), fontweight='bold', pad=12)    

# optional annotation on datapoints
#for i, txt in enumerate(usedForAnalysis):
#   ax2.annotate(txt, (normGeneFitness[i], tStat_abs[i])) 

# plot properties
fontProperties = {'family':'sans-serif','sans-serif':['Arial'],
            'weight' : 'normal', 'size' : 14}

ticks_font = font_manager.FontProperties(family='Arial', style='normal',
            size=14, weight='bold', stretch='normal')

rc('font',**fontProperties)

# Uncomment to save plot image
plt.savefig((workDir+'/fitness/'+condition+'_vs_'+baseline+'_volcanoPlot.png'),format='png', dpi=400, bbox_inches='tight', edgecolor='none')




###########################################################################
#### write list of ALL analyzed genes with tStat_abs and normGeneFit value ####
###########################################################################

# metrics for all analyzed genes
genesTable = np.loadtxt(genesDir+'/genesTable_forBarSeqProc_syk6.txt', skiprows=1, delimiter='\t', dtype=str)

fullGenes = []
fullNGF = []
fullTabs = []
y = []

for i in range(0,len(normGeneFitness)):
    if tStat_abs[i] >= 0:
        fullGenes.append(usedForAnalysis[i])
        fullNGF.append(normGeneFitness[i])
        fullTabs.append(tStat_abs[i])

for x in fullGenes:
    if (x in genesTable[:,0]):
        position = np.where(genesTable[:,0] == x)
        y.append(genesTable[position,(7,1)])

meta = np.concatenate(y, axis=0)       
    
# write file with everything
structuredArr = np.transpose(np.array([(fullGenes), (fullNGF), (fullTabs)]))
fullOutput = np.concatenate((structuredArr, meta), axis=1)
np.savetxt((workDir+'/fitness/'+condition+'_vs_'+baseline+'_allAnalyzedGenes.csv'), fullOutput, delimiter=',', fmt='%s', header='number,normGeneFit,tStat_abs,desc,gene', comments='')
