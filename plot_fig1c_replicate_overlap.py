#run after taking bedtools intersects with analyze_replicates.1.sh 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

direc = sys.argv[1]
infi = direc + '/macs2_results/' + sys.argv[2] + '_replicate_intersect.bed'
col='#89043d'
histox = []

with open(infi,'r') as peaks:
    for peak in peaks:
        histox.append(sum([1 for x in peak.split('\t')[3:] if int(x.strip()) > 0]))
    n_max = len(peak.split('\t'))-3

plt.figure(figsize=(4,2.5))
sns.set_style('white')
sns.distplot(histox,bins=np.linspace(0.5,n_max+0.5,n_max+1),kde=False,color=col)
plt.xlabel('replicate overlap')
plt.ylabel('# peaks')
plt.xlim([0.5,n_max+0.5])
plt.tight_layout()
plt.savefig('fig1/replicate_overlap_histogram_'+direc+'.pdf',dpi=500,bbox_inches='tight',transparent=True)

plt.clf()

plt.figure(figsize=(4,2.5))
perc_total = []
x_range = range(1,n_max+1)
for i in x_range:
    perc_total.append(sum([1 for h in histox if h >= i])*100./len(histox))

if perc_total[-1] == 0:
    x_range = x_range[:-1]
    perc_total = perc_total[:-1]
    n_max = n_max -1

print x_range
print perc_total

plt.bar(x_range,perc_total,color=col, align='center',width=1,alpha=0.7)
plt.ylabel('% of total peaks\nin >= r replicates')
plt.xlabel('r')
plt.xlim([0.5,n_max+0.5])
plt.ylim([0,100])
plt.tight_layout()
plt.savefig('fig1/fig1c_replicate_percent_total_'+direc+'.pdf',dpi=500,bbox_inches='tight',transparent=True)

