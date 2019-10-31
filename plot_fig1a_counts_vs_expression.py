import re
import sys
import numpy as np
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
import seaborn as sns
import matplotlib.pyplot as plt

def fsigmoid(x, a, b, c):
    return c*1.0 / (1.0 + np.exp(-a*(np.log10(x)-b)))

infi=sys.argv[1]
read_length=int(sys.argv[2])
cond=sys.argv[3]
pe=sys.argv[4]
fragment_length=sys.argv[5]
label=sys.argv[6]

print infi, read_length, cond, pe, fragment_length, label

if pe == "TRUE":
    read_length=min(read_length*2,fragment_length)
scale_by_length = True
direc='fig1'
ts_dict = {}
with open(infi,'r') as gtf:
    for line in gtf:
        tabbed = line.split('\t')
        transcript = re.search('(?<=transcript_id ")(.*?)(?=")',tabbed[8]).group(0)
        length = int(tabbed[-4])
        reads = float(tabbed[-3])
        tpm = float(tabbed[-2])
        count = int(tabbed[-1])

        if transcript not in ts_dict:
            ts_dict[transcript] = {}
            ts_dict[transcript]['counts'] = 0
            ts_dict[transcript]['reads'] = reads
            ts_dict[transcript]['length'] = length
        ts_dict[transcript]['counts'] += count

scale = [0,1,5,10,50,100,500,1000,5000,10000] 
if scale_by_length:
    scale = [0,0.05,0.1,0.5,1,5,10,50,100,500]
expression_levels_reads = {thresh:{x:[] for x in scale} for thresh in range(1,4)}
all_expression = []
all_counts = []
for transcript in ts_dict:
    reads = ts_dict[transcript]['reads']
    counts = ts_dict[transcript]['counts']
    if scale_by_length:
        reads = reads * read_length/ts_dict[transcript]['length'] 
    all_counts.append(counts)
    all_expression.append(reads)
    for thresh in expression_levels_reads:
        if counts >= thresh:
            c = 1
        else:
            c = 0
        for lev in scale:
            if reads <= lev:
                expression_levels_reads[thresh][lev].append(c)
                break

y = {thresh:[] for thresh in expression_levels_reads}
outfi = direc+'/expression_vs_m6A_detection_'+label+'_'+cond+'.txt'
with open(outfi,'w') as out:
    for thresh in expression_levels_reads:
        for lev in scale:
            if len(expression_levels_reads[thresh][lev]) > 0: 
                y[thresh].append(sum(expression_levels_reads[thresh][lev])*100./len(expression_levels_reads[thresh][lev]))
            else:
                y[thresh].append(0)
            out.write("{0}\t{1}\t{2}\t{3}\n".format(lev, str(sum(expression_levels_reads[thresh][lev])), str(len(expression_levels_reads[thresh][lev])), str(y[thresh][-1])))

colours = {1:'#89043d',2:'#726da8',3:'#a0d2db'}
sns.set_style('white')
fig = plt.figure(figsize=(4,4))
ax = plt.gca()
ax.set_xscale('log')
for thresh in y:
    # fit curve (cannot include x = 0 for log10 scale)
    popt, pcov = curve_fit(fsigmoid, scale[1:], y[thresh][1:], bounds=([0.0001,0,0],[np.inf,np.inf,100])) #, method='dogbox', bounds=([0.,max(scale)*2.],[0.,100.]))
    sigx = np.logspace(-3, np.log10(max(scale))*2, 500)
    sigy = fsigmoid(sigx, *popt)
    ax.plot(scale,y[thresh],'o',color=colours[thresh],label=thresh)
    ax.plot(sigx,sigy,color=colours[thresh])
plt.ylim([0,100])
plt.xlim([0.01,max(scale)*2])
plt.legend(title="n")
plt.ylabel('% transcripts with >= n m6A(m) peaks')
if scale_by_length:
    plt.xlabel('mean coverage\n(# reads*read length/transcript length)')
else:
    plt.xlabel('expression level (# reads)')
plt.tight_layout()
plt.savefig(direc+'/fig1a_expression_vs_m6A_detection_'+label+'_'+cond+'.pdf',bbox_inches='tight',transparent=True)

