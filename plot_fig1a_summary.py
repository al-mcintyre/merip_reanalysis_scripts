import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

print 'plotting summary for Fig 1a'
levels,fracs = [],[]
for fi in glob.glob('fig1/expression_vs_m6A_detection_*.txt'):
    with open(fi,'r') as infi:
        done = False
        for line in infi:
            level,meth,tot,frac = line.strip().split()
            levels.append(float(level))
            fracs.append(float(frac))
            if level == '0' and not done:
                done = True
            elif level == '0' and done:
                break
  
df = pd.DataFrame({'Expression':levels, '% methylated':fracs}) 

sns.boxplot(data=df,x='Expression',y='% methylated',color='#89043d')
plt.tight_layout()
plt.savefig('fig1/fig1a_expression_vs_m6A_detection_summary.pdf')
