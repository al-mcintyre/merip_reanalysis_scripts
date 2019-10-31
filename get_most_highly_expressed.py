import sys
import re
import numpy as np

genedict = {}
kallistout = sys.argv[1]
label = kallistout.split('_kallisto')[0]
read_length = float(sys.argv[2])
get_full_exons = sys.argv[3]
gtffi=sys.argv[4]
gtf=open(gtffi,'r')
expression_thresh=int(sys.argv[5])

transcript_dict = {}
transcript_gtf = {}
transcript_starts = {}

transcripts_over_thresh = set()

for line in gtf:
   tabbed = line.split('\t')
   if line[0] != '#' and len(tabbed) > 7:
      if tabbed[2] == 'transcript':
         try:
            gene = re.search('(?<=gene_name ")(.*?)(?=")',tabbed[8]).group(0)
            transcript = re.search('(?<=transcript_id ")(.*?)(?=")',tabbed[8]).group(0)
            if gene not in genedict:
                genedict[gene] = {}
                genedict[gene]['highest_t'] = ''
                genedict[gene]['highest_nreads'] = 0
                genedict[gene]['highest_tpm'] = 0
            ts = transcript.split('.')[0]
            transcript_dict[ts] = gene
         except AttributeError:
            print 'ERROR',line
      elif get_full_exons: 
          if tabbed[2] == 'exon':
              if ts not in transcript_gtf:
                  transcript_gtf[ts] = []
              transcript_gtf[ts].append(line)
      else:
         if tabbed[2] == 'start_codon':
            transcript_starts[ts] = int(tabbed[3])
         elif tabbed[2] in set(['CDS','UTR']):
            ts = re.search('(?<=transcript_id ")(.*?)(?=")',tabbed[8]).group(0).split('.')[0]
            if ts not in transcript_gtf:
                transcript_gtf[ts] = []
            transcript_gtf[ts].append(line)        

gtf.close()

highest_expressed = {}
not_found = 0
found = 0
with open(kallistout,'r') as infi:
    for line in infi:
      if line[:4] != 'targ':
        try:
            transcript,length,nreads,tpm = line.split()[0],int(line.split()[1]), float(line.split()[3]), float(line.strip().split()[4]) #or is TPM at -1?
        except ValueError:
            continue
        try:
            transcript = transcript.split('|')[0].split('.')[0]
            gene = transcript_dict[transcript]
            if tpm > genedict[gene]['highest_tpm']:
                genedict[gene]['highest_tpm'] = tpm
                genedict[gene]['length'] = length
                genedict[gene]['highest_t'] = transcript
                genedict[gene]['highest_nreads'] = nreads
            elif genedict[gene]['highest_tpm'] == 0 and nreads > genedict[gene]['highest_nreads']:
                #print nreads, line
                genedict[gene]['highest_t'] = transcript
                genedict[gene]['highest_nreads'] = nreads
            if nreads*read_length/length >= expression_thresh:
                transcripts_over_thresh.add(transcript)

            found +=1
        except KeyError:
            not_found +=1

print not_found, 'transcripts not found;', found, 'found'

gtfall = open(label+'_highest_expression_transcripts_all.gtf','w')
outfi = open(label+'_highest_expression_transcripts.txt','w')
gtffi = open(label+'_highest_expression_transcripts_over'+str(expression_thresh)+'.gtf','w')

with open(label+'_transcripts_over'+str(expression_thresh)+'.gtf','w') as outfi1:
    for ts in transcripts_over_thresh:
        if ts in transcript_gtf:
            for line in transcript_gtf[ts]:
                outfi1.write(line)

for gene in genedict:
    if genedict[gene]['highest_tpm'] > 0:
        ts = genedict[gene]['highest_t']
        outfi.write('{}\t{}\t{}\t{}\n'.format(gene,ts,genedict[gene]['highest_t'],genedict[gene]['highest_nreads'],genedict[gene]['highest_tpm']))
        #if genedict[gene]['highest_nreads'] > 500:
        try:
           for line in transcript_gtf[ts]:
               tabbed = line.split('\t')
               strand = tabbed[6]
               region = tabbed[2]
               start = int(tabbed[3])
               if region == 'UTR' and not get_full_exons:
                   if (strand == '-' and start < transcript_starts[ts]) or (strand == '+' and start > transcript_starts[ts]):
                        line = re.sub('UTR','3_UTR',line)
                   else:   
                        line = re.sub('UTR','5_UTR',line)
               line = line.strip()[:-1] + '\t.\t{}\t{}\t{}\n'.format(genedict[gene]['length'],genedict[gene]['highest_nreads'], genedict[gene]['highest_tpm'])
               if genedict[gene]['highest_nreads'] > 50:
                   gtffi.write(line)
               gtfall.write(line)
        except KeyError:
                pass

outfi.close()
gtffi.close()
gtfall.close()
