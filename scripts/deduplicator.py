from snakemake.utils import min_version
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import gzip
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
import numpy as np



#file = "/lustre/scratch117/cellgen/team218/gp7/Joe/Test/E6.5-1_i1-AAACTGTGCTACGAGC.fastq.gz"
UMI_flag =  snakemake.params["UMI_flag"]
file =  snakemake.input[0] 
out = gzip.open( snakemake.output[0] , "wt")
#out2 = open( snakemake.output[1] , "wt")

total_reads = []

UMI_dict = defaultdict(list)

with gzip.open(file, mode="rt") as fastq:
    
    for record in SeqIO.parse(fastq, "fastq"):

        UMI = ""
        for tag in record.id.split(";"):
            if tag.split(':')[0]==UMI_flag:
                UMI = tag.split(':')[1]
                
        min15_mean_phred = np.mean(sorted(record.letter_annotations["phred_quality"])[:15])
                
        total_reads.append( (UMI, record.id, record.seq, record.letter_annotations["phred_quality"], min15_mean_phred ) )
        
        
for read in sorted(total_reads, key=lambda x: x[4], reverse=True):
    
    UMI, record_id, record_seq, record_letter_annotations, min15_mean_phred = read
        
    info = (record_id, record_seq, record_letter_annotations)

    
    if UMI in UMI_dict:

        indexs = []
        new = True

        for i in UMI_dict[UMI]:

            identity_index = pairwise2.align.globalms(record_seq, i[1], 1, -1, -1, -0.5, score_only=True)/len(record_seq)
            if identity_index >= 0.95:

                new = False
                break

        if new:
            UMI_dict[UMI].append(info)        

    else:
        UMI_dict[UMI].append(info)
            

count = 0

for umi in UMI_dict:
    
    for info in UMI_dict[umi]:
        ID, seq, qual = info

        read = SeqRecord( seq , id = ID, description = "" )
        read.letter_annotations["phred_quality"] = qual   
        out.write(read.format("fastq"))
        
        count += 1
        
#out2.write(str(count))        
        
out.close()
#out2.close()
