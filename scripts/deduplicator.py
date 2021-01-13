#file = "/lustre/scratch117/cellgen/team218/gp7/Joe/Test/E6.5-1_i1-AAACTGTGCTACGAGC.fastq.gz"
UMI_flag = "RX"

file = "/nfs/team205/jdj1/AS/bam_noribo/fastq/E8.5-10_i22-AAACAAACAGGCAGTT.fastq.gz"        
out = gzip.open("/lustre/scratch117/cellgen/team218/gp7/Joe/Test/E8.5-10_i22-AAACAAACAGGCAGTT..deduplicated.fastq.gz", "wt")

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
            


for umi in UMI_dict:
    
    for info in UMI_dict[umi]:
        ID, seq, qual = info

        read = SeqRecord( seq , id = ID, description = "" )
        read.letter_annotations["phred_quality"] = qual   
        out.write(record.format("fastq"))
        
out.close()
