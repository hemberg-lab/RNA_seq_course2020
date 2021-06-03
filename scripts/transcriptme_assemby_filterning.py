import glob
import gzip
from collections import defaultdict 
import time
import numpy as np
import pybedtools

import sys
import csv

from Bio import SeqIO
from Bio.Seq import Seq


csv.field_size_limit(100000000)
csv.field_size_limit()




def Genomictabulator(fasta):



	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq



	f.close()


def main(genome_fasta, extended_ref_annotation, repeat_masker, black_list, out_gtf):

    new_transcript_exons_info = dict()
    new_transcript_exons = set()

    transcript_n_exons =  defaultdict(int)

    #with open("/lustre/scratch117/cellgen/team218/gp7/Joe/RNA_seq_snakepipes/gffcompare/extended_ref_annotation.gtf") as gtf:
    with open(extended_ref_annotation) as gtf:


        reader = csv.reader(gtf, delimiter="\t")

        exon_transcripts = defaultdict(list)



        annotated_exon = set()

        for row in reader:

            chrom = row[0]
            start = row[3]
            end = row[4]
            source = row[1]
            feature = row[2]
            strand = row[6]
            tags = row[8].strip(" ").split(";")

            info =  dict()

            info["chrom"] = chrom
            info["start"] =  start
            info["strand"] = strand
            info["end"] = end
            
            for t in tags:
                pair =  t.strip(" ").split(' "')
                if pair!=['']:
                    if len(pair)!=2:
                        print(row)
                    ID_type, ID  = pair

                    ID = ID.strip('"') 

                    info[ID_type] = ID
                    if ID_type == "transcript_id":
                        transcript = ID

            if feature == "exon":
                exon = (chrom, strand, start, end)
                exon_transcripts[exon].append(transcript)
                transcript_n_exons[transcript] += 1

                if source=="StringTie":
                    new_transcript_exons.add((chrom, start, end))
                    new_transcript_exons_info[exon] = info       
                else:
                    annotated_exon.add((chrom, start, end))

    annotated_exon_bed = pybedtools.BedTool(list(annotated_exon))
    new_transcript_exons_bed = pybedtools.BedTool(list(new_transcript_exons))

    not_intersecting_exons = set([tuple(x) for x in new_transcript_exons_bed.intersect(annotated_exon_bed, v=True)]) 

    #repeatmasker = pybedtools.BedTool("/lustre/scratch117/cellgen/team218/gp7/Genome/mm10/Tracks/Repeats/repeat_masker.mm10.bed")
    repeatmasker = pybedtools.BedTool(repeat_masker)

    repeat_overlap = new_transcript_exons_bed.window(repeatmasker, w=0).overlap(cols=[2,3,2,3])


    exon_repeat_range = defaultdict(set)

    for row in list(repeat_overlap):    


        exon = tuple(row[:3])
        exon_len = int(row[2])-int(row[1])



        e_start= int(row[1])
        e_end = int(row[2])
        r_start = int(row[4])
        r_end = int(row[5])

        ref_start = min(e_start, r_start)

        overlap = range(max(e_start-ref_start, r_start-ref_start), min(e_end-ref_start, r_end-ref_start))

        for i in overlap:
            exon_repeat_range[exon].add(i)


    exon_repeat_size = dict()

    for exon, overlap in exon_repeat_range.items():

        exon_repeat_size[exon] = len(overlap)





    Genome = dict()
    #Genomictabulator("/lustre/scratch117/cellgen/team218/gp7/Genome/GRCm38/GRCm38.fa")
    Genomictabulator(genome_fasta)


    #with open("/lustre/scratch117/cellgen/team218/gp7/Joe/RNA_seq_snakepipes/gffcompare/extended_ref_annotation.exon_black_list.tsv", "w") as out:
    with open(black_list) as out:
        writer = csv.writer(out, delimiter="\t")

        writer.writerow(["exon_coords", "exon_len",  "A_index", "strand", "first_exon", "last_exon", "max_poly(A)_len", "repeat_fraction", "exon_seq"])

        exon_black_list = set()

        for exon in new_transcript_exons_info:



            exon_info = new_transcript_exons_info[exon]

            exon_coords = " ".join((exon_info["chrom"], str(int(exon_info["start"])-1), exon_info["end"]))

            if (exon_info["chrom"], exon_info["start"], exon_info["end"]) in not_intersecting_exons:


                exon_seq =  Genome[exon_info["chrom"]][int(exon_info["start"])-1:int(exon_info["end"])]

                if exon_info["strand"] ==  "-":

                    exon_seq = exon_seq.reverse_complement()

                A_index = str(exon_seq).upper().count("A")/len(exon_seq)
                exon_seq = str(exon_seq).upper()

                x = "A"
                while x in exon_seq:
                    x +="A"
                    
                poliA_len = len(x)
                N_start_A = exon_seq[:10].count("A")
                n_repeats = exon_repeat_size.get((exon_info["chrom"], exon_info["start"], exon_info["end"]), 0)
                repeat_fraction =  n_repeats/len(exon_seq)

                first_exon = False
                last_exon = False


                if exon_info["strand"]=="-" and int(exon_info["exon_number"])==1:
                    last_exon = True
                elif exon_info["strand"]=="+" and int(exon_info["exon_number"])==transcript_n_exons[exon_info["transcript_id"]]:
                    last_exon = True
                elif exon_info["strand"]=="+" and int(exon_info["exon_number"])==1:
                    first_exon = True
                elif exon_info["strand"]=="-" and int(exon_info["exon_number"])==transcript_n_exons[exon_info["transcript_id"]]:
                    first_exon = True

                if len(exon_seq) <30:
                    exon_black_list.add(exon)

                if (first_exon or last_exon) and repeat_fraction>(1/3):
                    exon_black_list.add(exon)

                if last_exon and (A_index>0.5 or N_start_A>=8):
                    exon_black_list.add(exon)



                if exon in exon_black_list:
                    writer.writerow([exon_coords, len(exon_seq),  A_index, exon_info["strand"], first_exon, last_exon, exon_info["exon_number"], poliA_len, repeat_fraction, exon_seq])


    transcript_black_list = set([])

    for exon in exon_black_list:
        for t in exon_transcripts[exon]:
            transcript_black_list.add(t)


	#with open("/lustre/scratch117/cellgen/team218/gp7/Joe/RNA_seq_snakepipes/gffcompare/extended_ref_annotation.gtf") as gtf, \
	#open("/lustre/scratch117/cellgen/team218/gp7/Joe/RNA_seq_snakepipes/gffcompare/extended_ref_annotation.manual_filter.gtf", "w") as out:
		
    with open(extended_ref_annotation) as gtf, open(out_gtf, "w") as out:


	    reader = csv.reader(gtf, delimiter="\t")
	    writer = csv.writer(out, delimiter="\t")
	    exon_transcripts = defaultdict(list)
	    annotated_exon = set()

        for row in reader:

            chrom = row[0]
            start = row[3]
            end = row[4]
            source = row[1]
            feature = row[2]
            strand = row[6]
            tags = row[8].strip(" ").split(";")

            info =  dict()  

            for t in tags:
                pair =  t.strip(" ").split(' "')
                if pair!=['']:
                ID_type, ID  = pair

                ID = ID.strip('"') 

                info[ID_type] = ID
                if ID_type == "transcript_id":
                    transcript = ID

            #n_exons = transcript_n_exons[transcript]

            #if transcript not in transcript_black_list and n_exons>=3:
            if transcript not in transcript_black_list:
                writer.writerow(row)

#genome_fasta, extended_ref_annotation, repeat_masker, black_list, out_gtf
		
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
