import sys
import csv


def main(genome, gtf_ref, gtf_asembly):
	
	chromosomes = set([])

	f = open(genome)
	
	for chrfa in SeqIO.parse(f, "fasta"):
		chromosomes.add(chrfa.id)		
		
	for row in csv.reader(open(gtf_ref), delimiter = '\t'):
		if row[0] in chromosomes:
			print("\t".join(row))
			
	for row in csv.reader(open(gtf_asembly), delimiter = '\t'):
		if row[0] in chromosomes:
			print("\t".join(row))		
		

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
