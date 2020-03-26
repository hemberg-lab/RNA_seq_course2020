import csv, sys
from collections import defaultdict

def main(count_list):
	
	gene_count = defaultdict(list)
	
	for file in count_list:
		
		
		with open(file) as F:
			reader = csv.DictReader(filter(lambda row: row[0]!='#', F), delimiter="\t")
			
			for row in reader:
				
				samples = []
				
				for key in row.keys():
					if 'hisat2' in key:
						samples.append( key )
						
			sample_counts = [(x.split("/")[1].split(".")[0], row[x]) for x in samples]
			
			for i in sample_counts:
				 gene_count[row["Geneid"]].append(i)
			
			
	sample_IDs = [x[0] for x in gene_count[list(gene_count.keys())[0]]]
	
	header = ["gene"] + sample_IDs
	print( "\t".join(header) )
	
	for gene, sample_count_list in  gene_count.items():
		
		counts =  [x[1] for x in sample_count_list]
		out = [gene] + counts
		print( "\t".join(out) )
	

if __name__ == '__main__':
	main(sys.argv[1:])
