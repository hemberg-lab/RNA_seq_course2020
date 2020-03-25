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
							
						sample = key.split("/")[1].split(".")[0]
						samples.append( sample )
							

				
	
	
	
	
	



if __name__ == '__main__':
	main(sys.argv[1:])
